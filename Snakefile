import os.path
import dask
from pathlib import Path
import steps


configfile: "config.yaml"


shell.prefix(config["prefix"])

for k, v in config["dask"].items():
    dask.config.set({k: v})

data_dir = Path(config["data_dir"])


# These are quick steps that we don't want to submit jobs for
localrules:
    all,
    summary_table,


def ds_dir(wildcards):
    if hasattr(wildcards, "region_name"):
        chr = steps.parse_region(config["regions"][wildcards.region_name])[0]
    else:
        chr = wildcards.chrom_num
    return data_dir / "zarr_vcfs" / f"chr{chr}" / "data.zarr"


rule all:
    input:
        expand(
            data_dir / "{subset_name}-{filter_set}-region_summary_table.csv",
            subset_name=config["sample_subsets"].keys(),
            filter_set=config["filters"].keys(),
        ),
        expand(
            data_dir
            / "trees"
            / "{subset_name}-{region_name}-{filter_set}"
            / "{subset_name}-{region_name}-{filter_set}-truncate-{truncation}-mm{mismatch}-post-processed.trees",
            subset_name=config["sample_subsets"].keys(),
            region_name=config["regions"].keys(),
            filter_set=config["filters"].keys(),
            mismatch=config["mismatch_values"],
            truncation=[
                f"{c['lower']}-{c['upper']}-{c['multiplier']}"
                for c in config["truncate"]
            ],
        ),


rule bio2zarr_explode:
    input:
        vcf=lambda wildcards: config["vcf"].format(chrom=wildcards.chrom_num),
        tbi=lambda wildcards: config["vcf"].format(chrom=wildcards.chrom_num) + ".tbi",
    output:
        directory(data_dir / "exploded_vcfs" / "chr{chrom_num}"),
    threads: config["max_threads"]
    resources:
        mem_mb=config["max_mem"],
        time_min=config["max_time"],
        runtime=config["max_time"],
    run:
        from bio2zarr import vcf

        vcf.explode(
            output[0],
            [input.vcf],
            worker_processes=threads,
            column_chunk_size=config["bio2zarr"]["column_chunk_size"],
        )


rule bio2zarr_mkschema:
    input:
        data_dir / "exploded_vcfs" / "chr{chrom_num}",
    output:
        data_dir / "zarr_vcfs_schema" / "chr{chrom_num}" / "schema.json",
    threads: 1
    resources:
        mem_mb=16_000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        from bio2zarr import vcf

        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        with open(output[0], "w") as out:
            vcf.mkschema(
                input[0],
                out,
            )


rule bio2zarr_encode:
    input:
        data_dir / "exploded_vcfs" / "chr{chrom_num}",
        data_dir / "zarr_vcfs_schema" / "chr{chrom_num}" / "schema.json",
    output:
        data_dir / "zarr_vcfs" / "chr{chrom_num}" / "data.zarr" / ".vcf_done",
    threads: config["max_threads"]
    resources:
        mem_mb=config["max_mem"],
        time_min=config["max_time"],
        runtime=config["max_time"],
    run:
        from bio2zarr import vcf

        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        vcf.encode(
            input[0],
            output[0].replace(".vcf_done", ""),
            input[1],
            worker_processes=threads,
            max_memory=config["max_mem"],
        )
        Path(output[0]).touch()


rule load_ancestral_fasta:
    input:
        lambda wildcards: config["ancestral_fasta"].format(chrom=wildcards.chrom_num),
        data_dir / "zarr_vcfs" / "chr{chrom_num}" / "data.zarr" / ".vcf_done",
    output:
        directory(
            data_dir
            / "zarr_vcfs"
            / "chr{chrom_num}"
            / "data.zarr"
            / "variant_ancestral_allele"
        ),
        directory(
            data_dir
            / "zarr_vcfs"
            / "chr{chrom_num}"
            / "data.zarr"
            / "variant_low_quality_ancestral_allele_mask"
        ),
    threads: 1
    resources:
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        steps.load_ancestral_fasta(input, output, wildcards, config, params)


rule pre_subset_filters:
    input:
        data_dir / "zarr_vcfs" / "chr{chrom_num}" / "data.zarr" / ".vcf_done",
        data_dir
        / "zarr_vcfs"
        / "chr{chrom_num}"
        / "data.zarr"
        / "variant_ancestral_allele",
        data_dir
        / "zarr_vcfs"
        / "chr{chrom_num}"
        / "data.zarr"
        / "variant_low_quality_ancestral_allele_mask",
    output:
        directory(
            data_dir
            / "zarr_vcfs"
            / "chr{chrom_num}"
            / "data.zarr"
            / "variant_bad_ancestral_mask"
        ),
        directory(
            data_dir
            / "zarr_vcfs"
            / "chr{chrom_num}"
            / "data.zarr"
            / "variant_duplicate_position_mask"
        ),
        directory(
            data_dir
            / "zarr_vcfs"
            / "chr{chrom_num}"
            / "data.zarr"
            / "variant_not_biallelic_mask"
        ),
        directory(
            data_dir
            / "zarr_vcfs"
            / "chr{chrom_num}"
            / "data.zarr"
            / "variant_not_snps_mask"
        ),
        directory(
            data_dir
            / "zarr_vcfs"
            / "chr{chrom_num}"
            / "data.zarr"
            / "variant_no_ancestral_allele_mask"
        ),
    resources:
        dask_cluster=10,
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        from distributed import Client

        with Client(config["scheduler_address"]):
            steps.pre_subset_filters(input, output, wildcards, config, params)


rule region_mask:
    input:
        data_dir / "zarr_vcfs" / "chr{chrom_num}" / "data.zarr" / ".vcf_done",
    output:
        directory(
            data_dir
            / "zarr_vcfs"
            / "chr{chrom_num}"
            / "data.zarr"
            / "variant_{region_name}_region_mask"
        ),
    resources:
        dask_cluster=10,
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        from distributed import Client

        with Client(config["scheduler_address"]):
            steps.region_mask(input, output, wildcards, config, params)


rule sample_mask:
    input:
        data_dir / "zarr_vcfs" / "chr{chrom_num}" / "data.zarr" / ".vcf_done",
        lambda wildcards: config["sample_subsets"][wildcards.subset_name],
    output:
        directory(
            data_dir
            / "zarr_vcfs"
            / "chr{chrom_num}"
            / "data.zarr"
            / "sample_{subset_name}_subset_mask"
        ),
    resources:
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        steps.sample_mask(input, output, wildcards, config, params)


rule allele_counts:
    input:
        data_dir / "zarr_vcfs" / "chr{chrom_num}" / "data.zarr" / ".vcf_done",
        data_dir
        / "zarr_vcfs"
        / "chr{chrom_num}"
        / "data.zarr"
        / "sample_{subset_name}_subset_mask",
        data_dir
        / "zarr_vcfs"
        / "chr{chrom_num}"
        / "data.zarr"
        / "variant_ancestral_allele",
    output:
        directory(
            data_dir
            / "zarr_vcfs"
            / "chr{chrom_num}"
            / "data.zarr"
            / "variant_{subset_name}_subset_ref_count"
        ),
        directory(
            data_dir
            / "zarr_vcfs"
            / "chr{chrom_num}"
            / "data.zarr"
            / "variant_{subset_name}_subset_ancestral_count"
        ),
        directory(
            data_dir
            / "zarr_vcfs"
            / "chr{chrom_num}"
            / "data.zarr"
            / "variant_{subset_name}_subset_missing_count"
        ),
        directory(
            data_dir
            / "zarr_vcfs"
            / "chr{chrom_num}"
            / "data.zarr"
            / "variant_{subset_name}_subset_derived_count"
        ),
    resources:
        dask_cluster=10,
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        from distributed import Client

        with Client(config["scheduler_address"]):
            steps.allele_counts(input, output, wildcards, config, params)


rule subset_filters:
    input:
        lambda wildcards: ds_dir(wildcards) / ".vcf_done",
        lambda wildcards: expand(
            ds_dir(wildcards) / "{array_name}",
            array_name=[
                "variant_bad_ancestral_mask",
                "variant_no_ancestral_allele_mask",
                "variant_not_biallelic_mask",
                "variant_duplicate_position_mask",
                "variant_not_snps_mask",
                "variant_low_quality_ancestral_allele_mask",
                "variant_{subset_name}_subset_ref_count",
                "variant_{subset_name}_subset_ancestral_count",
                "variant_{subset_name}_subset_missing_count",
                "variant_{subset_name}_subset_derived_count",
                "sample_{subset_name}_subset_mask",
            ],
        ),
    output:
        data_dir
        / "zarr_vcfs"
        / "chr{chrom_num}"
        / "data.zarr"
        / ".{subset_name}_subset_{filter_set}_filters_done",
    resources:
        dask_cluster=10,
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        from distributed import Client

        with Client(config["scheduler_address"]):
            steps.subset_filters(input, output, wildcards, config, params)


rule site_density_mask:
    input:
        lambda wildcards: ds_dir(wildcards) / ".vcf_done",
        lambda wildcards: ds_dir(wildcards)
        / ".{subset_name}_subset_{filter_set}_filters_done",
        lambda wildcards: ds_dir(wildcards) / "variant_{region_name}_region_mask",
    output:
        data_dir
        / "zarr_vcfs"
        / "chr{chrom_num}"
        / "data.zarr"
        / ".{subset_name}_subset_{region_name}_region_{filter_set}_site_density_mask_done",
    resources:
        dask_cluster=10,
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        from distributed import Client

        with Client(config["scheduler_address"]):
            steps.site_density_mask(input, output, wildcards, config, params)


rule combined_mask:
    input:
        lambda wildcards: ds_dir(wildcards) / ".vcf_done",
        lambda wildcards: ds_dir(wildcards)
        / ".{subset_name}_subset_{filter_set}_filters_done",
        lambda wildcards: ds_dir(wildcards)
        / ".{subset_name}_subset_{region_name}_region_{filter_set}_site_density_mask_done",
        lambda wildcards: expand(
            ds_dir(wildcards) / "{array_name}",
            array_name=[
                "variant_bad_ancestral_mask",
                "variant_no_ancestral_allele_mask",
                "variant_not_biallelic_mask",
                "variant_duplicate_position_mask",
                "variant_not_snps_mask",
                "variant_low_quality_ancestral_allele_mask",
                "variant_{region_name}_region_mask",
                "sample_{subset_name}_subset_mask",
            ],
        ),
    output:
        directory(
            data_dir
            / "zarr_vcfs"
            / "chr{chrom_num}"
            / "data.zarr"
            / "variant_{subset_name}_subset_{region_name}_region_{filter_set}_mask"
        ),
    resources:
        dask_cluster=10,
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        from distributed import Client

        with Client(config["scheduler_address"]):
            steps.combined_mask(input, output, wildcards, config, params)


rule zarr_stats:
    input:
        lambda wildcards: ds_dir(wildcards) / ".vcf_done",
        lambda wildcards: ds_dir(wildcards)
        / "variant_{subset_name}_subset_{region_name}_region_{filter_set}_mask",
    output:
        data_dir
        / "zarr_stats"
        / "{subset_name}-{region_name}-{filter_set}"
        / "stats.json",
    resources:
        dask_cluster=10,
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        from distributed import Client

        with Client(config["scheduler_address"]):
            steps.zarr_stats(input, output, wildcards, config, params)


checkpoint summary_table:
    input:
        lambda wildcards: [
            data_dir
            / "zarr_stats"
            / f"{wildcards.subset_name}-{region_name}-{wildcards.filter_set}"
            / "stats.json"
            for region_name in config["regions"].keys()
        ],
    output:
        data_dir / "{subset_name}-{filter_set}-region_summary_table.csv",
    run:
        steps.summary_table(input, output, wildcards, config, params)


def get_ancestor_gen_memory(wildcards):
    import numpy
    import json

    # Use the checkpoint to check the stats file exists
    checkpoint_output = checkpoints.summary_table.get(
        subset_name=wildcards.subset_name, filter_set=wildcards.filter_set
    )
    region_stats = (
        data_dir
        / "zarr_stats"
        / f"{wildcards.subset_name}-{wildcards.region_name}-{wildcards.filter_set}"
        / "stats.json"
    )
    with open(region_stats, "r") as json_stats_f:
        with open(config["sample_subsets"][wildcards.subset_name], "r") as f:
            n_samples = len(numpy.genfromtxt(f, dtype=str))
        stats = json.load(json_stats_f)
        n_ploidy = stats["n_ploidy"]
        n_sites = stats["n_variants"]
        n_masked = stats["sites_masked"]
        mem = 16_000 + int(
            ((n_sites - n_masked) * n_samples * n_ploidy) / (8 * 1_048_576)
        )
        return mem


rule generate_ancestors:
    input:
        lambda wildcards: ds_dir(wildcards) / ".vcf_done",
        lambda wildcards: ds_dir(wildcards)
        / "variant_{subset_name}_subset_{region_name}_region_{filter_set}_mask",
        lambda wildcards: ds_dir(wildcards) / "sample_{subset_name}_subset_mask",
        data_dir
        / "zarr_stats"
        / "{subset_name}-{region_name}-{filter_set}"
        / "stats.json",
    output:
        data_dir
        / "ancestors"
        / "{subset_name}-{region_name}-{filter_set}"
        / "ancestors.zarr",
    threads: config["max_threads"]
    resources:
        mem_mb=get_ancestor_gen_memory,
        time_min=config["max_time"],
        runtime=config["max_time"],
    run:
        steps.generate_ancestors(input, output, wildcards, config, threads)


rule truncate_ancestors:
    input:
        data_dir
        / "ancestors"
        / "{subset_name}-{region_name}-{filter_set}"
        / "ancestors.zarr",
    output:
        data_dir
        / "ancestors"
        / "{subset_name}-{region_name}-{filter_set}"
        / "ancestors-truncate-{lower}-{upper}-{multiplier}.zarr",
    threads: 1
    resources:
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        import tsinfer
        import logging
        import shutil

        lower = float(wildcards.lower)
        upper = float(wildcards.upper)
        multiplier = float(wildcards.multiplier)

        logging.basicConfig(level=logging.INFO)
        if lower == 0 and upper == 0 and multiplier == 0:
            # No truncation, just copy the file using the OS
            shutil.copy(input[0], output[0])
        else:
            ancestors = tsinfer.AncestorData.load(input[0])
            truncated = ancestors.truncate_ancestors(
                lower, upper, length_multiplier=multiplier, path=output[0]
            )


rule match_ancestors:
    input:
        data_dir
        / "ancestors"
        / "{subset_name}-{region_name}-{filter_set}"
        / "ancestors-truncate-{lower}-{upper}-{multiplier}.zarr",
        lambda wildcards: ds_dir(wildcards) / ".vcf_done",
        lambda wildcards: ds_dir(wildcards)
        / "variant_{subset_name}_subset_{region_name}_region_{filter_set}_mask",
        lambda wildcards: ds_dir(wildcards) / "sample_{subset_name}_subset_mask",
    output:
        data_dir
        / "ancestors"
        / "{subset_name}-{region_name}-{filter_set}"
        / "ancestors-truncate-{lower}-{upper}-{multiplier}.trees",
        data_dir
        / "ancestors"
        / "{subset_name}-{region_name}-{filter_set}"
        / "ancestors-truncate-{lower}-{upper}-{multiplier}-performance_report.html",
    threads: config["max_threads"]
    resources:
        mem_mb=config["max_mem"],
        time_min=config["max_time"],
        runtime=config["max_time"],
    params:
        use_dask=config["match_ancestors"]["use_dask"],
    run:
        slug = f"{wildcards.subset_name}-{wildcards.region_name}-{wildcards.filter_set}-truncate-{wildcards.lower}-{wildcards.upper}-{wildcards.multiplier}"
        steps.match_ancestors(input, output, wildcards, config, threads, params, slug)


def get_sample_slices(subset_name):
    import numpy

    with open(config["sample_subsets"][subset_name], "r") as f:
        # FIXME! We need to know the ploidy here
        num_samples = len(numpy.genfromtxt(f, dtype=str)) * 2
        # Generate starts and ends for each chunk, of size config["match_samples"]["slice_size"]
        return [
            (i, min(i + config["match_samples"]["slice_size"] - 1, num_samples - 1))
            for i in range(0, num_samples, config["match_samples"]["slice_size"])
        ]


rule match_sample_paths:
    input:
        data_dir
        / "ancestors"
        / "{subset_name}-{region_name}-{filter_set}"
        / "ancestors-truncate-{lower}-{upper}-{multiplier}.trees",
        lambda wildcards: ds_dir(wildcards) / ".vcf_done",
        lambda wildcards: ds_dir(wildcards)
        / "variant_{subset_name}_subset_{region_name}_region_{filter_set}_mask",
        lambda wildcards: ds_dir(wildcards) / "sample_{subset_name}_subset_mask",
        lambda wildcards: config["recomb_map"].format(
            chrom=steps.parse_region(config["regions"][wildcards.region_name])[0]
        ),
    output:
        data_dir
        / "paths"
        / "{subset_name}-{region_name}-{filter_set}-truncate-{lower}-{upper}-{multiplier}-mm{mismatch}"
        / "sample-{sample_index_start}-{sample_index_end}.path",
    threads: 1
    resources:
        mem_mb=16000,
        time_min=60 * 4,
        runtime=60 * 4,
    run:
        steps.match_sample_path(input, output, wildcards, config, threads, params)


rule match_samples:
    input:
        data_dir
        / "ancestors"
        / "{subset_name}-{region_name}-{filter_set}"
        / "ancestors-truncate-{lower}-{upper}-{multiplier}.trees",
        lambda wildcards: ds_dir(wildcards) / ".vcf_done",
        lambda wildcards: ds_dir(wildcards)
        / "variant_{subset_name}_subset_{region_name}_region_{filter_set}_mask",
        lambda wildcards: ds_dir(wildcards) / "sample_{subset_name}_subset_mask",
        lambda wildcards: expand(
            data_dir
            / "paths"
            / "{subset_name}-{region_name}-{filter_set}-truncate-{lower}-{upper}-{multiplier}-mm{mismatch}"
            / "sample-{sample_slice[0]}-{sample_slice[1]}.path",
            sample_slice=get_sample_slices(wildcards.subset_name),
            allow_missing=True,
        ),
        lambda wildcards: config["recomb_map"].format(
            chrom=steps.parse_region(config["regions"][wildcards.region_name])[0]
        ),
    output:
        data_dir
        / "trees"
        / "{subset_name}-{region_name}-{filter_set}"
        / "{subset_name}-{region_name}-{filter_set}-truncate-{lower}-{upper}-{multiplier}-mm{mismatch}-raw.trees",
    # Minimal threads as we're using dask
    threads: 2
    resources:
        mem_mb=32000,
        time_min=config["max_time"],
        runtime=config["max_time"],
    run:
        slug = f"{wildcards.subset_name}-{wildcards.region_name}-{wildcards.filter_set}-truncate-{wildcards.lower}-{wildcards.upper}-{wildcards.multiplier}-mm-{wildcards.mismatch}"
        steps.match_samples(input, output, wildcards, config, threads, params, slug)


rule post_process:
    input:
        data_dir
        / "trees"
        / "{subset_name}-{region_name}-{filter_set}"
        / "{subset_name}-{region_name}-{filter_set}-truncate-{lower}-{upper}-{multiplier}-mm{mismatch}-raw.trees",
    output:
        data_dir
        / "trees"
        / "{subset_name}-{region_name}-{filter_set}"
        / "{subset_name}-{region_name}-{filter_set}-truncate-{lower}-{upper}-{multiplier}-mm{mismatch}-post-processed.trees",
    # Post process is currently not done in parallel
    threads: 2
    resources:
        mem_mb=64000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        import tsinfer
        import tskit

        ts = tskit.load(input[0])
        ts = tsinfer.post_process(ts)
        ts.dump(output[0])
