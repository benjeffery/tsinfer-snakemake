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


rule all:
    input:
        expand(
            data_dir / "{subset_name}-{filter}-region_summary_table.csv",
            subset_name=config["sample_subsets"].keys(),
            filter=config["filters"].keys(),
        ),
        expand(
            data_dir
            / "trees"
            / "{subset_name}-{region_name}-{filter}"
            / "{subset_name}-{region_name}-{filter}-truncate-{truncation}-mm{mismatch}-post-processed.trees",
            subset_name=config["sample_subsets"].keys(),
            region_name=config["regions"].keys(),
            filter=config["filters"].keys(),
            mismatch=config["mismatch_values"],
            truncation=[
                f"{c['lower']}-{c['upper']}-{c['multiplier']}"
                for c in config["truncate"]
            ],
        ),
        expand(
            data_dir
            / "ancestors"
            / "{subset_name}-{region_name}-{filter}"
            / "ancestors.png",
            subset_name=config["sample_subsets"].keys(),
            region_name=config["regions"].keys(),
            filter=config["filters"].keys(),
        ),


rule vcf_to_zarrs:
    input:
        vcf=lambda wildcards: config["vcf"].format(chrom=wildcards.chrom_num),
        tbi=lambda wildcards: config["vcf"].format(chrom=wildcards.chrom_num) + ".tbi",
    output:
        data_dir / "zarr_vcfs" / "chr{chrom_num}" / "data.zarr" / ".vcf_done",
        data_dir / "zarr_vcfs" / "chr{chrom_num}" / "performance_report.html",
    resources:
        dask_cluster=10,
        mem_mb=16000,
        time_min=24 * 60,
        runtime=24 * 60,
    params:
        target_part_size="5M",
        read_chunk_length=config["vcf_to_zarr"]["read_chunk_length"],
        temp_chunk_length=config["vcf_to_zarr"]["temp_chunk_length"],
        chunk_length=config["vcf_to_zarr"]["chunk_length"],
        chunk_width=config["vcf_to_zarr"]["chunk_width"],
        retain_temp_files=config["vcf_to_zarr"]["retain_temp_files"],
    run:
        steps.vcf_to_zarrs(input, output, wildcards, config, params)


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
    output:
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
            / "variant_bad_ancestral_mask"
        ),
        directory(
            data_dir
            / "zarr_vcfs"
            / "chr{chrom_num}"
            / "data.zarr"
            / "variant_no_ancestral_allele_mask"
        ),
    resources:
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        steps.pre_subset_filters(input, output, wildcards, config, params)


rule subset_zarr_vcf:
    input:
        lambda wildcards: [
            (
                data_dir
                / "zarr_vcfs"
                / f"chr{steps.parse_region(config['regions'][wildcards.region_name])[0]}"
                / "data.zarr"
                / suffix
            )
            for suffix in [
                ".vcf_done",
                "variant_ancestral_allele",
                "variant_duplicate_position_mask",
                "variant_bad_ancestral_mask",
                "variant_no_ancestral_allele_mask",
            ]
        ],
        lambda wildcards: config["sample_subsets"][wildcards.subset_name],
    output:
        data_dir
        / "zarr_vcfs_subsets"
        / "{subset_name}-{region_name}-{filter}"
        / "data.zarr"
        / ".subset_done",
    resources:
        dask_cluster=5,
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        steps.subset_zarr_vcf(input, output, wildcards, config, params)
        # We have to do this as snakemake doesn't like subsequent outputs that are children of this one
        Path(output[0]).touch()


rule post_subset_filters:
    input:
        data_dir
        / "zarr_vcfs_subsets"
        / "{subset_name}-{region_name}-{filter}"
        / "data.zarr"
        / ".subset_done",
    output:
        directory(
            data_dir
            / "zarr_vcfs_subsets"
            / "{subset_name}-{region_name}-{filter}"
            / "data.zarr"
            / "variant_mask"
        ),
    resources:
        dask_cluster=5,
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        steps.post_subset_filters(input, output, wildcards, config, params)


rule zarr_stats:
    input:
        data_dir
        / "zarr_vcfs_subsets"
        / "{subset_name}-{region_name}-{filter}"
        / "data.zarr"
        / "variant_mask",
    output:
        data_dir / "zarr_stats" / "{subset_name}-{region_name}-{filter}" / "stats.json",
    resources:
        dask_cluster=5,
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        steps.zarr_stats(input, output, wildcards, config, params)


checkpoint summary_table:
    input:
        lambda wildcards: [
            data_dir
            / "zarr_stats"
            / f"{wildcards.subset_name}-{region_name}-{wildcards.filter}"
            / "stats.json"
            for region_name in config["regions"].keys()
        ],
    output:
        data_dir / "{subset_name}-{filter}-region_summary_table.csv",
    run:
        steps.summary_table(input, output, wildcards, config, params)


def get_ancestor_gen_memory(wildcards):
    import numpy
    import json

    # Use the checkpoint to check the stats file exists
    checkpoint_output = checkpoints.summary_table.get(
        subset_name=wildcards.subset_name, filter=wildcards.filter
    )
    region_stats = (
        data_dir
        / "zarr_stats"
        / f"{wildcards.subset_name}-{wildcards.region_name}-{wildcards.filter}"
        / "stats.json"
    )
    with open(region_stats, "r") as json_stats_f:
        with open(config["sample_subsets"][wildcards.subset_name], "r") as f:
            n_samples = len(numpy.genfromtxt(f, dtype=str))
        stats = json.load(json_stats_f)
        n_ploidy = stats["n_ploidy"]
        n_sites = stats["n_variants"]
        n_masked = stats["sites_masked"]
        ac1 = sum([ac == 1 for ac in stats["allele_counts"]])
        mem = 16_000 + int(
            (((n_sites - n_masked) - ac1) * n_samples * n_ploidy) / (8 * 1_048_576)
        )
        return mem


rule generate_ancestors:
    input:
        data_dir
        / "zarr_vcfs_subsets"
        / "{subset_name}-{region_name}-{filter}"
        / "data.zarr"
        / "variant_mask",
        data_dir / "zarr_stats" / "{subset_name}-{region_name}-{filter}" / "stats.json",
    output:
        data_dir
        / "ancestors"
        / "{subset_name}-{region_name}-{filter}"
        / "ancestors.zarr",
    threads: config["max_threads"]
    resources:
        mem_mb=get_ancestor_gen_memory,
        time_min=config["max_time"],
        runtime=config["max_time"],
    run:
        steps.generate_ancestors(input, output, wildcards, config, threads)

rule ancestor_stats:
    input:
       data_dir
        / "ancestors"
        / "{subset_name}-{region_name}-{filter}"
        / "ancestors.zarr",
    output:
        data_dir
        / "ancestors"
        / "{subset_name}-{region_name}-{filter}"
        / "ancestors.png",
    threads: 1
    resources:
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        steps.ancestor_stats(input, output, wildcards, config)

rule truncate_ancestors:
    input:
        data_dir
        / "ancestors"
        / "{subset_name}-{region_name}-{filter}"
        / "ancestors.zarr",
    output:
        data_dir
        / "ancestors"
        / "{subset_name}-{region_name}-{filter}"
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
        / "{subset_name}-{region_name}-{filter}"
        / "ancestors-truncate-{lower}-{upper}-{multiplier}.zarr",
        data_dir
        / "zarr_vcfs_subsets"
        / "{subset_name}-{region_name}-{filter}"
        / "data.zarr"
        / "variant_mask",
    output:
        data_dir
        / "ancestors"
        / "{subset_name}-{region_name}-{filter}"
        / "ancestors-truncate-{lower}-{upper}-{multiplier}.trees",
        data_dir
        / "ancestors"
        / "{subset_name}-{region_name}-{filter}"
        / "ancestors-truncate-{lower}-{upper}-{multiplier}-performance_report.html",
    threads: config["max_threads"]
    resources:
        mem_mb=config["max_mem"],
        time_min=config["max_time"],
        runtime=config["max_time"],
    params:
        use_dask=config["match_ancestors"]["use_dask"],
    run:
        slug = f"{wildcards.subset_name}-{wildcards.region_name}-{wildcards.filter}-truncate-{wildcards.lower}-{wildcards.upper}-{wildcards.multiplier}"
        steps.match_ancestors(input, output, wildcards, config, threads, params, slug)


rule match_samples:
    input:
        data_dir
        / "ancestors"
        / "{subset_name}-{region_name}-{filter}"
        / "ancestors-truncate-{lower}-{upper}-{multiplier}.trees",
        data_dir
        / "zarr_vcfs_subsets"
        / "{subset_name}-{region_name}-{filter}"
        / "data.zarr"
        / "variant_mask",
        lambda wildcards: config["recomb_map"].format(
            chrom=steps.parse_region(config["regions"][wildcards.region_name])[0]
        ),
    output:
        data_dir
        / "trees"
        / "{subset_name}-{region_name}-{filter}"
        / "{subset_name}-{region_name}-{filter}-truncate-{lower}-{upper}-{multiplier}-mm{mismatch}-raw.trees",
        data_dir
        / "samples"
        / "{subset_name}-{region_name}-{filter}-truncate-{lower}-{upper}-{multiplier}-mm{mismatch}"
        / "performance_report.html",
    # Minimal threads as we're using dask
    threads: 2
    resources:
        mem_mb=32000,
        time_min=config["max_time"],
        runtime=config["max_time"],
    params:
        use_dask=config["match_samples"]["use_dask"],
    run:
        slug = f"{wildcards.subset_name}-{wildcards.region_name}-{wildcards.filter}-truncate-{wildcards.lower}-{wildcards.upper}-{wildcards.multiplier}-mm-{wildcards.mismatch}"
        steps.match_samples(input, output, wildcards, config, threads, params, slug)


rule post_process:
    input:
        data_dir
        / "trees"
        / "{subset_name}-{region_name}-{filter}"
        / "{subset_name}-{region_name}-{filter}-truncate-{lower}-{upper}-{multiplier}-mm{mismatch}-raw.trees",
    output:
        data_dir
        / "trees"
        / "{subset_name}-{region_name}-{filter}"
        / "{subset_name}-{region_name}-{filter}-truncate-{lower}-{upper}-{multiplier}-mm{mismatch}-post-processed.trees",
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
