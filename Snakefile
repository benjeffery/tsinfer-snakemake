import os.path
from pathlib import Path
import steps
import sys
import json


# The match ancestors step requires as much rescursion as
# there are ancestor groups
sys.setrecursionlimit(10000)


configfile: "config.yaml"


shell.prefix(config["prefix"])

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
            / "{subset_name}-{region_name}-{filter_set}-truncate-{truncation}-mm{mismatch}-post-processed-simplified-SDN-dated-{norm_int}.trees",
            subset_name=config["sample_subsets"].keys(),
            region_name=config["regions"].keys(),
            filter_set=config["filters"].keys(),
            mismatch=config["mismatch_values"],
            truncation=[
                f"{c['lower']}-{c['upper']}-{c['multiplier']}"
                for c in config["truncate"]
            ],
            norm_int=[f"{c['mutation_rate']}" for c in config["tsdate"]],
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

        # Snakemake creates the output dir - bio2zarr then complains that the output exists
        # so we need to remove it first
        try:
            os.rmdir(data_dir / "zarr_vcfs" / f"chr{wildcards.chrom_num}" / "data.zarr")
        except FileNotFoundError:
            pass
        vcf.encode(
            input[0],
            output[0].replace(".vcf_done", ""),
            input[1],
            worker_processes=threads,
            max_memory=config["max_mem"],
        )
        # Remove the consolidated metadata file
        os.remove(output[0].replace(".vcf_done", ".zmetadata"))
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
        lambda wildcards: ds_dir(wildcards) / ".vcf_done",
        lambda wildcards: ds_dir(wildcards) / "variant_ancestral_allele",
        lambda wildcards: ds_dir(wildcards)
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
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        steps.pre_subset_filters(input, output, wildcards, config, params)


rule region_mask:
    input:
        lambda wildcards: ds_dir(wildcards) / ".vcf_done",
    output:
        directory(
            data_dir
            / "zarr_vcfs"
            / "chr{chrom_num}"
            / "data.zarr"
            / "variant_{region_name}_region_mask"
        ),
    resources:
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        steps.region_mask(input, output, wildcards, config, params)


rule sample_mask:
    input:
        lambda wildcards: ds_dir(wildcards) / ".vcf_done",
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
        lambda wildcards: ds_dir(wildcards) / ".vcf_done",
        lambda wildcards: ds_dir(wildcards) / "sample_{subset_name}_subset_mask",
        lambda wildcards: ds_dir(wildcards) / "variant_ancestral_allele",
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
        threads=config["max_threads"],
        mem_mb=config["max_mem"],
        time_min=4 * 60,
        runtime=4 * 60,
    run:
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
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
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
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
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
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
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


checkpoint match_ancestors_init:
    input:
        lambda wildcards: ds_dir(wildcards) / ".vcf_done",
        data_dir
        / "ancestors"
        / "{subset_name}-{region_name}-{filter_set}"
        / "ancestors-truncate-{lower}-{upper}-{multiplier}.zarr",
    output:
        data_dir
        / "ancestors_working"
        / "{subset_name}-{region_name}-{filter_set}"
        / "ancestors-truncate-{lower}-{upper}-{multiplier}"
        / "metadata.json",
    threads: 1
    resources:
        mem_mb=32000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        import tsinfer
        import logging
        import matplotlib.pyplot as plt

        logging.basicConfig(level=logging.INFO)
        tsinfer.match_ancestors_batch_init(
            working_dir=Path(output[0]).parent,
            sample_data_path=input[0].replace(".vcf_done", ""),
            ancestor_data_path=input[1],
            min_work_per_job=config["match_ancestors"]["min_work_per_job"],
            sgkit_samples_mask_name=f"sample_{wildcards.subset_name}_subset_mask",
            sites_mask_name=f"variant_{wildcards.subset_name}_subset_{wildcards.region_name}_region_{wildcards.filter_set}_mask",
            path_compression=True,
            precision=15,
        )
        with open(output[0], "r") as f:
            md = json.load(f)
            plt.plot(
                [
                    len(g["partitions"]) if g["partitions"] is not None else 1
                    for g in md["ancestor_grouping"]
                ]
            )
            plt.xlabel("Group")
            plt.ylabel("Number of partitions")
            plt.title(
                f"Number of partitions {wildcards.subset_name}-{wildcards.region_name}-{wildcards.filter_set}"
            )
            plt.yscale("log")
            # on a second axis plot num ancestors
            plt.twinx()
            plt.plot(
                [len(g["ancestors"]) for g in md["ancestor_grouping"]],
                color="red",
            )
            plt.yscale("log")
            plt.ylabel("Number of ancestors")
            plt.savefig(output[0].replace("metadata.json", "partitions.png"))


# Cluster filesystems are slow so cache the ancestor groupings
ancestor_groupings_cache = {}


def ancestor_groupings(checkpoint_output):
    output_path = checkpoint_output.output[0]
    if output_path not in ancestor_groupings_cache:
        with open(output_path, "r") as f:
            md = json.load(f)
            ancestor_groupings_cache[output_path] = md["ancestor_grouping"]
    return ancestor_groupings_cache[output_path]


def num_partitions(checkpoint_output):
    return ancestor_groupings(checkpoint_output)[int(wildcards.group)]["partitions"]


# This function decides if a group should be processed in a single job or partitioned
def match_ancestor_group_input(wildcards):
    checkpoint_output = checkpoints.match_ancestors_init.get(**wildcards)
    match_dir = Path(checkpoint_output.output[0]).parent
    groupings = ancestor_groupings(checkpoint_output)
    group_index = int(wildcards.group)
    partitions = groupings[group_index]["partitions"]
    # Don't use partitions if we are in the first half of the groups, or there are none
    if (
        partitions is not None
        and len(partitions) > config["max_threads"]
        and group_index > len(groupings) // 2
    ):
        return expand(
            f"{match_dir}/group_{group_index}/partition_" + "{partition}.pkl",
            partition=range(len(partitions)),
            allow_missing=True,
        )

    # This group is small enough to do locally
    # search back until we find a group that requires partitioning, or we reach the start, or we have enough groups
    max_groups = config["match_ancestors"]["max_groups"]
    for i in range(group_index, max(group_index - max_groups, 0), -1):
        if (
            groupings[i]["partitions"] is not None
            and len(groupings[i]["partitions"]) > config["max_threads"]
            and i > len(groupings) // 2
        ):
            return match_dir / f"ancestors_{i}.trees"
    if group_index - max_groups > 0:
        return match_dir / f"ancestors_{group_index-max_groups}.trees"
    else:
        return checkpoint_output.output[0]


def match_ancestors_group_num_threads(wildcards):
    input_ = match_ancestor_group_input(wildcards)
    if isinstance(input_, list):
        input_ = input_[0]
    input_ = str(input_)
    # If the input is a pickle file, we are finalising a split group
    if ".pkl" in input_:
        return 1
    # The first group is trivial
    if "metadata.json" in input_:
        first_group = 0
    else:
        first_group = int(re.search(r"ancestors_(\d+).trees", input_).group(1))
    last_group = int(wildcards.group)
    groupings = ancestor_groupings(checkpoints.match_ancestors_init.get(**wildcards))
    most_seen_ancestors = 1
    for i in range(first_group, last_group + 1):
        most_seen_ancestors = max(most_seen_ancestors, len(groupings[i]["ancestors"]))
    return min(most_seen_ancestors, config["max_threads"])


def match_ancestors_group_ram(wildcards):
    input_ = match_ancestor_group_input(wildcards)
    if isinstance(input_, list):
        input_ = input_[0]
    if ".pkl" in str(input_):
        return 16000
    else:
        return config["max_mem"]


rule match_ancestors_group:
    input:
        match_ancestor_group_input,
    output:
        data_dir
        / "ancestors_working"
        / "{subset_name}-{region_name}-{filter_set}"
        / "ancestors-truncate-{lower}-{upper}-{multiplier}"
        / "ancestors_{group}.trees",
    threads: match_ancestors_group_num_threads
    resources:
        mem_mb=match_ancestors_group_ram,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        import tsinfer
        import logging

        logging.basicConfig(level=logging.INFO)
        if ".pkl" in input[0]:
            tsinfer.match_ancestors_batch_group_finalise(
                Path(input[0]).parent.parent,
                group_index=int(wildcards.group),
            )
        else:
            output_group = int(wildcards.group)
            if "metadata.json" in input[0]:
                input_group = -1
            else:
                input_group = int(
                    re.search(r"ancestors_(\d+).trees", input[0]).group(1)
                )
            tsinfer.match_ancestors_batch_groups(
                Path(input[0]).parent,
                input_group + 1,
                output_group + 1,
                threads,
            )


rule match_ancestors_group_partition:
    input:
        lambda wildcards: (
            data_dir
            / "ancestors_working"
            / "{subset_name}-{region_name}-{filter_set}"
            / "ancestors-truncate-{lower}-{upper}-{multiplier}"
            / f"ancestors_{int(wildcards.group)-1}.trees"
        ),
    output:
        data_dir
        / "ancestors_working"
        / "{subset_name}-{region_name}-{filter_set}"
        / "ancestors-truncate-{lower}-{upper}-{multiplier}"
        / "group_{group}"
        / "partition_{partition}.pkl",
    threads: 1
    resources:
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        import tsinfer
        import logging

        logging.basicConfig(level=logging.INFO)
        tsinfer.match_ancestors_batch_group_partition(
            Path(input[0]).parent,
            group_index=int(wildcards.group),
            partition_index=int(wildcards.partition),
        )


def last_ancestor_group(wildcards):
    checkpoint_output = checkpoints.match_ancestors_init.get(**wildcards)
    groupings = ancestor_groupings(checkpoint_output)
    return len(groupings) - 1


rule match_ancestors_final:
    input:
        lambda wildcards: (
            data_dir
            / "ancestors_working"
            / "{subset_name}-{region_name}-{filter_set}"
            / "ancestors-truncate-{lower}-{upper}-{multiplier}"
            / f"ancestors_{last_ancestor_group(wildcards)}.trees"
        ),
    output:
        data_dir
        / "ancestors"
        / "{subset_name}-{region_name}-{filter_set}"
        / "ancestors-truncate-{lower}-{upper}-{multiplier}.trees",
    threads: 1
    resources:
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        import tsinfer
        import logging

        logging.basicConfig(level=logging.INFO)
        ts = tsinfer.match_ancestors_batch_finalise(
            Path(input[0]).parent,
        )
        ts.dump(output[0])


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
    threads: 1
    resources:
        mem_mb=64000,
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
    threads: 1
    resources:
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        import tsinfer
        import tskit

        ts = tskit.load(input[0])
        ts = tsinfer.post_process(ts)
        ts.dump(output[0])


rule full_simplify:
    input:
        data_dir
        / "trees"
        / "{subset_name}-{region_name}-{filter_set}"
        / "{subset_name}-{region_name}-{filter_set}-truncate-{lower}-{upper}-{multiplier}-mm{mismatch}-post-processed.trees",
    output:
        data_dir
        / "trees"
        / "{subset_name}-{region_name}-{filter_set}"
        / "{subset_name}-{region_name}-{filter_set}-truncate-{lower}-{upper}-{multiplier}-mm{mismatch}-post-processed-simplified.trees",
    threads: 1
    resources:
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        import tskit

        ts = tskit.load(input[0])
        ts = ts.simplify()
        ts.dump(output[0])


rule split_disjoint_nodes:
    input:
        data_dir
        / "trees"
        / "{subset_name}-{region_name}-{filter_set}"
        / "{subset_name}-{region_name}-{filter_set}-truncate-{lower}-{upper}-{multiplier}-mm{mismatch}-post-processed-simplified.trees",
    output:
        data_dir
        / "trees"
        / "{subset_name}-{region_name}-{filter_set}"
        / "{subset_name}-{region_name}-{filter_set}-truncate-{lower}-{upper}-{multiplier}-mm{mismatch}-post-processed-simplified-SDN.trees",
    threads: 1
    resources:
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        import tskit
        import tsdate

        ts = tskit.load(input[0])
        ts = tsdate.util.split_disjoint_nodes(ts)
        ts.dump(output[0])


rule tsdate:
    input:
        data_dir
        / "trees"
        / "{subset_name}-{region_name}-{filter_set}"
        / "{subset_name}-{region_name}-{filter_set}-truncate-{lower}-{upper}-{multiplier}-mm{mismatch}-post-processed-simplified-SDN.trees",
    output:
        data_dir
        / "trees"
        / "{subset_name}-{region_name}-{filter_set}"
        / "{subset_name}-{region_name}-{filter_set}-truncate-{lower}-{upper}-{multiplier}-mm{mismatch}-post-processed-simplified-SDN-dated-{mut_rate}.trees",
    threads: 1
    resources:
        mem_mb=16000,
        time_min=4 * 60,
        runtime=4 * 60,
    run:
        import tskit
        import tsdate

        ts = tskit.load(input[0])
        ts = tsdate.date(
            ts,
            mutation_rate=float(wildcards.mut_rate),
            progress=True,
            method="variational_gamma",
        )
        ts.dump(output[0])
