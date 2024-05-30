import yaml
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

with open("resources.yaml", "r") as resource_file:
    resources_config = yaml.safe_load(resource_file)


def get_resource(rule_name, resource_type):
    if rule_name in resources_config and resource_type in resources_config[rule_name]:
        return resources_config[rule_name][resource_type]
    else:
        return resources_config["default"][resource_type]


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
            / "{subset_name}-{region_name}-{filter_set}-truncate-{truncation}-mm{mismatch}-post-processed-simplified-SDN-dated-{mut_rate}.trees",
            subset_name=config["sample_subsets"].keys(),
            region_name=config["regions"].keys(),
            filter_set=config["filters"].keys(),
            mismatch=config["mismatch_values"],
            truncation=[
                f"{c['lower']}-{c['upper']}-{c['multiplier']}"
                for c in config["truncate"]
            ],
            mut_rate=[f"{c['mutation_rate']}" for c in config["tsdate"]],
        ),


checkpoint bio2zarr_dexplode_init:
    input:
        vcf=lambda wildcards: config["vcf"].format(chrom=wildcards.chrom_num),
        tbi=lambda wildcards: config["vcf"].format(chrom=wildcards.chrom_num) + ".tbi",
    output:
        data_dir / "exploded_vcfs_md" / "chr{chrom_num}.metadata.json",
    threads: get_resource("bio2zarr_dexplode_init", "threads")
    resources:
        mem_mb=get_resource("bio2zarr_dexplode_init", "mem_mb"),
        time_min=get_resource("bio2zarr_dexplode_init", "time_min"),
        runtime=get_resource("bio2zarr_dexplode_init", "time_min"),
    run:
        (data_dir / "exploded_vcfs" / f"chr{wildcards.chrom_num}").mkdir(
            parents=True, exist_ok=True
        )
        shell(
            f"python -m bio2zarr vcf2zarr dexplode-init --json --force --num-partitions {config['bio2zarr']['num_partitions']} {input.vcf} {data_dir}/exploded_vcfs/chr{wildcards.chrom_num} > {output[0]}"
        )


rule bio2zarr_dexplode_partition:
    input:
        data_dir / "exploded_vcfs_md" / "chr{chrom_num}.metadata.json",
    output:
        data_dir / "exploded_vcfs" / "chr{chrom_num}" / ".done_{partition}",
    threads: get_resource("bio2zarr_dexplode_partition", "threads")
    resources:
        mem_mb=get_resource("bio2zarr_dexplode_partition", "mem_mb"),
        time_min=get_resource("bio2zarr_dexplode_partition", "time_min"),
        runtime=get_resource("bio2zarr_dexplode_partition", "time_min"),
    run:
        shell(
            f"python -m bio2zarr vcf2zarr dexplode-partition {Path(output[0]).parent} {wildcards.partition} && touch {output[0]}"
        )


def dexplode_partitions(wildcards):
    import json

    checkpoint_output = checkpoints.bio2zarr_dexplode_init.get(
        chrom_num=wildcards.chrom_num
    )
    with open(checkpoint_output.output[0], "r") as f:
        md = json.load(f)
    return [
        data_dir / "exploded_vcfs" / f"chr{wildcards.chrom_num}" / f".done_{p}"
        for p in range(md["num_partitions"])
    ]


rule bio2zarr_dexplode_finalise:
    input:
        dexplode_partitions,
    output:
        data_dir / "exploded_vcfs" / "chr{chrom_num}" / ".done",
    threads: get_resource("bio2zarr_dexplode_finalise", "threads")
    resources:
        mem_mb=get_resource("bio2zarr_dexplode_finalise", "mem_mb"),
        time_min=get_resource("bio2zarr_dexplode_finalise", "time_min"),
        runtime=get_resource("bio2zarr_dexplode_finalise", "time_min"),
    run:
        shell(
            f"python -m bio2zarr vcf2zarr dexplode-finalise {Path(output[0]).parent} && touch {output[0]}"
        )


rule bio2zarr_mkschema:
    input:
        data_dir / "exploded_vcfs" / "chr{chrom_num}" / ".done",
    output:
        data_dir / "zarr_vcfs_schema" / "chr{chrom_num}" / "schema.json",
    threads: get_resource("bio2zarr_mkschema", "threads")
    resources:
        mem_mb=get_resource("bio2zarr_mkschema", "mem_mb"),
        time_min=get_resource("bio2zarr_mkschema", "time_min"),
        runtime=get_resource("bio2zarr_mkschema", "time_min"),
    run:
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        shell(
            f"python -m bio2zarr vcf2zarr mkschema {Path(input[0]).parent} > {output[0]}"
        )


checkpoint bio2zarr_dencode_init:
    input:
        data_dir / "exploded_vcfs" / "chr{chrom_num}" / ".done",
        data_dir / "zarr_vcfs_schema" / "chr{chrom_num}" / "schema.json",
    output:
        data_dir / "zarr_vcfs_md" / "chr{chrom_num}.metadata.json",
    threads: get_resource("bio2zarr_dencode_init", "threads")
    resources:
        mem_mb=get_resource("bio2zarr_dencode_init", "mem_mb"),
        time_min=get_resource("bio2zarr_dencode_init", "time_min"),
        runtime=get_resource("bio2zarr_dencode_init", "time_min"),
    run:
        (data_dir / "zarr_vcfs" / f"chr{wildcards.chrom_num}" / "data.zarr").mkdir(
            parents=True, exist_ok=True
        )
        shell(
            f"python -m bio2zarr vcf2zarr dencode-init --json --force --num-partitions {config['bio2zarr']['num_partitions']} --variants-chunk-size {config['bio2zarr']['variants_chunk_size']} {Path(input[0]).parent} {data_dir}/zarr_vcfs/chr{wildcards.chrom_num}/data.zarr > {output[0]}"
        )


def dencode_partitions(wildcards):
    import json

    checkpoint_output = checkpoints.bio2zarr_dencode_init.get(
        chrom_num=wildcards.chrom_num
    )
    with open(checkpoint_output.output[0], "r") as f:
        md = json.load(f)
    return [
        ds_dir(wildcards) /
        / "data.zarr"
        / f".done_{p}"
        for p in range(md["num_partitions"])
    ]


rule bio2zarr_dencode_partition:
    input:
        data_dir / "zarr_vcfs_md" / "chr{chrom_num}.metadata.json",
    output:
        data_dir / "zarr_vcfs" / "chr{chrom_num}" / "data.zarr" / ".done_{partition}",
    threads: get_resource("bio2zarr_dencode_partition", "threads")
    resources:
        mem_mb=get_resource("bio2zarr_dencode_partition", "mem_mb"),
        time_min=get_resource("bio2zarr_dencode_partition", "time_min"),
        runtime=get_resource("bio2zarr_dencode_partition", "time_min"),
    run:
        shell(
            f"python -m bio2zarr vcf2zarr dencode-partition {Path(output[0]).parent} {wildcards.partition} && touch {output[0]}"
        )


rule bio2zarr_dencode_finalise:
    input:
        dencode_partitions,
    output:
        data_dir / "zarr_vcfs" / "chr{chrom_num}" / "data.zarr" / ".vcf_done",
    threads: get_resource("bio2zarr_dencode_finalise", "threads")
    resources:
        mem_mb=get_resource("bio2zarr_dencode_finalise", "mem_mb"),
        time_min=get_resource("bio2zarr_dencode_finalise", "time_min"),
        runtime=get_resource("bio2zarr_dencode_finalise", "time_min"),
    run:
        shell(
            f"python -m bio2zarr vcf2zarr dencode-finalise {Path(output[0]).parent} && touch {output[0]}"
        )
        # Remove the metadata files so thaat we can write extra arrays later without updating it.
        shell(f"rm -rf {Path(output[0]).parent}/.zmetadata*")
        shell(f"touch {Path(output[0])}")


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
    threads: get_resource("load_ancestral_fasta", "threads")
    resources:
        mem_mb=get_resource("load_ancestral_fasta", "mem_mb"),
        time_min=get_resource("load_ancestral_fasta", "time_min"),
        runtime=get_resource("load_ancestral_fasta", "time_min"),
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
    threads: get_resource("pre_subset_filters", "threads")
    resources:
        mem_mb=get_resource("pre_subset_filters", "mem_mb"),
        time_min=get_resource("pre_subset_filters", "time_min"),
        runtime=get_resource("pre_subset_filters", "time_min"),
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
    threads: get_resource("region_mask", "threads")
    resources:
        mem_mb=get_resource("region_mask", "mem_mb"),
        time_min=get_resource("region_mask", "time_min"),
        runtime=get_resource("region_mask", "time_min"),
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
    threads: get_resource("sample_mask", "threads")
    resources:
        mem_mb=get_resource("sample_mask", "mem_mb"),
        time_min=get_resource("sample_mask", "time_min"),
        runtime=get_resource("sample_mask", "time_min"),
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
    threads: get_resource("allele_counts", "threads")
    resources:
        mem_mb=get_resource("allele_counts", "mem_mb"),
        time_min=get_resource("allele_counts", "time_min"),
        runtime=get_resource("allele_counts", "time_min"),
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
    threads: get_resource("subset_filters", "threads")
    resources:
        mem_mb=get_resource("subset_filters", "mem_mb"),
        time_min=get_resource("subset_filters", "time_min"),
        runtime=get_resource("subset_filters", "time_min"),
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
    threads: get_resource("site_density_mask", "threads")
    resources:
        mem_mb=get_resource("site_density_mask", "mem_mb"),
        time_min=get_resource("site_density_mask", "time_min"),
        runtime=get_resource("site_density_mask", "time_min"),
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
    threads: get_resource("combined_mask", "threads")
    resources:
        mem_mb=get_resource("combined_mask", "mem_mb"),
        time_min=get_resource("combined_mask", "time_min"),
        runtime=get_resource("combined_mask", "time_min"),
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
    threads: get_resource("zarr_stats", "threads")
    resources:
        mem_mb=get_resource("zarr_stats", "mem_mb"),
        time_min=get_resource("zarr_stats", "time_min"),
        runtime=get_resource("zarr_stats", "time_min"),
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
    threads: get_resource("generate_ancestors", "threads")
    resources:
        mem_mb=get_resource("generate_ancestors", "mem_mb"),
        time_min=get_resource("generate_ancestors", "time_min"),
        runtime=get_resource("generate_ancestors", "time_min"),
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
    threads: get_resource("truncate_ancestors", "threads")
    resources:
        mem_mb=get_resource("truncate_ancestors", "mem_mb"),
        time_min=get_resource("truncate_ancestors", "time_min"),
        runtime=get_resource("truncate_ancestors", "time_min"),
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
    threads: get_resource("match_ancestors_init","threads")
    resources:
        mem_mb=get_resource("match_ancestors_init","mem_mb"),
        time_min=get_resource("match_ancestors_init","time_min"),
        runtime=get_resource("match_ancestors_init", "time_min"),
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


def ancestor_groupings(checkpoint_output):
    with open(checkpoint_output.output[0], "r") as f:
        md = json.load(f)
    return md["ancestor_grouping"]


def num_partitions(checkpoint_output):
    return ancestor_groupings(checkpoint_output)[int(wildcards.group)]["partitions"]


# This function decides if a group should be processed in a single job or partitioned
def match_ancestor_group_input(wildcards):
    checkpoint_output = checkpoints.match_ancestors_init.get(**wildcards)
    match_dir = Path(checkpoint_output.output[0]).parent
    groupings = ancestor_groupings(checkpoint_output)
    group_index = int(wildcards.group)
    partitions = groupings[group_index]["partitions"]
    if partitions is not None:
        return expand(
            f"{match_dir}/group_{group_index}/partition_" + "{partition}.pkl",
            partition=range(len(partitions)),
            allow_missing=True,
        )

    # This group is small enough to do locally
    # search back until we find a group that requires partitioning, or we reach the start, or we have enough groups
    max_groups = config["match_ancestors"]["max_groups"]
    for i in range(group_index, max(group_index - max_groups, 0), -1):
        if groupings[i]["partitions"] is not None:
            return match_dir / f"ancestors_{i}.trees"
    if group_index - max_groups > 0:
        return match_dir / f"ancestors_{group_index-max_groups}.trees"
    else:
        return checkpoint_output.output[0]


def match_ancestors_group_num_threads(wildcards):
    input_ = match_ancestor_group_input(wildcards)
    if isinstance(input_, list):
        input_ = input_[0]
    if ".pkl" in str(input_):
        return 1
    else:
        return config["max_threads"]


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
            for group in range(input_group + 1, output_group + 1):
                tsinfer.match_ancestors_batch_group(
                    Path(input[0]).parent,
                    group,
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
    threads: get_resource("match_ancestors_group_partition", "threads")
    resources:
        get_resource("match_ancestors_group_partition", "mem_mb"),
        get_resource("match_ancestors_group_partition", "time_min"),
        get_resource("match_ancestors_group_partition", "time_min"),
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
    threads: get_resource("match_sample_paths", "threads")
    resources:
        mem_mb=get_resource("match_sample_paths", "mem_mb"),
        time_min=get_resource("match_sample_paths", "time_min"),
        runtime=get_resource("match_sample_paths", "time_min"),
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
    threads: get_resource("match_samples", "threads")
    resources:
        mem_mb=get_resource("match_samples", "mem_mb"),
        time_min=get_resource("match_samples", "time_min"),
        runtime=get_resource("match_samples", "time_min"),
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
    threads: get_resource("post_process", "threads")
    resources:
        mem_mb=get_resource("post_process", "mem_mb"),
        time_min=get_resource("post_process", "time_min"),
        runtime=get_resource("post_process", "time_min"),
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
    threads: get_resource("full_simplify", "threads")
    resources:
        mem_mb=get_resource("full_simplify", "mem_mb"),
        time_min=get_resource("full_simplify", "time_min"),
        runtime=get_resource("full_simplify", "time_min"),
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
    threads: get_resource("split_disjoint_nodes", "threads")
    resources:
        mem_mb=get_resource("split_disjoint_nodes", "mem_mb"),
        time_min=get_resource("split_disjoint_nodes", "time_min"),
        runtime=get_resource("split_disjoint_nodes", "time_min"),
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
    threads: get_resource("tsdate", "threads")
    resources:
        mem_mb=get_resource("tsdate", "mem_mb"),
        time_min=get_resource("tsdate", "time_min"),
        runtime=get_resource("tsdate", "time_min"),
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
