import json
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import tsinfer

logging.basicConfig(level=logging.INFO)
tsinfer.match_ancestors_batch_init(
    work_dir=Path(snakemake.output[0]).parent,
    variant_data_path=snakemake.input[0].replace(".vcf_done", ""),
    ancestor_data_path=snakemake.input[1],
    ancestral_state="variant_ancestral_allele",
    min_work_per_job=snakemake.config["match_ancestors"]["min_work_per_job"],
    sample_mask=f"sample_{snakemake.wildcards.subset_name}_subset_mask",
    site_mask=f"variant_{snakemake.wildcards.subset_name}_subset_{snakemake.wildcards.region_name}_region_{snakemake.wildcards.filter_set}_mask",
    path_compression=True,
    precision=15,
)
with open(snakemake.output[0]) as f:
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
        f"Number of partitions {snakemake.wildcards.subset_name}-"
        + f"{snakemake.wildcards.region_name}-{snakemake.wildcards.filter_set}"
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
    plt.savefig(snakemake.output[0].replace("metadata.json", "partitions.png"))
