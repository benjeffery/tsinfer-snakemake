import logging
import os
from pathlib import Path

import tsinfer

logging.basicConfig(level=logging.INFO)
data_dir = Path(snakemake.config["data_dir"])
sample_data = tsinfer.VariantData(
    snakemake.input[0].replace(".vcf_done", ""),
    sample_mask=f"sample_{snakemake.wildcards.subset_name}_subset_mask",
    site_mask=f"variant_{snakemake.wildcards.subset_name}_subset_{snakemake.wildcards.region_name}_region_{snakemake.wildcards.filter_set}_mask",
    ancestral_state="variant_ancestral_allele",
)
logdir = data_dir / "progress" / "generate_ancestors"
os.makedirs(logdir, exist_ok=True)
logname = (
    f"{snakemake.wildcards.subset_name}-"
    + f"{snakemake.wildcards.region_name}-{snakemake.wildcards.filter_set}.log"
)
with open(logdir / logname, "w") as log_f:
    ancestors = tsinfer.generate_ancestors(
        sample_data,
        path=snakemake.output[0],
        genotype_encoding=1,
        num_threads=snakemake.threads,
        progress_monitor=tsinfer.progress.ProgressMonitor(
            tqdm_kwargs={"file": log_f, "mininterval": 30}
        ),
    )
if ancestors.num_ancestors == 0:
    raise ValueError("No ancestors generated")
if ancestors.num_sites == 0:
    raise ValueError("No sites generated")
