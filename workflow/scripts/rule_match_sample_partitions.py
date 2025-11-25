import logging
from pathlib import Path

import tsinfer

logging.basicConfig(level=logging.INFO)
work_dir = Path(snakemake.output[0]).parent
tsinfer.match_samples_batch_partition(
    work_dir,
    partition_index=int(snakemake.wildcards.partition),
)
