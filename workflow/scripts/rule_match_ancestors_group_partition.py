import logging
from pathlib import Path

import tsinfer

logging.basicConfig(level=logging.INFO)
tsinfer.match_ancestors_batch_group_partition(
    Path(snakemake.input[0]).parent,
    group_index=int(snakemake.wildcards.group),
    partition_index=int(snakemake.wildcards.partition),
)
