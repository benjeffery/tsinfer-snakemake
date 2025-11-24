import logging
import re
from pathlib import Path

import tsinfer

logging.basicConfig(level=logging.INFO)
if ".pkl" in snakemake.input[0]:
    tsinfer.match_ancestors_batch_group_finalise(
        Path(snakemake.input[0]).parent.parent,
        group_index=int(snakemake.wildcards.group),
    )
else:
    output_group = int(snakemake.wildcards.group)
    if "metadata.json" in snakemake.input[0]:
        input_group = -1
    else:
        input_group = int(
            re.search(r"ancestors_(\d+).trees", snakemake.input[0]).group(1)
        )
    tsinfer.match_ancestors_batch_groups(
        Path(snakemake.input[0]).parent,
        input_group + 1,
        output_group + 1,
        snakemake.threads,
    )
