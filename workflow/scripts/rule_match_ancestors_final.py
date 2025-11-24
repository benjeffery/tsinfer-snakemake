import logging
from pathlib import Path

import tsinfer

logging.basicConfig(level=logging.INFO)
ts = tsinfer.match_ancestors_batch_finalise(
    Path(snakemake.input[0]).parent,
)
ts.dump(snakemake.output[0])
