import logging
from pathlib import Path

import tsinfer

logging.basicConfig(level=logging.INFO)
work_dir = Path(snakemake.input[0]).parent
ts = tsinfer.match_samples_batch_finalise(work_dir)
ts.dump(snakemake.output[0])
