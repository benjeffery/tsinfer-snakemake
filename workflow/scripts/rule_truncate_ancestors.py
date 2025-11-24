import logging
import shutil

import tsinfer

lower = float(snakemake.wildcards.lower)
upper = float(snakemake.wildcards.upper)
multiplier = float(snakemake.wildcards.multiplier)
logging.basicConfig(level=logging.INFO)
if lower == 0 and upper == 0 and multiplier == 0:
    # No truncation, just copy the file using the OS
    shutil.copy(snakemake.input[0], snakemake.output[0])
else:
    ancestors = tsinfer.AncestorData.load(snakemake.input[0])
    truncated = ancestors.truncate_ancestors(
        lower, upper, length_multiplier=multiplier, path=snakemake.output[0]
    )
