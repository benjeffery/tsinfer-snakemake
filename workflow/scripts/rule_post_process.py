import tsinfer
import tskit

ts = tskit.load(snakemake.input[0])
ts = tsinfer.post_process(ts)
ts.dump(snakemake.output[0])
