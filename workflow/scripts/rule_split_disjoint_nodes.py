import tsdate
import tskit

ts = tskit.load(snakemake.input[0])
ts = tsdate.util.split_disjoint_nodes(ts)
ts.dump(snakemake.output[0])
