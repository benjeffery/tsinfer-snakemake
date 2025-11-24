import tskit

ts = tskit.load(snakemake.input[0])
ts = ts.simplify()
ts.dump(snakemake.output[0])
