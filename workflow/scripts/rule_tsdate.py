import tsdate
import tskit

ts = tskit.load(snakemake.input[0])
ts = tsdate.date(
    ts,
    mutation_rate=float(snakemake.wildcards.mut_rate),
    progress=True,
    method="variational_gamma",
)
ts.dump(snakemake.output[0])
