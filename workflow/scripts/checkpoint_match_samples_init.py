import logging
from pathlib import Path

import msprime
import numpy
import tsinfer
import tskit


def build_maps(anc_ts, mismatch, recomb_map):
    mismatch = float(mismatch)
    if mismatch > 0:
        inference_pos = anc_ts.tables.sites.position
        rate_map = msprime.RateMap.read_hapmap(recomb_map, position_col=1, rate_col=2)
        genetic_dists = tsinfer.Matcher.recombination_rate_to_dist(
            rate_map, inference_pos
        )
        recombination_map = tsinfer.Matcher.recombination_dist_to_prob(genetic_dists)
        # Set 0 probabilities to a small value
        recombination_map[recombination_map == 0] = 1e-19
        mismatch_ratio = mismatch
        num_alleles = 2
        mismatch_map = numpy.full(
            len(inference_pos),
            tsinfer.Matcher.mismatch_ratio_to_prob(
                mismatch_ratio, numpy.median(genetic_dists), num_alleles
            ),
        )
    else:
        recombination_map = None
        mismatch_map = None
    return recombination_map, mismatch_map


logging.basicConfig(level=logging.INFO)
data_dir = Path(snakemake.config["data_dir"])
anc_ts = tskit.load(snakemake.input[0])
recomb_map = snakemake.input[-1]
(data_dir / "progress" / "match_samples").mkdir(parents=True, exist_ok=True)
Path(snakemake.output[0]).parent.mkdir(parents=True, exist_ok=True)
recombination_map, mismatch_map = build_maps(
    anc_ts, snakemake.wildcards.mismatch, recomb_map
)
tsinfer.match_samples_batch_init(
    Path(snakemake.output[0]).parent,
    snakemake.input[1].replace(".vcf_done", ""),
    "variant_ancestral_allele",
    snakemake.input[0],
    snakemake.config["match_samples"]["min_work_per_job"],
    sample_mask=f"sample_{snakemake.wildcards.subset_name}_subset_mask",
    site_mask=f"variant_{snakemake.wildcards.subset_name}_subset_{snakemake.wildcards.region_name}_region_{snakemake.wildcards.filter_set}_mask",
    path_compression=True,
    recombination=recombination_map,
    mismatch=mismatch_map,
    precision=15,
    post_process=False,
)
