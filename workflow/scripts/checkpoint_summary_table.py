import csv
import json
import os

from common_functions import parse_region


def number_to_SI(number):
    """Convert a number to a string with SI units, unless it is a string already."""
    units = ["", "K", "M", "G", "T", "P", "E", "Z", "Y"]
    unit = 0
    if isinstance(number, str):
        return number
    while number > 1000:
        number /= 1000
        unit += 1
    return f"{number:.2f}{units[unit]}"


header = [
    "region_name",
    "vcf_size",
    "zarr_size",
    "n_variants",
    "n_samples",
    "num_sites_triallelic",
    "sites_bad_ancestral",
    "sites_no_ancestral",
    "sites_duplicate_pos",
    "sites_masked",
    "inference_nbytes",
    "inference_bitpack_nbytes",
]
with open(snakemake.output[0], "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=header)
    writer.writeheader()

    for vcf_stats in snakemake.input:
        with open(vcf_stats) as json_stats_f:
            stats = json.load(json_stats_f)
            n_sites = stats["n_variants"]
            n_samples = stats["n_samples"]
            n_ploidy = stats["n_ploidy"]
            n_masked = stats["sites_masked"]
            chrom, start, stop = parse_region(
                snakemake.config["regions"][stats["name"]]
            )
            row_dict = {
                "region_name": stats["name"],
                "vcf_size": os.path.getsize(
                    snakemake.config["vcf"].format(chrom=chrom)
                ),
                "zarr_size": stats["size"],
                "n_variants": n_sites,
                "n_samples": n_samples,
                "sites_masked": n_masked,
                "inference_nbytes": (n_sites - n_masked) * n_samples * n_ploidy,
                "inference_bitpack_nbytes": (
                    (n_sites - n_masked) * n_samples * n_ploidy
                )
                / 8,
            }
            row_dict = {k: number_to_SI(v) for k, v in row_dict.items()}
            writer.writerow(row_dict)
