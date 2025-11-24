import json
import os

import matplotlib.pyplot as plt
import numpy
import pandas as pd
import sgkit
from common_functions import make_filter_key
from common_functions import sliding_window_density

ds_dir = snakemake.input[0].replace(".vcf_done", "")
ds = sgkit.load_dataset(ds_dir)
ds = ds.sel(
    samples=~ds[f"sample_{snakemake.wildcards.subset_name}_subset_mask"].values,
    variants=~ds[f"variant_{snakemake.wildcards.region_name}_region_mask"].values,
)
out = {}
out["dataset_summary"] = str(ds)
out["name"] = snakemake.wildcards.region_name
out["n_samples"] = ds.dims["samples"]
out["n_variants"] = ds.dims["variants"]
out["n_ploidy"] = ds.dims["ploidy"]
for filter_name in snakemake.config["filters"][snakemake.wildcards.filter_set]:
    filter_key = make_filter_key(
        snakemake.wildcards.subset_name,
        filter_name,
        snakemake.config["filters"][snakemake.wildcards.filter_set][filter_name],
        snakemake.wildcards.region_name,
        snakemake.wildcards.filter_set,
    )
    out[filter_key] = int((ds[filter_key]).sum())
out["sites_masked"] = int(
    (
        ds[
            f"variant_{snakemake.wildcards.subset_name}_subset_{snakemake.wildcards.region_name}_region_{snakemake.wildcards.filter_set}_mask"
        ]
    ).sum()
)
total_size = 0
for dirpath, _, filenames in os.walk(ds_dir):
    for f in filenames:
        fp = os.path.join(dirpath, f)
        total_size += os.path.getsize(fp)
out["size"] = total_size
with open(snakemake.output[0], "w") as f:
    f.write(json.dumps(out))

counts = [
    (
        "Ref allele",
        ds[f"variant_{snakemake.wildcards.subset_name}_subset_ref_count"].values,
    ),
    (
        "Ancestral allele",
        ds[f"variant_{snakemake.wildcards.subset_name}_subset_ancestral_count"].values,
    ),
    (
        "Missing allele",
        ds[f"variant_{snakemake.wildcards.subset_name}_subset_missing_count"].values,
    ),
    (
        "Derived allele",
        ds[f"variant_{snakemake.wildcards.subset_name}_subset_derived_count"].values,
    ),
]

# Plot the allele count spectrum
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111)

num_samples = ds.dims["samples"] * ds.dims["ploidy"]
for i, (label, count) in enumerate(counts):
    ax.hist(count, bins=200, log=True, histtype="step", label=label)
    text = ", ".join(
        f"{i}: {(count == i).sum()}"
        for i in [0, 1, 2, num_samples - 2, num_samples - 1, num_samples]
    )
    ax.text(
        0.75,
        0.95 - (i * 0.05),
        f"{label} {text}",
        fontsize=8,
        horizontalalignment="right",
        verticalalignment="top",
        transform=ax.transAxes,
    )
ax.set_title(
    f"Raw allele count spectrum - "
    f"{snakemake.wildcards.subset_name}-{snakemake.wildcards.region_name}-{snakemake.wildcards.filter_set}"
)
ax.set_xlabel("Allele count")
ax.set_ylabel("Number of sites")
ax.legend(loc="upper right")
fig.savefig(
    f"{snakemake.config['data_dir']}/zarr_stats/{snakemake.wildcards.subset_name}-{snakemake.wildcards.region_name}-{snakemake.wildcards.filter_set}/ac-raw.png"
)

fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111)

final_mask = ds[
    f"variant_{snakemake.wildcards.subset_name}_subset_{snakemake.wildcards.region_name}_region_{snakemake.wildcards.filter_set}_mask"
].values
for i, (label, count) in enumerate(counts):
    ax.hist(count[final_mask], bins=200, log=True, histtype="step", label=label)
    text = ", ".join(
        f"{i}: {(count[final_mask] == i).sum()}"
        for i in [0, 1, 2, num_samples - 2, num_samples - 1, num_samples]
    )
    ax.text(
        0.75,
        0.95 - (i * 0.05),
        f"{label} {text}",
        fontsize=8,
        horizontalalignment="right",
        verticalalignment="top",
        transform=ax.transAxes,
    )
ax.set_title(
    f"Filtered allele count spectrum - "
    f"{snakemake.wildcards.subset_name}-{snakemake.wildcards.region_name}-{snakemake.wildcards.filter_set}"
)
ax.set_xlabel("Allele count")
ax.set_ylabel("Number of sites")
ax.legend(loc="upper right")

fig.savefig(
    f"{snakemake.config['data_dir']}/zarr_stats/"
    f"{snakemake.wildcards.subset_name}-{snakemake.wildcards.region_name}-{snakemake.wildcards.filter_set}/"
    f"ac-filtered.png"
)

# Plot site density
fig = plt.figure(figsize=(20, 12))
ax = fig.add_subplot(111)

# First, get the total counts in each bin for all sites
counts_all, bins = numpy.histogram(ds.variant_position, bins=200)
counts_all = counts_all.astype(float)
counts_all[counts_all == 0] = numpy.nan
bin_centers = (bins[:-1] + bins[1:]) / 2

# Plot for each filter
for filter_name in snakemake.config["filters"][snakemake.wildcards.filter_set]:
    filter_key = make_filter_key(
        snakemake.wildcards.subset_name,
        filter_name,
        snakemake.config["filters"][snakemake.wildcards.filter_set][filter_name],
        snakemake.wildcards.region_name,
        snakemake.wildcards.filter_set,
    )
    mask = ds[filter_key].values
    counts_filter, _ = numpy.histogram(ds.variant_position[~mask], bins=bins)
    fraction = counts_filter / counts_all
    ax.plot(
        bin_centers,
        fraction,
        label=f"{filter_name} - {mask.sum() / ds.dims['variants']:.2f}",
    )

# Also plot for the final variant mask
counts_final_mask, _ = numpy.histogram(ds.variant_position[~final_mask], bins=bins)
fraction_final_mask = counts_final_mask / counts_all
ax.plot(
    bin_centers,
    fraction_final_mask,
    label=f"variant_mask - " f"{(numpy.sum(final_mask) / ds.dims['variants']):.2f}",
)

ax.set_title(
    f"Filters passing fraction - "
    f"{snakemake.wildcards.subset_name}-{snakemake.wildcards.region_name}-{snakemake.wildcards.filter_set}"
)
ax.set_xlabel("Position")
ax.set_ylabel("Fraction of sites passing")
# Put the legend outside the plot
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.85])
ax.legend(
    loc="upper center",
    bbox_to_anchor=(0.5, -0.1),
    fancybox=True,
    shadow=True,
    ncol=4,
)
# fig.tight_layout()
fig.savefig(
    f"{snakemake.config['data_dir']}/zarr_stats/"
    f"{snakemake.wildcards.subset_name}-{snakemake.wildcards.region_name}-{snakemake.wildcards.filter_set}/"
    f"filter-fractions.png"
)

filter_config = snakemake.config["filters"][snakemake.wildcards.filter_set]
# Plot raw site density after all filters, but before the site_density filter
# to inform threshold choices
other_filters_mask = numpy.full_like(ds["variant_position"], False, dtype=bool)
for filter_name in set(filter_config) - {"site_density"}:
    filter_key = make_filter_key(
        snakemake.wildcards.subset_name, filter_name, filter_config[filter_name]
    )
    other_filters_mask |= ds[filter_key].values

all_sites = ds.variant_position.values
fig = plt.figure(figsize=(20, 12))
ax = fig.add_subplot(111)

try:
    window_size = filter_config["site_density"]["window_size"]
except KeyError:
    # Default to a sensible value
    window_size = 1000
all_sites_count = sliding_window_density(
    numpy.full_like(all_sites, True, dtype=bool), all_sites, window_size
)
normalised_all_sites_count = (all_sites_count / window_size) * 1000
filtered_sites_count = sliding_window_density(
    other_filters_mask, all_sites, window_size
)
normalised_filtered_sites_count = (filtered_sites_count / window_size) * 1000

df = pd.DataFrame(
    {
        "position": numpy.arange(min(all_sites), max(all_sites) + 2 - window_size),
        "all_sites_count": normalised_all_sites_count,
        "filtered_sites_count": normalised_filtered_sites_count,
    }
)
summary_window_size = len(df) // 10000

rolling_all_sites = df["all_sites_count"].rolling(
    window=summary_window_size, min_periods=1
)
all_sites_stats = rolling_all_sites.agg(["max"])
rolling_filtered_sites = df["filtered_sites_count"].rolling(
    window=summary_window_size, min_periods=1
)
filtered_sites_stats = rolling_filtered_sites.agg(["max"])
summary_df = pd.concat([df["position"], all_sites_stats, filtered_sites_stats], axis=1)

# Optionally rename columns for clarity
summary_df.columns = ["position", "all_sites_max", "filtered_sites_max"]

ax.plot(summary_df["position"], summary_df["all_sites_max"], label="All sites")
ax.plot(summary_df["position"], summary_df["filtered_sites_max"], label="Passing sites")

# If site_density is in the filter snakemake.config, plot the threshold
if "site_density" in filter_config:
    # Plot the threshold
    threshold = filter_config["site_density"]["threshold_sites_per_kbp"]
    ax.axhline(
        threshold,
        color="red",
        label=f"Site density threshold ({threshold:.2f} sites/kbp)",
    )
ax.set_title(
    f"Site density ({window_size} bp sliding window) - "
    f"{snakemake.wildcards.subset_name}-{snakemake.wildcards.region_name}-{snakemake.wildcards.filter_set}"
)
ax.set_xlabel("Position")
ax.set_ylabel("Number of sites per kb")
# Put the legend outside the plot
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.85])
ax.legend(
    loc="upper center",
    bbox_to_anchor=(0.5, -0.1),
    fancybox=True,
    shadow=True,
)
fig.savefig(
    f"{snakemake.config['data_dir']}/zarr_stats/"
    f"{snakemake.wildcards.subset_name}-{snakemake.wildcards.region_name}-{snakemake.wildcards.filter_set}/"
    f"site-density.png"
)
