from pathlib import Path

import numpy
import pandas as pd
import sgkit
import xarray
from common_functions import make_filter_key
from common_functions import sliding_window_density

ds_dir = snakemake.input[0].replace(".vcf_done", "")
ds = sgkit.load_dataset(ds_dir)
chunks = ds.variant_position.chunks
filter_config = snakemake.config["filters"][snakemake.wildcards.filter_set]

# Site density needs to be run after all other filters
if "site_density" in filter_config:
    # First create a site mask based on all the other filters
    all_positions = ds["variant_position"]
    all_filters_mask = xarray.full_like(all_positions, False, dtype=bool)
    for filter_name in set(filter_config) - {"site_density"}:
        all_filters_mask |= ds[
            make_filter_key(
                snakemake.wildcards.subset_name, filter_name, filter_config[filter_name]
            )
        ]
    all_filters_mask |= ds[f"variant_{snakemake.wildcards.region_name}_region_mask"]
    all_filters_mask = all_filters_mask.values
    all_positions = all_positions.values
    # Retrieve some snakemake.config
    site_density_config = filter_config["site_density"]
    window_size = site_density_config["window_size"]
    count_threshold = (
        site_density_config["threshold_sites_per_kbp"] / 1000
    ) * window_size

    used_sites_count = sliding_window_density(
        all_filters_mask, all_positions, window_size
    )

    site_density_mask = numpy.full_like(ds["variant_position"], False, dtype=bool)
    # If none of the sites are above the threshold, mask all sites
    # we have to do this as argmax returns 0 if there are no True values
    if not numpy.any(used_sites_count >= count_threshold):
        site_density_mask[:] = True
    else:
        # Find the start of the window where the count goes over the threshold
        # from the left
        start = numpy.argmax(used_sites_count >= count_threshold)

        # Find the start of the window where the count goes over the threshold
        # from the right
        end = (
            len(used_sites_count)
            - 1
            - numpy.argmax(used_sites_count[::-1] >= count_threshold)
        )

        # These values are relative to the start of the sites, so we need
        # to add the first site position to get the absolute position
        first_site = all_positions[0]
        start += first_site
        end += first_site

        # Find the index of start and end positions
        start = numpy.argmax(all_positions >= start)
        end = numpy.argmax(all_positions >= end)
        site_density_mask[:start] = True
        site_density_mask[end:] = True
    site_density_mask_key = make_filter_key(
        snakemake.wildcards.subset_name,
        "site_density",
        site_density_config,
        snakemake.wildcards.region_name,
        snakemake.wildcards.filter_set,
    )
    site_density_mask = xarray.DataArray(
        site_density_mask, dims=["variants"], name=site_density_mask_key
    )
    site_density_mask = site_density_mask.chunk(chunks).compute()
    ds.update({site_density_mask_key: site_density_mask})
    sgkit.save_dataset(
        ds.drop_vars(set(ds.data_vars) - {site_density_mask_key}),
        ds_dir,
        mode="a",
        consolidated=False,
    )

    # Find regions where the density is below the threshold
    below_threshold = used_sites_count < count_threshold
    start_indices = numpy.where(
        numpy.diff(numpy.concatenate(([False], below_threshold, [False])))
    )[0][::2]
    end_indices = numpy.where(
        numpy.diff(numpy.concatenate(([False], below_threshold, [False])))
    )[0][1::2]
    actual_starts = start_indices + first_site
    actual_ends = end_indices + first_site
    lengths = actual_ends - actual_starts
    low_density_data = pd.DataFrame(
        {
            "Start": actual_starts,
            "End": actual_ends,
            "Length": lengths,
        }
    )
    low_density_data.to_csv(
        f"{ds_dir}/{site_density_mask_key}_low_density_regions.csv", index=False
    )
Path(snakemake.output[0]).touch()
