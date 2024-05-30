# Imports are function-local as imports can be very expensive on an NFS
# mount that has high latency.
import filters


def parse_region(region):
    chrom, rest = region.split(":")
    start, stop = map(int, rest.split("-"))
    return chrom, start, stop


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


def load_ancestral_fasta(input, output, wildcards, config, params):  # noqa: A002
    import pyfaidx
    import numpy
    import sgkit
    import xarray

    fasta = pyfaidx.Fasta(input[0])
    seq_name = list(fasta.keys())[0]
    ancestral_sequence = numpy.asarray(fasta[seq_name], dtype="U1")

    ds_dir = input[1].replace(".vcf_done", "")
    ds = sgkit.load_dataset(ds_dir)
    ancestral_states = ancestral_sequence[ds["variant_position"].values - 1]
    ancestral_states_upper = numpy.char.upper(ancestral_states)
    low_quality_mask = ancestral_states != ancestral_states_upper

    ancestral_states = xarray.DataArray(
        data=ancestral_states_upper, dims=["variants"], name="variant_ancestral_allele"
    )
    # Store this now as we won't have lowercase alleles in the zarr
    # as tsinfer sees them as different alleles
    low_quality_mask = xarray.DataArray(
        data=low_quality_mask,
        dims=["variants"],
        name="variant_low_quality_ancestral_allele_mask",
    )
    ds.update(
        {
            "variant_ancestral_allele": ancestral_states,
            "variant_low_quality_ancestral_allele_mask": low_quality_mask,
        }
    )
    sgkit.save_dataset(
        ds.drop_vars(
            set(ds.data_vars)
            - {"variant_ancestral_allele", "variant_low_quality_ancestral_allele_mask"}
        ),
        ds_dir,
        mode="a",
        consolidated=False,
    )


def make_filter_key(
    subset_name, filter_name, filter_kwargs=None, region_name=None, filter_set=None
):
    filter_kwargs = filter_kwargs or {}
    ret = "variant_"
    if filter_name not in filters.SUBSET_INDEPENDENT_FILTERS:
        ret += f"{subset_name}_subset_"
    if filter_name in filters.REGION_DEPENDENT_FILTERS:
        assert region_name is not None
        assert filter_set is not None
        ret += f"{region_name}_region_"
        ret += f"{filter_set}_"
    ret += f"{filter_name}_"
    if filter_kwargs is not None and len(filter_kwargs) > 0:
        ret += "_".join([f"{k}_{v}" for k, v in sorted(filter_kwargs.items())]) + "_"
    return ret + "mask"


def pre_subset_filters(input, output, wildcards, config, params):  # noqa: A002
    import sgkit
    import xarray

    ds_dir = input[0].replace(".vcf_done", "")
    ds = sgkit.load_dataset(ds_dir)
    chunks = ds.variant_position.chunks

    for filter_name in filters.SUBSET_INDEPENDENT_FILTERS:
        mask = getattr(filters, filter_name)(ds, None)
        # Rename to match sgkit convention
        filter_key = make_filter_key(None, filter_name)
        mask = mask.rename(filter_key).chunk(chunks)
        # Materialize the mask to avoid blosc decompression error
        mask = mask.values
        mask = xarray.DataArray(mask, dims=["variants"], name=filter_key)
        ds.update({filter_key: mask})
        sgkit.save_dataset(
            ds.drop_vars(set(ds.data_vars) - {filter_key}),
            ds_dir,
            mode="a",
            consolidated=False,
        )


def sliding_window_density(mask, positions, window_size):
    import numpy

    # Make an array of positions of used sites
    used_sites_positions = positions[mask]
    # Create a boolean array covering each base between the start and end of all
    # the sites, marking bases with sites as True
    first_site = positions[0]
    last_site = positions[-1]
    bool_array = numpy.zeros(last_site - first_site + 1, dtype=bool)
    bool_array[used_sites_positions - first_site] = True
    # Convolve the boolean array with a window of size window_size to get the
    # number of sites in each sliding window
    window = numpy.ones(window_size)
    used_sites_count = numpy.convolve(bool_array, window, mode="valid")
    return used_sites_count


def region_mask(input, output, wildcards, config, params):  # noqa: A002
    import sgkit

    ds_dir = input[0].replace(".vcf_done", "")
    ds = sgkit.load_dataset(ds_dir)
    chrom, start, end = parse_region(config["regions"][wildcards.region_name])
    mask_name = f"variant_{wildcards.region_name}_region_mask"
    # get the index of the contig
    contig_index = (
        ds["contig_id"].values.tolist().index(config["contig_name"].format(chrom=chrom))
    )
    mask = (
        (ds["variant_contig"] != contig_index)
        | (ds["variant_position"] < start)
        | (ds["variant_position"] >= end)
    )
    mask = mask.rename(mask_name)
    ds.update({mask_name: mask})
    sgkit.save_dataset(
        ds.drop_vars(set(ds.data_vars) - {mask_name}),
        ds_dir,
        mode="a",
        consolidated=False,
    )


def sample_mask(input, output, wildcards, config, params):  # noqa: A002
    import numpy
    import sgkit
    import xarray

    ds_dir = input[0].replace(".vcf_done", "")
    ds = sgkit.load_dataset(ds_dir)
    with open(input[-1]) as f:
        sample_ids = numpy.genfromtxt(f, dtype=str)
    sample_ids = xarray.DataArray(sample_ids, dims="sample")
    sample_mask = ds.sample_id.isin(sample_ids)
    if sample_mask.sum() != len(sample_ids):
        raise ValueError(
            f"Could not find all samples in dataset. "
            f"Failed to find {sample_ids[~sample_mask].values}"
        )
    sgkit_samples_mask_name = f"sample_{wildcards.subset_name}_subset_mask"
    array = xarray.DataArray(
        ~sample_mask, dims=["samples"], name=sgkit_samples_mask_name
    )
    ds.update({sgkit_samples_mask_name: array})
    sgkit.save_dataset(
        ds.drop_vars(set(ds.data_vars) - {sgkit_samples_mask_name}),
        ds_dir,
        mode="a",
        consolidated=False,
    )


def allele_counts(input, output, wildcards, config, params):  # noqa: A002
    import sgkit
    import numpy
    import xarray

    ds_dir = input[0].replace(".vcf_done", "")
    ds = sgkit.load_dataset(ds_dir)
    sample_mask = ds[f"sample_{wildcards.subset_name}_subset_mask"].values
    subset_ds = ds.sel(samples=~sample_mask)
    ac = (
        sgkit.count_call_alleles(subset_ds)["call_allele_count"]
        .sum(dim="samples")
        .values
    )

    num_samples = subset_ds.dims["samples"] * subset_ds.dims["ploidy"]
    ref_count = ac[:, 0]

    # Calculate the ancestral allele index
    allele_matches = (
        subset_ds["variant_allele"] == subset_ds["variant_ancestral_allele"]
    ).values
    ancestral_indices = numpy.argmax(allele_matches, axis=1)
    # Mark unknown ancestral alleles as REF. This is just for plots not for inference
    ancestral_indices[numpy.sum(allele_matches, axis=1) == 0] = 0
    # Use the index to get the ancestral allele count
    ancestral_count = ac[numpy.arange(len(ancestral_indices)), ancestral_indices]

    total_ac = ac.sum(axis=1)
    missing_count = num_samples - total_ac
    derived_count = total_ac - missing_count - ancestral_count

    for var_name in [
        "ref_count",
        "ancestral_count",
        "missing_count",
        "derived_count",
    ]:
        array = locals()[var_name]
        # Convert to an xarray DataArray
        array_name = f"variant_{wildcards.subset_name}_subset_{var_name}"
        array = xarray.DataArray(array, dims=["variants"], name=array_name)
        ds.update({array_name: array})
        sgkit.save_dataset(
            ds.drop_vars(set(ds.data_vars) - {array_name}),
            ds_dir,
            mode="a",
            consolidated=False,
        )


def subset_filters(input, output, wildcards, config, params):  # noqa: A002
    import sgkit
    import filters
    from pathlib import Path

    ds_dir = input[0].replace(".vcf_done", "")
    ds = sgkit.load_dataset(ds_dir)
    # We don't need to subset here as the filters are using allele counts that
    # are already subset
    chunks = ds.variant_position.chunks
    filter_config = config["filters"][wildcards.filter_set]
    for filter_name in (set(filter_config) - {"site_density"}) - set(
        filters.SUBSET_INDEPENDENT_FILTERS
    ):
        filter_kwargs = filter_config[filter_name]
        filter_key = make_filter_key(wildcards.subset_name, filter_name, filter_kwargs)
        mask = getattr(filters, filter_name)(
            ds, wildcards.subset_name, **(filter_kwargs or {})
        )
        mask = mask.rename(filter_key).chunk(chunks)
        ds.update({filter_key: mask})
        sgkit.save_dataset(
            ds.drop_vars(set(ds.data_vars) - {filter_key}),
            ds_dir,
            mode="a",
            consolidated=False,
        )
    Path(output[0]).touch()


def site_density_mask(input, output, wildcards, config, params):  # noqa: A002
    import xarray
    import numpy
    import sgkit
    import pandas as pd
    from pathlib import Path

    ds_dir = input[0].replace(".vcf_done", "")
    ds = sgkit.load_dataset(ds_dir)
    chunks = ds.variant_position.chunks
    filter_config = config["filters"][wildcards.filter_set]

    # Site density needs to be run after all other filters
    if "site_density" in filter_config:
        # First create a site mask based on all the other filters
        all_positions = ds["variant_position"]
        all_filters_mask = xarray.full_like(all_positions, False, dtype=bool)
        for filter_name in set(filter_config) - {"site_density"}:
            all_filters_mask |= ds[
                make_filter_key(
                    wildcards.subset_name, filter_name, filter_config[filter_name]
                )
            ]
        all_filters_mask |= ds[f"variant_{wildcards.region_name}_region_mask"]
        all_filters_mask = all_filters_mask.values
        all_positions = all_positions.values
        # Retrieve some config
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
            wildcards.subset_name,
            "site_density",
            site_density_config,
            wildcards.region_name,
            wildcards.filter_set,
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
    Path(output[0]).touch()


def combined_mask(input, output, wildcards, config, params):  # noqa: A002
    import sgkit
    import xarray

    ds_dir = input[0].replace(".vcf_done", "")
    ds = sgkit.load_dataset(ds_dir)
    chunks = ds.variant_position.chunks
    filter_config = config["filters"][wildcards.filter_set]

    final_mask = xarray.full_like(ds["variant_position"], False, dtype=bool)
    for filter_name in set(filter_config.keys()) - {"site_density"}:
        final_mask |= ds[
            make_filter_key(
                wildcards.subset_name,
                filter_name,
                filter_config[filter_name],
                wildcards.region_name,
            )
        ].values
    final_mask |= ds[f"variant_{wildcards.region_name}_region_mask"].values
    final_mask_key = (
        f"variant_{wildcards.subset_name}_subset_{wildcards.region_name}_"
        f"region_{wildcards.filter_set}_mask"
    )
    final_mask = final_mask.rename(final_mask_key).chunk(chunks).compute()
    ds.update({final_mask_key: final_mask})

    sgkit.save_dataset(
        ds.drop_vars(set(ds.data_vars) - {final_mask_key}),
        ds_dir,
        mode="a",
        consolidated=False,
    )


def zarr_stats(input, output, wildcards, config, params):  # noqa: A002
    import sgkit
    import json
    import os
    import numpy
    import matplotlib.pyplot as plt
    import pandas as pd

    ds_dir = input[0].replace(".vcf_done", "")
    ds = sgkit.load_dataset(ds_dir)
    ds = ds.sel(
        samples=~ds[f"sample_{wildcards.subset_name}_subset_mask"].values,
        variants=~ds[f"variant_{wildcards.region_name}_region_mask"].values,
    )
    out = {}
    out["dataset_summary"] = str(ds)
    out["name"] = wildcards.region_name
    out["n_samples"] = ds.dims["samples"]
    out["n_variants"] = ds.dims["variants"]
    out["n_ploidy"] = ds.dims["ploidy"]
    for filter_name in config["filters"][wildcards.filter_set]:
        filter_key = make_filter_key(
            wildcards.subset_name,
            filter_name,
            config["filters"][wildcards.filter_set][filter_name],
            wildcards.region_name,
            wildcards.filter_set,
        )
        out[filter_key] = int((ds[filter_key]).sum())
    out["sites_masked"] = int(
        (
            ds[
                f"variant_{wildcards.subset_name}_subset_{wildcards.region_name}_region_{wildcards.filter_set}_mask"
            ]
        ).sum()
    )
    total_size = 0
    for dirpath, _, filenames in os.walk(ds_dir):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            total_size += os.path.getsize(fp)
    out["size"] = total_size
    with open(output[0], "w") as f:
        f.write(json.dumps(out))

    counts = [
        ("Ref allele", ds[f"variant_{wildcards.subset_name}_subset_ref_count"].values),
        (
            "Ancestral allele",
            ds[f"variant_{wildcards.subset_name}_subset_ancestral_count"].values,
        ),
        (
            "Missing allele",
            ds[f"variant_{wildcards.subset_name}_subset_missing_count"].values,
        ),
        (
            "Derived allele",
            ds[f"variant_{wildcards.subset_name}_subset_derived_count"].values,
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
        f"{wildcards.subset_name}-{wildcards.region_name}-{wildcards.filter_set}"
    )
    ax.set_xlabel("Allele count")
    ax.set_ylabel("Number of sites")
    ax.legend(loc="upper right")
    fig.savefig(
        f"{config['data_dir']}/zarr_stats/{wildcards.subset_name}-{wildcards.region_name}-{wildcards.filter_set}/ac-raw.png"
    )

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)

    final_mask = ds[
        f"variant_{wildcards.subset_name}_subset_{wildcards.region_name}_region_{wildcards.filter_set}_mask"
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
        f"{wildcards.subset_name}-{wildcards.region_name}-{wildcards.filter_set}"
    )
    ax.set_xlabel("Allele count")
    ax.set_ylabel("Number of sites")
    ax.legend(loc="upper right")

    fig.savefig(
        f"{config['data_dir']}/zarr_stats/"
        f"{wildcards.subset_name}-{wildcards.region_name}-{wildcards.filter_set}/"
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
    for filter_name in config["filters"][wildcards.filter_set]:
        filter_key = make_filter_key(
            wildcards.subset_name,
            filter_name,
            config["filters"][wildcards.filter_set][filter_name],
            wildcards.region_name,
            wildcards.filter_set,
        )
        mask = ds[filter_key].values
        counts_filter, _ = numpy.histogram(ds.variant_position[~mask], bins=bins)
        fraction = counts_filter / counts_all
        ax.plot(
            bin_centers,
            fraction,
            label=f"{filter_name} - {mask.sum()/ds.dims['variants']:.2f}",
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
        f"{wildcards.subset_name}-{wildcards.region_name}-{wildcards.filter_set}"
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
        f"{config['data_dir']}/zarr_stats/"
        f"{wildcards.subset_name}-{wildcards.region_name}-{wildcards.filter_set}/"
        f"filter-fractions.png"
    )

    filter_config = config["filters"][wildcards.filter_set]
    # Plot raw site density after all filters, but before the site_density filter
    # to inform threshold choices
    other_filters_mask = numpy.full_like(ds["variant_position"], False, dtype=bool)
    for filter_name in set(filter_config) - {"site_density"}:
        filter_key = make_filter_key(
            wildcards.subset_name, filter_name, filter_config[filter_name]
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
    summary_df = pd.concat(
        [df["position"], all_sites_stats, filtered_sites_stats], axis=1
    )

    # Optionally rename columns for clarity
    summary_df.columns = ["position", "all_sites_max", "filtered_sites_max"]

    ax.plot(summary_df["position"], summary_df["all_sites_max"], label="All sites")
    ax.plot(
        summary_df["position"], summary_df["filtered_sites_max"], label="Passing sites"
    )

    # If site_density is in the filter config, plot the threshold
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
        f"{wildcards.subset_name}-{wildcards.region_name}-{wildcards.filter_set}"
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
        f"{config['data_dir']}/zarr_stats/"
        f"{wildcards.subset_name}-{wildcards.region_name}-{wildcards.filter_set}/"
        f"site-density.png"
    )


def summary_table(input, output, wildcards, config, params):  # noqa: A002
    import csv
    import json
    import os

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
    with open(output[0], "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()

        for vcf_stats in input:
            with open(vcf_stats) as json_stats_f:
                stats = json.load(json_stats_f)
                n_sites = stats["n_variants"]
                n_samples = stats["n_samples"]
                n_ploidy = stats["n_ploidy"]
                n_masked = stats["sites_masked"]
                chrom, start, stop = parse_region(config["regions"][stats["name"]])
                row_dict = {
                    "region_name": stats["name"],
                    "vcf_size": os.path.getsize(config["vcf"].format(chrom=chrom)),
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


def generate_ancestors(input, output, wildcards, config, threads):  # noqa: A002
    import tsinfer
    import logging
    import os
    from pathlib import Path

    logging.basicConfig(level=logging.INFO)
    data_dir = Path(config["data_dir"])
    sample_data = tsinfer.SgkitSampleData(
        input[0].replace(".vcf_done", ""),
        sgkit_samples_mask_name=f"sample_{wildcards.subset_name}_subset_mask",
        sites_mask_name=f"variant_{wildcards.subset_name}_subset_{wildcards.region_name}_region_{wildcards.filter_set}_mask",
    )
    os.makedirs(data_dir / "progress" / "generate_ancestors", exist_ok=True)
    with open(
        data_dir
        / "progress"
        / "generate_ancestors"
        / f"{wildcards.subset_name}-{wildcards.region_name}-{wildcards.filter_set}.log",
        "w",
    ) as log_f:
        ancestors = tsinfer.generate_ancestors(
            sample_data,
            path=output[0],
            genotype_encoding=1,
            num_threads=threads,
            progress_monitor=tsinfer.progress.ProgressMonitor(
                tqdm_kwargs={"file": log_f, "mininterval": 30}
            ),
        )
    if ancestors.num_ancestors == 0:
        raise ValueError("No ancestors generated")
    if ancestors.num_sites == 0:
        raise ValueError("No sites generated")


def build_maps(anc_ts, mismatch, recomb_map):
    import msprime
    import tsinfer
    import numpy

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


def match_sample_path(input, output, wildcards, config, threads, params):  # noqa: A002
    import tsinfer
    import tskit
    import os
    import logging

    logging.basicConfig(level=logging.INFO)
    anc_ts = tskit.load(input[0])
    sample_data = tsinfer.SgkitSampleData(
        input[1].replace(".vcf_done", ""),
        sgkit_samples_mask_name=f"sample_{wildcards.subset_name}_subset_mask",
        sites_mask_name=f"variant_{wildcards.subset_name}_subset_{wildcards.region_name}_region_{wildcards.filter_set}_mask",
    )
    recomb_map = input[-1]
    os.makedirs(os.path.dirname(output[0]), exist_ok=True)
    recombination_map, mismatch_map = build_maps(anc_ts, wildcards.mismatch, recomb_map)
    tsinfer.match_samples_slice_to_disk(
        sample_data,
        anc_ts,
        (int(wildcards.sample_index_start), int(wildcards.sample_index_end) + 1),
        output[0],
        path_compression=True,
        recombination=recombination_map,
        mismatch=mismatch_map,
        precision=15,
    )


def match_samples(
    input, output, wildcards, config, threads, params, slug  # noqa: A002
):
    import tsinfer
    import logging
    import tskit
    from pathlib import Path
    import os

    logging.basicConfig(level=logging.INFO)
    data_dir = Path(config["data_dir"])
    anc_ts = tskit.load(input[0])
    recomb_map = input[-1]
    sample_data = tsinfer.SgkitSampleData(
        input[1].replace(".vcf_done", ""),
        sgkit_samples_mask_name=f"sample_{wildcards.subset_name}_subset_mask",
        sites_mask_name=f"variant_{wildcards.subset_name}_subset_{wildcards.region_name}_region_{wildcards.filter_set}_mask",
    )
    print(sample_data.num_samples)
    os.makedirs(data_dir / "progress" / "match_samples", exist_ok=True)
    os.makedirs(os.path.dirname(output[0]), exist_ok=True)
    recombination_map, mismatch_map = build_maps(anc_ts, wildcards.mismatch, recomb_map)
    with open(data_dir / "progress" / "match_samples" / f"{slug}.log", "w") as log_f:
        ts = tsinfer.match_samples(
            sample_data,
            anc_ts,
            match_data_dir=os.path.dirname(input[-2]),
            path_compression=True,
            num_threads=threads,
            recombination=recombination_map,
            mismatch=mismatch_map,
            precision=15,
            progress_monitor=tsinfer.progress.ProgressMonitor(
                tqdm_kwargs={"file": log_f, "mininterval": 30}
            ),
            post_process=False,
        )
    ts.dump(output[0])
