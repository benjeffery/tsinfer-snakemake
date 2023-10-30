# Imports are function-local as imports can be very expensive on an NFS
# mount that has high latency.


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


def vcf_to_zarrs(input, output, wildcards, config, params):  # noqa: A002
    import sgkit.io.vcf
    from pathlib import Path
    from dask.distributed import Client
    from dask.distributed import performance_report

    with Client(config["scheduler_address"]), performance_report(filename=output[1]):
        sgkit.io.vcf.vcf_reader.vcf_to_zarr(
            input[0],
            output[0].replace(".vcf_done", ""),
            target_part_size=params.target_part_size,
            read_chunk_length=params.read_chunk_length,
            temp_chunk_length=params.temp_chunk_length,
            chunk_length=params.chunk_length,
            chunk_width=params.chunk_width,
            tempdir=config["temp_dir"],
            retain_temp_files=params.retain_temp_files,
        )
    Path(str(output[0])).touch()


def load_ancestral_fasta(input, output, wildcards, config, params):  # noqa: A002
    import pyfaidx
    import numpy
    import sgkit
    import xarray

    fasta = pyfaidx.Fasta(input[0])
    seq_name = list(fasta.keys())[0]
    print("Ancestral sequence:", seq_name)
    ancestral_sequence = numpy.asarray(fasta[seq_name], dtype="U1")

    ds_dir = input[1].replace(".vcf_done", "")
    ds = sgkit.load_dataset(ds_dir)
    ancestral_states = ancestral_sequence[ds["variant_position"].values - 1]
    ancestral_states_upper = numpy.char.upper(ancestral_states)
    low_quality_mask = ancestral_states == ancestral_states_upper

    ancestral_states = xarray.DataArray(
        data=ancestral_states_upper, dims=["variants"], name="variant_ancestral_allele"
    )
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
    )


def pre_subset_filters(input, output, wildcards, config, params):  # noqa: A002
    import sgkit
    import filters
    from distributed import Client

    with Client(config["scheduler_address"]):
        ds_dir = input[0].replace(".vcf_done", "")
        ds = sgkit.load_dataset(ds_dir)
        chunks = ds.variant_position.chunks

        for filter_name in [
            "not_snps",
            "bad_ancestral",
            "no_ancestral_allele",
            "not_biallelic",
            "duplicate_position",
        ]:
            mask = getattr(filters, filter_name)(ds)
            # Rename to match sgkit convention
            filter_name = f"variant_{filter_name}_mask"
            mask = mask.rename(filter_name).chunk(chunks).compute()
            ds.update({filter_name: mask})
            sgkit.save_dataset(
                ds.drop_vars(set(ds.data_vars) - {filter_name}), ds_dir, mode="a"
            )


def subset_zarr_vcf(input, output, wildcards, config, params):  # noqa: A002
    import xarray
    import numpy
    import sgkit
    from pathlib import Path
    from distributed import Client

    with Client(config["scheduler_address"]):
        chrom, start, end = parse_region(config["regions"][wildcards.region_name])
        ds = sgkit.load_dataset(input[0].replace(".vcf_done", ""))
        with open(input[-1]) as f:
            sample_ids = numpy.genfromtxt(f, dtype=str)
        sample_ids = xarray.DataArray(sample_ids, dims="sample")
        sample_mask = ds.sample_id.isin(sample_ids)
        variant_mask = (ds["variant_position"] >= start) & (
            ds["variant_position"] < end
        )
        ds = ds.sel(samples=sample_mask.values, variants=variant_mask.values)
        ds = ds.unify_chunks()
        sgkit.save_dataset(ds, output[0].replace(".subset_done", ""), auto_rechunk=True)
        Path(str(output[0])).touch()


def post_subset_filters(input, output, wildcards, config, params):  # noqa: A002
    import numpy
    import sgkit
    import filters
    import xarray
    from distributed import Client

    with Client(config["scheduler_address"]) as client:
        ds_dir = input[0].replace(".subset_done", "")
        ds = sgkit.load_dataset(ds_dir)

        # Some filters need allele counts
        ac = (
            sgkit.count_call_alleles(ds, merge=True)["call_allele_count"]
            .sum(dim="samples")
            .values
        )  # .values here so we can use numpy indexing below
        sgkit.save_dataset(
            ds.drop_vars(set(ds.data_vars) - {"variant_allele_count"}), ds_dir, mode="a"
        )

        num_samples = ds.dims["samples"] * ds.dims["ploidy"]

        # Normal non-ref allele count:
        ref_count = ac[:, 0]

        # Calculate the ancestral allele index
        allele_matches = (ds["variant_allele"] == ds["variant_ancestral_allele"]).values
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
            array = xarray.DataArray(
                array, dims=["variants"], name=f"variant_{var_name}"
            )
            ds.update({f"variant_{var_name}": array})
            sgkit.save_dataset(
                ds.drop_vars(set(ds.data_vars) - {f"variant_{var_name}"}),
                ds_dir,
                mode="a",
            )

        # Calculate the missing filters
        # (note could use dask better here by combining filters into one call)
        chunks = ds.variant_position.chunks
        filter_config = config["filters"][wildcards.filter]
        for filter_name, filter_kwargs in filter_config.items():
            if filter_name not in ds.keys():
                mask = getattr(filters, filter_name)(ds, **(filter_kwargs or {}))
                # Rename to match sgkit convention
                filter_name = f"variant_{filter_name}_mask"
                mask = mask.rename(filter_name).chunk(chunks).compute()
                ds.update({filter_name: mask})
                sgkit.save_dataset(
                    ds.drop_vars(set(ds.data_vars) - {filter_name}), ds_dir, mode="a"
                )

        mask = xarray.full_like(ds["variant_position"], True, dtype=bool)
        for filter_name in filter_config:
            filter_name = f"variant_{filter_name}_mask"
            mask &= ds[filter_name]
        mask = mask.rename("variant_mask").chunk(chunks).compute()
        ds.update({"variant_mask": mask})
        sgkit.save_dataset(
            ds.drop_vars(set(ds.data_vars) - {"variant_mask"}), ds_dir, mode="a"
        )


def zarr_stats(input, output, wildcards, config, params):  # noqa: A002
    import sgkit
    import json
    import os
    from dask.distributed import Client
    import matplotlib.pyplot as plt

    ds_dir = input[0].replace("variant_mask", "")
    with Client(config["scheduler_address"]):
        ds = sgkit.load_dataset(ds_dir)
        out = {}
        out["dataset_summary"] = str(ds)
        out["name"] = wildcards.region_name
        out["n_samples"] = ds.dims["samples"]
        out["n_variants"] = ds.dims["variants"]
        out["n_ploidy"] = ds.dims["ploidy"]
        # Flatten the ploidy dimension as tsinfer sees each phased haplotype as a sample
        gt = ds.call_genotype.stack(haplotype=("samples", "ploidy"))
        out["allele_counts"] = (gt > 0).sum(dim=["haplotype"]).to_numpy().tolist()
        out["missing_count"] = int((gt == -1).sum())
        for filter_name in config["filters"][wildcards.filter]:
            filter_name = f"variant_{filter_name}_mask"
            out[filter_name] = int((~ds[filter_name]).sum())
        out["sites_masked"] = int((~ds.variant_mask).sum())
        total_size = 0
        for dirpath, _, filenames in os.walk(ds_dir):
            for f in filenames:
                fp = os.path.join(dirpath, f)
                total_size += os.path.getsize(fp)
        out["size"] = total_size
        with open(output[0], "w") as f:
            f.write(json.dumps(out))

        counts = [
            ("Ref allele", ds.variant_ref_count),
            ("Ancestral allele", ds.variant_ancestral_count),
            ("Derived allele", ds.variant_derived_count),
            ("Missing allele", ds.variant_missing_count),
        ]

        # Plot the allele count spectrum
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111)

        num_samples = ds.dims["samples"] * ds.dims["ploidy"]
        for i, (label, count) in enumerate(counts):
            ax.hist(count, bins=200, log=True, histtype="step", label=label)
            text = ", ".join(
                f"{i}: {(count == i).sum().values}"
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
            f"{wildcards.subset_name}-{wildcards.region_name}-{wildcards.filter}"
        )
        ax.set_xlabel("Allele count")
        ax.set_ylabel("Number of sites")
        ax.legend(loc="upper right")
        fig.savefig(
            f"{config['data_dir']}/zarr_stats/{wildcards.subset_name}-{wildcards.region_name}-{wildcards.filter}/ac-raw.png"
        )

        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111)

        mask = ds.variant_mask.values
        for i, (label, count) in enumerate(counts):
            ax.hist(count[mask], bins=200, log=True, histtype="step", label=label)
            text = ", ".join(
                f"{i}: {(count[mask] == i).sum().values}"
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
            f"{wildcards.subset_name}-{wildcards.region_name}-{wildcards.filter}"
        )
        ax.set_xlabel("Allele count")
        ax.set_ylabel("Number of sites")
        ax.legend(loc="upper right")
        # Add a text box with the number of sites that have an allele count
        # of 0, 1, 2, n-2, n-1, n where n is the number of samples

        fig.savefig(
            f"{config['data_dir']}/zarr_stats/"
            f"{wildcards.subset_name}-{wildcards.region_name}-{wildcards.filter}/"
            f"ac-filtered.png"
        )

        # Plot site density
        fig = plt.figure(figsize=(20, 12))
        ax = fig.add_subplot(111)
        ax.hist(
            ds.variant_position, bins=200, log=True, histtype="step", label="All Sites"
        )
        for filter_name in config["filters"][wildcards.filter]:
            mask = ds[f"variant_{filter_name}_mask"].values
            ax.hist(
                ds.variant_position[mask],
                bins=200,
                log=True,
                histtype="step",
                label=f"{filter_name} - {mask.sum()/ds.dims['variants']:.2f}",
            )

        ax.hist(
            ds.variant_position[ds.variant_mask.values],
            bins=200,
            log=True,
            histtype="step",
            label=f"variant_mask - "
            f"{(ds.variant_mask.sum() / ds.dims['variants']).values:.2f}",
        )
        ax.set_title(
            f"Site density - "
            f"{wildcards.subset_name}-{wildcards.region_name}-{wildcards.filter}"
        )
        ax.set_xlabel("Position")
        ax.set_ylabel("Number of sites passing")
        box = ax.get_position()
        ax.set_position(
            [box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.85]
        )
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
            f"{wildcards.subset_name}-{wildcards.region_name}-{wildcards.filter}/"
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
        "ac==0",
        "ac==1",
        "ac==2",
        "missing_genotypes",
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
                ac0 = sum([ac == 0 for ac in stats["allele_counts"]])
                ac1 = sum([ac == 1 for ac in stats["allele_counts"]])
                ac2 = sum([ac == 2 for ac in stats["allele_counts"]])
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
                    "ac==0": ac0,
                    "ac==1": ac1,
                    "ac==2": ac2,
                    "missing_genotypes": stats["missing_count"],
                    "sites_masked": n_masked,
                    "inference_nbytes": ((n_sites - n_masked) - ac1)
                    * n_samples
                    * n_ploidy,
                    "inference_bitpack_nbytes": (
                        ((n_sites - n_masked) - ac1) * n_samples * n_ploidy
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
    sample_data = tsinfer.SgkitSampleData(input[0].replace("variant_mask", ""))
    os.makedirs(data_dir / "progress" / "generate_ancestors", exist_ok=True)
    with open(
        data_dir
        / "progress"
        / "generate_ancestors"
        / f"{wildcards.subset_name}-{wildcards.region_name}.log",
        "w",
    ) as log_f:
        tsinfer.generate_ancestors(
            sample_data,
            path=output[0],
            genotype_encoding=1,
            num_threads=threads,
            progress_monitor=tsinfer.progress.ProgressMonitor(
                tqdm_kwargs={"file": log_f, "mininterval": 30}
            ),
        )

def ancestor_stats(input, output, wildcards, config):  # noqa: A002
    import tsinfer
    from pathlib import Path
    import matplotlib.pyplot as plt

    data_dir = Path(config["data_dir"])
    ancestors = tsinfer.AncestorData.load(input[0])

    time = ancestors.ancestors_time[:]
    focal_sites = ancestors.ancestors_focal_sites[:]
    start = ancestors.ancestors_start[:]
    end = ancestors.ancestors_end[:]

    dpi = 3840 // 16
    fig, ax = plt.subplots(figsize=(16, 9), dpi=dpi)
    ax.hlines(time, start, end, color='blue', alpha=0.5, linewidth=.5)
    # Flatten focal sites and times for scatter plot
    all_focal_sites = [site for sublist in focal_sites for site in sublist]
    all_times = [time for time, focal_sites in zip(time, focal_sites) for _
                 in focal_sites]
    # Plot the focal sites using rasterization for speed
    ax.scatter(all_focal_sites, all_times, color='red', s=.5, rasterized=True)

    ax.set_xlabel("Genetic Position")
    ax.set_ylabel("Time")
    ax.set_ylim(0, 1)  # Setting y-axis limits to [0, 1]
    ax.set_title("Ancestral Genetic Sequences")

    # Save the image with 4K resolution
    plt.tight_layout()
    plt.savefig(output[0], dpi=dpi)
    plt.close()

def match_ancestors(
    input, output, wildcards, config, threads, params, slug  # noqa: A002
):
    import tsinfer
    import logging
    import os
    import asyncio
    import uuid
    import subprocess
    from pathlib import Path
    from dask.distributed import Client
    from dask.distributed import performance_report

    logging.basicConfig(level=logging.INFO)
    data_dir = Path(config["data_dir"])
    ancestors = tsinfer.AncestorData.load(input[0])
    sample_data = tsinfer.SgkitSampleData(input[1].replace("variant_mask", ""))
    os.makedirs(data_dir / "progress" / "match_ancestors", exist_ok=True)
    os.makedirs(data_dir / "resume" / "match_ancestors", exist_ok=True)
    os.makedirs(os.path.dirname(output[0]), exist_ok=True)

    async def run_match_ancestors_with_workers():
        # We want to be able to use local CPUs for the small epochs, and then scale
        # out for the large ones. When scaling out this would leave the local CPUs idle
        # so we launch local workers too
        os.makedirs(data_dir / "local_worker_logs", exist_ok=True)
        with open(
            f"{data_dir}/local_worker_logs/{slug}-{str(uuid.uuid4())}.log", "w"
        ) as log_file:
            worker_process = subprocess.Popen(
                ["dask-worker", config["scheduler_address"]],
                stdout=log_file,
                stderr=subprocess.STDOUT,
            )
            try:
                with open(
                    data_dir / "progress" / "match_ancestors" / f"{slug}.log", "w"
                ) as log_f, Client(
                    config["scheduler_address"]
                ) as client, performance_report(
                    filename=output[1]
                ):
                    # Disable dask's aggressive memory management
                    client.amm.stop()
                    ts = tsinfer.match_ancestors(
                        sample_data,
                        ancestors,
                        path_compression=True,
                        num_threads=threads,
                        precision=15,
                        progress_monitor=tsinfer.progress.ProgressMonitor(
                            tqdm_kwargs={"file": log_f, "mininterval": 30}
                        ),
                        resume_lmdb_file=str(
                            data_dir / "resume" / "match_ancestors" / f"{slug}.lmdb"
                        ),
                        use_dask=params.use_dask,
                    )
                ts.dump(output[0])
            finally:
                worker_process.terminate()
                worker_process.wait()

    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    loop.run_until_complete(run_match_ancestors_with_workers())


def match_samples(
    input, output, wildcards, config, threads, params, slug  # noqa: A002
):
    import tsinfer
    import logging
    import tskit
    import numpy
    import msprime
    from pathlib import Path
    import os
    from dask.distributed import Client
    from dask.distributed import performance_report

    logging.basicConfig(level=logging.INFO)
    data_dir = Path(config["data_dir"])
    anc_ts = tskit.load(input[0])
    sample_data = tsinfer.SgkitSampleData(input[1].replace("variant_mask", ""))
    os.makedirs(data_dir / "progress" / "match_samples", exist_ok=True)
    os.makedirs(data_dir / "resume" / "match_samples", exist_ok=True)
    os.makedirs(os.path.dirname(output[0]), exist_ok=True)
    mismatch = float(wildcards.mismatch)
    if mismatch > 0:
        inference_pos = anc_ts.tables.sites.position
        rate_map = msprime.RateMap.read_hapmap(input[-1], position_col=1, rate_col=2)
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
    with open(
        data_dir / "progress" / "match_samples" / f"{slug}.log", "w"
    ) as log_f, Client(config["scheduler_address"]) as client, performance_report(
        filename=output[1]
    ):
        # Disable dask's aggressive memory management
        client.amm.stop()
        ts = tsinfer.match_samples(
            sample_data,
            anc_ts,
            path_compression=True,
            num_threads=threads,
            recombination=recombination_map,
            mismatch=mismatch_map,
            precision=15,
            progress_monitor=tsinfer.progress.ProgressMonitor(
                tqdm_kwargs={"file": log_f, "mininterval": 30}
            ),
            resume_lmdb_file=str(
                data_dir / "resume" / "match_samples" / f"{slug}.lmdb"
            ),
            use_dask=params.use_dask,
            post_process=False,
        )
    ts.dump(output[0])
