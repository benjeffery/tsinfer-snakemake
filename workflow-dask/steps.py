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


def vcf_to_zarrs(input, output, wildcards, config, params):
    import sgkit.io.vcf
    from pathlib import Path
    from dask.distributed import Client
    from dask.distributed import performance_report

    with Client(config["scheduler_address"]) as client, performance_report(
            filename=output[1]):
        sgkit.io.vcf.vcf_reader.vcf_to_zarr(
            input[0],
            output[0].replace(".vcf_done", ""),
            target_part_size=params.target_part_size,
            read_chunk_length=params.read_chunk_length,
            temp_chunk_length=params.temp_chunk_length,
            chunk_length=params.chunk_length,
            chunk_width=params.chunk_width,
            tempdir=config["temp_dir"],
            retain_temp_files=params.retain_temp_files
        )
    Path(str(output[0])).touch()

def load_ancestral_fasta(input, output, wildcards, config, params):
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
    ancestral_states = numpy.char.upper(ancestral_sequence[ds['variant_position'].values-1])
    ancestral_states = xarray.DataArray(data=ancestral_states, dims=["variants"], name="variant_ancestral_allele")
    ds.update({"variant_ancestral_allele": ancestral_states})
    sgkit.save_dataset(ds.drop_vars(set(ds.data_vars) - {"variant_ancestral_allele"}), ds_dir, mode="a")

def mask_sites(input, output, wildcards, config, params):
    from xarray import full_like
    import numpy
    import sgkit

    ds_dir = input[0].replace(".vcf_done", "")
    ds = sgkit.load_dataset(ds_dir)
    chunks = ds.variant_position.chunks

    # Mask bad ancestral sites
    wanted_variants = (ds['variant_ancestral_allele'] != '-') & (
                ds['variant_ancestral_allele'] != '.') & (
                                  ds['variant_ancestral_allele'] != 'N')
    wanted_variants = wanted_variants.rename("variant_bad_ancestral_mask").chunk(
        chunks).compute()
    assert set(numpy.unique(ds['variant_ancestral_allele'][wanted_variants])) == {'A',
                                                                                  'C',
                                                                                  'G',
                                                                                  'T'}
    ds.update({"variant_bad_ancestral_mask": wanted_variants})
    sgkit.save_dataset(ds.drop_vars(set(ds.data_vars) - {"variant_bad_ancestral_mask"}),
                       ds_dir, mode="a")

    # Mask sites where the ancestral state is not an allele
    wanted_variants = ((ds['variant_ancestral_allele'] == ds['variant_allele'][:, 0]) |
                       (ds['variant_ancestral_allele'] == ds['variant_allele'][:, 1]))
    wanted_variants = wanted_variants.rename("variant_no_ancestral_allele_mask").chunk(
        chunks).compute()
    ds.update({"variant_no_ancestral_allele_mask": wanted_variants})
    sgkit.save_dataset(
        ds.drop_vars(set(ds.data_vars) - {"variant_no_ancestral_allele_mask"}), ds_dir,
        mode="a")

    # Mask duplicate sites with duplicate position
    pos = ds['variant_position']
    pos_shift_left = full_like(pos, -1)
    pos_shift_left[0:-1] = pos[1:]
    pos_shift_right = full_like(pos, -1)
    pos_shift_right[1:] = pos[:-1]
    wanted_variants = (pos != pos_shift_left) & (pos != pos_shift_right)
    wanted_variants = wanted_variants.rename("variant_duplicate_position_mask").chunk(
        chunks).compute()
    ds.update({"variant_duplicate_position_mask": wanted_variants})
    sgkit.save_dataset(
        ds.drop_vars(set(ds.data_vars) - {"variant_duplicate_position_mask"}), ds_dir,
        mode="a")

    ## Create the combined mask - tsinfer now copes with missing and bad ancestral sites so don't use those masks
    wanted_variants = ds['variant_duplicate_position_mask']
    wanted_variants = wanted_variants.rename("variant_mask").chunk(chunks).compute()
    ds.update({"variant_mask": wanted_variants})
    sgkit.save_dataset(ds.drop_vars(set(ds.data_vars) - {"variant_mask"}), ds_dir,
                       mode="a")



def subset_zarr_vcf(input, output, wildcards, config, params):
    import numpy
    import sgkit
    from pathlib import Path
    from dask.distributed import Client

    with Client(config["scheduler_address"]) as client:
        chrom, start, end = parse_region(config['regions'][wildcards.region_name])
        ds = sgkit.load_dataset(input[0].replace(".vcf_done", ""))
        with open(input[-1], 'r') as f:
            sample_ids = numpy.genfromtxt(f, dtype=str)
        sample_mask = numpy.isin(ds.sample_id.values, sample_ids)
        variant_position = ds['variant_position'].values
        variant_mask = numpy.logical_and(variant_position >= start, variant_position < end)
        ds = ds.sel(samples=sample_mask, variants=variant_mask)
        ds = ds.unify_chunks()
        sgkit.save_dataset(ds, output[0].replace(".vcf_done", ""), auto_rechunk=True)
        Path(str(output[0])).touch()

def zarr_stats(input, output, wildcards, config, params):
    import sgkit
    import json
    import numpy
    import os
    from dask.distributed import Client
    import matplotlib.pyplot as plt
    import xarray as xr

    ds_dir = input[0].replace(".vcf_done", "")
    with Client(config["scheduler_address"]) as client:
        ds = sgkit.load_dataset(ds_dir)
        out = {}
        out['dataset_summary'] = str(ds)
        out['name'] = wildcards.region_name
        out['n_samples'] = ds.dims['samples']
        out['n_variants'] = ds.dims['variants']
        out['n_ploidy'] = ds.dims['ploidy']
        # Flatten the ploidy dimension as tsinfer sees each phased haplotype as a sample
        gt = (ds.call_genotype.stack(haplotype=("samples", "ploidy")))
        out['allele_counts'] = (gt > 0).sum(dim=['haplotype']).to_numpy().tolist()
        out['missing_count'] = int((gt == -1).sum())
        out['num_sites_triallelic'] = int(((gt > 1).sum(dim=["haplotype"]) > 0).sum())
        out['sites_bad_ancestral'] = int((~ds.variant_bad_ancestral_mask).sum())
        out['sites_no_ancestral'] = int((~ds.variant_no_ancestral_allele_mask).sum())
        out['sites_duplicate_pos'] = int((~ds.variant_duplicate_position_mask).sum())
        out['sites_masked'] = int((~ds.variant_mask).sum())
        total_size = 0
        for dirpath, dirnames, filenames in os.walk(ds_dir):
            for f in filenames:
                fp = os.path.join(dirpath, f)
                total_size += os.path.getsize(fp)
        out['size'] = total_size
        with open(output[0], "w") as f:
            f.write(json.dumps(out))

        # Store the allele counts for later filtering
        # Note that no missingness is assumed
        ac = sgkit.count_call_alleles(ds, merge=False)['call_allele_count'].sum(dim="samples").values
        num_samples = ds.dims['samples'] * ds.dims['ploidy']

        # Normal non-ref allele count:
        ref_ac = ac[:, 0]

        # Calculate the ancestral allele index
        allele_matches = (ds['variant_allele'] == ds['variant_ancestral_allele']).values
        ancestral_indices = numpy.argmax(allele_matches, axis=1)
        # Mark unknown ancestral alleles as REF
        ancestral_indices[numpy.sum(allele_matches, axis=1) == 0] = 0

        # Calculate the allele count for the ancestral allele
        aa_ac = ac[numpy.arange(len(ancestral_indices)), ancestral_indices]

        #Plot the allele count spectrum
        plt.hist(ref_ac, bins=200, log=True, histtype='step', label="Ref allele count")
        plt.hist(aa_ac, bins=200, log=True, histtype='step', label="Ancestral allele count")
        plt.title('Allele count spectrum')
        plt.xlabel('Allele count')
        plt.ylabel('Number of sites')
        plt.legend()
        plt.savefig(f"{config['data_dir']}/zarr_stats/{wildcards.subset_name}-{wildcards.region_name}/ac.png")

        # Save the ancestral allele counts to the dataset
        aa_ac = xr.DataArray(aa_ac, dims=["variants"], name="variant_ancestral_allele_counts")
        ds.update({"variant_ancestral_allele_counts": aa_ac})
        sgkit.save_dataset(
             ds.drop_vars(set(ds.data_vars) - {"variant_ancestral_allele_counts"}),
         ds_dir, mode="a")


def summary_table(input, output, wildcards, config, params):
    import csv
    import json
    import os

    header = ["region_name", "vcf_size", "zarr_size", "n_variants", "n_samples", "ac==0",
              "ac==1",
              "ac==2", "missing_genotypes", "num_sites_triallelic",
              "sites_bad_ancestral", "sites_no_ancestral",
              "sites_duplicate_pos", "sites_masked", "inference_nbytes",
              "inference_bitpack_nbytes"]
    with open(output[0], "w", newline='') as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()

        for vcf_stats in input:
            with open(vcf_stats, "r") as json_stats_f:
                stats = json.load(json_stats_f)
                ac0 = sum([ac == 0 for ac in stats['allele_counts']])
                ac1 = sum([ac == 1 for ac in stats['allele_counts']])
                ac2 = sum([ac == 2 for ac in stats['allele_counts']])
                n_sites = stats['n_variants']
                n_samples = stats['n_samples']
                n_ploidy = stats['n_ploidy']
                n_masked = stats['sites_masked']
                chrom, start, stop = parse_region(config["regions"][stats['name']])
                row_dict = {
                    "region_name": stats['name'],
                    "vcf_size": os.path.getsize(config['vcf'].format(chrom=chrom)),
                    "zarr_size": stats['size'],
                    "n_variants": n_sites,
                    "n_samples": n_samples,
                    "ac==0": ac0,
                    "ac==1": ac1,
                    "ac==2": ac2,
                    "missing_genotypes": stats['missing_count'],
                    "num_sites_triallelic": stats['num_sites_triallelic'],
                    "sites_bad_ancestral": stats['sites_bad_ancestral'],
                    "sites_no_ancestral": stats['sites_no_ancestral'],
                    "sites_duplicate_pos": stats['sites_duplicate_pos'],
                    "sites_masked": n_masked,
                    "inference_nbytes": ((
                                                     n_sites - n_masked) - ac1) * n_samples * n_ploidy,
                    "inference_bitpack_nbytes": (((
                                                              n_sites - n_masked) - ac1) * n_samples * n_ploidy) / 8
                }
                row_dict = {k: number_to_SI(v) for k, v in row_dict.items()}
                writer.writerow(row_dict)


def match_ancestors(input, output, wildcards, config, threads):
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
    data_dir = Path(config['data_dir'])
    ancestors = tsinfer.AncestorData.load(input[0])
    sample_data = tsinfer.SgkitSampleData(input[1])
    subset_name = wildcards.subset_name
    region_name = wildcards.region_name
    os.makedirs(data_dir / "progress" / "match_ancestors", exist_ok=True)
    os.makedirs(data_dir / "resume" / "match_ancestors", exist_ok=True)
    os.makedirs(os.path.dirname(output[0]), exist_ok=True)
    slug = f"{subset_name}-{region_name}-truncate-{wildcards.lower}-{wildcards.upper}-{wildcards.multiplier}"
    async def run_match_ancestors_with_workers():
        # We want to be able to use local CPUs for the small epochs, and then scale out for the large ones.
        # When scaling out this would leave the local CPUs idle so we launch local workers too
        os.makedirs(data_dir / "local_worker_logs", exist_ok=True)
        with open(f"{data_dir}/local_worker_logs/{slug}-{str(uuid.uuid4())}.log",
                  "w") as log_file:
            worker_process = subprocess.Popen(
                ["dask-worker", config['scheduler_address']],
                stdout=log_file,
                stderr=subprocess.STDOUT,
            )
            try:
                with open(
                        data_dir / "progress" / "match_ancestors" / f"{slug}.log",
                        "w") as log_f, Client(config["scheduler_address"]) as client, performance_report(filename=output[1]):
                    # Disable dask's aggressive memory management
                    client.amm.stop()
                    ts = tsinfer.match_ancestors(
                        sample_data,
                        ancestors,
                        path_compression=True,
                        num_threads=threads,
                        precision=15,
                        progress_monitor=tsinfer.progress.ProgressMonitor(
                            tqdm_kwargs={'file': log_f, 'mininterval': 30}),
                        resume_lmdb_file=str(
                            data_dir / "resume" / "match_ancestors" / f"{slug}.lmdb"),
                        use_dask=True,
                    )
                ts.dump(output[0])
            finally:
                worker_process.terminate()
                worker_process.wait()

    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    loop.run_until_complete(run_match_ancestors_with_workers())

def match_samples(input, output, wildcards, config, threads):
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
    data_dir = Path(config['data_dir'])
    anc_ts = tskit.load(input[0])
    sample_data = tsinfer.SgkitSampleData(input[1])
    os.makedirs(data_dir / "progress" / "match_samples", exist_ok=True)
    os.makedirs(data_dir / "resume" / "match_samples", exist_ok=True)
    os.makedirs(os.path.dirname(output[0]), exist_ok=True)
    mismatch = float(wildcards.mismatch)
    if mismatch > 0:
        inference_pos = anc_ts.tables.sites.position
        rate_map = msprime.RateMap.read_hapmap(
            input[-1],
            position_col=1,
            rate_col=2
        )
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
                mismatch_ratio,
                numpy.median(genetic_dists),
                num_alleles
            )
        )
    else:
        recombination_map = None
        mismatch_map = None
    slug= f"{wildcards.subset_name}-{wildcards.region_name}-mm-{wildcards.mismatch}-truncate-{wildcards.lower}-{wildcards.upper}-{wildcards.multiplier}"
    with open(data_dir / "progress" / "match_samples" / f"{slug}.log","w") as log_f,        Client(config["scheduler_address"]) as client,        performance_report(filename=output[1]):
        #Disable dask's aggressive memory management
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
                tqdm_kwargs={'file': log_f, 'mininterval': 30}),
            resume_lmdb_file=str(
                data_dir / "resume" / "match_samples" / f"{slug}.lmdb"),
            use_dask=True,
        )
    ts.dump(output[0])
