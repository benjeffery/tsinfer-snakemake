import numpy
import sgkit
import xarray

ds_dir = snakemake.input[0].replace(".vcf_done", "")
ds = sgkit.load_dataset(ds_dir)
sample_mask = ds[f"sample_{snakemake.wildcards.subset_name}_subset_mask"].values
subset_ds = ds.sel(samples=~sample_mask)
ac = sgkit.count_call_alleles(subset_ds)["call_allele_count"].sum(dim="samples").values

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
    array_name = f"variant_{snakemake.wildcards.subset_name}_subset_{var_name}"
    array = xarray.DataArray(array, dims=["variants"], name=array_name)
    ds.update({array_name: array})
    sgkit.save_dataset(
        ds.drop_vars(set(ds.data_vars) - {array_name}),
        ds_dir,
        mode="a",
        consolidated=False,
    )
