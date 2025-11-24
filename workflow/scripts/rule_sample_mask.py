import numpy
import sgkit
import xarray

ds_dir = snakemake.input[0].replace(".vcf_done", "")
ds = sgkit.load_dataset(ds_dir)
with open(snakemake.input[-1]) as f:
    sample_ids = numpy.genfromtxt(f, dtype=str)
sample_ids = xarray.DataArray(sample_ids, dims="sample")
sample_mask = ds.sample_id.isin(sample_ids)
if sample_mask.sum() != len(sample_ids):
    raise ValueError(
        f"Could not find all samples in dataset. "
        f"Failed to find {sample_ids[~sample_mask].values}"
    )
sample_mask_name = f"sample_{snakemake.wildcards.subset_name}_subset_mask"
array = xarray.DataArray(~sample_mask, dims=["samples"], name=sample_mask_name)
ds.update({sample_mask_name: array})
sgkit.save_dataset(
    ds.drop_vars(set(ds.data_vars) - {sample_mask_name}),
    ds_dir,
    mode="a",
    consolidated=False,
)
