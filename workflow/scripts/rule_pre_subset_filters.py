import sgkit
import xarray
from common_functions import filters
from common_functions import make_filter_key

ds_dir = snakemake.input[0].replace(".vcf_done", "")
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
