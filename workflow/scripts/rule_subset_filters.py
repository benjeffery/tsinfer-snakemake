from pathlib import Path

import filters
import sgkit
from common_functions import make_filter_key

ds_dir = snakemake.input[0].replace(".vcf_done", "")
ds = sgkit.load_dataset(ds_dir)
# We don't need to subset here as the filters are using allele counts that
# are already subset
chunks = ds.variant_position.chunks
filter_config = snakemake.config["filters"][snakemake.wildcards.filter_set]
for filter_name in (set(filter_config) - {"site_density"}) - set(
    filters.SUBSET_INDEPENDENT_FILTERS
):
    filter_kwargs = filter_config[filter_name]
    filter_key = make_filter_key(
        snakemake.wildcards.subset_name, filter_name, filter_kwargs
    )
    mask = getattr(filters, filter_name)(
        ds, snakemake.wildcards.subset_name, **(filter_kwargs or {})
    )
    mask = mask.rename(filter_key).chunk(chunks)
    ds.update({filter_key: mask})
    sgkit.save_dataset(
        ds.drop_vars(set(ds.data_vars) - {filter_key}),
        ds_dir,
        mode="a",
        consolidated=False,
    )
Path(snakemake.output[0]).touch()
