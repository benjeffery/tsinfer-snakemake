import sgkit
import xarray
from common_functions import make_filter_key

ds_dir = snakemake.input[0].replace(".vcf_done", "")
ds = sgkit.load_dataset(ds_dir)
chunks = ds.variant_position.chunks
filter_config = snakemake.config["filters"][snakemake.wildcards.filter_set]

final_mask = xarray.full_like(ds["variant_position"], False, dtype=bool)
for filter_name in set(filter_config.keys()) - {"site_density"}:
    final_mask |= ds[
        make_filter_key(
            snakemake.wildcards.subset_name,
            filter_name,
            filter_config[filter_name],
            snakemake.wildcards.region_name,
        )
    ].values
final_mask |= ds[f"variant_{snakemake.wildcards.region_name}_region_mask"].values
final_mask_key = (
    f"variant_{snakemake.wildcards.subset_name}_subset_{snakemake.wildcards.region_name}_"
    f"region_{snakemake.wildcards.filter_set}_mask"
)
final_mask = final_mask.rename(final_mask_key).chunk(chunks).compute()
ds.update({final_mask_key: final_mask})

sgkit.save_dataset(
    ds.drop_vars(set(ds.data_vars) - {final_mask_key}),
    ds_dir,
    mode="a",
    consolidated=False,
)
