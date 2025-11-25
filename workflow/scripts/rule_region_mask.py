import sgkit
from common_functions import parse_region

ds_dir = snakemake.input[0].replace(".vcf_done", "")
ds = sgkit.load_dataset(ds_dir)
chrom, start, end = parse_region(
    snakemake.config["regions"][snakemake.wildcards.region_name]
)
mask_name = f"variant_{snakemake.wildcards.region_name}_region_mask"
# get the index of the contig
contig_index = (
    ds["contig_id"]
    .values.tolist()
    .index(snakemake.config["contig_name"].format(chrom=chrom))
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
