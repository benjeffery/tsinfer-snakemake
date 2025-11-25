import numpy
import pyfaidx
import sgkit
import xarray

fasta = pyfaidx.Fasta(snakemake.input[0])
seq_name = list(fasta.keys())[0]
ancestral_sequence = numpy.asarray(fasta[seq_name], dtype="U1")

ds_dir = snakemake.input[1].replace(".vcf_done", "")
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
