import numpy
import xarray

SUBSET_INDEPENDENT_FILTERS = [
    "bad_ancestral",
    "no_ancestral_allele",
    "not_biallelic",
    "duplicate_position",
    "not_snps",
    "low_quality_ancestral_allele",
]

def bad_ancestral(ds, subset):
    # Filter sites that have a bad ancestral state
    aa = ds["variant_ancestral_allele"]
    wanted_variants = (aa != "-") & (aa != ".") & (aa != "N")
    # ".values" here as we can't index dask arrays with a dask boolean array
    assert set(numpy.unique(aa[wanted_variants.values])) == {"A", "C", "G", "T"}
    return ~wanted_variants


def no_ancestral_allele(ds, subset):
    # Filter sites where the ancestral state is not a seen allele
    aa = ds["variant_ancestral_allele"]
    alleles = ds["variant_allele"]
    mask = False
    for i in range(ds.dims["alleles"]):
        mask |= aa == alleles[:, i]
    return ~mask


def not_biallelic(ds, subset):
    # Mask sites that are not biallelic
    allele = ds["variant_allele"]
    return ~((allele != "").sum(dim="alleles") == 2)


def duplicate_position(ds, subset):
    # Mask sites that have duplicate positions
    pos = ds["variant_position"]
    pos_shift_left = xarray.full_like(pos, -1)
    pos_shift_left[0:-1] = pos[1:]
    pos_shift_right = xarray.full_like(pos, -1)
    pos_shift_right[1:] = pos[:-1]
    return ~((pos != pos_shift_left) & (pos != pos_shift_right))


def not_snps(ds, subset):
    # Mask sites that are not single character SNPs
    def check_lengths(arr):
        return numpy.all([len(str(s)) <= 1 for s in arr])

    return ~(xarray.apply_ufunc(
        check_lengths,
        ds["variant_allele"],
        input_core_dims=[["alleles"]],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[bool],
        dask_gufunc_kwargs={"allow_rechunk": True},
    ))

def low_quality_ancestral_allele(ds, subset):
    # Mask sites that have a lowercase ancestral allele - for
    # tsinfer we store all alleles as uppercase, so we use the
    # array stored earlier when the fasta was read
    return ds["variant_low_quality_ancestral_allele_mask"]

def low_allele_count(ds, subset, min_derived, min_ancestral):
    # Mask sites that have too few derived or ancestral alleles
    return ~((ds[f"variant_{subset}_subset_ancestral_count"] >= min_ancestral) & (
        ds[f"variant_{subset}_subset_derived_count"] >= min_derived
    ))



