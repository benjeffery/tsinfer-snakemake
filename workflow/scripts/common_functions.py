# Imports are function-local as imports can be very expensive on an NFS
# mount that has high latency.
try:
    import filters
except ImportError:
    from scripts import filters


def parse_region(region):
    chrom, rest = region.split(":")
    start, stop = map(int, rest.split("-"))
    return chrom, start, stop


def make_filter_key(
    subset_name, filter_name, filter_kwargs=None, region_name=None, filter_set=None
):
    filter_kwargs = filter_kwargs or {}
    ret = "variant_"
    if filter_name not in filters.SUBSET_INDEPENDENT_FILTERS:
        ret += f"{subset_name}_subset_"
    if filter_name in filters.REGION_DEPENDENT_FILTERS:
        assert region_name is not None
        assert filter_set is not None
        ret += f"{region_name}_region_"
        ret += f"{filter_set}_"
    ret += f"{filter_name}_"
    if filter_kwargs is not None and len(filter_kwargs) > 0:
        ret += "_".join([f"{k}_{v}" for k, v in sorted(filter_kwargs.items())]) + "_"
    return ret + "mask"


def sliding_window_density(mask, positions, window_size):
    import numpy

    # Make an array of positions of used sites
    used_sites_positions = positions[mask]
    # Create a boolean array covering each base between the start and end of all
    # the sites, marking bases with sites as True
    first_site = positions[0]
    last_site = positions[-1]
    bool_array = numpy.zeros(last_site - first_site + 1, dtype=bool)
    bool_array[used_sites_positions - first_site] = True
    # Convolve the boolean array with a window of size window_size to get the
    # number of sites in each sliding window
    window = numpy.ones(window_size)
    used_sites_count = numpy.convolve(bool_array, window, mode="valid")
    return used_sites_count
