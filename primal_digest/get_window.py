
"""
These are binary searches implementions
"""

def get_r_window_FAST2(kmers, start: int, end: int):
    """
    This will perform a binary search on the list of kmers.
    The Kmer start positions will be used as the search value
    """
    included_kmers = []
    n_kmers = len(kmers)
    high = n_kmers - 1
    low = 0
    mid = 0
    while low <= high:
        mid = (high + low) // 2
        # If the midpoint is inside the range
        if start <= kmers[mid].start <= end:
            while True:
                # Walk back until first value, or the first position
                if mid == 0:
                    break
                elif kmers[mid - 1].start >= start:
                    mid -= 1
                else:
                    break
            # Mid is now the first value so walk forwards
            while True:
                if mid < n_kmers and kmers[mid].start <= end:
                    included_kmers.append(kmers[mid])
                    mid += 1
                else:
                    return included_kmers
        # If start is greater ignore the left half
        elif kmers[mid].start < start:
            low = mid + 1
        # If start is smaller ignore the right half
        elif kmers[mid].start > end:
            high = mid - 1

    # If the code reaches here there are no KMERS within the list inside the range
    ## Return an empty list for continuity
    return []

    # These are modifed get_window_fast for the new kmers classes


def get_f_window_FAST2(kmers, start: int, end: int) :
    """
    This will perform a binary search on the list of kmers.
    The Kmer end position will be used as the search value
    """
    included_kmers = []
    n_kmers = len(kmers)
    high = n_kmers - 1
    low = 0
    mid = 0
    while low <= high:
        mid = (high + low) // 2
        # If the midpoint is inside the range
        if start <= kmers[mid].end <= end:
            while True:
                # Walk back until first value, or the first position
                if mid == 0:
                    break
                elif kmers[mid - 1].end >= start:
                    mid -= 1
                else:
                    break
            # Mid is now the first value so walk forwards
            while True:
                if mid < n_kmers and kmers[mid].end <= end:
                    included_kmers.append(kmers[mid])
                    mid += 1
                else:
                    return included_kmers
        # If start is greater ignore the left half
        elif kmers[mid].end < start:
            low = mid + 1
        # If start is smaller ignore the right half
        elif kmers[mid].end > end:
            high = mid - 1

    # If the code reaches here there are no KMERS within the list inside the range
    ## Return an empty list for continuity
    return []


def get_pp_window(pp, fp_start: int, fp_end: int, rp_end: int):
    """
    This will perform a semi-binary search on the list of primerpairs.
    The primerpair.fprimer.end position will be used as the search value
    """
    included_pp = []
    n_pp = len(pp)
    high = n_pp - 1
    low = 0
    mid = 0
    while low <= high:
        mid = (high + low) // 2
        # If the midpoint is inside the range
        if fp_start <= pp[mid].fprimer.end <= fp_end:
            while True:
                # Walk back until first value, or the first position
                if mid == 0:
                    break
                elif pp[mid - 1].fprimer.end >= fp_start:
                    mid -= 1
                else:
                    break
            # Mid is now the first value so walk forwards
            while mid < n_pp and min(pp[mid].fprimer.starts()) < fp_end:
                if pp[mid].rprimer.start > rp_end:
                    included_pp.append(pp[mid])
                    mid += 1
            return included_pp
        # If start is greater ignore the left half
        elif pp[mid].fprimer.end < fp_start:
            low = mid + 1
        # If start is smaller ignore the right half
        elif pp[mid].fprimer.end > fp_end:
            high = mid - 1

    # If the code reaches here there are no KMERS within the list inside the range
    ## Return an empty list for continuity
    return []

    # These are modifed get_window_fast for the new kmers classes