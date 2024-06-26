from math import sqrt

from primalscheme3.core.config import Config


def ol_pp_score(
    pp_r_start: int, pp_n_p: int, leading_edge: int, config: Config
) -> float:
    """
    Higher score is better
    """
    dist_extend = pp_r_start - config.min_overlap - leading_edge
    prop_extended = dist_extend / config.amplicon_size_max

    return prop_extended**2 / sqrt(pp_n_p)


def bt_ol_pp_score(
    pp_r_start: int, pp_n_p: int, leading_edge: int, config: Config
) -> float:
    """
    Higher score is better
    """
    dist_extend = pp_r_start - config.min_overlap - leading_edge
    prop_extended = dist_extend / config.amplicon_size_max

    return prop_extended**2 / (sqrt(pp_n_p) / 2)


def walk_pp_score(pp_f_end: int, pp_n_p: int, leading_edge: int) -> float:
    """Higher score is better"""

    # If this is this is positive the primerpair is before the leading edge (HALF COVERAGE!!!)
    dist_from_leading_edge = leading_edge - pp_f_end

    # If the primerpair is before the leading edge favor it massively
    if dist_from_leading_edge > 0:
        return dist_from_leading_edge / sqrt(pp_n_p)
    else:
        return dist_from_leading_edge * sqrt(pp_n_p)
