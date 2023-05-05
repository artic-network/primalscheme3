from math import sqrt

def ol_pp_score(pp_r_start: int, pp_n_p: int,leading_edge: int, min_overlap: int) -> int:
    return -(pp_r_start - min_overlap - leading_edge) ** 2 / sqrt(pp_n_p)

def walk_pp_score(pp_f_end: int, pp_n_p: int, leading_edge: int) -> int:
    return (pp_f_end- leading_edge) * sqrt(pp_n_p)

    



