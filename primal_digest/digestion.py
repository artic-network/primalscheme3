# Modules
from primal_digest.thermo import calc_tm
from primal_digest.config import ALL_DNA
from primal_digest.classes import RKmer, FKmer
from primal_digest.thermo import *
from primal_digest.seq_functions import expand_ambs, get_most_common_base

# Externals
from typing import Iterable
from itertools import product
import numpy as np
from collections import Counter
from multiprocessing import Pool


def walk_right(
    array: np.ndarray,
    col_index_right: int,
    col_index_left: int,
    row_index: int,
    seq_str: str,
    cfg: dict,
) -> set[str] | None:
    passing_str = set()

    if calc_tm(seq_str, cfg) < cfg["primer_tm_min"]:
        # Check the primer cannot walk out of array size
        if col_index_right < array.shape[1] - 1:
            new_base = array[row_index, col_index_right]
        else:
            return None

        # Fix incomplete ends
        if new_base == "":
            new_base = get_most_common_base(array, col_index_right + 1)
        new_string = (seq_str + new_base).replace("-", "")

        # If Sequence contains N prevent it from being added
        if new_string.__contains__("N"):
            return None
        else:
            exp_new_string: set[str] = expand_ambs([new_string])

        passing_str = set()
        for exp_str in exp_new_string:
            # Try/expect is to catch max Recursion depth Error fairly cleanly
            try:
                results = walk_right(
                    array,
                    col_index_right + 1,
                    col_index_left,
                    row_index,
                    exp_str,
                    cfg,
                )
            except RecursionError:
                return None

            if results is not None:
                [passing_str.add(x) for x in results]
            else:
                return None
    else:
        return {seq_str}

    if len(passing_str) >= 1:
        return passing_str


def walk_left(
    array: np.ndarray,
    col_index_right: int,
    col_index_left: int,
    row_index: int,
    seq_str: str,
    cfg: dict,
) -> set[str] | None:
    """
    This will take a string and its indexes and will recurvisly walk left
    until either tm is reached
    """
    passing_str = set()
    if calc_tm(seq_str, cfg) < cfg["primer_tm_min"]:
        # Prevents walking out of array size
        if col_index_left > 0:
            new_base = array[row_index, col_index_left - 1]
        else:
            return None

        # Ensure it can repair truncated regions
        if new_base == "":
            new_base = get_most_common_base(array, col_index_left - 1)
        new_string = (new_base + seq_str).replace("-", "")

        # Prevents seqs with an N
        if new_string.__contains__("N"):
            return None
        else:
            exp_new_string: set[str] = expand_ambs([new_string])

        passing_str = set()
        for exp_str in exp_new_string:
            try:
                results = walk_left(
                    array,
                    col_index_right,
                    col_index_left - 1,
                    row_index,
                    exp_str,
                    cfg,
                )
            except RecursionError:
                return None

            if results is not None:
                [passing_str.add(x) for x in results]
            else:
                return None
    else:
        return {seq_str}

    if len(passing_str) >= 1:
        return passing_str


def mp_r_digest(data: tuple[np.ndarray, dict, int, int]) -> RKmer | None:
    align_array: np.ndarray = data[0]
    cfg: dict = data[1]
    start_col = data[2]
    offset = data[3]

    # Check that there are no gaps at the first base
    if "-" in align_array[:, start_col]:
        return None

    total_col_seqs = set()
    for row_index in range(0, align_array.shape[0]):
        start_array = align_array[row_index, start_col : start_col + 14]
        start_seq = "".join(start_array).replace("-", "")

        if (
            start_array[0] == ""
        ):  # If the kmer starts on an invalid base, skip this start position
            total_col_seqs.add(None)
            break

        if not start_seq:  # If the start seq is empty go to the next row
            continue

        seqs = walk_right(
            array=align_array,
            col_index_right=start_col + 14,
            col_index_left=start_col,
            row_index=row_index,
            seq_str=start_seq,
            cfg=cfg,
        )
        # If the seq contains N, expand_amps will return {} which is parsed up the stack as None
        if seqs == None:
            total_col_seqs.add(None)
            break

        # Get a union of the seqs
        if seqs:
            total_col_seqs = total_col_seqs | seqs

    if not total_col_seqs.__contains__(None):
        # Thermo check the kmers
        if thermo_check_kmers(total_col_seqs, cfg) and not forms_hairpin(
            total_col_seqs, cfg=cfg
        ):
            tmp_kmer = RKmer(start=start_col + offset, seqs=total_col_seqs)
            tmp_kmer = tmp_kmer.reverse_complement()
            return tmp_kmer
    else:
        return None


def mp_f_digest(data: tuple[np.ndarray, dict, int, int]) -> FKmer | None:
    align_array: np.ndarray = data[0]
    cfg: dict = data[1]
    end_col = data[2]
    offset = data[3]

    # Check that there are no gaps at the first base
    if "-" in align_array[:, end_col]:
        return None

    total_col_seqs = set()
    for row_index in range(0, align_array.shape[0]):
        start_seq = "".join(align_array[row_index, end_col - 14 : end_col]).replace(
            "-", ""
        )

        if not start_seq:  # If the start seq is empty go to the next row
            continue

        seqs = walk_left(
            array=align_array,
            col_index_right=end_col,
            col_index_left=end_col - 14,
            row_index=row_index,
            seq_str=start_seq,
            cfg=cfg,
        )
        # If the seq contains N, expand_amps will return {} which is parsed up the stack as None
        if seqs == None:
            total_col_seqs.add(None)
            break

            # Get a union of the seqs
        if seqs:
            total_col_seqs = total_col_seqs | seqs

    if not total_col_seqs.__contains__(None):
        # Thermo check the kmers
        if thermo_check_kmers(total_col_seqs, cfg) and not forms_hairpin(
            total_col_seqs, cfg=cfg
        ):
            return FKmer(end=end_col + offset, seqs=total_col_seqs)
        else:
            return None


def digest(
    msa_array, cfg, thermo_cfg, offset=0
) -> tuple[list[FKmer | None], list[RKmer | None]]:
    with Pool(cfg["n_cores"]) as p:
        fprimer_mp = p.map(
            mp_f_digest,
            [
                (msa_array, thermo_cfg, end_col, offset)
                for end_col in range(14, msa_array.shape[1])
            ],
        )
    mp_thermo_pass_fkmers = [x for x in fprimer_mp if x is not None]
    mp_thermo_pass_fkmers.sort(key=lambda fkmer: fkmer.end)

    # RPrimers digestion via MP
    with Pool(cfg["n_cores"]) as p:
        rprimer_mp = p.map(
            mp_r_digest,
            [
                (msa_array, thermo_cfg, start_col, offset)
                for start_col in range(msa_array.shape[1] - 14)
            ],
        )
    mp_thermo_pass_rkmers = [x for x in rprimer_mp if x is not None]
    mp_thermo_pass_rkmers.sort(key=lambda rkmer: rkmer.start)

    return (mp_thermo_pass_fkmers, mp_thermo_pass_rkmers)
