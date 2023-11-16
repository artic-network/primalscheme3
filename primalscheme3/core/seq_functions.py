from itertools import product
from typing import Iterable
import numpy as np
from collections import Counter
from math import log2

# Module Imports
from primalscheme3.core.config import (
    ALL_DNA,
    ALL_BASES,
    AMB_BASES,
    AMBIGUOUS_DNA_COMPLEMENT,
)


def reverse_complement(kmer_seq: str) -> str:
    rev_seq = kmer_seq[::-1]
    return "".join(AMBIGUOUS_DNA_COMPLEMENT[base.upper()] for base in rev_seq)


def get_most_common_base(array: np.ndarray, col_index: int) -> str:
    col_slice = array[:, col_index]
    counter = Counter(col_slice)
    del counter[""]
    return counter.most_common(1)[0][0]


def remove_end_insertion(msa_array: np.ndarray) -> np.ndarray:
    """
    Removes leading and trailing "-" from an msa
    """
    tmp_array = msa_array
    ncols = tmp_array.shape[1]
    for row_index in range(0, tmp_array.shape[0]):
        # Solves the 5' end
        for col_index in range(0, ncols):
            if tmp_array[row_index, col_index] == "-":
                tmp_array[row_index, col_index] = ""
            else:
                break
        for rev_col_index in range(ncols - 1, 0, -1):
            if tmp_array[row_index, rev_col_index] == "-":
                tmp_array[row_index, rev_col_index] = ""
            else:
                break
    return tmp_array


def expand_ambs(seqs: Iterable[str]) -> set[str] | None:
    """
    Takes a list / set of strings and returns a set with all ambs expanded
    Return None on invalid bases (Including N)
    """
    returned_seq = set()

    for seq in seqs:
        bases = {base for base in seq}

        # if invalid bases are in sequence return None
        if not bases.issubset(ALL_BASES):
            return None

        # If there is any amb_bases in the seq
        if bases & AMB_BASES:
            expanded_seqs = set(map("".join, product(*map(ALL_DNA.get, seq))))  # type: ignore
            for exp_seq in expanded_seqs:
                returned_seq.add(exp_seq)
        else:
            returned_seq.add(seq)
    return returned_seq


def extend_ambiguous_base(base: str) -> list[str]:
    """Return list of all possible sequences given an ambiguous DNA input"""
    return [x for x in ALL_DNA.get(base, "N")]


def calc_entropy(probs: list[float]) -> float:
    return -sum([p * log2(p) for p in probs if p])


def calc_probs(bases: list[str]) -> list[float]:
    """Calculate the probability/proportion of each base in each column"""
    all_bases = [
        y for sublist in (extend_ambiguous_base(x) for x in bases) for y in sublist
    ]
    counter = Counter(all_bases)
    num_invalids = counter.pop("N", 0)

    return [v / (len(all_bases) - num_invalids) for _, v in counter.items()]


def entropy_score_array(msa: np.ndarray) -> list[float]:
    """
    Creates an list with the entropy at each index
    """
    score_array: list[float] = [0] * msa.shape[1]
    # Iterate over colums
    for col in range(msa.shape[1]):
        proportions = calc_probs(list(msa[:, col]))
        score_array[col] = calc_entropy(proportions)
    return score_array
