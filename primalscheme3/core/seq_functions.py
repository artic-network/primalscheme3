from collections import Counter
from itertools import product
from math import log2

import numpy as np

# Module Imports
from primalscheme3.core.config import (
    ALL_BASES,
    ALL_BASES_WITH_N,
    ALL_DNA,
    ALL_DNA_WITH_N,
    AMB_BASES,
    AMBIGUOUS_DNA_COMPLEMENT,
    SIMPLE_BASES,
)


def reverse_complement(kmer_seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    rev_seq = kmer_seq[::-1]
    return complement(rev_seq)


def complement(seq: str) -> str:
    """Return the complement of a DNA sequence."""
    return "".join(AMBIGUOUS_DNA_COMPLEMENT[base.upper()] for base in seq)


def get_most_common_base(array: np.ndarray, col_index: int) -> str:
    values, counts = np.unique(array[:, col_index], return_counts=True)
    counter = dict(zip(values, counts))
    counter.pop("", None)
    return max(counter, key=counter.get)  # type: ignore


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


def expand_ambs(seqs: list[str] | set[str]) -> set[str] | None:
    """
    Takes a list / set of strings and returns a set with all ambs expanded
    Return None on invalid bases (Including N)
    """
    returned_seq = set()

    for seq in seqs:
        bases = {*seq}

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


def expand_all_ambs(seqs: list[str] | set[str]) -> set[str] | None:
    """
    Expands all ambiguous bases in a sequence. Including N
    """
    returned_seq = set()

    for seq in seqs:
        bases = {*seq}

        # if invalid bases are in sequence return None
        if not bases.issubset(ALL_BASES_WITH_N):
            return None

        # If there is any amb_bases in the seq

        expanded_seqs = set(map("".join, product(*map(ALL_DNA_WITH_N.get, seq))))  # type: ignore
        for exp_seq in expanded_seqs:
            returned_seq.add(exp_seq)

    return returned_seq


def extend_ambiguous_base(base: str) -> list[str]:
    """Return list of all possible sequences given an ambiguous DNA input"""
    return [*ALL_DNA.get(base, "N")]


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
        value, counts = np.unique(msa[:, col], return_counts=True)
        count_dict = dict(zip(value, counts))

        # Remove non DNA bases
        count_dict = {k: v for k, v in count_dict.items() if k in ALL_BASES}

        parsed_counts = {base: 0 for base in "ACGT"}
        # Expand ambiguous bases
        for base in count_dict:
            if base in AMB_BASES:
                amb_count = count_dict.get(base, 0)  # Remove the ambiguous base
                for expanded_base in extend_ambiguous_base(base):
                    if expanded_base != "N":
                        parsed_counts[expanded_base] += amb_count
            else:
                parsed_counts[base] += count_dict.get(base, 0)
        # Remove Invalid bases
        parsed_counts = {k: v for k, v in count_dict.items() if k in SIMPLE_BASES}

        # Calculate the proportions (probabilities) of each base
        proportions = [v / sum(parsed_counts.values()) for v in parsed_counts.values()]

        score_array[col] = calc_entropy(proportions)
    return score_array
