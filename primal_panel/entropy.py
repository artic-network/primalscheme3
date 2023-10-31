# Submodules imports
from primal_digest.config import ALL_DNA


from collections import Counter
import numpy as np
from math import log2, sqrt
from itertools import product


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
    num_invalids = counter.pop("N", None)

    return [v / (len(all_bases) - num_invalids) for _, v in counter.items()]


def entropy_score_array(msa: np.ndarray) -> list[int]:
    """
    Creates an list with the entropy at each index
    """
    score_array = [0] * msa.shape[1]
    # Iterate over colums
    for col in range(msa.shape[1]):
        proportions = calc_probs(list(msa[:, col]))
        score_array[col] = calc_entropy(proportions)
    return score_array
