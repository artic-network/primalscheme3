# Core imports
from primalscheme3.core.bedfiles import BedLine, read_bedfile
from primalscheme3.core.config import ALL_DNA

from primaldimer_py import (
    calc_at_offset_py,
)  # calc_at_offset_py(seq1: &str, seq2: &str, offset: i32) -> f64


import pathlib

MATCHES: dict[tuple, bool] = {
    ("A", "T"): True,
    ("T", "A"): True,
    ("G", "C"): True,
    ("C", "G"): True,
}

SIMPLE_DNA = {"A", "C", "G", "T"}


def create_cigar(seq1: str, seq2: str) -> str:
    cigar = []
    for seq1base, seq2base in zip(seq1, seq2):
        # Guard against non DNA bases
        if seq1base not in SIMPLE_DNA or seq2base not in SIMPLE_DNA:
            cigar.append(" ")
            continue
        # Check if the bases match
        if MATCHES.get((seq1base, seq2base), False):
            cigar.append("|")
        else:
            cigar.append(".")
    return "".join(cigar)


def create_str(seq1: str, seq2: str, offset: int, score: float) -> str:
    """
    Returns a string representing the interaction between two sequences.
    Example:
        offset = -1
        seq1:   AGCATCATGCTAGCT
                ..|.|..|.|.|.|..
        seq2:    AGCATCATGCTAGCT

    :param seq1: The first sequence  in 5'-3' orientation
    :param seq2: The second sequence in 5'-3' orientation
    :param offset: The offset of the sequences
    :return: A string representing the interaction
    """

    seq2 = seq2[::-1]  # Reverse the sequence

    # Add some direction annotation
    seq1 = f"5'-{seq1}-3' >"
    seq2 = f"3'-{seq2}-5'"

    # Put the sequences in the correct offset
    if offset < 0:
        seq2 = " " * abs(offset) + seq2
    elif offset > 0:
        seq1 = " " * offset + seq1

    # Pad the sequences so they are the same length
    max_length = max(len(seq1), len(seq2))
    seq1 = seq1.ljust(max_length)
    seq2 = seq2.ljust(max_length)
    # Create the cigar string
    cigar = create_cigar(seq1, seq2)

    return f"score: {round(score,2)}\n{seq1}\n{cigar}\n{seq2}\n"


def interaction(seq1: str, seq2: str, threshold: float) -> list[str]:
    """
    Find interactions between two sequences based on a given threshold.

    :param seq1 (str): The first sequence.
    :param seq2 (str): The second sequence.
    :param threshold (float): The threshold for the interaction score.
    :return  list[str]: A list of interactions between the two sequences.
    """
    interactions = []

    for offset in range(-(len(seq1) - 2), len(seq2) - len(seq1)):
        score = calc_at_offset_py(seq1, seq2, offset)
        if score <= threshold:
            interactions.append(create_str(seq1, seq2, offset, score))

    return interactions


def visulise_interactions(args) -> None:
    """
    Calculate the interaction score between two sequences.
    If the score is less than the threshold
    :param seq1: Bedfile path
    :param threshold: The threshold for the interaction score
    :return: None
    """
    bedpath = args.bedfile
    threshold = args.threshold

    # Read in the bedfile
    bedfile = read_bedfile(bedpath)

    # Split the bedfile into pools
    pools = [[] for _ in {bedline.pool for bedline in bedfile}]
    for bedline in bedfile:
        pools[bedline.pool].append(bedline)

    # Loop over the pools
    for pool in pools:
        # Loop over the bedlines in the pool
        for bedline1 in pool:
            for bedline2 in pool:
                for l in interaction(bedline1.sequence, bedline2.sequence, threshold):
                    print(bedline1.primername, bedline2.primername)
                    print(l)
