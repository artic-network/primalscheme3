# Core imports
import pathlib

from primalbedtools.bedfiles import BedLineParser
from primalschemers._core import (
    calc_at_offset_py,  # type: ignore
)

from primalscheme3.core.seq_functions import expand_all_ambs

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

    return f"score: {round(score, 2)}\n{seq1}\n{cigar}\n{seq2}\n"


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


def visualise_interactions(bedpath: pathlib.Path, threshold: float) -> None:
    """
    Calculate the interaction score between two sequences.
    If the score is less than the threshold
    :param seq1: Bedfile path
    :param threshold: The threshold for the interaction score
    :return: None
    """

    # Read in the bedfile
    _header, bedlines = BedLineParser.from_file(bedpath)

    # Split the bedfile into pools
    pools = [[] for _ in {bedline.pool for bedline in bedlines}]
    for bedline in bedlines:
        pools[bedline.ipool].append(bedline)

    tested = 0
    interactions = 0

    # Loop over the pools
    for pool in pools:
        # Loop over the bedlines in the pool
        for bedline1 in pool:
            for bedline2 in pool:
                tested += 1

                bedline1_seqs = list(expand_all_ambs([bedline1.sequence]))  # type: ignore
                bedline2_seqs = list(expand_all_ambs([bedline2.sequence]))  # type: ignore

                if bedline1_seqs is None or bedline2_seqs is None:
                    continue

                for bedline1_seq in bedline1_seqs:
                    for bedline2_seq in bedline2_seqs:
                        for line in interaction(bedline1_seq, bedline2_seq, threshold):
                            interactions += 1
                            print(bedline1.primername, bedline2.primername)
                            print(line)

    print(
        f"Tested {tested} possible combinations and found {interactions} interactions"
    )
