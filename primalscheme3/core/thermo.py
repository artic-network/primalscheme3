from primer3 import calc_tm as p3_calc_tm, calc_hairpin as p3_calc_hairpin
from typing import Iterable
from itertools import groupby

from enum import Enum


class THERMORESULT(Enum):
    # THERMORESULT.value == 0 is a pass
    PASS = 0
    HIGH_GC = 1
    LOW_GC = 2
    HIGH_TM = 3
    LOW_TM = 4
    MAX_HOMOPOLY = 5


def calc_tm(kmer_seq, cfg: dict) -> float:
    """Return Tm for the kmer sequence."""
    return p3_calc_tm(
        kmer_seq,
        mv_conc=cfg["mv_conc"],
        dv_conc=cfg["dv_conc"],
        dntp_conc=cfg["dntp_conc"],
        dna_conc=cfg["dna_conc"],
    )


def calc_hairpin_tm(seq: str, config: dict) -> float:
    """
    Calculate the hairpin formation thermodynamics of a DNA sequence.
    Returns tm.
    """
    return p3_calc_hairpin(
        seq,
        mv_conc=config["mv_conc"],
        dv_conc=config["dv_conc"],
        dntp_conc=config["dntp_conc"],
        dna_conc=config["dna_conc"],
    ).tm


def calc_hairpin_struct(seq: str, config: dict) -> float:
    """
    Calculate the hairpin formation thermodynamics of a DNA sequence.
    Returns tm.
    """
    return p3_calc_hairpin(
        seq,
        mv_conc=config["mv_conc"],
        dv_conc=config["dv_conc"],
        dntp_conc=config["dntp_conc"],
        dna_conc=config["dna_conc"],
        output_structure=True,
    ).ascii_structure_lines


def forms_hairpin(seqs: Iterable[str], cfg: dict) -> bool:
    """
    Given a iterable of strings it will check the hairpin tm of all seqs
    If any form haripins it will return True
    If all clear it will return False
    """
    for seq in seqs:
        if calc_hairpin_tm(seq, cfg) > cfg["primer_hairpin_th_max"]:
            return True
    return False


def gc(kmer_seq: str) -> float:
    return round(100.0 * (kmer_seq.count("G") + kmer_seq.count("C")) / len(kmer_seq), 1)


def max_homo(kmer_seq) -> int:
    """Return max homopolymer length for the kmer sequence."""
    if not kmer_seq:
        return 0
    return max(sum(1 for _ in group) for _, group in groupby(kmer_seq))


def passes_thermo_checks(kmer_seq: str, cfg: dict) -> THERMORESULT:
    """Are all kmer thermo values below threshold?.

    Evaluation order.
    GC CHECK
    TM CHECK
    HOMOPOLY CHECK
    PASS

    Args:
        kmer_seq (str): The kmer sequence to be checked.
        cfg (dict): The configuration dictionary containing threshold values.

    Returns:
        THERMORESULT: The result of the thermo checks.
    """

    # Check for gc in range
    kmer_gc = gc(kmer_seq)
    if kmer_gc > cfg["primer_gc_max"]:
        return THERMORESULT.HIGH_GC
    elif kmer_gc < cfg["primer_gc_min"]:
        return THERMORESULT.LOW_GC

    # Check for tm in range
    kmer_tm = calc_tm(kmer_seq, cfg)
    if kmer_tm > cfg["primer_tm_max"]:
        return THERMORESULT.HIGH_TM
    elif kmer_tm < cfg["primer_tm_min"]:
        return THERMORESULT.LOW_TM

    # Check for maxhomopolymer
    if max_homo(kmer_seq) > cfg["primer_homopolymer_max"]:
        return THERMORESULT.MAX_HOMOPOLY

    return THERMORESULT.PASS


def thermo_check_kmers(kmers: Iterable[str], cfg: dict) -> THERMORESULT:
    """
    Will call passes_thermo_checks on each kmer sequence in the kmers list
    Will stop evaluating on first error

    Args:
        kmers (Iterable[str]): A list of kmer sequences to be evaluated.
        cfg (dict): A dictionary containing configuration settings.

    Returns:
        THERMORESULT: The result of the thermo checks. THERMORESULT.PASS if all kmers pass the checks, otherwise the first encountered error.

    """
    for kmer in kmers:
        result = passes_thermo_checks(kmer, cfg)
        if result == THERMORESULT.PASS:
            continue
        else:
            return result

    return THERMORESULT.PASS
