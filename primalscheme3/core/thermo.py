from enum import Enum
from itertools import groupby
from typing import Iterable

from primer3 import calc_hairpin as p3_calc_hairpin
from primer3 import calc_tm as p3_calc_tm

from primalscheme3.core.config import Config


class THERMORESULT(Enum):
    # THERMORESULT.value == 0 is a pass
    PASS = 0
    HIGH_GC = 1
    LOW_GC = 2
    HIGH_TM = 3
    LOW_TM = 4
    MAX_HOMOPOLY = 5
    HAIRPIN = 6


def calc_tm(kmer_seq, mv_conc, dv_conc, dntp_conc, dna_conc) -> float:
    """Return Tm for the kmer sequence."""
    return p3_calc_tm(
        kmer_seq,
        mv_conc=mv_conc,
        dv_conc=dv_conc,
        dntp_conc=dntp_conc,
        dna_conc=dna_conc,
    )


def calc_hairpin_tm(seq: str, mv_conc, dv_conc, dntp_conc, dna_conc) -> float:
    """
    Calculate the hairpin formation thermodynamics of a DNA sequence.
    Returns tm.
    """
    return p3_calc_hairpin(
        seq,
        mv_conc=mv_conc,
        dv_conc=dv_conc,
        dntp_conc=dntp_conc,
        dna_conc=dna_conc,
    ).tm


def calc_hairpin_struct(seq: str, mv_conc, dv_conc, dntp_conc, dna_conc) -> float:
    """
    Calculate the hairpin formation thermodynamics of a DNA sequence.
    Returns tm.
    """
    return p3_calc_hairpin(
        seq,
        mv_conc=mv_conc,
        dv_conc=dv_conc,
        dntp_conc=dntp_conc,
        dna_conc=dna_conc,
        output_structure=True,
    ).ascii_structure_lines


def forms_hairpin(seqs: Iterable[str], config: Config) -> bool:
    """
    Given a iterable of strings it will check the hairpin tm of all seqs
    If any form hairpins it will return True
    If all clear it will return False
    """
    for seq in seqs:
        if (
            calc_hairpin_tm(
                seq,
                mv_conc=config.mv_conc,
                dv_conc=config.dv_conc,
                dntp_conc=config.dntp_conc,
                dna_conc=config.dna_conc,
            )
            > config.primer_hairpin_th_max
        ):
            return True
    return False


def gc(kmer_seq: str) -> float:
    return round(100.0 * (kmer_seq.count("G") + kmer_seq.count("C")) / len(kmer_seq), 1)


def max_homo(kmer_seq) -> int:
    """Return max homopolymer length for the kmer sequence."""
    if not kmer_seq:
        return 0
    return max(sum(1 for _ in group) for _, group in groupby(kmer_seq))


def thermo_check(kmer_seq: str, config: Config) -> THERMORESULT:
    """Are all kmer thermo values below threshold?.

    Evaluation order.
    GC CHECK
    TM CHECK
    HOMOPOLY CHECK
    HAIRPIN
    PASS

    Args:
        kmer_seq (str): The kmer sequence to be checked.
        cfg (dict): The configuration dictionary containing threshold values.

    Returns:
        THERMORESULT: The result of the thermo checks.
    """
    # Check for gc in range
    kmer_gc = gc(kmer_seq)
    if kmer_gc > config.primer_gc_max:
        return THERMORESULT.HIGH_GC
    elif kmer_gc < config.primer_gc_min:
        return THERMORESULT.LOW_GC

    # Check for tm in range
    kmer_tm = calc_tm(
        kmer_seq,
        mv_conc=config.mv_conc,
        dv_conc=config.dv_conc,
        dna_conc=config.dna_conc,
        dntp_conc=config.dntp_conc,
    )
    if kmer_tm > config.primer_tm_max:
        return THERMORESULT.HIGH_TM
    elif kmer_tm < config.primer_tm_min:
        return THERMORESULT.LOW_TM

    # Check for maxhomopolymer
    if max_homo(kmer_seq) > config.primer_homopolymer_max:
        return THERMORESULT.MAX_HOMOPOLY

    # Check for hairpin
    if (
        calc_hairpin_tm(
            kmer_seq,
            mv_conc=config.mv_conc,
            dv_conc=config.dv_conc,
            dntp_conc=config.dntp_conc,
            dna_conc=config.dna_conc,
        )
        > config.primer_hairpin_th_max
    ):
        return THERMORESULT.HAIRPIN

    return THERMORESULT.PASS


def thermo_check_all_kmers(
    kmers: Iterable[str], config: Config
) -> dict[str, THERMORESULT]:
    return {kmer: thermo_check(kmer, config) for kmer in kmers}


def thermo_check_kmers(kmers: Iterable[str], config: Config) -> THERMORESULT:
    """
    Will call thermo_check on each kmer sequence in the kmers list
    Will stop evaluating on first error

    Args:
        kmers (Iterable[str]): A list of kmer sequences to be evaluated.
        cfg (dict): A dictionary containing configuration settings.

    Returns:
        THERMORESULT: The result of the thermo checks. THERMORESULT.PASS if all kmers pass the checks, otherwise the first encountered error.

    """
    for kmer in kmers:
        result = thermo_check(kmer, config)
        if result == THERMORESULT.PASS:
            continue
        else:
            return result

    return THERMORESULT.PASS
