from primer3 import calc_tm as p3_calc_tm, calc_hairpin as p3_calc_hairpin
from typing import Iterable
from itertools import groupby


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


def passes_thermo_checks(kmer_seq, cfg) -> bool:
    """Are all kmer thermo values below threshold?"""
    return (
        (cfg["primer_gc_min"] <= gc(kmer_seq) <= cfg["primer_gc_max"])
        and (cfg["primer_tm_min"] <= calc_tm(kmer_seq, cfg) <= cfg["primer_tm_max"])
        and (max_homo(kmer_seq) <= cfg["primer_homopolymer_max"])
    )


def thermo_check_kmers(kmers: Iterable[str], cfg) -> bool:
    """
    Will call passes_thermo_checks on each kmer in the kmers list
    Will stop evaluating on false

    False means fail
    """
    for kmer in kmers:
        if not passes_thermo_checks(kmer, cfg):
            return False
    return True
