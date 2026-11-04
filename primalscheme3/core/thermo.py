from collections.abc import Iterable
from enum import Enum
from itertools import groupby

import primalschemers as ps

from primalscheme3.core.config import Config
from primalscheme3.core.seq_functions import reverse_complement

# Thermodynamics come from primalschemers: thal, as primer3 2.3.0, for
# hairpins and dimers, and the same nearest-neighbour model for an oligo's Tm
# and annealing, so these agree with the Rust digest.


class THERMO_RESULT(Enum):
    # THERMO_RESULT.value == 0 is a pass
    PASS = 0
    HIGH_GC = 1
    LOW_GC = 2
    HIGH_TM = 3
    LOW_TM = 4
    MAX_HOMOPOLY = 5
    HAIRPIN = 6
    TO_LONG = 7
    # Failed one of the Rust digest's checks (see the kmer's statuses)
    FAIL = 8
    HIGH_ANNEALING = 9
    LOW_ANNEALING = 10


ANNEALING_DIFF = 5


def calc_annealing(kmer_seq, mv_conc, dv_conc, dntp_conc, dna_conc, temp_c) -> float:
    """Percentage of the kmer bound to its complement at temp_c."""
    return ps.calc_annealing(
        kmer_seq, dna_conc, mv_conc, dv_conc, dntp_conc, 0.0, 0.6, 0.0, temp_c
    )


def calc_binding(kmer_seq: str, template_seq: str, config: Config):
    return calc_binding_at(
        kmer_seq,
        template_seq,
        mv_conc=config.mv_conc,
        dv_conc=config.dv_conc,
        dntp_conc=config.dntp_conc,
        dna_conc=config.dna_conc,
        temp_c=config.primer_annealing_tempc,
    )


def calc_binding_at(
    kmer_seq: str,
    template_seq: str,
    mv_conc: float,
    dv_conc: float,
    dntp_conc: float,
    dna_conc: float,
    temp_c: float,
):
    """calc_binding at explicit conditions rather than a Config's."""
    return ps.calc_heterodimer(
        kmer_seq,
        reverse_complement(template_seq),
        mv_conc=mv_conc,
        dv_conc=dv_conc,
        dntp_conc=dntp_conc,
        dna_conc=dna_conc,
        temp_c=temp_c,
    )


def calc_annealing_hetro(kmer_seq: str, template_seq: str, config: Config) -> float:
    """Percentage of the kmer bound to the template site at the annealing
    temperature (0 if no structure forms)."""
    return ps.calc_annealing_hetero(
        kmer_seq,
        reverse_complement(template_seq),
        config.primer_annealing_tempc,
        mv_conc=config.mv_conc,
        dv_conc=config.dv_conc,
        dntp_conc=config.dntp_conc,
        dna_conc=config.dna_conc,
    )


def calc_annealing_profile(
    kmer_seq, mv_conc, dv_conc, dntp_conc, dna_conc, min_temp=40, max_temp=80
) -> dict[int, float]:
    temps = list(range(min_temp, max_temp + 1))
    profile = ps.calc_annealing_profile(
        kmer_seq,
        [float(t) for t in temps],
        mv_conc=mv_conc,
        dv_conc=dv_conc,
        dntp_conc=dntp_conc,
        dna_conc=dna_conc,
    )
    return dict(zip(temps, profile, strict=True))


def calc_tm(
    kmer_seq, mv_conc, dv_conc, dntp_conc, dna_conc, formamide_conc: float = 0.8
) -> float:
    """Return Tm for the kmer sequence, as the digest checks it."""
    return ps.calc_tm(
        kmer_seq,
        mv_conc=mv_conc,
        dv_conc=dv_conc,
        dntp_conc=dntp_conc,
        dna_conc=dna_conc,
        formamide_conc=formamide_conc,
    )


def calc_hairpin_tm(seq: str, mv_conc, dv_conc, dntp_conc, dna_conc) -> float:
    """
    Calculate the hairpin formation thermodynamics of a DNA sequence.
    Returns tm.
    """
    return ps.calc_hairpin(
        seq,
        mv_conc=mv_conc,
        dv_conc=dv_conc,
        dntp_conc=dntp_conc,
        dna_conc=dna_conc,
    ).tm


def calc_hairpin_struct(
    seq: str, mv_conc, dv_conc, dntp_conc, dna_conc
) -> list[str] | None:
    """
    Calculate the hairpin formation thermodynamics of a DNA sequence.
    Returns the ASCII structure lines.
    """
    return ps.calc_hairpin(
        seq,
        mv_conc=mv_conc,
        dv_conc=dv_conc,
        dntp_conc=dntp_conc,
        dna_conc=dna_conc,
        output_structure=True,
    ).ascii_structure_lines


def forms_hairpin(seqs: list[str], config: Config) -> bool:
    """
    Given a iterable of strings it will check the hairpin tm of all seqs
    If any form hairpins it will return True
    If all clear it will return False
    """
    results = ps.calc_hairpin_batch(
        list(seqs),
        mv_conc=config.mv_conc,
        dv_conc=config.dv_conc,
        dntp_conc=config.dntp_conc,
        dna_conc=config.dna_conc,
    )
    return any(r.tm > config.primer_hairpin_th_max for r in results)


def gc(kmer_seq: str) -> float:
    return round(100.0 * (kmer_seq.count("G") + kmer_seq.count("C")) / len(kmer_seq), 1)


def max_homo(kmer_seq) -> int:
    """Return max homopolymer length for the kmer sequence."""
    if not kmer_seq:
        return 0
    return max(sum(1 for _ in group) for _, group in groupby(kmer_seq))


def thermo_check(kmer_seq: str, config: Config) -> THERMO_RESULT:
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
        THERMO_RESULT: The result of the thermo checks.
    """
    # Check for gc in range
    kmer_gc = gc(kmer_seq)
    if kmer_gc > config.primer_gc_max:
        return THERMO_RESULT.HIGH_GC
    elif kmer_gc < config.primer_gc_min:
        return THERMO_RESULT.LOW_GC

    # Check length
    if len(kmer_seq) > config.primer_size_max:
        return THERMO_RESULT.TO_LONG

    if config.use_annealing and config.primer_annealing_prop is not None:
        kmer_annealing = calc_annealing(
            kmer_seq,
            mv_conc=config.mv_conc,
            dv_conc=config.dv_conc,
            dna_conc=config.dna_conc,
            dntp_conc=config.dntp_conc,
            temp_c=config.primer_annealing_tempc,
        )
        if kmer_annealing > config.primer_annealing_prop + ANNEALING_DIFF:
            return THERMO_RESULT.HIGH_ANNEALING
        elif kmer_annealing < config.primer_annealing_prop - ANNEALING_DIFF:
            return THERMO_RESULT.LOW_ANNEALING
    else:
        # Check for tm in range
        kmer_tm = calc_tm(
            kmer_seq,
            mv_conc=config.mv_conc,
            dv_conc=config.dv_conc,
            dna_conc=config.dna_conc,
            dntp_conc=config.dntp_conc,
            formamide_conc=config.formamide_conc,
        )
        # As the Rust digest: fail at or above the maximum
        if kmer_tm >= config.primer_tm_max:
            return THERMO_RESULT.HIGH_TM
        elif kmer_tm < config.primer_tm_min:
            return THERMO_RESULT.LOW_TM

    # Check for maxhomopolymer
    if max_homo(kmer_seq) > config.primer_homopolymer_max:
        return THERMO_RESULT.MAX_HOMOPOLY

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
        return THERMO_RESULT.HAIRPIN

    return THERMO_RESULT.PASS


def thermo_check_all_kmers(
    kmers: Iterable[str], config: Config
) -> dict[str, THERMO_RESULT]:
    return {kmer: thermo_check(kmer, config) for kmer in kmers}


def thermo_check_kmers(kmers: Iterable[str], config: Config) -> THERMO_RESULT:
    """
    Will call thermo_check on each kmer sequence in the kmers list
    Will stop evaluating on first error

    Args:
        kmers (Iterable[str]): A list of kmer sequences to be evaluated.
        cfg (dict): A dictionary containing configuration settings.

    Returns:
        THERMO_RESULT: The result of the thermo checks. THERMO_RESULT.PASS if all kmers pass the checks, otherwise the first encountered error.

    """
    for kmer in kmers:
        result = thermo_check(kmer, config)
        if result == THERMO_RESULT.PASS:
            continue
        else:
            return result

    return THERMO_RESULT.PASS
