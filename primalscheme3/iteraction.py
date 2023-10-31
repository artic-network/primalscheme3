import parasail
from primal_digest.thermo import *
import re

CIGAR_REGEX = re.compile(r"(?P<len>\d+)(?P<op>\D+)")


def cig_check(regex_matches: list[tuple[str, str]], query: str, cfg: dict) -> bool:
    if regex_matches[0][1] == "=":  # Nothing to extend
        return False
    overlap = ""
    indels = 0
    matches = 0
    mismatches = 0
    for group in regex_matches:  # Iterate over the sequences
        if group[1] == "I" or group[1] == "D":  # Skip I or D ops
            indels += int(group[0])
            continue
        if group[1] == "X":  # Mismatch
            if not overlap:  # Must be 3' anchored
                return False
            mismatches += int(group[0])
        if group[1] == "=":  # Match
            i = sum((indels, matches, mismatches))
            overlap += query[i : i + int(group[0])]
            matches += int(group[0])
        identity = matches / (matches + mismatches)
        ol_tm = calc_tm(overlap, cfg)
        return ol_tm > cfg["dimer_max_tm"] and identity > cfg["dimer_min_identity"]
    return False  # no interaction


def parasail_align(seq1: str, seq2: str) -> parasail.Traceback:
    OPEN = 10
    EXTEND = 5
    # Semi-Global, do not penalize gaps at beginning and end of both sequences
    trace = parasail.sg_trace_scan(seq1, seq2, OPEN, EXTEND, parasail.dnafull)
    return trace


def seqs_may_interact(seq1: str, seq2: str, cfg, verbose: bool = False) -> bool:
    trace = parasail_align(seq1, seq2)
    traceback = trace.get_traceback()
    matches = re.findall(CIGAR_REGEX, trace.cigar.decode.decode())
    interaction_predicted = cig_check(matches, traceback.query, cfg) or cig_check(
        matches[::-1], traceback.ref[::-1], cfg
    )
    return interaction_predicted


def interact_predict(seqs: Iterable[str], cfg: dict) -> bool:
    """
    Given a list of strings, it look for interactions between all strings
    Will only check unique combinations of strings
    True means there is an interaciton between the strings
    False means there is no inter
    """
    set_of_interactions: set[set[str]] = set()
    for sequence1 in seqs:
        for sequence2 in seqs:
            sub = frozenset((sequence1, sequence2))
            if len(sub) > 1:
                set_of_interactions.add(sub)
    for seq_interaction in set_of_interactions:
        tmp = list(seq_interaction)
        if seqs_may_interact(seq1=tmp[0], seq2=tmp[1], cfg=cfg):
            return True
    return False


def all_inter_checker(seq1s: Iterable[str], seq2s: Iterable[str], cfg) -> bool:
    """
    Give two iterables for strings, it will check all seq1s vs seqs2
    True means interaction
    """
    inter_results = set()
    for seq1 in seq1s:
        for seq2 in seq2s:
            inter_results.add(seqs_may_interact(seq1=seq1, seq2=seq2, cfg=cfg))
    return inter_results.__contains__(True)
