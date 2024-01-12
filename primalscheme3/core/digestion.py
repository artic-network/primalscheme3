# Modules
from primalscheme3.core.thermo import calc_tm, thermo_check_kmers, forms_hairpin
from primalscheme3.core.seq_functions import expand_ambs, get_most_common_base
from primalscheme3.core.classes import RKmer, FKmer, PrimerPair
from primalscheme3.core.get_window import get_r_window_FAST2
from primalscheme3.core.errors import (
    WalksOut,
    GapOnSetBase,
    ContainsInvalidBase,
    CustomErrors,
    ERROR_SET,
    CustomRecursionError,
    WalksTooFar,
)

# Submodules
from primaldimer_py import do_pools_interact_py  # type: ignore

# Externals
from tqdm import tqdm
import numpy as np
from multiprocessing import Pool
import itertools
import networkx as nx
from collections import Counter
from typing import Callable, Set, Union
from enum import Enum


class DIGESTION_ERROR(Enum):
    """
    Enum for the different types of errors that can occur during digestion
    """

    WALKS_OUT = "WalksOut"
    CONTAINS_INVALID_BASE = "ContainsInvalidBase"
    CUSTOM_RECURSION_ERROR = "CustomRecursionError"
    CUSTOM_ERRORS = "CustomErrors"
    GAP_ON_SET_BASE = "GapOnSetBase"
    THERMO_FAIL = "ThermoFail"
    HAIRPIN_FAIL = "HairpinFail"
    DIMER_FAIL = "DimerFail"  # Interaction within the kmer
    WALK_TO_FAR = "WalkToFar"  # When indels causes the walk to go to far
    AMB_FAIL = "AmbFail"  # Generic error for when the error is unknown


def parse_error(results: set) -> DIGESTION_ERROR:
    """
    Parses the error set for the error that occured
    As only one error is returned, there is an arbitrary heirarchy of errors
    - CONTAINS_INVALID_BASE > GAP_ON_SET_BASE > WALKS_OUT > CUSTOM_RECURSION_ERROR > WALK_TO_FAR > CUSTOM_ERRORS
    """
    if ContainsInvalidBase() in results:
        return DIGESTION_ERROR.CONTAINS_INVALID_BASE
    elif GapOnSetBase() in results:
        return DIGESTION_ERROR.GAP_ON_SET_BASE
    elif WalksOut() in results:
        return DIGESTION_ERROR.WALKS_OUT
    elif CustomRecursionError() in results:
        return DIGESTION_ERROR.CUSTOM_RECURSION_ERROR
    elif CustomErrors() in results:
        return DIGESTION_ERROR.CUSTOM_ERRORS
    elif WalksTooFar() in results:
        return DIGESTION_ERROR.WALK_TO_FAR
    else:  # Return a generic error
        return DIGESTION_ERROR.AMB_FAIL


def _mp_pp_inter_free(data: tuple[PrimerPair, dict]) -> bool:
    """
    True means interaction
    """
    pp = data[0]
    cfg = data[1]
    # If they do interact return None
    return do_pools_interact_py(
        [*pp.fprimer.seqs], [*pp.rprimer.seqs], cfg["dimerscore"]
    )


def generate_valid_primerpairs(
    fkmers: list[FKmer],
    rkmers: list[RKmer],
    cfg: dict,
    msa_index: int,
    disable_progress_bar: bool = False,
) -> list[PrimerPair]:
    """Generates valid primer pairs for a given set of forward and reverse kmers.

    Args:
        fkmers: A list of forward kmers.
        rkmers: A list of reverse kmers.
        cfg: A dictionary containing configuration parameters.
        msa_index: An integer representing the index of the multiple sequence alignment.
        disable_progress_bar: A boolean indicating whether to disable the progress bar.

    Returns:
        A list of valid primer pairs.
    """
    ## Generate all primerpairs without checking
    non_checked_pp = []
    for fkmer in fkmers:
        fkmer_start = min(fkmer.starts())
        # Get all rkmers that would make a valid amplicon
        pos_rkmer = get_r_window_FAST2(
            kmers=rkmers,
            start=fkmer_start + cfg["amplicon_size_min"],
            end=fkmer_start + cfg["amplicon_size_max"],
        )
        for rkmer in pos_rkmer:
            non_checked_pp.append(PrimerPair(fkmer, rkmer, msa_index))

    ## Interaction check all the primerpairs
    iter_free_primer_pairs = []
    if cfg["n_cores"] == 1:
        for pp in tqdm(
            non_checked_pp,
            desc="Generating PrimerPairs",
            disable=disable_progress_bar,
        ):
            if not pp.inter_free(cfg):
                iter_free_primer_pairs.append(pp)
    else:
        with Pool(cfg["n_cores"]) as p:
            mp_pp_bool = tqdm(
                p.imap_unordered(
                    _mp_pp_inter_free, ((pp, cfg) for pp in non_checked_pp)
                ),
                total=len(non_checked_pp),
                desc="Generating PrimerPairs",
                disable=disable_progress_bar,
            )
            for bool, pp in zip(mp_pp_bool, non_checked_pp):
                if not bool:
                    iter_free_primer_pairs.append(pp)

    iter_free_primer_pairs.sort(key=lambda pp: (pp.fprimer.end, -pp.rprimer.start))
    return iter_free_primer_pairs


def walk_right(
    array: np.ndarray,
    col_index_right: int,
    col_index_left: int,
    row_index: int,
    seq_str: str,
    cfg: dict,
) -> Union[set[str], Exception]:
    """
    Walks to the right of the array and returns a set of valid sequences.

    Args:
        array: A numpy array of DNA sequences.
        col_index_right: The current column index to the right.
        col_index_left: The current column index to the left.
        row_index: The current row index.
        seq_str: The current sequence string.
        cfg: A dictionary of configuration parameters.

    Returns:
        A set of valid DNA sequences or an exception if an error occurs.

    Raises:
        WalksOut: If the function walks out of the array size.
        ContainsInvalidBase: If the sequence contains an invalid base.
    """
    # Guard for correct tm
    if calc_tm(seq_str, cfg) >= cfg["primer_tm_min"]:
        return {seq_str}

    # Guard prevents walking out of array size
    if col_index_right >= array.shape[1] - 1 or col_index_left >= array.shape[1] - 1:
        raise WalksOut()

    # Guard for walking too far
    if col_index_right - col_index_left >= cfg["primer_max_walk"]:
        raise WalksTooFar()

    new_base = array[row_index, col_index_right]

    # Fix incomplete ends
    if new_base == "":
        new_base = get_most_common_base(array, col_index_right + 1)
    new_string = (seq_str + new_base).replace("-", "")

    # Prevent Ns from being added
    if "N" in new_string:
        raise ContainsInvalidBase()

    # Guard for invalid bases in the sequence
    exp_new_string: set[str] | None = expand_ambs([new_string])
    if exp_new_string is None:
        raise ContainsInvalidBase()

    passing_str = []
    for exp_str in exp_new_string:
        results = wrap_walk(
            walk_right,
            array,
            col_index_right + 1,
            col_index_left,
            row_index,
            exp_str,
            cfg,
        )
        passing_str.extend(results)

    return passing_str  # type: ignore


def walk_left(
    array: np.ndarray,
    col_index_right: int,
    col_index_left: int,
    row_index: int,
    seq_str: str,
    cfg: dict,
) -> set[str] | Exception:
    """
    Recursively walks left from a given starting position in a 2D numpy array of DNA bases,
    constructing a set of valid DNA sequences that meet certain criteria.

    Args:
        array: A 2D numpy array of DNA bases.
        col_index_right: The rightmost column index of the region of interest.
        col_index_left: The current leftmost column index of the region of interest.
        row_index: The current row index of the region of interest.
        seq_str: The current DNA sequence being constructed.
        cfg: A dictionary of configuration parameters.

    Returns:
        A set of valid DNA sequences that meet the criteria specified in the function body.

    Raises:
        WalksOut: If the function attempts to walk out of the array.
        ContainsInvalidBase: If the constructed sequence contains an invalid DNA base.
    """

    # Guard prevents walking out of array size
    if col_index_left <= 0 or col_index_right <= 0:
        raise WalksOut()

    # Guard for correct tm
    if calc_tm(seq_str, cfg) >= cfg["primer_tm_min"]:
        return {seq_str}

    # Guard for walking too far
    if col_index_right - col_index_left >= cfg["primer_max_walk"]:
        raise WalksTooFar()

    new_base = array[row_index, col_index_left - 1]

    # Ensure it can repair truncated regions
    if new_base == "":
        new_base = get_most_common_base(array, col_index_left - 1)
    new_string = (new_base + seq_str).replace("-", "")

    # Guard prevents seqs with an N
    if "N" in new_string:
        raise ContainsInvalidBase()

    # If invalid bases return None
    exp_new_string: set[str] | None = expand_ambs([new_string])
    if exp_new_string is None:
        raise ContainsInvalidBase()

    passing_str = []
    for exp_str in exp_new_string:
        results = wrap_walk(
            walk_left,
            array=array,
            col_index_right=col_index_right,
            col_index_left=col_index_left - 1,
            row_index=row_index,
            seq_str=exp_str,
            cfg=cfg,
        )
        passing_str.extend(results)

    return passing_str  # type: ignore


def wrap_walk(
    walkfunction: Callable,
    array: np.ndarray,
    col_index_right: int,
    col_index_left: int,
    row_index: int,
    seq_str: str,
    cfg: dict,
) -> list[str | Exception]:
    return_list = []
    try:
        seqs = walkfunction(
            array=array,
            col_index_right=col_index_right,
            col_index_left=col_index_left,
            row_index=row_index,
            seq_str=seq_str,
            cfg=cfg,
        )
    except CustomErrors as e:
        return_list.append(e)
    except Exception as e:
        raise e
    else:
        return_list.extend(seqs)

    return return_list


def mp_r_digest(
    data: tuple[np.ndarray, dict, int, float]
) -> RKmer | tuple[int, DIGESTION_ERROR]:
    """
    This will try and create a RKmer started at the given index
    :data: A tuple of (align_array, cfg, start_col, min_freq)
    :return: A RKmer object or a tuple of (start_col, error)
    """
    align_array: np.ndarray = data[0]
    cfg: dict = data[1]
    start_col: int = data[2]
    min_freq: float = data[3]

    # Check for gap frequency on first base
    first_base_counter = Counter(align_array[:, start_col])
    del first_base_counter[""]
    num_seqs = sum(first_base_counter.values())
    first_base_freq = {k: v / num_seqs for k, v in first_base_counter.items()}

    # If the freq of gap is above minfreq
    if first_base_freq.get("-", 0) > min_freq:
        return (start_col, DIGESTION_ERROR.GAP_ON_SET_BASE)

    # If the initial slice is outside the range of the array
    if start_col + cfg["primer_size_min"] >= align_array.shape[1]:
        return (start_col, DIGESTION_ERROR.WALKS_OUT)

    # Create a counter
    total_col_seqs = Counter()
    for row_index in range(0, align_array.shape[0]):
        # Check if this row starts on a gap, and if so update the counter and skip
        if align_array[row_index, start_col] == "-":
            total_col_seqs.update([GapOnSetBase()])
            continue

        start_array = align_array[
            row_index, start_col : start_col + cfg["primer_size_min"]
        ]
        start_seq = "".join(start_array).replace("-", "")

        if not start_seq:  # If the start seq is empty go to the next row
            continue

        # Get all sequences
        results = wrap_walk(
            walk_right,
            array=align_array,
            col_index_right=start_col + cfg["primer_size_min"],
            col_index_left=start_col,
            row_index=row_index,
            seq_str=start_seq,
            cfg=cfg,
        )
        # If all mutations matter, return on any Error
        if min_freq == 0 and set(results) & ERROR_SET:
            return (start_col, parse_error(set(results)))

        # Add the results to the Counter
        total_col_seqs.update(results)

    total_values = sum(total_col_seqs.values())
    # Filter out values below the threshold freq
    wanted_seqs = {k for k, v in total_col_seqs.items() if v / total_values > min_freq}

    # Guard: If wanted_seqs contains errors return None
    if ERROR_SET & wanted_seqs:
        return (start_col, parse_error(wanted_seqs))

    # Create the Kmer
    tmp_kmer = RKmer(start=start_col, seqs=wanted_seqs)
    tmp_kmer.seqs = tmp_kmer.reverse_complement()
    # Downsample the seqs if asked
    if cfg["reducekmers"]:
        tmp_kmer.seqs = reduce_kmers(
            seqs=tmp_kmer.seqs,
            max_edit_dist=cfg["editdist_max"],
            end_3p=cfg["editdist_end3p"],
        )
    # Thermo check the kmers
    if not thermo_check_kmers(tmp_kmer.seqs, cfg):
        return (start_col, DIGESTION_ERROR.THERMO_FAIL)
    # Check for hairpins
    if forms_hairpin(tmp_kmer.seqs, cfg=cfg):
        return (start_col, DIGESTION_ERROR.HAIRPIN_FAIL)
    # Check for dimer
    if do_pools_interact_py([*tmp_kmer.seqs], [*tmp_kmer.seqs], cfg["dimerscore"]):
        return (start_col, DIGESTION_ERROR.DIMER_FAIL)
    # All checks pass return the kmer
    return tmp_kmer


def mp_f_digest(
    data: tuple[np.ndarray, dict, int, float]
) -> FKmer | tuple[int, DIGESTION_ERROR]:
    """
    This will try and create a FKmer ended at the given index
    :data: A tuple of (align_array, cfg, end_col, min_freq)
    :return: A FKmer object or a tuple of (end_col, error)
    """
    align_array: np.ndarray = data[0]
    cfg: dict = data[1]
    end_col: int = data[2]
    min_freq: float = data[3]

    # Check for gap frequency on first base
    first_base_counter = Counter(align_array[:, end_col])
    del first_base_counter[""]
    num_seqs = sum(first_base_counter.values())
    first_base_freq = {k: v / num_seqs for k, v in first_base_counter.items()}

    # If the freq of gap is above minfreq
    if first_base_freq.get("-", 0) > min_freq:
        return (end_col, DIGESTION_ERROR.GAP_ON_SET_BASE)

    # If the initial slice is outside the range of the array
    if end_col - cfg["primer_size_min"] < 0:
        return (end_col, DIGESTION_ERROR.WALKS_OUT)

    total_col_seqs = Counter()
    for row_index in range(0, align_array.shape[0]):
        # Check if this row starts on a gap, and if so update the counter and skip
        if align_array[row_index, end_col] == "-":
            total_col_seqs.update([GapOnSetBase()])
            # Skip to next row
            continue

        start_seq = "".join(
            align_array[row_index, end_col - cfg["primer_size_min"] : end_col]
        ).replace("-", "")

        if not start_seq:  # If the start seq is empty go to the next row
            continue

        results = wrap_walk(
            walk_left,
            array=align_array,
            col_index_right=end_col,
            col_index_left=end_col - cfg["primer_size_min"],
            row_index=row_index,
            seq_str=start_seq,
            cfg=cfg,
        )
        if min_freq == 0 and set(results) & ERROR_SET:
            return (end_col, parse_error(set(results)))

        total_col_seqs.update(results)

    total_values = sum(total_col_seqs.values())
    # Filter out values below the threshold freq
    wanted_seqs = {k for k, v in total_col_seqs.items() if v / total_values > min_freq}

    # Guard: If wanted_seqs contains errors return None
    if ERROR_SET & wanted_seqs:
        return (end_col, parse_error(wanted_seqs))

    # DownSample the seqs if asked
    if cfg["reducekmers"]:
        wanted_seqs = reduce_kmers(
            seqs=wanted_seqs,
            max_edit_dist=cfg["editdist_max"],
            end_3p=cfg["editdist_end3p"],
        )
    # Thermo check the kmers
    if not thermo_check_kmers(wanted_seqs, cfg):
        return (end_col, DIGESTION_ERROR.THERMO_FAIL)
    if forms_hairpin(wanted_seqs, cfg=cfg):
        return (end_col, DIGESTION_ERROR.HAIRPIN_FAIL)

    if do_pools_interact_py(list(wanted_seqs), list(wanted_seqs), cfg["dimerscore"]):
        return (end_col, DIGESTION_ERROR.DIMER_FAIL)

    return FKmer(end=end_col, seqs=wanted_seqs)


def hamming_dist(s1, s2) -> int:
    """
    Return the number of subsitutions, starting from the 3p end
    """
    return sum((x != y for x, y in zip(s1[::-1], s2[::-1])))


def reduce_kmers(seqs: set[str], max_edit_dist: int = 1, end_3p: int = 6) -> set[str]:
    """
    Reduces a set of DNA sequences by clustering them based on their 3' end, and then minimizing the edit distance between
    all tails within the same 3' cluster. The resulting set of sequences will have at most `max_edit_dist` differences
    between any two sequences, and will all have a common 3' end of length `end_3p`.

    Args:
        seqs: A set of DNA sequences to be reduced.
        max_edit_dist: The maximum edit distance allowed between any two sequences in the same 3' cluster. Defaults to 1.
        end_3p: The length of the 3' end to use for clustering. Defaults to 6.

    Returns:
        A set of reduced DNA sequences, where each sequence has a common 3' end of length `end_3p`, and at most
        `max_edit_dist` differences between any two sequences.
    """
    ## Cluster sequences via the 3p end
    p3_end_dict: dict[str, set[str]] = {}
    for sequence in seqs:
        p3_end = sequence[-end_3p:]
        p5_tail = sequence[:-end_3p]
        if p3_end in p3_end_dict:
            p3_end_dict[p3_end].add(p5_tail)
        else:
            p3_end_dict[p3_end] = {p5_tail}

    ## Minimise edit distance between all tails within the same p3 cluster
    for p3_end, p5_tails in p3_end_dict.items():
        # If only one sequence skip
        if len(p5_tails) <= 1:
            continue

        # Create a linkage graph
        G = nx.Graph()
        G.add_nodes_from(p5_tails)
        for s1, s2 in itertools.combinations(p5_tails, 2):
            if hamming_dist(s1, s2) <= max_edit_dist:
                # Add edges if the distance is <= hamming dist max
                G.add_edge(s1, s2)

        # Find the most connected sequence
        sorted_sequences = sorted(
            p5_tails, key=lambda seq: (len(list(G.neighbors(seq))), seq), reverse=True
        )

        # Seqs which are included in the scheme
        included_seqs = set()
        # Seqs which have a closely related sequence included
        accounted_seqs = set()

        for sequence in sorted_sequences:
            # If the sequence is not accounted for and not included
            if sequence not in accounted_seqs and sequence not in included_seqs:
                included_seqs.add(sequence)
                # Add all the neighbor into accounted seqs
                for neighbor in G.neighbors(sequence):
                    accounted_seqs.add(neighbor)

        # Update the p3_end_dict to contain the downsampled tails
        p3_end_dict[p3_end] = included_seqs

    seqs = set()
    ## Regenerate all the sequences from p3_end_dict
    for k, v in p3_end_dict.items():
        for seq in v:
            seqs.add(f"{seq}{k}")
    return seqs


def digest(
    msa_array: np.ndarray,
    cfg: dict,
    indexes: tuple[list[int], list[int]] | None = None,
    logger: None = None,
) -> tuple[list[FKmer], list[RKmer]]:
    """
    Digest the given MSA array and return the FKmers and RKmers.

    :param msa_array: The input MSA array.
    :param cfg: A dictionary containing configuration parameters.
    :param indexes: A tuple of MSA indexes for (FKmers, RKmers), or None to use all indexes.
    :param logger: None or the Logguru logger object.
    :return: A tuple containing lists of sorted FKmers and RKmers.
    """
    # Guard for invalid indexes
    if indexes is not None:
        if min(indexes[0]) < 0 or max(indexes[1]) >= msa_array.shape[1]:
            raise IndexError("FIndexes are out of range")
        if min(indexes[1]) < 0 or max(indexes[1]) >= msa_array.shape[1]:
            raise IndexError("RIndexes are out of range")

    # Get the indexes to digest
    findexes = (
        indexes[0]
        if indexes is not None
        else range(cfg["primer_size_min"], msa_array.shape[1])
    )
    rindexes = (
        indexes[1]
        if indexes is not None
        else range(msa_array.shape[1] - cfg["primer_size_min"])
    )

    # Create the MP Pool
    with Pool(cfg["n_cores"]) as p:
        # Generate the FKmers via MP
        fprimer_mp = p.map(
            mp_f_digest,
            ((msa_array, cfg, end_col, cfg["minbasefreq"]) for end_col in findexes),
        )

        pass_fprimer_mp = [x for x in fprimer_mp if type(x) is FKmer and x.seqs]
        pass_fprimer_mp.sort(key=lambda fkmer: fkmer.end)

        # Generate the FKmers via MP
        rprimer_mp = p.map(
            mp_r_digest,
            ((msa_array, cfg, start_col, cfg["minbasefreq"]) for start_col in rindexes),
        )
        pass_rprimer_mp = [x for x in rprimer_mp if type(x) is RKmer and x.seqs]
        # mp_thermo_pass_rkmers = [x for x in rprimer_mp if x is not None]
        pass_rprimer_mp.sort(key=lambda rkmer: rkmer.start)

        # If a logger has been provided dumb the error stats
        if logger is not None:
            # Log the fkmer errors
            for fkmer_result in fprimer_mp:
                if type(fkmer_result) is tuple:
                    logger.debug(
                        "FKmer: <red>{end_col}\t{error}</>",
                        end_col=fkmer_result[0],
                        error=fkmer_result[1].value,
                    )
                else:
                    logger.debug(
                        "FKmer: <green>{end_col}</>: AllPass", end_col=fkmer_result.end  # type: ignore
                    )
            # log the rkmer errors
            for rkmer_result in rprimer_mp:
                if type(rkmer_result) is tuple:
                    logger.debug(
                        "RKmer: <red>{start_col}\t{error}</>",
                        start_col=rkmer_result[0],
                        error=rkmer_result[1].value,
                    )
                else:
                    logger.debug(
                        "RKmer: <green>{start_col}</>: AllPass",
                        start_col=rkmer_result.start,  # type: ignore
                    )

    return (pass_fprimer_mp, pass_rprimer_mp)
