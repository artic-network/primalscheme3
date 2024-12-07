# Modules
from collections import Counter
from enum import Enum
from typing import Callable, Union

import numpy as np

# Submodules
from primaldimer_py import (
    do_pools_interact_py,  # type: ignore
    which_kmers_pools_interact,  # type: ignore
)

from primalscheme3.core.classes import FKmer, PrimerPair, RKmer
from primalscheme3.core.config import Config
from primalscheme3.core.errors import (
    ERROR_SET,
    ContainsInvalidBase,
    CustomErrors,
    CustomRecursionError,
    GapOnSetBase,
    WalksOut,
    WalksTooFar,
)
from primalscheme3.core.get_window import get_r_window_FAST2
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.core.seq_functions import (
    expand_ambs,
    get_most_common_base,
    reverse_complement,
)
from primalscheme3.core.thermo import (
    THERMORESULT,
    calc_tm,
    thermo_check_kmers,
)


class DIGESTION_ERROR(Enum):
    """
    Enum for the different types of errors that can occur during digestion
    """

    WALKS_OUT = "WalksOut"
    CONTAINS_INVALID_BASE = "ContainsInvalidBase"
    CUSTOM_RECURSION_ERROR = "CustomRecursionError"
    CUSTOM_ERRORS = "CustomErrors"
    GAP_ON_SET_BASE = "GapOnSetBase"
    HAIRPIN_FAIL = "HairpinFail"
    DIMER_FAIL = "DimerFail"  # Interaction within the kmer
    WALK_TO_FAR = "WalkToFar"  # When indels causes the walk to go to far
    AMB_FAIL = "AmbFail"  # Generic error for when the error is unknown

    # Thermo errors
    THERMO_HIGH_GC = "HighGC"
    THERMO_LOW_GC = "LowGC"
    THERMO_HIGH_TM = "HighTM"
    THERMO_LOW_TM = "LowTM"
    THERMO_MAX_HOMOPOLY = "MaxHomopoly"
    THERMO_HAIRPIN = "Hairpin"
    NO_SEQUENCES = "NoSequences"


def parse_error(results: set[CustomErrors | str]) -> DIGESTION_ERROR:
    """
    Parses the error set for the error that occurred
    As only one error is returned, there is an arbitrary hierarchy of errors
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


def parse_thermo_error(result: THERMORESULT) -> DIGESTION_ERROR:
    """
    Parses the THERMORESULT for the error that occurred
    """
    match result:
        case THERMORESULT.HIGH_GC:
            return DIGESTION_ERROR.THERMO_HIGH_GC
        case THERMORESULT.LOW_GC:
            return DIGESTION_ERROR.THERMO_LOW_GC
        case THERMORESULT.HIGH_TM:
            return DIGESTION_ERROR.THERMO_HIGH_TM
        case THERMORESULT.LOW_TM:
            return DIGESTION_ERROR.THERMO_LOW_TM
        case THERMORESULT.MAX_HOMOPOLY:
            return DIGESTION_ERROR.THERMO_MAX_HOMOPOLY
        case THERMORESULT.HAIRPIN:
            return DIGESTION_ERROR.THERMO_HAIRPIN
        case _:
            raise ValueError("Unknown error occurred")


def parse_error_list(
    error_list: list[str | CustomErrors],
) -> list[str | DIGESTION_ERROR]:
    """
    Parses a list of errors and returns a list of DIGESTION_ERROR
    """
    return_list = []
    for result in error_list:
        if isinstance(result, str):
            return_list.append(result)
        elif isinstance(result, CustomErrors):
            return_list.append(parse_error({result}))
    return return_list


def generate_valid_primerpairs(
    fkmers: list[FKmer],
    rkmers: list[RKmer],
    amplicon_size_min: int,
    amplicon_size_max: int,
    dimerscore: float,
    msa_index: int,
    progress_manager: ProgressManager,
    chrom: str = "",
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
    checked_pp = []
    pt = progress_manager.create_sub_progress(
        iter=fkmers, process="Generating primer pairs", chrom=chrom
    )
    for fkmer in pt:
        fkmer_start = min(fkmer.starts())
        # Get all rkmers that would make a valid amplicon
        pos_rkmer = get_r_window_FAST2(
            kmers=rkmers,
            start=fkmer_start + amplicon_size_min,
            end=fkmer_start + amplicon_size_max,
        )
        for rkmer in pos_rkmer:
            # Check for interactions
            if not which_kmers_pools_interact(
                [fkmer], [rkmer], dimerscore, calc_all=False
            ):
                checked_pp.append(PrimerPair(fkmer, rkmer, msa_index))

        # Update the count
        pt.manual_update(count=len(checked_pp))

    checked_pp.sort(key=lambda pp: (pp.fprimer.end, -pp.rprimer.start))
    return checked_pp


def walk_right(
    array: np.ndarray,
    col_index_right: int,
    col_index_left: int,
    row_index: int,
    seq_str: str,
    config: Config,
) -> Union[set[str], Exception]:
    """
    Walks to the right of the array and returns a set of valid sequences.

    Args:
        array: A numpy array of DNA sequences.
        col_index_right: The current column index to the right.
        col_index_left: The current column index to the left.
        row_index: The current row index.
        seq_str: The current sequence string.
        config: The configuration object.

    Returns:
        A set of valid DNA sequences or an exception if an error occurs.

    Raises:
        WalksOut: If the function walks out of the array size.
        ContainsInvalidBase: If the sequence contains an invalid base.
    """
    # Guard for correct tm
    if (
        calc_tm(
            seq_str,
            mv_conc=config.mv_conc,
            dv_conc=config.dv_conc,
            dna_conc=config.dna_conc,
            dntp_conc=config.dntp_conc,
        )
        >= config.primer_tm_min
    ):
        return {seq_str}

    # Guard prevents walking out of array size
    if col_index_right >= array.shape[1] - 1 or col_index_left >= array.shape[1] - 1:
        raise WalksOut()

    # Guard for walking too far
    if col_index_right - col_index_left >= config.primer_max_walk:
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
            config,
        )
        passing_str.extend(results)

    return passing_str  # type: ignore


def walk_left(
    array: np.ndarray,
    col_index_right: int,
    col_index_left: int,
    row_index: int,
    seq_str: str,
    config: Config,
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
    if (
        calc_tm(
            seq_str,
            mv_conc=config.mv_conc,
            dv_conc=config.dv_conc,
            dna_conc=config.dna_conc,
            dntp_conc=config.dntp_conc,
        )
        >= config.primer_tm_min
    ):
        return {seq_str}

    # Guard for walking too far
    if col_index_right - col_index_left >= config.primer_max_walk:
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
            config=config,
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
    config: Config,
) -> list[str | CustomErrors]:
    return_list = []
    try:
        seqs = walkfunction(
            array=array,
            col_index_right=col_index_right,
            col_index_left=col_index_left,
            row_index=row_index,
            seq_str=seq_str,
            config=config,
        )
    except CustomErrors as e:
        return_list.append(e)
    except Exception as e:
        raise e
    else:
        return_list.extend(seqs)

    return return_list


def r_digest_to_count(
    align_array: np.ndarray, config: Config, start_col: int, min_freq: float
) -> tuple[int, dict[str | DIGESTION_ERROR, int]]:
    """
    Returns the count of each sequence / error at a given index
    A value of -1 in the return dict means the function returned early, and not all seqs were counted. Only used for WALKS_OUT and GAP_ON_SET_BASE
    """

    ### Process early return conditions
    # If the initial slice is outside the range of the array
    if start_col + config.primer_size_min >= align_array.shape[1]:
        return (start_col, {DIGESTION_ERROR.WALKS_OUT: -1})

    # Check for gap frequency on first base
    base, counts = np.unique(align_array[:, start_col], return_counts=True)
    first_base_counter = dict(zip(base, counts))
    first_base_counter.pop("", None)

    num_seqs = np.sum(counts)
    first_base_freq = {k: v / num_seqs for k, v in first_base_counter.items()}

    # If the freq of gap is above minfreq
    if first_base_freq.get("-", 0) > min_freq:
        return (start_col, {DIGESTION_ERROR.GAP_ON_SET_BASE: -1})

    ### Calculate the total number of sequences
    # Create a counter
    total_col_seqs: Counter[str | DIGESTION_ERROR] = Counter()
    for row_index in range(0, align_array.shape[0]):
        # Check if this row starts on a gap, and if so update the counter and skip
        if align_array[row_index, start_col] == "-":
            total_col_seqs.update([DIGESTION_ERROR.GAP_ON_SET_BASE])
            continue

        start_array = align_array[
            row_index, start_col : start_col + config.primer_size_min
        ]
        start_seq = "".join(start_array).replace("-", "")

        # Prevent Ns from being added
        if "N" in start_seq:
            total_col_seqs.update([DIGESTION_ERROR.CONTAINS_INVALID_BASE])
            continue

        if not start_seq:  # If the start seq is empty go to the next row
            continue

        # Get all sequences
        results = wrap_walk(
            walk_right,
            array=align_array,
            col_index_right=start_col + config.primer_size_min,
            col_index_left=start_col,
            row_index=row_index,
            seq_str=start_seq,
            config=config,
        )
        # If all mutations matter, return on any Error
        if min_freq == 0 and set(results) & ERROR_SET:
            return (start_col, {parse_error(set(results)): -1})

        # Add the results to the Counter
        total_col_seqs.update(
            {seq: 1 / len(results) for seq in parse_error_list(results)}
        )

    return (start_col, dict(total_col_seqs))


def process_seqs(
    seq_counts: dict[str | DIGESTION_ERROR, int], min_freq, ignore_n: bool = False
) -> DIGESTION_ERROR | dict[str, float]:
    """
    Takes the output from *_digest_to_count and returns a set of valid sequences. Or the error that occurred.

    Args:
        col (int): The column number.
        seq_counts (dict[str | DIGESTION_ERROR, int]): A dictionary containing sequence counts.
        min_freq: The minimum frequency threshold.

    Returns:
        DIGESTION_ERROR | dict[str, float]: either an error or a dictionary of parsed sequences.
    """
    # Check for early return conditions
    for error, count in seq_counts.items():
        if count == -1 and isinstance(error, DIGESTION_ERROR):
            return error

    # Remove Ns if asked
    if ignore_n:
        seq_counts.pop(DIGESTION_ERROR.CONTAINS_INVALID_BASE, None)

    total_values = sum(seq_counts.values())
    # Filter out values below the threshold freq
    above_freq_seqs: dict[str | DIGESTION_ERROR, float] = {
        k: v / total_values
        for k, v in seq_counts.items()
        if v / total_values > min_freq
    }
    parsed_seqs: dict[str, float] = {}
    # Guard: If wanted_seqs contains errors return None
    for seq in above_freq_seqs.keys():
        if isinstance(seq, DIGESTION_ERROR):
            return seq
        elif isinstance(seq, str):
            parsed_seqs[seq] = above_freq_seqs[seq]

    return parsed_seqs


def r_digest_index(
    align_array: np.ndarray, config: Config, start_col: int, min_freq: float
) -> RKmer | tuple[int, DIGESTION_ERROR]:
    """
    This will try and create a RKmer started at the given index
    :align_array: The alignment array
    :config: The configuration object
    :start_col: The column index to start the RKmer
    :min_freq: The minimum frequency threshold

    :return: A RKmer object or a tuple of (start_col, error)
    """
    # Count how many times each sequence / error occurs
    _start_col, seq_counts = r_digest_to_count(align_array, config, start_col, min_freq)
    tmp_parsed_seqs = process_seqs(seq_counts, min_freq, ignore_n=config.ignore_n)
    if isinstance(tmp_parsed_seqs, DIGESTION_ERROR):
        return (start_col, tmp_parsed_seqs)
    elif isinstance(tmp_parsed_seqs, dict):
        parsed_seqs = tmp_parsed_seqs
    else:
        raise ValueError("Unknown error occurred")

    if not parsed_seqs:
        return (start_col, DIGESTION_ERROR.NO_SEQUENCES)

    # Create the Kmer
    rc_seqs = [reverse_complement(seq) for seq in parsed_seqs.keys()]
    tmp_kmer = RKmer(start_col, rc_seqs)

    # Thermo check the kmers
    thermo_result = thermo_check_kmers(tmp_kmer.seqs, config)
    match thermo_result:
        case THERMORESULT.PASS:
            pass
        case _:
            return (start_col, parse_thermo_error(thermo_result))

    # Check for dimer
    if do_pools_interact_py([*tmp_kmer.seqs], [*tmp_kmer.seqs], config.dimer_score):
        return (start_col, DIGESTION_ERROR.DIMER_FAIL)
    # All checks pass return the kmer
    return tmp_kmer


def f_digest_to_count(
    align_array: np.ndarray, config: Config, end_col: int, min_freq: float
) -> tuple[int, dict[str | DIGESTION_ERROR, int]]:
    """
    This will try and create a FKmer ended at the given index
    :return: A FKmer object or a tuple of (end_col, error)
    """

    # Check for gap frequency on first base
    base, counts = np.unique(
        align_array[:, end_col], return_counts=True
    )  # -1 for non-inclusive end
    first_base_counter = dict(zip(base, counts))
    first_base_counter.pop("", None)

    num_seqs = np.sum(counts)

    first_base_freq = {k: v / num_seqs for k, v in first_base_counter.items()}

    # If the freq of gap is above minfreq
    if first_base_freq.get("-", 0) > min_freq:
        return (end_col, {DIGESTION_ERROR.GAP_ON_SET_BASE: -1})

    # If the initial slice is outside the range of the array
    if end_col - config.primer_size_min < 0:
        return (end_col, {DIGESTION_ERROR.WALKS_OUT: -1})

    total_col_seqs: Counter[str | DIGESTION_ERROR] = Counter()
    for row_index in range(0, align_array.shape[0]):
        # Check if this row starts on a gap, and if so update the counter and skip
        if align_array[row_index, end_col] == "-":
            total_col_seqs.update([DIGESTION_ERROR.GAP_ON_SET_BASE])
            # Skip to next row
            continue

        start_seq = "".join(
            align_array[row_index, end_col - config.primer_size_min : end_col]
        ).replace("-", "")

        # Prevent Ns from being added
        if "N" in start_seq:
            total_col_seqs.update([DIGESTION_ERROR.CONTAINS_INVALID_BASE])
            continue

        if not start_seq:  # If the start seq is empty go to the next row
            continue

        results = wrap_walk(
            walk_left,
            array=align_array,
            col_index_right=end_col,
            col_index_left=end_col - config.primer_size_min,
            row_index=row_index,
            seq_str=start_seq,
            config=config,
        )
        if min_freq == 0 and set(results) & ERROR_SET:
            return (end_col, {parse_error(set(results)): -1})

        # Add the results to the Counter
        total_col_seqs.update(
            {seq: 1 / len(results) for seq in parse_error_list(results)}
        )

    return (end_col, dict(total_col_seqs))


def f_digest_index(
    align_array: np.ndarray, config: Config, end_col: int, min_freq: float
) -> FKmer | tuple[int, DIGESTION_ERROR]:
    """
    This will try and create a FKmer ended at the given index
    :align_array: The alignment array
    :config: The configuration object
    :end_col: The column index to end the FKmer
    :min_freq: The minimum frequency threshold

    :return: A FKmer object or a tuple of (end_col, error)
    """

    # Count how many times each sequence / error occurs
    _end_col, seq_counts = f_digest_to_count(align_array, config, end_col, min_freq)
    tmp_parsed_seqs = process_seqs(seq_counts, min_freq, ignore_n=config.ignore_n)
    if isinstance(tmp_parsed_seqs, DIGESTION_ERROR):
        return (end_col, tmp_parsed_seqs)
    elif isinstance(tmp_parsed_seqs, dict):
        parsed_seqs = tmp_parsed_seqs
    else:
        raise ValueError("Unknown error occurred")

    # Thermo check the kmers
    thermo_result = thermo_check_kmers({*parsed_seqs.keys()}, config)
    match thermo_result:
        case THERMORESULT.PASS:
            pass
        case _:
            return (end_col, parse_thermo_error(thermo_result))

    if do_pools_interact_py(
        [*parsed_seqs.keys()], [*parsed_seqs.keys()], config.dimer_score
    ):
        return (end_col, DIGESTION_ERROR.DIMER_FAIL)

    if not parsed_seqs:
        return (end_col, DIGESTION_ERROR.NO_SEQUENCES)

    return FKmer(end_col, list(parsed_seqs.keys()))


def hamming_dist(s1, s2) -> int:
    """
    Return the number of substitutions, starting from the 3p end
    """
    return sum((x != y for x, y in zip(s1[::-1], s2[::-1])))


def f_digest(
    msa_array: np.ndarray, config: Config, findexes: list[int], logger
) -> list[FKmer]:
    fkmers = []
    for findex in findexes:
        fkmer = f_digest_index(msa_array, config, findex, config.min_base_freq)

        # Append valid FKmers
        if isinstance(fkmer, FKmer) and fkmer.seqs:
            fkmers.append(fkmer)

        # Log the Digestion
        if logger is not None:
            if isinstance(fkmer, tuple):
                logger.debug(f"FKmer: [red]{fkmer[0]}[/red]\t{fkmer[1].value}")
            else:
                logger.debug(f"FKmer: [green]{fkmer.end}[/green]: AllPass")

    return fkmers


def r_digest(
    msa_array: np.ndarray, config: Config, rindexes: list[int], logger
) -> list[RKmer]:
    rkmers = []
    for rindex in rindexes:
        rkmer = r_digest_index(msa_array, config, rindex, config.min_base_freq)

        # Append valid RKmers
        if isinstance(rkmer, RKmer) and rkmer.seqs:
            rkmers.append(rkmer)

        # Log the Digestion
        if logger is not None:
            if isinstance(rkmer, tuple):
                logger.debug(f"RKmer: [red]{rkmer[0]}[/red]\t{rkmer[1].value}")
            else:
                logger.debug(f"RKmer: [green]{rkmer.start}[/green]: AllPass")
    return rkmers


def digest(
    msa_array: np.ndarray,
    config: Config,
    progress_manager: ProgressManager,
    indexes: tuple[list[int], list[int]] | None = None,
    logger=None,
    chrom: str = "",
) -> tuple[list[FKmer], list[RKmer]]:
    """
    Digest the given MSA array and return the FKmers and RKmers.

    :param msa_array: The input MSA array.
    :param cfg: A dictionary containing configuration parameters.
    :param indexes: A tuple of MSA indexes for (FKmers, RKmers), or None to use all indexes.
    :param logger: None or the logger object.
    :return: A tuple containing lists of sorted FKmers and RKmers.
    """
    # Guard for invalid indexes
    if indexes is not None:
        if min(indexes[0]) < 0 or max(indexes[0]) >= msa_array.shape[1]:
            raise IndexError("FIndexes are out of range")
        if min(indexes[1]) < 0 or max(indexes[1]) >= msa_array.shape[1]:
            raise IndexError("RIndexes are out of range")

    # Get the indexes to digest
    findexes = (
        indexes[0]
        if indexes is not None
        else range(config.primer_size_min, msa_array.shape[1])
    )
    rindexes = (
        indexes[1]
        if indexes is not None
        else range(msa_array.shape[1] - config.primer_size_min)
    )

    # Digest the findexes
    fkmers = []
    pt = progress_manager.create_sub_progress(
        iter=findexes, process="Creating forward primers", chrom=chrom
    )
    for findex in pt:
        fkmer = f_digest_index(msa_array, config, findex, config.min_base_freq)

        if logger is not None:
            if isinstance(fkmer, tuple):
                logger.debug(f"{chrom}:FKmer: {fkmer[0]}\t{fkmer[1].value}")
            else:
                logger.debug(f"{chrom}:FKmer: {fkmer.end}\tPass")

        # Append valid FKmers
        if isinstance(fkmer, FKmer) and fkmer.seqs:
            fkmers.append(fkmer)

        # Update the count
        pt.manual_update(count=len(fkmers))

    # Digest the rindexes
    rkmers = []
    pt = progress_manager.create_sub_progress(
        iter=rindexes, process="Creating reverse primers", chrom=chrom
    )
    for rindex in pt:
        rkmer = r_digest_index(msa_array, config, rindex, config.min_base_freq)

        if logger is not None:
            if isinstance(rkmer, tuple):
                logger.debug(f"{chrom}:RKmer: {rkmer[0]}\t{rkmer[1].value}")
            else:
                logger.debug(f"{chrom}:RKmer: {rkmer.start}\tPass")

        # Append valid RKmers
        if isinstance(rkmer, RKmer) and rkmer.seqs:
            rkmers.append(rkmer)

        # Update the count
        pt.manual_update(count=len(rkmers))

    return (fkmers, rkmers)
