# Modules
from primal_digest.thermo import calc_tm
from primal_digest.classes import RKmer, FKmer, PrimerPair
from primal_digest.thermo import calc_tm, thermo_check_kmers, forms_hairpin
from primal_digest.seq_functions import expand_ambs, get_most_common_base
from primal_digest.get_window import get_r_window_FAST2
from primal_digest.errors import (
    WalksOut,
    GapOnSetBase,
    ContainsInvalidBase,
    CustomErrors,
    ERROR_SET,
    CustomRecursionError,
)

# Submodules
from primaldimer_py import do_pools_interact_py


# Externals
from tqdm import tqdm
import numpy as np
from multiprocessing import Pool
import itertools
import networkx as nx
from collections import Counter


def _mp_pp_inter_free(data: tuple[PrimerPair, dict]) -> bool:
    """
    True means interaction
    """
    pp = data[0]
    cfg = data[1]
    # If they do interact return None
    return do_pools_interact_py(
        list(pp.fprimer.seqs), list(pp.rprimer.seqs), cfg["dimerscore"]
    )


def generate_valid_primerpairs(
    fkmers: list[FKmer],
    rkmers: list[RKmer],
    cfg: dict,
    msa_index: int,
    disable_progress_bar: bool = False,
) -> list[PrimerPair]:
    """
    Generates all valid primers pairs, and then runs interaction checker and returns passing PP
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
    with Pool(cfg["n_cores"]) as p:
        mp_pp_bool = tqdm(
            p.imap_unordered(_mp_pp_inter_free, ((pp, cfg) for pp in non_checked_pp)),
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
) -> set[str] | Exception:
    # Guard for correct tm
    if calc_tm(seq_str, cfg) >= cfg["primer_tm_min"]:
        return {seq_str}

    # Guard prevents walking out of array size
    if col_index_right >= array.shape[1] - 1:
        raise WalksOut()

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

    return passing_str


def walk_left(
    array: np.ndarray,
    col_index_right: int,
    col_index_left: int,
    row_index: int,
    seq_str: str,
    cfg: dict,
) -> set[str] | Exception:
    """
    This will take a string and its indexes and will recurvisly walk left
    until either tm is reached
    """

    # Guard for correct tm
    if calc_tm(seq_str, cfg) >= cfg["primer_tm_min"]:
        return {seq_str}

    # Guard prevents walking out of array size
    if col_index_left <= 0:
        raise WalksOut()

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

    return passing_str


def wrap_walk(
    walkfunction,
    array,
    col_index_right,
    col_index_left,
    row_index,
    seq_str,
    cfg,
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
    except RecursionError:
        return_list.append(CustomRecursionError())
    else:
        return_list.extend(seqs)

    return return_list


def mp_r_digest(data: tuple[np.ndarray, dict, int, float]) -> RKmer | None:
    align_array: np.ndarray = data[0]
    cfg: dict = data[1]
    start_col = data[2]
    min_freq = data[3]

    # Check for gap frequency on first base
    first_base_counter = Counter(align_array[:, start_col])
    del first_base_counter[""]
    num_seqs = sum(first_base_counter.values())
    first_base_freq = {k: v / num_seqs for k, v in first_base_counter.items()}

    # If the freq of gap is above minfreq
    if first_base_freq.get("-", 0) > min_freq:
        return None

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
        # If all mutations matter, return None on any Error
        if min_freq == 0 and set(results) & ERROR_SET:
            return None

        # Add the results to the Counter
        total_col_seqs.update(results)

    total_values = sum(total_col_seqs.values())
    # Filter out values below the threshold freq
    wanted_seqs = {k for k, v in total_col_seqs.items() if v / total_values > min_freq}

    # Guard: If wanted_seqs contains errors return None
    if ERROR_SET & wanted_seqs:
        return None

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
    if (
        thermo_check_kmers(tmp_kmer.seqs, cfg)
        and not forms_hairpin(tmp_kmer.seqs, cfg=cfg)
        and not do_pools_interact_py(
            list(tmp_kmer.seqs), list(tmp_kmer.seqs), cfg["dimerscore"]
        )
    ):
        return tmp_kmer
    else:
        return None


def mp_f_digest(data: tuple[np.ndarray, dict, int, int]) -> FKmer | None:
    align_array: np.ndarray = data[0]
    cfg: dict = data[1]
    end_col = data[2]
    min_freq = data[3]

    # Check for gap frequency on first base
    first_base_counter = Counter(align_array[:, end_col])
    del first_base_counter[""]
    num_seqs = sum(first_base_counter.values())
    first_base_freq = {k: v / num_seqs for k, v in first_base_counter.items()}

    # If the freq of gap is above minfreq
    if first_base_freq.get("-", 0) > min_freq:
        return None

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
            return None

        total_col_seqs.update(results)

    total_values = sum(total_col_seqs.values())
    # Filter out values below the threshold freq
    wanted_seqs = {k for k, v in total_col_seqs.items() if v / total_values > min_freq}

    # Guard: If wanted_seqs contains errors return None
    if ERROR_SET & wanted_seqs:
        return None

    # DownSample the seqs if asked
    if cfg["reducekmers"]:
        wanted_seqs = reduce_kmers(
            seqs=wanted_seqs,
            max_edit_dist=cfg["editdist_max"],
            end_3p=cfg["editdist_end3p"],
        )
    # Thermo check the kmers
    if (
        thermo_check_kmers(wanted_seqs, cfg)
        and not forms_hairpin(wanted_seqs, cfg=cfg)
        and not do_pools_interact_py(
            list(wanted_seqs), list(wanted_seqs), cfg["dimerscore"]
        )
    ):
        return FKmer(end=end_col, seqs=wanted_seqs)
    else:
        return None


def hamming_dist(s1, s2) -> int:
    """
    Return the number of subsitutions, starting from the 3p end
    """
    return sum((x != y for x, y in zip(s1[::-1], s2[::-1])))


def reduce_kmers(seqs: set[str], max_edit_dist: int = 1, end_3p: int = 6) -> set[str]:
    ## Cluster sequences via the 3p end
    p3_end_dict: dict[str : set[str]] = {}
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
    msa_array,
    cfg,
    indexes: tuple[list[int], list[int]] | bool = False,
) -> tuple[list[FKmer], list[RKmer]]:
    """
    Digest the given MSA array and return the FKmers and RKmers.

    :param msa_array: The input MSA array.
    :param cfg: A dictionary containing configuration parameters.
    :param indexes: A tuple of MSA indexes for (FKmers, RKmers), or False to use all indexes.
    :param disable_progress_bar: True to disable to the progress bar.
    :return: A tuple containing lists of sorted FKmers and RKmers.
    """
    # Get the indexes to digest
    findexes = (
        indexes[0] if indexes else range(cfg["primer_size_min"], msa_array.shape[1])
    )
    rindexes = (
        indexes[1] if indexes else range(msa_array.shape[1] - cfg["primer_size_min"])
    )

    # Create the MP Pool
    with Pool(cfg["n_cores"]) as p:
        # Generate the FKmers via MP
        fprimer_mp = p.map(
            mp_f_digest,
            ((msa_array, cfg, end_col, cfg["minbasefreq"]) for end_col in findexes),
        )
        pass_fprimer_mp = [x for x in fprimer_mp if x is not None and x.seqs]
        pass_fprimer_mp.sort(key=lambda fkmer: fkmer.end)

        # Generate the FKmers via MP
        rprimer_mp = p.map(
            mp_r_digest,
            ((msa_array, cfg, start_col, cfg["minbasefreq"]) for start_col in rindexes),
        )
        pass_rprimer_mp = [x for x in rprimer_mp if x is not None and x.seqs]
        # mp_thermo_pass_rkmers = [x for x in rprimer_mp if x is not None]
        pass_rprimer_mp.sort(key=lambda rkmer: rkmer.start)

    return (pass_fprimer_mp, pass_rprimer_mp)
