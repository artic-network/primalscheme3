from primal_digest.thermo import *
from primal_digest.cli import cli
from primal_digest.config import (
    AMBIGUOUS_DNA_COMPLEMENT,
    config_dict,
    thermo_config,
    ALL_DNA,
)
import numpy as np
from Bio import SeqIO
from collections import Counter
from primal_digest.iteraction import all_inter_checker
from Bio import Seq


from multiprocessing import Pool
import kmertools
import json

# Added
from diskcache import Cache
import hashlib
import sys
import pathlib
from itertools import product

"""
This is a test of a new dynamic digestion algo
"""


def remove_end_insertion(msa_array):
    tmp_array = msa_array
    ncols = tmp_array.shape[1]
    for row_index in range(0, tmp_array.shape[0]):
        # Solves the 5' end
        for col_index in range(0, ncols):
            if tmp_array[row_index, col_index] == "-":
                tmp_array[row_index, col_index] = ""
            else:
                break
        for rev_col_index in range(ncols - 1, 0, -1):
            if tmp_array[row_index, rev_col_index] == "-":
                tmp_array[row_index, rev_col_index] = ""
            else:
                break
    return tmp_array


def expand_ambs(seqs: Iterable[str]) -> set[set]:
    """
    Takes a list / set of strings and returns a set wuth all ambs expanded
    It can return an empty set
    """
    returned_seq = set()
    amb_bases = {"Y", "W", "R", "B", "H", "V", "D", "K", "M", "S"}
    for seq in seqs:
        if seq.__contains__("N"):
            break
        # If there is any amb_bases in the seq
        if {base for base in seq} & amb_bases:
            expanded_seqs = set(map("".join, product(*map(ALL_DNA.get, seq))))
            for exp_seq in expanded_seqs:
                returned_seq.add(exp_seq)
        else:
            returned_seq.add(seq)
    return returned_seq


class FKmer:
    end: int
    seqs: set[str]

    def __init__(self, end, seqs) -> None:
        self.end = end
        self.seqs = seqs

    def len(self) -> set[int]:
        return {len(x) for x in self.seqs}

    def starts(self) -> set[int]:
        return {self.end - len(x) for x in self.seqs}

    def __str__(self, referance, amplicon_prefix, pool) -> str:
        string_list = []
        counter = 0
        seqs = list(self.seqs)
        seqs.sort()
        for seq in seqs:
            string_list.append(
                f"{referance}\t{self.end-len(seq)}\t{self.end}\t{amplicon_prefix}_LEFT_{counter}\t{pool}\t+\t{seq}\n"
            )
            counter += 1
        return "".join(string_list)

    def __hash__(self) -> int:
        seqs = list(self.seqs)
        seqs.sort()
        return hash(f"{self.end}{self.seqs}")

    def __eq__(self, other):
        if isinstance(other, FKmer):
            return self.__hash__() == other.__hash__()


class RKmer:
    start: int
    seqs: set[str]

    def __init__(self, start, seqs) -> None:
        self.start = start
        self.seqs = seqs

    def len(self) -> set[int]:
        return {len(x) for x in self.seqs}

    def ends(self) -> set[int]:
        return {len(x) + self.start for x in self.seqs}

    def __str__(self, referance, amplicon_prefix, pool) -> str:
        string_list = []
        counter = 0
        seqs = list(self.seqs)
        seqs.sort()
        for seq in seqs:
            string_list.append(
                f"{referance}\t{self.start}\t{self.start+len(seq)}\t{amplicon_prefix}_RIGHT_{counter}\t{pool}\t-\t{seq}\n"
            )
            counter += 1
        return "".join(string_list)

    def reverse_complement(self):
        rev_seqs = {x[::-1] for x in self.seqs}
        self.seqs = {
            "".join(AMBIGUOUS_DNA_COMPLEMENT[base.upper()] for base in seq)
            for seq in rev_seqs
        }
        return self

    def __hash__(self) -> int:
        seqs = list(self.seqs)
        seqs.sort()
        return hash(f"{self.start}{self.seqs}")

    def __eq__(self, other):
        if isinstance(other, RKmer):
            return self.__hash__() == other.__hash__()


class PrimerPair:
    fprimer: FKmer
    rprimer: RKmer
    amplicon_number: int
    pool: int
    msa_index: int

    def __init__(self, fprimer, rprimer, amplicon_number=-1, pool=-1, msa_index=-1):
        self.fprimer = fprimer
        self.rprimer = rprimer
        self.amplicon_number = amplicon_number
        self.pool = pool
        self.msa_index = msa_index

    def set_amplicon_number(self, amplicon_number) -> None:
        self.amplicon_number = amplicon_number

    def set_pool_number(self, pool_number) -> None:
        self.amplicon_number = pool_number

    def start(self) -> int:
        return max(self.fprimer.starts())

    def end(self) -> int:
        return min(self.rprimer.ends())

    def inter_free(self, cfg) -> bool:
        """
        True means interaction
        """
        return all_inter_checker(self.fprimer.seqs, self.rprimer.seqs, cfg=cfg)

    def all_seqs(self) -> set[str]:
        return [x for x in self.fprimer.seqs] + [x for x in self.rprimer.seqs]

    def __hash__(self) -> int:
        return hash(f"{self.start}{self.end}{self.all_seqs()}")

    def __str__(self, ref_name, amplicon_prefix):
        return self.fprimer.__str__(
            referance=f"{ref_name}",
            amplicon_prefix=f"{amplicon_prefix}_{self.amplicon_number}",
            pool=self.pool + 1,
        ) + self.rprimer.__str__(
            referance=f"{ref_name}",
            amplicon_prefix=f"{amplicon_prefix}_{self.amplicon_number}",
            pool=self.pool + 1,
        )


def mp_pp_inter_free(data: tuple[PrimerPair, dict]) -> bool:
    """
    True means interaction
    """
    pp = data[0]
    cfg = data[1]
    return all_inter_checker(pp.fprimer.seqs, pp.rprimer.seqs, cfg=cfg)


# These are modifed get_window_fast for the new kmers classes
def get_r_window_FAST2(kmers: list[RKmer], start: int, end: int) -> list[RKmer]:
    """
    This will perform a binary search on the list of kmers.
    The Kmer start positions will be used as the search value
    """
    included_kmers: list[RKmer] = []
    n_kmers = len(kmers)
    high = n_kmers - 1
    low = 0
    mid = 0
    while low <= high:
        mid = (high + low) // 2
        # If the midpoint is inside the range
        if start <= kmers[mid].start <= end:
            while True:
                # Walk back until first value, or the first position
                if mid == 0:
                    break
                elif kmers[mid - 1].start >= start:
                    mid -= 1
                else:
                    break
            # Mid is now the first value so walk forwards
            while True:
                if mid < n_kmers and kmers[mid].start <= end:
                    included_kmers.append(kmers[mid])
                    mid += 1
                else:
                    return included_kmers
        # If start is greater ignore the left half
        elif kmers[mid].start < start:
            low = mid + 1
        # If start is smaller ignore the right half
        elif kmers[mid].start > end:
            high = mid - 1

    # If the code reaches here there are no KMERS within the list inside the range
    ## Return an empty list for continuity
    return []

    # These are modifed get_window_fast for the new kmers classes


def get_f_window_FAST2(kmers: list[FKmer], start: int, end: int) -> list[FKmer]:
    """
    This will perform a binary search on the list of kmers.
    The Kmer end position will be used as the search value
    """
    included_kmers: list[FKmer] = []
    n_kmers = len(kmers)
    high = n_kmers - 1
    low = 0
    mid = 0
    while low <= high:
        mid = (high + low) // 2
        # If the midpoint is inside the range
        if start <= kmers[mid].end <= end:
            while True:
                # Walk back until first value, or the first position
                if mid == 0:
                    break
                elif kmers[mid - 1].end >= start:
                    mid -= 1
                else:
                    break
            # Mid is now the first value so walk forwards
            while True:
                if mid < n_kmers and kmers[mid].end <= end:
                    included_kmers.append(kmers[mid])
                    mid += 1
                else:
                    return included_kmers
        # If start is greater ignore the left half
        elif kmers[mid].end < start:
            low = mid + 1
        # If start is smaller ignore the right half
        elif kmers[mid].end > end:
            high = mid - 1

    # If the code reaches here there are no KMERS within the list inside the range
    ## Return an empty list for continuity
    return []


def get_most_common_base(array: np.ndarray, col_index: int) -> str:
    col_slice = array[:, col_index]
    counter = Counter(col_slice)
    del counter[""]
    return counter.most_common(1)[0][0]


def is_unique(seqs: Iterable[str], start: int, end: int, msa: np.ndarray) -> bool:
    """
    This function should see if the exact Kmer sequence or the RC, is found in any region of any sequences other
    than the expected binding site.


    False means not unquie
    True means exact kmer seqs are not found anywhere other than expected binding site
    """
    # This creates a tmp_array which has the primer binding region removed
    if start == 0:  # If the kmer bind at first pos
        tmp_array = msa[:, end:]
    elif end == msa.shape[1]:  # If the kmers bind at the last positions
        tmp_array = msa[:, :start]
    else:
        left_array = msa[:, :start]
        right_array = msa[:, end:]
        tmp_array = np.concatenate((left_array, right_array), axis=1)

    # Iterate over rows in the tmp msa
    for row in tmp_array:
        seq = "".join(row)

        # for each kmer seq
        for kmer_seq in seqs:
            rc_seq = str(
                Seq.Seq(kmer_seq).reverse_complement()
            )  # TODO FIX THIS to not need import

            if seq.__contains__(kmer_seq) or seq.__contains__(
                rc_seq
            ):  # Only looks for exact matches
                return False

    return True  # If no early break


def walk_right(
    array: np.ndarray,
    col_index_right: int,
    col_index_left: int,
    row_index: int,
    seq_str: str,
    cfg: dict,
) -> set[str]:

    passing_str = set()

    if calc_tm(seq_str, cfg) < cfg["primer_tm_min"]:
        if col_index_right < array.shape[1] - 1:
            new_base = array[row_index, col_index_right]
        else:
            return None

        if new_base == "":
            new_base = get_most_common_base(array, col_index_right + 1)
        new_string = (seq_str + new_base).replace("-", "")

        if new_string.__contains__("N"):
            return None
        else:
            exp_new_string: set[str] = expand_ambs([new_string])

        passing_str = set()
        for exp_str in exp_new_string:
            results = walk_right(
                array,
                col_index_right + 1,
                col_index_left,
                row_index,
                exp_str,
                cfg,
            )
            if results is not None:
                [passing_str.add(x) for x in results]
            else:
                return None
    else:
        return {seq_str}

    if len(passing_str) >= 1:
        return passing_str


def walk_left(
    array: np.ndarray,
    col_index_right: int,
    col_index_left: int,
    row_index: int,
    seq_str: str,
    cfg: dict,
) -> set[str]:
    """
    This will take a string and its indexes and will recurvisly walk left
    until either tm is reached
    """
    passing_str = set()
    if calc_tm(seq_str, cfg) < cfg["primer_tm_min"]:
        if col_index_left > 0:
            new_base = array[row_index, col_index_left - 1]
        else:
            return None

        # Ensure it can repair truncated regions
        if new_base == "":
            new_base = get_most_common_base(array, col_index_left - 1)
        new_string = (new_base + seq_str).replace("-", "")

        if new_string.__contains__("N"):
            return None
        else:
            exp_new_string: set[str] = expand_ambs([new_string])

        passing_str = set()
        for exp_str in exp_new_string:
            results = walk_left(
                array,
                col_index_right,
                col_index_left - 1,
                row_index,
                exp_str,
                cfg,
            )
            if results is not None:
                [passing_str.add(x) for x in results]
            else:
                return None
    else:
        return {seq_str}

    if len(passing_str) >= 1:
        return passing_str


def mp_r_digest(data: tuple[np.ndarray, dict, int]) -> RKmer:
    align_array: np.ndarray = data[0]
    cfg: dict = data[1]
    start_col = data[2]

    total_col_seqs = set()
    for row_index in range(0, align_array.shape[0]):
        start_array = align_array[row_index, start_col : start_col + 14]
        start_seq = "".join(start_array).replace("-", "")

        if (
            start_array[0] == ""
        ):  # If the kmer starts on an invalid base, skip this start position
            total_col_seqs.add(None)
            break

        if not start_seq:  # If the start seq is empty go to the next row
            continue

        seqs = walk_right(
            array=align_array,
            col_index_right=start_col + 14,
            col_index_left=start_col,
            row_index=row_index,
            seq_str=start_seq,
            cfg=cfg,
        )
        # If the seq contains N, expand_amps will return {} which is parsed up the stack as None
        if seqs == None:
            total_col_seqs.add(None)
            break

        # Get a union of the seqs
        if seqs:
            total_col_seqs = total_col_seqs | seqs

    if not total_col_seqs.__contains__(None):
        # Thermo check the kmers
        if thermo_check_kmers(total_col_seqs, cfg) and not forms_hairpin(
            total_col_seqs, cfg=cfg
        ):
            tmp_kmer = RKmer(start=start_col, seqs=total_col_seqs)
            tmp_kmer = tmp_kmer.reverse_complement()
            return tmp_kmer
    else:
        return None


def mp_f_digest(data: tuple[np.ndarray, dict, int]) -> FKmer:
    align_array: np.ndarray = data[0]
    cfg: dict = data[1]
    end_col = data[2]

    total_col_seqs = set()
    for row_index in range(0, align_array.shape[0]):
        start_seq = "".join(align_array[row_index, end_col - 14 : end_col]).replace(
            "-", ""
        )

        if not start_seq:  # If the start seq is empty go to the next row
            continue

        seqs = walk_left(
            array=align_array,
            col_index_right=end_col,
            col_index_left=end_col - 14,
            row_index=row_index,
            seq_str=start_seq,
            cfg=cfg,
        )
        # If the seq contains N, expand_amps will return {} which is parsed up the stack as None
        if seqs == None:
            total_col_seqs.add(None)
            break

            # Get a union of the seqs
        if seqs:
            total_col_seqs = total_col_seqs | seqs

    if not total_col_seqs.__contains__(None):
        # Thermo check the kmers
        if thermo_check_kmers(total_col_seqs, cfg) and not forms_hairpin(
            total_col_seqs, cfg=cfg
        ):
            return FKmer(end=end_col, seqs=total_col_seqs)
        else:
            return None


def main():
    args = cli()
    ARG_MSA = args.msa
    OUTPUT_DIR = pathlib.Path(args.output).absolute()

    cfg = config_dict
    thermo_cfg = thermo_config

    # Primer Digestion settings
    thermo_cfg["kmer_size_max"] = 36
    thermo_cfg["kmer_size_min"] = 14
    thermo_cfg["primer_gc_min"] = args.primer_gc_min
    thermo_cfg["primer_gc_max"] = args.primer_gc_max
    thermo_cfg["primer_tm_min"] = args.primer_tm_min
    thermo_cfg["primer_tm_max"] = args.primer_tm_max

    # Run Settings
    cfg["refname"] = args.refnames
    cfg["n_cores"] = args.cores
    cfg["output_prefix"] = args.prefix
    cfg["output_dir"] = str(OUTPUT_DIR)
    cfg["msa_paths"] = ARG_MSA
    cfg["mismatches_self"] = args.mismatches_self
    cfg["mismatches_alt"] = args.mismatches_alt
    cfg["amplicon_size_max"] = args.ampliconsizemax
    cfg["amplicon_size_min"] = args.ampliconsizemin
    cfg["min_overlap"] = args.minoverlap
    cfg["use_cache"] = args.use_cache
    cfg["force"] = args.force

    msa_index_to_ref_name = {
        index: msa_name for index, msa_name in enumerate(cfg["refname"])
    }

    cfg["msa_index_to_ref_name"] = msa_index_to_ref_name

    # See if the output dir already exsits
    if OUTPUT_DIR.is_dir():
        if not args.force:
            sys.exit(
                f"ERROR: {OUTPUT_DIR} already exists, please use --force to override"
            )

    # Read in the MSAs
    msa_list: list[np.ndarray] = []
    for msa_path in ARG_MSA:
        records = SeqIO.parse(msa_path, "fasta")
        align_array = np.array([record.seq.upper() for record in records], dtype=str)
        align_array = remove_end_insertion(align_array)
        msa_list.append(align_array)

    raw_msa_list: list[np.ndarray] = []
    for msa_path in ARG_MSA:
        records = SeqIO.parse(msa_path, "fasta")
        align_array = np.array([record.seq.upper() for record in records], dtype=str)
        raw_msa_list.append(align_array)

    # Generate the Kmers for each array in msa_list
    unique_f_r_msa: list[list[list[FKmer], list[RKmer]]] = []

    for msa_index, msa_array in enumerate(msa_list):
        # msa_kmers = []

        # Check if the MSA has already been digested. (But not unique checked)
        cache = Cache("./.primal-digest-cache")
        msa_checksum = hashlib.md5(bytes(msa_array.tobytes())).hexdigest()
        tmp_data = cfg["msa_checksums"]
        tmp_data.append((msa_index, msa_checksum))
        cfg["msa_checksums"] = tmp_data  # Add the checksum to the cfg

        thermo_checksum = hashlib.md5(
            json.dumps(thermo_cfg, sort_keys=True).encode()
        ).hexdigest()

        total_checksum = hashlib.md5(
            bytes((msa_checksum + thermo_checksum).encode())
        ).hexdigest()

        if total_checksum in cache and cfg["use_cache"]:
            f_r_kmers = cache[total_checksum]
            print(f"Found cache")
        else:
            ## Mutlicore the digestion for f and r primers
            # FPrimers digestion via MP
            with Pool(cfg["n_cores"]) as p:
                fprimer_mp = p.map(
                    mp_f_digest,
                    [
                        (msa_array, thermo_cfg, end_col)
                        for end_col in range(14, msa_array.shape[1])
                    ],
                )
            mp_thermo_pass_fkmers = [x for x in fprimer_mp if x is not None]
            mp_thermo_pass_fkmers.sort(key=lambda fkmer: fkmer.end)

            # RPrimers digestion via MP
            with Pool(cfg["n_cores"]) as p:
                rprimer_mp = p.map(
                    mp_r_digest,
                    [
                        (msa_array, thermo_cfg, start_col)
                        for start_col in range(msa_array.shape[1] - 14)
                    ],
                )
            mp_thermo_pass_rkmers = [x for x in rprimer_mp if x is not None]
            mp_thermo_pass_rkmers.sort(key=lambda rkmer: rkmer.start)

            f_r_kmers = (mp_thermo_pass_fkmers, mp_thermo_pass_rkmers)
            cache[total_checksum] = f_r_kmers

        # msa_f_r_kmers.append(f_r_kmers)

        mp_thermo_pass_fkmers, mp_thermo_pass_rkmers = f_r_kmers
        # Use the custom Rust unique checker

        f_kmer_bools = []
        r_kmer_bools = []

        for index, raw_msa_array in enumerate(raw_msa_list):
            referance_seqs = ["".join(x) for x in raw_msa_array]
            if msa_index == index:
                r_results = kmertools.rkmer_is_unique(
                    rkmers=[(x.start, list(x.seqs)) for x in mp_thermo_pass_rkmers],
                    referance_seqs=referance_seqs,
                    n_cores=cfg["n_cores"],
                    mismatches=cfg["mismatches_self"],
                    detect_expected=False,
                )  # Should be false
                f_results = kmertools.fkmer_is_unique(
                    fkmers=[(x.end, list(x.seqs)) for x in mp_thermo_pass_fkmers],
                    referance_seqs=referance_seqs,
                    n_cores=cfg["n_cores"],
                    mismatches=cfg["mismatches_self"],
                    detect_expected=False,
                )
                f_kmer_bools.append(f_results)
                r_kmer_bools.append(r_results)

            else:
                r_results = kmertools.rkmer_is_unique(
                    rkmers=[(x.start, list(x.seqs)) for x in mp_thermo_pass_rkmers],
                    referance_seqs=referance_seqs,
                    n_cores=cfg["n_cores"],
                    mismatches=cfg["mismatches_alt"],
                    detect_expected=True,
                )
                f_results = kmertools.fkmer_is_unique(
                    fkmers=[(x.end, list(x.seqs)) for x in mp_thermo_pass_fkmers],
                    referance_seqs=referance_seqs,
                    n_cores=cfg["n_cores"],
                    mismatches=cfg["mismatches_alt"],
                    detect_expected=True,
                )
                f_kmer_bools.append(f_results)
                r_kmer_bools.append(r_results)

        ## Combine the data for each kmer for each msa
        unique_f_bools = [all(x) for x in zip(*f_kmer_bools)]
        unique_r_bools = [all(x) for x in zip(*r_kmer_bools)]

        ## Get the unique kmers
        unique_fkmers = [
            fkmer
            for (bool, fkmer) in zip(unique_f_bools, mp_thermo_pass_fkmers)
            if bool
        ]
        unique_rkmers = [
            rkmer
            for (bool, rkmer) in zip(unique_r_bools, mp_thermo_pass_rkmers)
            if bool
        ]

        # Append them to the msa
        unique_f_r_msa.append((unique_fkmers, unique_rkmers))

    msa_primerpairs = []
    # Generate all valid primerpairs for each msa
    for msa_index, unique_fr_kmers in enumerate(unique_f_r_msa):
        wanted_fkmers = unique_fr_kmers[0]
        wanted_rkmers = unique_fr_kmers[1]

        # Generate all primer pairs
        primer_pairs: list[PrimerPair] = []
        for f in wanted_fkmers:
            pos_r = get_r_window_FAST2(
                kmers=wanted_rkmers,
                start=min(f.starts()) + cfg["amplicon_size_min"],
                end=min(f.starts()) + cfg["amplicon_size_max"],
            )
            for r in pos_r:
                primer_pairs.append(PrimerPair(f, r))

        # Filter out primerpairs if f primer seqs interact with r primer seqs
        with Pool(cfg["n_cores"]) as p:
            mp_pp_bool = p.map(
                mp_pp_inter_free, [(pp, thermo_cfg) for pp in primer_pairs]
            )

        iter_free_primer_pairs: list[PrimerPair] = [
            pp for (bool, pp) in zip(mp_pp_bool, primer_pairs) if not bool
        ]
        iter_free_primer_pairs.sort(key=lambda pp: (pp.start(), -pp.end()))

        msa_primerpairs.append(iter_free_primer_pairs)

    # TODO Start creating the scheme for each msa
    pools: list[list[PrimerPair]] = [[], []]

    current_pool = 0
    last_primer_pair = 0
    for msa_index, primerpairs_in_msa in enumerate(msa_primerpairs):
        counter = 1
        keep_looping = True

        while True and keep_looping:
            new_primerpair = False

            all_seqs_in_pool = [
                y
                for sublist in (x.all_seqs() for x in pools[current_pool])
                for y in sublist
            ]
            if last_primer_pair == 0:  # If picking first amplicon for the pool
                pp = primerpairs_in_msa[0]
                pp.pool = 0
                pp.amplicon_number = counter
                pp.msa_index = msa_index
                pools[pp.pool].append(pp)

                last_primer_pair = pp
                current_pool = (current_pool + 1) % 2
                counter += 1
            elif counter == 1 and last_primer_pair != 0:
                # When the first primer of a new msa is addded to an exsiting pool
                all_seqs_in_other_pool = [
                    y
                    for sublist in (x.all_seqs() for x in pools[(current_pool + 1) % 2])
                    for y in sublist
                ]
                for pp in primerpairs_in_msa:
                    # It doesnt matter which pool this primer goes into, so check both and pick the best.
                    # Try current pool

                    if not all_inter_checker(
                        pp.all_seqs(), all_seqs_in_pool, thermo_cfg
                    ):
                        # If there are no interactions add the primer to the pool
                        pp.pool = current_pool
                        pp.amplicon_number = counter
                        pp.msa_index = msa_index
                        pools[pp.pool].append(pp)

                        last_primer_pair = pp
                        current_pool += 1
                        current_pool = current_pool % 2
                        counter += 1
                        break
                    elif not all_inter_checker(
                        pp.all_seqs(), all_seqs_in_other_pool, thermo_cfg
                    ):
                        # If there are no interactions add the primer to the other pool
                        pp.pool = (current_pool + 1) % 2
                        pp.amplicon_number = counter
                        pp.msa_index = msa_index
                        pools[pp.pool].append(pp)

                        last_primer_pair = pp
                        current_pool = (current_pool + 1) % 2
                        counter += 1
                        break
                    else:
                        # If this primer has interactions with both pools move to next primerpair
                        continue
            elif pools[current_pool] and counter == 2:
                # This if for the second pp in a new msa.
                # There are primers already in the pool (other msa), but we can ignore the last primer in the same pool's end pos,
                # as we do not care about overlaps between msa
                pos_pp = [
                    x
                    for x in primerpairs_in_msa
                    if x.fprimer.end
                    < last_primer_pair.rprimer.start - cfg["min_overlap"]
                    and x.start() > last_primer_pair.fprimer.end
                    and x.rprimer.start > max(last_primer_pair.rprimer.ends())
                ]
                if pos_pp:
                    # If there are primerpairs in the wanted region
                    for pp in pos_pp:
                        if not all_inter_checker(
                            pp.all_seqs(), all_seqs_in_pool, thermo_cfg
                        ):
                            # If there is no interaction, use the pp
                            pp.pool = current_pool
                            pp.amplicon_number = counter
                            pp.msa_index = msa_index
                            pools[pp.pool].append(pp)

                            last_primer_pair = pp
                            current_pool += 1
                            current_pool = current_pool % 2
                            counter += 1

                            new_primerpair = True
                            break

                if new_primerpair == False:
                    # If there are no primerpairs in the wanted region or no interfree pp was found in the region start walking
                    next_pos_pp = [
                        x
                        for x in primerpairs_in_msa
                        if x.fprimer.end
                        > last_primer_pair.rprimer.start - cfg["min_overlap"]
                    ]
                    if not next_pos_pp:
                        # If there are no primerpairs left, stop generation for this msa
                        keep_looping = False
                        break
                    else:
                        # If there are pp left
                        next_pos_pp.sort(
                            key=lambda pp: (pp.fprimer.end, -pp.rprimer.start)
                        )
                        for pp in next_pos_pp:
                            if not all_inter_checker(
                                pp.all_seqs(), all_seqs_in_pool, thermo_cfg
                            ):
                                # If there is no interaction, use the pp
                                pp.pool = current_pool
                                pp.amplicon_number = counter
                                pp.msa_index = msa_index
                                pools[pp.pool].append(pp)

                                last_primer_pair = pp
                                current_pool += 1
                                current_pool = current_pool % 2
                                counter += 1

                                new_primerpair = True
                                break

            else:
                # Do the normal pp forward walk
                print(last_primer_pair.end())

                if pools[current_pool]:  # If there are primers in the same pool
                    last_primer_in_same_pool = pools[current_pool][-1]
                    pos_pp = [
                        x
                        for x in primerpairs_in_msa
                        if x.fprimer.end
                        < last_primer_pair.rprimer.start - cfg["min_overlap"]
                        and x.start() > max(last_primer_in_same_pool.rprimer.ends())
                        and x.rprimer.start > max(last_primer_pair.rprimer.ends())
                    ]
                else:
                    pos_pp = [
                        x
                        for x in primerpairs_in_msa
                        if x.fprimer.end
                        < last_primer_pair.rprimer.start - cfg["min_overlap"]
                        and x.start() > last_primer_pair.fprimer.end
                        and x.rprimer.start > max(last_primer_pair.rprimer.ends())
                    ]  # Stops the pp from using the same rprimer

                if not pos_pp:  # If there are no primerpairs left in the wanted region
                    # Find the next valid primerpair
                    pos_pp = [
                        x
                        for x in primerpairs_in_msa
                        if x.fprimer.end
                        > last_primer_pair.rprimer.start - cfg["min_overlap"]
                    ]
                    # If there are no more primerpairs left
                    if not pos_pp:
                        keep_looping = False
                        break
                    pos_pp.sort(key=lambda pp: (pp.fprimer.end, -pp.rprimer.start))

                    for pp in pos_pp:
                        if all_inter_checker(
                            pp.all_seqs(), all_seqs_in_pool, thermo_cfg
                        ):
                            continue  # If interaction skip
                        else:
                            pp.pool = current_pool
                            pp.amplicon_number = counter
                            pp.msa_index = msa_index
                            pools[pp.pool].append(pp)
                            last_primer_pair = pp
                            current_pool = (current_pool + 1) % 2
                            counter += 1

                            new_primerpair = True
                            break

                    if (
                        new_primerpair
                    ):  # If a gapped primerpair has been found, then jump to the next iteration
                        continue

                else:
                    pos_pp.sort(
                        key=lambda pp: (
                            (last_primer_pair.rprimer.start - pp.fprimer.end)
                            * len(pp.all_seqs()),
                            -pp.rprimer.start,
                        )
                    )

                if not pools[current_pool]:  # if the current pool is empty and
                    pp = pos_pp[0]
                    pp.pool = current_pool
                    pp.amplicon_number = counter
                    pp.msa_index = msa_index
                    pools[pp.pool].append(pp)

                    last_primer_pair = pp
                    current_pool = (current_pool + 1) % 2
                    counter += 1
                else:  # If the current pool is not empty
                    for pp in pos_pp:
                        if all_inter_checker(
                            pp.all_seqs(), all_seqs_in_pool, thermo_cfg
                        ):
                            continue  # If interaction skip
                        else:
                            pp.pool = current_pool
                            pp.amplicon_number = counter
                            pp.msa_index = msa_index
                            pools[pp.pool].append(pp)

                            last_primer_pair = pp
                            current_pool += 1
                            current_pool = current_pool % 2
                            counter += 1
                            new_primerpair = True
                            break

                    if not new_primerpair:  # If all pp in the window have interactions
                        pos_pp = [
                            x
                            for x in primerpairs_in_msa
                            if x.fprimer.end
                            > last_primer_pair.rprimer.start - cfg["min_overlap"]
                        ]
                        # If there are no more primerpairs left
                        if not pos_pp:
                            keep_looping = False
                            break
                        pos_pp.sort(key=lambda pp: (pp.fprimer.end, -pp.rprimer.start))

                        for pp in pos_pp:
                            if all_inter_checker(
                                pp.all_seqs(), all_seqs_in_pool, thermo_cfg
                            ):
                                continue  # If interaction skip
                            else:
                                pp.pool = current_pool
                                pp.amplicon_number = counter
                                pp.msa_index = msa_index
                                pools[pp.pool].append(pp)
                                last_primer_pair = pp
                                current_pool = (current_pool + 1) % 2
                                counter += 1

                                new_primerpair = True
                                break

    all_pp = [y for sublist in (x for x in pools) for y in sublist]
    all_pp.sort(key=lambda pp: (pp.msa_index, pp.amplicon_number))

    # Prints the primer.bed to stout (testing only)
    # for pp in all_pp:
    #    st = pp.__str__("test_seq", "test_scheme")
    #    print(st, end="")

    # Create the output
    pathlib.Path.mkdir(OUTPUT_DIR, exist_ok=True)

    # Write bed file
    with open(OUTPUT_DIR / f"{cfg['output_prefix']}.primer.bed", "w") as outfile:
        for pp in all_pp:
            st = pp.__str__(
                msa_index_to_ref_name.get(pp.msa_index, "NA"),
                msa_index_to_ref_name.get(pp.msa_index, "NA"),
            )
            outfile.write(st)

    # Write the config file, combining the cfg and thermo_cfg
    with open(OUTPUT_DIR / f"config.json", "w") as outfile:
        comb = thermo_cfg | cfg
        outfile.write(json.dumps(comb, sort_keys=True))


if __name__ == "__main__":
    main()
