import random
import numpy as np
from Bio import SeqIO
from typing import Iterable
from itertools import product

import dbm.ndbm

# Module imports
from primal_digest.classes import FKmer, RKmer, AMBIGUOUS_DNA_COMPLEMENT
from primal_digest.config import AMB_BASES, ALL_DNA

MUTATIONS = {
    "A": "CGT",
    "C": "AGT",
    "G": "CAT",
    "T": "CGA",
}


class MatchDB:
    """
    Match encoding
    delim between matches:  b'*'
    delim between values:   b';'
    """

    def __init__(self, path) -> None:
        self.db = dbm.ndbm.open(path, "n")

    def write_unique(self, sequence, match: tuple[int, int]):
        """This will only write unqiue values to the db"""
        match_bstr = f"{match[0]};{match[1]}".encode()

        if db_bstr := self.db.get(sequence):
            bmatches = {x for x in db_bstr.split(b"*")}

            # If the new match bstr is not already in db
            if match_bstr not in bmatches:
                self.db[sequence] = db_bstr + b"*" + match_bstr

        else:
            self.db[sequence] = match_bstr

    def read_matches(self, sequence) -> list[list]:
        parsed_matches = []
        if db_string := self.db.get(sequence):
            matches = [bmatch for bmatch in db_string.split(b"*")]

            for match in matches:
                d = match.split(b";")
                parsed_matches.append([int(d[0]), int(d[1])])

        return parsed_matches

    def find_match(self, sequence) -> set[tuple]:
        """
        Find all matches for a single sequence.

        :param sequence: The sequence to search for.
        :return: A list of matches, where each match is a list containing sequence details and orientation indicator.
        """
        matches = []

        # If the sequence is found
        if db_f_matches := self.read_matches(sequence):
            for match in db_f_matches:
                matches.append(match + ["+"])

        if db_r_matches := self.read_matches(reverse_complement(sequence)):
            for match in db_r_matches:
                matches.append(match + ["-"])

        return matches

    def find_matches(self, seqs: Iterable[str], fuzzy: bool = False) -> set[tuple]:
        """
        Find all matches for a collection of sequences.

        :param seqs: An iterable of sequences to search for.
        :param fuzzy: If True, consider fuzzy matches with single mismatches.
        :return: A set of matches, each represented as a tuple.
        """
        if fuzzy:
            search_seqs = {
                fseq
                for fseq in (generate_single_mismatches(seq) for seq in seqs)
                for fseq in fseq
            }
        else:
            search_seqs = seqs

        # Find all matches
        matches = {
            tuple(m) for m in (self.find_match(st) for st in search_seqs) for m in m
        }
        return matches

    def find_fkmer(
        self,
        fkmer: FKmer,
        fuzzy: bool,
        kmersize: int,
        msaindex: int,
        remove_expected: bool,
    ):
        """
        Find matches for the given FKmer.

        :param fkmer: The FKmer to search for.
        :param fuzzy: If True, consider fuzzy matches with single mismatches.
        :param kmersize: The size of the kmer.
        :param msaindex: The index of the MSA.
        :param remove_expected: If True, remove expected matches.
        :return: A set of unexpected matches for the FKmer.
        """
        kmer_seqs = {x[-kmersize:] for x in fkmer.seqs}
        matches = self.find_matches(kmer_seqs, fuzzy)

        # Filter out expected matches
        if remove_expected:
            return {
                match
                for match in matches
                if match[1] != fkmer.end - kmersize and match[0] == msaindex
            }
        else:
            return matches

    def find_rkmer(
        self,
        rkmer: FKmer,
        fuzzy: bool,
        kmersize: int,
        msaindex: int,
        remove_expected: bool,
    ):
        """
        Find unexpected matches for the given RKmer.

        :param rkmer: The RKmer to search for.
        :param fuzzy: If True, consider fuzzy matches with single mismatches.
        :param kmersize: The size of the kmer.
        :param msaindex: The index of the MSA.
        :param remove_expected: If True, remove expected matches.
        :return: A set of unexpected matches for the RKmer.
        """
        kmer_seqs = {x[:kmersize] for x in rkmer.seqs}
        matches = self.find_matches(kmer_seqs, fuzzy)

        # Filter out expected matches
        if remove_expected:
            return {
                match
                for match in matches
                if match[1] != rkmer.end and match[0] == msaindex
            }
        else:
            return matches

    def keys(self):
        return self.db.keys()

    def get(self, sequence, default=None):
        return self.db.get(sequence, default)


def hamming_dist(seq1, seq2) -> int:
    return sum([x != y for x, y in zip(seq1, seq2)])


def expand_ambs(seqs: Iterable[str]) -> set[set]:
    """
    Takes a list / set of strings and returns a set with all ambs expanded
    It can return an empty set
    """
    returned_seq = set()

    for seq in seqs:
        # If there is any amb_bases in the seq
        if {base for base in seq} & AMB_BASES:
            expanded_seqs = set(map("".join, product(*map(ALL_DNA.get, seq))))
            for exp_seq in expanded_seqs:
                returned_seq.add(exp_seq)
        else:
            returned_seq.add(seq)
    return returned_seq


def generate_single_mismatches(base_seq: str) -> set[str]:
    return_seqs = set([base_seq])
    base_seq = [x for x in base_seq]
    for mut_index, base in enumerate(base_seq):
        for alt_base in MUTATIONS.get(base):
            return_seqs.add(
                "".join(base_seq[0:mut_index] + [alt_base] + base_seq[mut_index + 1 :])
            )

    return return_seqs


def reverse_complement(kmer_seq: str):
    rev_seq = kmer_seq[::-1]
    return "".join(AMBIGUOUS_DNA_COMPLEMENT[base.upper()] for base in rev_seq)


def random_dna(size, seed) -> str:
    random.seed(seed)
    return "".join(random.choice("CGTA") for _ in range(size))


def digest_kmers(
    seq: str, kmer_size: int, msa_index
) -> dict[str : set[tuple[int, int]]]:
    """Return is (msa_index, start_index)"""

    return_dict = {}
    for i in range(len(seq) + 1 - kmer_size):
        # Prevent a kmer starting on invalid base
        if seq[i] in {"", "-"}:
            continue

        kmer = "".join(seq[i : i + kmer_size]).replace("-", "")

        if len(kmer) != kmer_size:
            counter = 1
            # Keep walking right until the kmer is the correct size
            while counter + kmer_size + i < len(seq) and len(kmer) < kmer_size:
                kmer += seq[i + kmer_size + counter]
                kmer.replace("-", "")
                counter += 1

        current_result = return_dict.get(kmer, set())
        current_result.add((msa_index, i))
        return_dict[kmer] = current_result

    return return_dict


def digest_kmers_db(
    db: MatchDB, seq: str, kmer_size: int, msa_index
) -> dict[str : set[tuple[int, int]]]:
    """Return is (msa_index, start_index)"""

    for i in range(len(seq) + 1 - kmer_size):
        # Prevent a kmer starting on invalid base
        if seq[i] in {"", "-"}:
            continue

        kmer = "".join(seq[i : i + kmer_size]).replace("-", "")

        if len(kmer) != kmer_size:
            counter = 1

            # Keep walking right until the kmer is the correct size or walks out of index
            while counter + kmer_size + i < len(seq) and len(kmer) < kmer_size:
                new_base = seq[i + kmer_size + counter]

                # If the new base in valid
                if new_base != "-" and new_base != "":
                    kmer += new_base

                counter += 1

        # Check Kmer is correct size, and doesn't contain N
        if len(kmer) == kmer_size and "N" not in kmer:
            db.write_unique(kmer, (msa_index, i))


def digest_msa(msa_index: int, msa_array: np.ndarray, kmer_size=20) -> dict[str:list]:
    return_dict: dict[str : list[int]] = {}

    for row_index in range(msa_array.shape[0]):
        row_dict = digest_kmers(
            msa_array[row_index],
            kmer_size=kmer_size,
            msa_index=msa_index,
        )
        # merge the dicts
        for k, v in row_dict.items():
            if k in return_dict:
                return_dict[k] = return_dict[k] | v
            else:
                return_dict[k] = v

    return return_dict


def detect_new_products(
    new_matches: set[tuple],
    old_fmatches: set[tuple],
    old_rmatches: set[tuple],
    product_size: int = 2000,
) -> bool:
    """
    Detects if adding the new matches will interact with the old matches
    """
    # Split the new matches in forward and reverse
    fmatches = set()
    rmatches = set()
    for newmatch in new_matches:
        if newmatch[2] == "+":
            fmatches.add(newmatch)
        elif newmatch[2] == "-":
            rmatches.add(newmatch)

    # Check the new fmatches against the old rmatches
    for fmatch in fmatches:
        for old_rmatch in old_rmatches:
            # If from same msa
            if fmatch[0] == old_rmatch[0]:
                # If within product distance
                if 0 < old_rmatch[1] - fmatch[1] < product_size:
                    return True

    # Check the new rmatches against the old fmatches
    for rmatch in rmatches:
        for old_fmatch in old_fmatches:
            # If from same msa
            if rmatch[0] == old_fmatch[0]:
                # If within product distance
                if 0 < rmatch[1] - old_fmatch[1] < product_size:
                    return True

    return False


def detect_products(
    matches: set[tuple[tuple[int, int], str]],
    product_size=2000,
) -> bool:
    """
    False means no product
    matches should be union of all matches (matches1 | matches2)
    """

    # TODO write a hash function that colides when a product is formed.
    ## Will make O(N) rather than than O(N^2)

    # Split the matches in forward and reverse
    fmatches = set()
    rmatches = set()
    for match in matches:
        if match[2] == "+":
            fmatches.add(match)
        elif match[2] == "-":
            rmatches.add(match)

    # If all matches are in the same direction
    if not fmatches or not rmatches:
        return False

    for fmatch in fmatches:
        for rmatch in rmatches:
            # If from same msa
            if fmatch[0][0] == rmatch[0][0]:
                # If within product distance
                if 0 < rmatch[0][1] - fmatch[0][1] < product_size:
                    return True

    return False


def main():
    ARG_MSA = [
        "/Users/kentcg/primal-pannel/nibsc_viral_mutliplex/rhinovirus_a39--185907.fasta",
        "/Users/kentcg/primal-pannel/nibsc_viral_mutliplex/human_astrovirus--1868658.fasta",
        "/Users/kentcg/primal-pannel/nibsc_viral_mutliplex/human_respirovirus_1--12730.fasta",
        "/Users/kentcg/primal-pannel/nibsc_viral_mutliplex/human_orthorubulavirus_2--2560525.fasta",
        "/Users/kentcg/primal-pannel/nibsc_viral_mutliplex/RSVA_208893_final.fasta",
        "/Users/kentcg/primal-pannel/nibsc_viral_mutliplex/human_adenovirus_2--10515.fasta",
        "/Users/kentcg/primal-pannel/nibsc_viral_mutliplex/human_orthorubulavirus_4--2560526.fasta",
        "/Users/kentcg/primal-pannel/nibsc_viral_mutliplex/human_respirovirus_3--11216.fasta",
        "/Users/kentcg/primal-pannel/nibsc_viral_mutliplex/align_coxsackie_virus_b4--12073.fasta",
        "/Users/kentcg/primal-pannel/nibsc_viral_mutliplex/human_metapneumovirus--162145.fasta",
        "/Users/kentcg/primal-pannel/nibsc_viral_mutliplex/align_human_parechovirus--195055.fasta",
        "/Users/kentcg/primal-pannel/nibsc_viral_mutliplex/human_adenovirus_41--10524.fasta",
        "/Users/kentcg/primal-pannel/nibsc_viral_mutliplex/align_coronavirus-229E.fasta",
        "/Users/kentcg/primal-pannel/nibsc_viral_mutliplex/align_human_alphaherpesvirus_2.fasta",
        "/Users/kentcg/primal-pannel/nibsc_viral_mutliplex/align_human_gammaherpesvirus_4.fasta",
        "/Users/kentcg/primal-pannel/nibsc_viral_mutliplex/human_alphaherpesvirus_1.fasta",
        "/Users/kentcg/primal-pannel/nibsc_viral_mutliplex/align_human_betaherpesvirus_5.fasta",
        "/Users/kentcg/primal-pannel/nibsc_viral_mutliplex/align_human_alphaherpesvirus_3.fasta",
    ]
    matchdb = MatchDB("GRCh38_test")

    # Read in and digest each MSA
    for msa_index, msa_path in enumerate(ARG_MSA):
        records = SeqIO.parse(msa_path, "fasta")
        # For each sequence in the MSA
        for record in records:
            print(record.id)
            digest_kmers_db(matchdb, str(record.seq.upper()), 20, msa_index)

    fkmer = FKmer(
        200,
        seqs={
            "AAAACGGCAACTGCTGATATGCTC",
            "CGGTAACTGCTGATATGCTC",
            "CGGTAACTGCTGATATGCTC",
            "CGGCAACTGCTGATATGCTC",
            "CGGCAACTGCTGATATGCTC",
            "CGGTAACTGCTGATATGCTC",
        },
    )
    Rkmer = RKmer(
        200,
        seqs={
            "AAAACGGCAACTGCTGATATGCTC",
            "CGGTAACTGCTGATATGCTC",
            "CGGTAACTGCTGATATGCTC",
            "CGGCAACTGCTGATATGCTC",
            "CGGCAACTGCTGATATGCTC",
            "CGGTAACTGCTGATATGCTC",
        },
    )

    # Match Kmers
    print(
        matchdb.find_fkmer(
            fkmer=fkmer, fuzzy=True, remove_expected=False, msaindex=0, kmersize=20
        )
    )

    print(detect_new_products({(0, 100, "+")}, {(0, 100, "+")}, {(0, 200, "-")}))


if __name__ == "__main__":
    main()


## It need to find the reverse complement of the primer in the forward sequence.
kmerrc = "GTACGTCGATAG"
kmer = "CTATCGACGTAC"
"ACGATCGACTATCGACGTACGACATCGGACAGCAGATGTCGTACGTGATAGCTGCATGGTACGTCGATAG"
"TGCTAGCTGATAGCTGCATGCTGTAGCCTGTCGTCTACAGCATGCACTATCGACGTACCATGCAGCTATC"
