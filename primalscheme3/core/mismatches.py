# Does this work?
try:
    import dbm.ndbm as db
except ImportError:
    import dbm.dumb as db

from typing import Iterable

import dnaio

# Module imports
from primalscheme3.core.seq_functions import expand_ambs, reverse_complement

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

    def __init__(self, path, msas_paths: list[str], kmer_size: int) -> None:
        self.db = db.open(str(path), "n")

        # Read in and digest each MSA
        for msa_index, msa_path in enumerate(msas_paths):
            with dnaio.open(msa_path) as file:
                # For each sequence in the MSA
                for entry in file:
                    if entry.sequence is not None:
                        self._digest_kmers_into_db(
                            entry.sequence.upper(), kmer_size, msa_index
                        )

    def _write_unique(self, sequence, match: tuple[int, int]):
        """This will only write unqiue values to the db"""
        match_bstr = f"{match[0]};{match[1]}".encode()

        if db_bstr := self.db.get(sequence):
            bmatches = {*db_bstr.split(b"*")}

            # If the new match bstr is not already in db
            if match_bstr not in bmatches:
                self.db[sequence] = db_bstr + b"*" + match_bstr

        else:
            self.db[sequence] = match_bstr

    def read_matches(self, sequence) -> list[list]:
        parsed_matches = []
        if db_string := self.db.get(sequence):
            matches = [*db_string.split(b"*")]

            for match in matches:
                d = match.split(b";")
                parsed_matches.append([int(d[0]), int(d[1])])

        return parsed_matches

    def find_match(self, sequence) -> list[list]:
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
        fkmer,
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
        kmer_seqs = {x[-kmersize:] for x in fkmer.seqs()}
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
        rkmer,
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
        kmer_seqs = {x[:kmersize] for x in rkmer.seqs()}
        matches = self.find_matches(kmer_seqs, fuzzy)

        # Filter out expected matches
        if remove_expected:
            return {
                match
                for match in matches
                if match[1] != rkmer.start and match[0] == msaindex
            }
        else:
            return matches

    def keys(self):
        return self.db.keys()

    def get(self, sequence, default=None):
        return self.db.get(sequence, default)

    def _digest_kmers_into_db(self, seq: str, kmer_size: int, msa_index):
        """
        Digests a sequence into kmers and adds them to the db

        :param seq: The sequence to digest.
        :param kmer_size: The Kmer size for digestion.
        :param kmersize: The size of the kmer.
        :param msa_index: The index of the MSA.
        """

        for i in range(len(seq) + 1 - kmer_size):
            # Prevent a kmer starting on invalid base
            if seq[i] in {"", "-"}:
                continue

            kmer = "".join(seq[i : i + kmer_size]).replace("-", "")

            if len(kmer) != kmer_size:
                counter = 0

                # Keep walking right until the kmer is the correct size or walks out of index
                while counter + kmer_size + i < len(seq) and len(kmer) < kmer_size:
                    new_base = seq[i + kmer_size + counter]

                    # If the new base in valid
                    if new_base != "-" and new_base != "":
                        kmer += new_base

                    counter += 1

            # Guard Check Kmer is correct size, and doesn't contain N
            if len(kmer) != kmer_size or "N" in kmer:
                continue
            # Expand any ambiguous bases
            if exp_kmers := expand_ambs([kmer]):
                # Write each sequence into the db
                for exp_kmer in exp_kmers:
                    self._write_unique(exp_kmer, (msa_index, i))
            else:
                continue


def generate_single_mismatches(base_seq: str) -> set[str]:
    """
    Generates a set of sequences with all single-base mismatches for the given base sequence.

    :param base_seq: The input base sequence.
    :return: A set containing base sequences with single-base mismatches
    """
    return_seqs = set([base_seq])
    base_seq_list = [x for x in base_seq]  # Split the seq into bases
    for mut_index, base in enumerate(base_seq_list):
        # handle invalid bases
        if alt_bases := MUTATIONS.get(base):
            for alt_base in alt_bases:
                return_seqs.add(
                    "".join(
                        base_seq_list[0:mut_index]
                        + [alt_base]
                        + base_seq_list[mut_index + 1 :]
                    )
                )
        else:
            raise ValueError(f"Invalid base '{base}' in sequence '{base_seq_list}'")
    return return_seqs


def detect_new_products(
    new_matches: set[tuple],
    old_matches: set[tuple],
    product_size: int = 2000,
) -> bool:
    """
    Detects if adding the new matches will result in interactions with the old matches.

    :param new_matches: A set of new matched sequences represented as tuples containing
                        the match position, MSA index, and orientation indicator.
    :param old_matches: A set of old matched sequences used for comparison.
    :param product_size: The maximum allowable product size.
    :return: True if interactions between new and old matches are detected, False otherwise.
    """
    # Split the new matches in forward and reverse
    fmatches = set()
    rmatches = set()
    for newmatch in new_matches:
        if newmatch[2] == "+":
            fmatches.add(newmatch)
        elif newmatch[2] == "-":
            rmatches.add(newmatch)

    # Split the old matches in forward and reverse
    old_fmatches = set()
    old_rmatches = set()
    for oldmatch in old_matches:
        if oldmatch[2] == "+":
            old_fmatches.add(oldmatch)
        elif oldmatch[2] == "-":
            old_rmatches.add(oldmatch)

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


def detect_products(matches: set[tuple[int, int, str]], product_size=2000) -> bool:
    """
    Detect the presence of potential product formations based on matched sequences.

    :param matches: A set of matched sequences represented as tuples containing the match position,
                    MSA index, and orientation indicator.
    :param product_size: The maximum allowable product size.
    :return: True if potential product formations are detected, False otherwise.
    """

    # TODO write a hash function that colides when a product is formed.
    ## Will make O(N) rather than than O(N^2)

    # Split the mew matches in forward and reverse
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
            if fmatch[0] == rmatch[0]:
                # If within product distance
                if 0 < rmatch[1] - fmatch[1] < product_size:
                    return True

    return False


## It need to find the reverse complement of the primer in the forward sequence.
# kmerrc = "GTACGTCGATAG"
# kmer = "CTATCGACGTAC"
# "ACGATCGACTATCGACGTACGACATCGGACAGCAGATGTCGTACGTGATAGCTGCATGGTACGTCGATAG"
# "TGCTAGCTGATAGCTGCATGCTGTAGCCTGTCGTCTACAGCATGCACTATCGACGTACCATGCAGCTATC"
