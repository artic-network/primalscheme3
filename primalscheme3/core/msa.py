import pathlib
from uuid import uuid4

import dnaio
import numpy as np
from primalschemers._core import (
    Digester,  # type: ignore
    FKmer,  # type: ignore
    RKmer,  # type: ignore
)

from primalscheme3.core.classes import PrimerPair
from primalscheme3.core.config import IUPAC_ALL_ALLOWED_DNA, Config, MappingType
from primalscheme3.core.digestion import digest, generate_valid_primerpairs
from primalscheme3.core.errors import (
    MSAFileInvalid,
    MSAFileInvalidBase,
    MSAFileInvalidLength,
)
from primalscheme3.core.mapping import create_mapping, ref_index_to_msa
from primalscheme3.core.seq_functions import remove_end_insertion
from primalscheme3.core.thermo import forms_hairpin


def parse_msa(msa_path: pathlib.Path) -> tuple[np.ndarray, dict]:
    """
    Parses a multiple sequence alignment (MSA) file in FASTA format.

    This function reads an MSA file, validates its format and content, and returns a numpy array of the sequences
    and a dictionary with additional information. It checks for sequences of different lengths, empty columns,
    and non-DNA characters. It also removes end insertions from the sequences.

    Args:
        msa_path (pathlib.Path): The path to the MSA file to be parsed.

    Returns:
        tuple: A tuple containing two elements:
            - np.ndarray: A 2D numpy array where each row represents a sequence in the MSA and each column represents a position in the alignment.
            - dict: A dictionary with additional information about the MSA (currently not implemented, returns an empty dict).

    Raises:
        MSAFileInvalidLength: If the MSA contains sequences of different lengths.
        MSAFileInvalid: If the MSA file is empty or not in FASTA format.
        ValueError: If the MSA contains empty columns.
        MSAFileInvalidBase: If the MSA contains non-DNA characters.
    """
    try:
        records_index = {}
        with dnaio.open(msa_path) as input_msa:
            for record in input_msa:
                parsed_record = record.id.replace("/", "_")
                # Deal with duplicates
                if parsed_record in records_index:
                    raise ValueError(f"Duplicate ID ({parsed_record})")
                records_index[parsed_record] = record.sequence.upper()

    except (ValueError, dnaio.exceptions.FastaFormatError) as e:
        raise MSAFileInvalid(f"{msa_path.name}: {e}") from e

    try:
        array = np.array(
            [list(record) for record in records_index.values()],
            dtype="U1",
            ndmin=2,  # Enforce 2D array even if one genome
        )
    except ValueError as e:
        raise MSAFileInvalidLength(
            f"MSA ({msa_path.name}): contains sequences of different lengths. Please ensure MSA is aligned!"
        ) from e

    # Check for empty MSA, caused by no records being parsed
    if array.size == 0:
        raise MSAFileInvalid(
            f"No sequences in MSA ({msa_path.name}). Please ensure the MSA uses .fasta format."
        )

    empty_set = {"", "-"}

    empty_col_indexes = []
    # Check for empty columns and non DNA characters
    for col_index in range(0, array.shape[1]):
        slice: set[str] = set(array[:, col_index])
        # Check for empty columns
        if slice.issubset(empty_set):
            empty_col_indexes.append(col_index)
        # Check for non DNA characters
        if slice.difference(IUPAC_ALL_ALLOWED_DNA):
            base_str = ", ".join(slice.difference(IUPAC_ALL_ALLOWED_DNA))
            raise MSAFileInvalidBase(
                f"MSA ({msa_path.name}) contains non DNA characters ({base_str}) at column: {col_index}"
            )
    # Remove empty columns
    # array = np.delete(array, empty_col_indexes, axis=1)

    # Remove end insertions
    array = remove_end_insertion(array)

    return array, records_index


class MSA:
    # Provided
    name: str
    path: str
    msa_index: int

    # Calculated on init
    array: np.ndarray
    _uuid: str
    _chrom_name: str  # only used in the primer.bed file and html report
    _mapping_array: np.ndarray
    _ref_to_msa: dict[int, int]
    _seq_dict: dict

    # Calculated on evaluation
    fkmers: list[FKmer]
    rkmers: list[RKmer]
    primerpairs: list[PrimerPair]

    # Rust Nonsense
    _digester: Digester

    def __init__(
        self,
        name: str,
        path: pathlib.Path,
        msa_index: int,
        mapping: str,
        progress_manager,
        config: Config,
        logger=None,
    ) -> None:
        self.name = name
        self.path = str(path)
        self.msa_index = msa_index
        self.logger = logger
        self.progress_manager = progress_manager

        # Add empty lists for the primerpairs, fkmers and rkmers
        self.primerpairs = []
        self.fkmers = []
        self.rkmers = []

        # digester
        self._digester = Digester(
            msa_path=str(path.absolute()),
            ncores=config.ncores,
            remap=config.mapping == MappingType.FIRST,
        )

        # Read in the MSA
        try:
            self.array, self._seq_dict = parse_msa(path)
        except Exception as e:
            # Log the error and raise it
            if self.logger:
                self.logger.error(f"MSA: {self.name} failed QC: {e}")
            raise e

        # Create the mapping array
        # Goes from msa idx -> ref idx
        if mapping == "consensus":
            self._chrom_name = self.name + "_consensus"
            self._mapping_array = np.array([*range(len(self.array[0]))])
        elif mapping == "first":
            self._mapping_array, self.array = create_mapping(self.array, 0)
            self._chrom_name = list(self._seq_dict)[0]
        else:
            raise ValueError(f"Mapping method: {mapping} not recognised")

        # Goes from ref idx -> msa idx
        self._ref_to_msa = ref_index_to_msa(self._mapping_array)

        # Assign a UUID
        self._uuid = str(uuid4())[:8]

        if "/" in self._chrom_name:
            new_chromname = self._chrom_name.replace("/", "_")
            warning_str = (
                f"Replacing '/' with '_'. '{self._chrom_name}' -> '{new_chromname}'"
            )
            if self.logger:
                self.logger.warning(warning_str)
            else:
                print(warning_str)
            self._chrom_name = new_chromname

        # Check length
        if len(self._chrom_name) > 200:  # limit is 255
            new_chromname = self._chrom_name[:200]
            if self.logger:
                self.logger.warning(
                    f"Chromname '{self._chrom_name}' is too long, "
                    f"limit is 100 characters. Truncating to '{new_chromname}'"
                )
            self._chrom_name = new_chromname

    def digest_rs(
        self,
        config: Config,
        indexes: tuple[list[int], list[int]] | None = None,  # type: ignore
    ) -> None:
        if indexes is None:
            indexes: tuple[None, None] = (None, None)

        self.fkmers, self.rkmers, logs = self._digester.digest(
            findexes=indexes[0],
            rindexes=indexes[1],
            primer_len_min=config.primer_size_min,
            primer_len_max=config.primer_size_max,
            primer_gc_max=config.primer_gc_max / 100,
            primer_gc_min=config.primer_gc_min / 100,
            primer_tm_max=config.primer_tm_max,
            primer_tm_min=config.primer_tm_min,
            primer_annealing_prop=config.primer_annealing_prop,  # if None will use TM
            annealing_temp_c=config.primer_annealing_tempc,
            max_walk=config.primer_max_walk,
            max_homopolymers=config.primer_homopolymer_max,
            min_freq=config.min_base_freq,
            ignore_n=config.ignore_n,
            dimerscore=config.dimer_score,
        )

        # Log
        if self.logger:
            for s in logs:
                self.logger.debug(s)

            self.logger.info("Starting Primer Hairpin Check")

        ## Hairpin check the kmers
        self.fkmers = [x for x in self.fkmers if not forms_hairpin(x.seqs(), config)]
        self.rkmers = [x for x in self.rkmers if not forms_hairpin(x.seqs(), config)]

    #    def digest_f_to_count(self, config: Config, findexes: list[int] | None = None):
    #        counts = self._digester.digest_f_to_count(
    #            findexes=findexes,
    #            primer_len_min=config.primer_size_min,
    #            primer_len_max=config.primer_size_max,
    #            primer_gc_max=config.primer_gc_max / 100,
    #            primer_gc_min=config.primer_gc_min / 100,
    #            primer_tm_max=config.primer_tm_max,
    #            primer_tm_min=config.primer_tm_min,
    #            primer_annealing_prop=config.primer_annealing_prop,  # if None will use TM
    #            annealing_temp_c=config.primer_annealing_tempc,
    #            max_walk=config.primer_max_walk,
    #            max_homopolymers=config.primer_homopolymer_max,
    #            min_freq=config.min_base_freq,
    #            ignore_n=config.ignore_n,
    #            dimerscore=config.dimer_score,
    #        )
    #        return counts

    def digest(
        self,
        config: Config,
        indexes: tuple[list[int], list[int]] | None = None,
    ) -> None:
        """
        Digest the given MSA array and return the FKmers and RKmers.

        :param cfg: A dictionary containing configuration parameters.
        :param indexes: A tuple of MSA indexes for (FKmers, RKmers), or False to use all indexes.
        :return: None (Class is updated inplace)
        """
        # Create all the kmers
        self.fkmers, self.rkmers = digest(
            msa_array=self.array,
            config=config,
            indexes=indexes,
            logger=self.logger,
            progress_manager=self.progress_manager,
            chrom=self.name,
        )
        # remap the fkmer and rkmers if needed
        if self._mapping_array is not None:
            mapping_set = set(self._mapping_array)

            remaped_fkmers = [fkmer.remap(self._mapping_array) for fkmer in self.fkmers]  # type: ignore
            self.fkmers = [
                x
                for x in remaped_fkmers
                if x is not None and x.end in mapping_set and min(x.starts()) >= 0
            ]
            remaped_rkmers = [rkmer.remap(self._mapping_array) for rkmer in self.rkmers]  # type: ignore
            self.rkmers = [
                x
                for x in remaped_rkmers
                if x is not None
                and x.start in mapping_set
                and max(x.ends()) < self.array.shape[1]
            ]

    def generate_primerpairs(
        self, amplicon_size_min: int, amplicon_size_max: int, dimerscore: float
    ) -> None:
        self.primerpairs = generate_valid_primerpairs(
            fkmers=self.fkmers,
            rkmers=self.rkmers,
            amplicon_size_min=amplicon_size_min,
            amplicon_size_max=amplicon_size_max,
            dimerscore=dimerscore,
            msa_index=self.msa_index,
            progress_manager=self.progress_manager,
            chrom=self.name,
        )
        # Update primerpairs to include the chrom_name and amplicon_prefix
        for primerpair in self.primerpairs:
            primerpair.chrom_name = self._chrom_name
            primerpair.amplicon_prefix = self._uuid

    def write_msa_to_file(self, path: pathlib.Path):
        # Write all the consensus sequences to a single file
        with dnaio.FastaWriter(path, line_length=60) as msa_outfile:
            for id, seq in self._seq_dict.items():
                msa_outfile.write(dnaio.SequenceRecord(name=id, sequence=seq))
