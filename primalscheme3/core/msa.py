import os
import pathlib
import re
from uuid import uuid4

import dnaio
import numpy as np
import tqdm
from primalbedtools.bedfiles import CHROM_REGEX
from primalschemers import (
    Digester,
    FKmer,
    RKmer,
    do_pool_interact,
)

from primalscheme3.core.classes import PrimerPair
from primalscheme3.core.config import IUPAC_ALL_ALLOWED_DNA, Config, MappingType
from primalscheme3.core.downsample import (
    _ds_pool_init,
    _ds_pool_worker,
    delta_ct_downsample_seqs,
    downsample_cloud,
    forms_dimer,
)
from primalscheme3.core.errors import (
    MSAFileInvalid,
    MSAFileInvalidBase,
    MSAFileInvalidLength,
)
from primalscheme3.core.get_window import get_r_window_FAST2
from primalscheme3.core.mapping import create_mapping, ref_index_to_msa
from primalscheme3.core.primerscorer import PrimerScorer
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.core.seq_functions import remove_end_insertion
from primalscheme3.core.thermo import forms_hairpin

_INVALID_CHROM_CHARS = re.compile(r"[^a-zA-Z0-9_.]")


def parse_chrom_name(old_chrom: str) -> str:
    """
    Sanitize a raw FASTA/genome id into a valid bed chrom name, replacing every
    character not allowed by primalbedtools' CHROM_REGEX (e.g. '|', '/', '-', spaces)
    with '_'. Applied to every genome id at MSA parse time, so anything read back out
    of ``_seq_dict`` (representative-subset headers, write_msa_to_file, etc.) is
    already a valid chrom name.
    """
    return _INVALID_CHROM_CHARS.sub("_", old_chrom)


def parse_msa(msa_path: pathlib.Path) -> tuple[np.ndarray, dict[str, str]]:
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
        records_index: dict[str, str] = {}
        with dnaio.open(msa_path) as input_msa:
            for record in input_msa:
                parsed_record_id = parse_chrom_name(record.id)
                # Deal with duplicates
                if parsed_record_id in records_index:
                    raise ValueError(f"Duplicate ID ({parsed_record_id})")
                records_index[parsed_record_id] = record.sequence.upper()

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
            if not do_pool_interact(fkmer.seqs_bytes(), rkmer.seqs_bytes(), dimerscore):
                checked_pp.append(PrimerPair(fkmer, rkmer, msa_index))

        # Update the count
        pt.manual_update(count=len(checked_pp))

    checked_pp.sort(key=lambda pp: (pp.fprimer.end, -pp.rprimer.start))
    return checked_pp


def digest_kwargs(config: Config) -> dict:
    """The primer design and reaction settings Digester.digest takes from a
    Config: everything except the indexes, thermo_check and return_all."""
    return {
        "primer_len_min": config.primer_size_min,
        "primer_len_max": config.primer_size_max,
        "primer_gc_max": config.primer_gc_max / 100,
        "primer_gc_min": config.primer_gc_min / 100,
        "primer_tm_max": config.primer_tm_max,
        "primer_tm_min": config.primer_tm_min,
        "primer_annealing_prop": config.primer_annealing_prop,  # None: use Tm
        "primer_annealing_diff": config.primer_annealing_diff,
        "annealing_temp_c": config.primer_annealing_tempc,
        "max_walk": config.primer_max_walk,
        "max_homopolymers": config.primer_homopolymer_max,
        "min_freq": config.min_base_freq,
        "ignore_n": config.ignore_n,
        "dimerscore": config.dimer_score,
        "primer_hairpin_th_max": config.primer_hairpin_th_max,
        "mv_conc": config.mv_conc,
        "dv_conc": config.dv_conc,
        "dntp_conc": config.dntp_conc,
        "dna_conc": config.dna_conc,
        "formamide_conc": config.formamide_conc,
    }


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
    _seq_dict: dict[str, str]

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
        self.create_digester(
            path=path,
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

        # Create the mapping arrays + set chrom name
        self.set_reference_genome(config.mapping)

        # Assign a UUID
        self._uuid = str(uuid4())[:8]

        # Check chromname
        if not re.match(CHROM_REGEX, self._chrom_name):
            raise ValueError(
                f"chrom must match '{CHROM_REGEX}'. Got (`{self._chrom_name}`)"
            )

        # Check length
        if len(self._chrom_name) > 200:  # limit is 255
            new_chromname = self._chrom_name[:200]
            if self.logger:
                self.logger.warning(
                    f"Chromname '{self._chrom_name}' is too long, "
                    f"limit is 100 characters. Truncating to '{new_chromname}'"
                )
            self._chrom_name = new_chromname

    def create_digester(self, path: pathlib.Path, ncores: int, remap: bool):
        """
        Creates an instance of the Rust Digester Object.
        """
        # digester
        self._digester = Digester(
            msa_path=str(path.absolute()),
            ncores=ncores,
            remap=remap,
        )

    def digest_rs(
        self,
        config: Config,
        indexes: tuple[list[int], list[int]] | None = None,  # type: ignore
        rs_thermo=True,
        py_hairpin=True,
    ) -> None:
        """
        Digest the MSA using the Rust Digester and optionally downsample and filter kmers.

        This method performs the following steps:
        1. Uses the Rust Digester to generate FKmers and RKmers from the MSA, with configurable parameters.
        2. Optionally disables Rust-side thermo checking if downsampling is enabled.
        3. If downsampling is enabled, uses Python logic to select and filter kmers in parallel (using multiprocessing),
           based on sequence abundance and thermodynamic properties.
        4. Optionally applies a Python-side hairpin filter to remove kmers that form hairpins.
        5. Updates the object's FKmers and RKmers in place.

        Args:
            config (Config): Configuration object with all primer design parameters.
            indexes (tuple[list[int], list[int]] | None): Optional tuple of indexes for FKmers and RKmers to digest.
            rs_thermo (bool): If True, enables the Rust thermo checks (including hairpins). Always on when downsampling, which digests with return_all and uses the checks' per-sequence statuses.
            py_hairpin (bool): If True and rs_thermo is False, applies the Python hairpin filter.

        Returns:
            None. Updates self.fkmers and self.rkmers in place.

        Raises:
            ValueError: If chrom name does not match expected regex or is too long.
        """
        if indexes is None:
            indexes: tuple[None, None] = (None, None)

        # Downsampling needs every primer in each cloud, failing ones too, with
        # the Rust checks' verdict on each: return_all gives them with statuses
        if config.downsample:
            rs_thermo = True

        self.fkmers, self.rkmers, logs = self._digester.digest(
            findexes=indexes[0],
            rindexes=indexes[1],
            thermo_check=rs_thermo,
            return_all=config.downsample,
            **digest_kwargs(config),
        )
        # Log
        if self.logger:
            for s in logs:
                self.logger.debug(s)

        # Start the downsample
        if config.downsample:
            model_path = (
                os.getenv("PRIMALSCHEME3_PRIMER_SCORER_MODEL")
                or config.primer_scorer_model
            )
            # Load once here to fail early, and record exactly which model ran
            loaded = PrimerScorer(model_path)
            config.primer_scorer_model = str(model_path)
            config.primer_scorer_model_sha256 = loaded.sha256
            # Optionally dump a downsample graph SVG per (non-trivial) cloud plus a
            # downsample_summary.csv (pre -> post counts + coverage prop). Gated by
            # config.downsample_plot as it is slow and verbose.
            plot_dir = None
            if config.downsample_plot:
                plot_dir = config.output / "work" / "downsample_plots"
                os.makedirs(plot_dir, exist_ok=True)
            plot_outdir = str(plot_dir) if plot_dir is not None else None

            # One cloud per kmer: the sequences a normal digest would consider,
            # their counts, which fail the Rust checks, and dimer pairs
            f_clouds = [(k, downsample_cloud(k)) for k in self.fkmers]
            r_clouds = [(k, downsample_cloud(k)) for k in self.rkmers]
            f_clouds = [(k, c) for k, c in f_clouds if c is not None]
            r_clouds = [(k, c) for k, c in r_clouds if c is not None]

            def rebuild(cloud, new_seqs, kmer_cls, pos):
                if new_seqs is None:
                    return None
                # The selection avoids dimers; this guards against a regression
                assert not forms_dimer(cloud, new_seqs), "downsampled cloud dimerises"
                counts = dict(zip(cloud.seqs, cloud.counts, strict=True))
                return kmer_cls(
                    [s.encode() for s in new_seqs], pos, [counts[s] for s in new_seqs]
                )

            # Plotting writes SVGs via matplotlib and does not parallelise cleanly,
            # so the plot/debug path stays serial. ncores <= 1 also runs serially.
            if config.downsample_plot or config.ncores <= 1:

                def run(clouds, label, prefix, pos_of):
                    kept = []
                    for kmer, cloud in tqdm.tqdm(clouds, desc=label):
                        new_seqs = delta_ct_downsample_seqs(
                            seq_counts=dict(zip(cloud.seqs, cloud.counts, strict=True)),
                            config=config,
                            ps=loaded,
                            visualise=config.downsample_plot and len(cloud.seqs) > 1,
                            outdir=plot_outdir,
                            name=f"{self._chrom_name}_{prefix}_{pos_of(kmer)}",
                            failed_seqs=cloud.failed,
                            dimers=cloud.dimers,
                        )
                        kept.append(new_seqs)
                    return kept

                f_res = run(f_clouds, "downsampling FPrimers", "F", lambda k: k.end)
                r_res = run(r_clouds, "downsampling RPrimers", "R", lambda k: k.start)
            else:
                from multiprocessing import Pool

                f_tasks = [
                    (i, c.seqs, c.counts, c.failed, c.dimers)
                    for i, (_, c) in enumerate(f_clouds)
                ]
                r_tasks = [
                    (i, c.seqs, c.counts, c.failed, c.dimers)
                    for i, (_, c) in enumerate(r_clouds)
                ]
                chunk = max(1, (len(f_tasks) + len(r_tasks)) // (config.ncores * 4))
                with Pool(
                    config.ncores,
                    initializer=_ds_pool_init,
                    initargs=(model_path, config),
                ) as pool:
                    f_out = list(
                        tqdm.tqdm(
                            pool.imap_unordered(
                                _ds_pool_worker, f_tasks, chunksize=chunk
                            ),
                            total=len(f_tasks),
                            desc="downsampling FPrimers",
                        )
                    )
                    r_out = list(
                        tqdm.tqdm(
                            pool.imap_unordered(
                                _ds_pool_worker, r_tasks, chunksize=chunk
                            ),
                            total=len(r_tasks),
                            desc="downsampling RPrimers",
                        )
                    )
                f_res = [kept for _, kept in sorted(f_out, key=lambda x: x[0])]
                r_res = [kept for _, kept in sorted(r_out, key=lambda x: x[0])]

            new_fkmers = [
                rebuild(c, kept, FKmer, k.end)
                for (k, c), kept in zip(f_clouds, f_res, strict=True)
            ]
            new_rkmers = [
                rebuild(c, kept, RKmer, k.start)
                for (k, c), kept in zip(r_clouds, r_res, strict=True)
            ]
            self.fkmers = [k for k in new_fkmers if k is not None]
            self.rkmers = [k for k in new_rkmers if k is not None]

        # The Rust thermo check includes the hairpin check (and downsampling
        # uses it), so the Python filter is only for unchecked digests
        if py_hairpin and not rs_thermo:
            ## Hairpin check the kmers
            if self.logger:
                self.logger.info("Starting Primer Hairpin Check")
            self.fkmers = [
                x for x in self.fkmers if not forms_hairpin(x.seqs(), config)
            ]
            self.rkmers = [
                x for x in self.rkmers if not forms_hairpin(x.seqs(), config)
            ]

    def set_reference_genome(self, mapping: MappingType | int):
        """
        Given the mapping type, or the index changes the MSA to use that genome as the reference.
        - Updates/create the mapping arrays
        - Updates _chromname
        """

        # Create the mapping array
        # Goes from msa idx -> ref idx
        if mapping == MappingType.CONSENSUS:
            self._chrom_name = self.name + "_consensus"
            self._mapping_array = np.array([*range(len(self.array[0]))])
        elif mapping == MappingType.FIRST:
            self._mapping_array, self.array = create_mapping(self.array, 0)
            self._chrom_name = list(self._seq_dict)[0]
        elif isinstance(mapping, int):
            self._mapping_array, self.array = create_mapping(self.array, mapping)
            self._chrom_name = list(self._seq_dict)[mapping]
        else:
            raise ValueError(f"Mapping method: {mapping} not recognised")

        # Goes from ref idx -> msa idx
        self._ref_to_msa = ref_index_to_msa(self._mapping_array)

        # Parse the chrom name
        new_chromname = parse_chrom_name(self._chrom_name)
        if new_chromname != self._chrom_name:
            warning_str = (
                f"Sanitizing chrom name to match '{CHROM_REGEX}'. "
                f"'{self._chrom_name}' -> '{new_chromname}'"
            )
            if self.logger:
                self.logger.warning(warning_str)
            else:
                print(warning_str)
            self._chrom_name = new_chromname

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
                msa_outfile.write(
                    dnaio.SequenceRecord(name=parse_chrom_name(id), sequence=seq)
                )
