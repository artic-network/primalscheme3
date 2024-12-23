import hashlib
import json
import pathlib
import shutil
from enum import Enum

from click import UsageError
from primaldimer_py import do_pools_interact_py  # type: ignore

# Core imports
from primalscheme3.core.bedfiles import read_in_bedprimerpairs
from primalscheme3.core.config import Config
from primalscheme3.core.digestion import (
    DIGESTION_ERROR,
    f_digest_to_count,
    r_digest_to_count,
)
from primalscheme3.core.logger import setup_rich_logger
from primalscheme3.core.msa import MSA
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.core.seq_functions import reverse_complement
from primalscheme3.core.thermo import THERMORESULT, thermo_check


class NewPrimerStatus(Enum):
    VALID = "valid"
    PRESENT = "present"
    FAILED = "failed"


class SeqStatus:
    seq: str | None
    count: int
    thermo_status: THERMORESULT | DIGESTION_ERROR

    def __init__(
        self, seq: str | None, count: int, thermo_status: THERMORESULT | DIGESTION_ERROR
    ):
        self.seq = seq
        self.count = count
        self.thermo_status = thermo_status

    def __str__(self) -> str:
        return f"{self.seq}\t{self.count}\t{self.thermo_status}"


def detect_early_return(seq_counts: dict[str | DIGESTION_ERROR, int]) -> bool:
    """
    Checks for an early return condition, will return True condition is met
    """
    # Check for early return conditions
    for error, count in seq_counts.items():
        if count == -1 and type(error) == DIGESTION_ERROR:
            return True
    return False


def report_check(
    seqstatus: SeqStatus,
    current_primer_seqs: set[str],
    seqs_in_pools: list[list[str]],
    pool: int,
    dimerscore: float,
    logger,
    config: Config,
) -> bool:
    """
    Will carry out the checks and report the results via the logger. Will return False if the seq should not be added
    """

    report_seq = seqstatus.seq if seqstatus.seq is not None else "DIGESTION_ERROR"
    report_seq = report_seq.rjust(config.primer_size_max + 5, " ")

    # Check it passed thermo
    if seqstatus.thermo_status != THERMORESULT.PASS or seqstatus.seq is None:
        logger.warning(
            f"{report_seq}\t{seqstatus.count}\t[red]{NewPrimerStatus.FAILED.value}[/red]: {seqstatus.thermo_status}",
        )
        return False

    # Check it is a new seq
    if seqstatus.seq in current_primer_seqs:
        logger.info(
            f"{report_seq}\t{seqstatus.count}\t[blue]{NewPrimerStatus.PRESENT.value}[/blue]: In scheme",
        )
        return False

    # Check for minor allele
    if seqstatus.count < 0:
        logger.warning(
            f"{report_seq}\t{seqstatus.count}\t[red]{NewPrimerStatus.FAILED.value}[/red]: Minor allele",
        )
        return False

    # Check for dimer with pool
    if do_pools_interact_py([seqstatus.seq], seqs_in_pools[pool], dimerscore):
        logger.warning(
            f"{report_seq}\t{seqstatus.count}\t[red]{NewPrimerStatus.FAILED.value}[/red]: Interaction with pool",
        )
        return False

    # Log the seq
    logger.info(
        f"{report_seq}\t{seqstatus.count}\t[green]{NewPrimerStatus.VALID.value}[/green]: Can be added",
    )

    return True


def repair(
    config_path: pathlib.Path,
    msa_path: pathlib.Path,
    bedfile_path: pathlib.Path,
    output_dir: pathlib.Path,
    force: bool,
    pm: ProgressManager | None,
):
    OUTPUT_DIR = pathlib.Path(output_dir).absolute()  # Keep absolute path

    # Read in the config file
    with open(config_path) as f:
        base_cfg = json.load(f)

    msa_data = base_cfg["msa_data"]

    # Parse params from the config
    config = Config(**base_cfg)
    base_cfg = config.to_dict()
    base_cfg["msa_data"] = msa_data

    config.min_base_freq = 0.01

    # See if the output dir already exists
    if OUTPUT_DIR.is_dir() and not force:
        raise UsageError(f"{OUTPUT_DIR} already exists, please use --force to override")

    # Create the output dir and a work subdir
    pathlib.Path.mkdir(OUTPUT_DIR, exist_ok=True)
    pathlib.Path.mkdir(OUTPUT_DIR / "work", exist_ok=True)

    ## Set up the logger
    logger = setup_rich_logger(str(OUTPUT_DIR / "work" / "file.log"))

    ## Set up the progress manager
    if pm is None:
        pm = ProgressManager()

    # Read in the MSA file
    msa_obj = MSA(
        name=msa_path.stem,
        path=msa_path,
        msa_index=0,
        mapping=base_cfg["mapping"],
        logger=logger,
        progress_manager=pm,
    )
    logger.info(
        f"Read in MSA: [blue]{msa_path.name}[/blue] ({msa_obj._chrom_name})\t"
        f"seqs:[green]{msa_obj.array.shape[0]}[/green]\t"
        f"cols:[green]{msa_obj.array.shape[1]}[/green]"
    )
    # Check for a '/' in the chromname
    if "/" in msa_obj._chrom_name:
        new_chromname = msa_obj._chrom_name.split("/")[0]
        logger.warning(
            f"Having a '/' in the chromname {msa_obj._chrom_name} "
            f"will cause issues with figure generation bedfile output. "
            f"Parsing chromname [yellow]{msa_obj._chrom_name}[/yellow] -> [green]{new_chromname}[/green]"
        )
        msa_obj._chrom_name = new_chromname

    # Update the base_cfg with the new msa
    # Create MSA checksum
    with open(msa_path, "rb") as f:
        msa_checksum = hashlib.file_digest(f, "md5").hexdigest()

    current_msa_index = max([int(x) for x in base_cfg["msa_data"].keys()])
    base_cfg["msa_data"][str(current_msa_index + 1)] = {
        "msa_name": msa_obj.name,
        "msa_path": str("work/" + msa_path.name),
        "msa_chromname": msa_obj._chrom_name,
        "msa_uuid": msa_obj._uuid,
        "msa_checksum": msa_checksum,
    }
    # Copy the MSA file to the work dir
    local_msa_path = OUTPUT_DIR / "work" / msa_path.name
    shutil.copy(msa_path, local_msa_path)

    # Read in the bedfile
    all_primerpairs, _header = read_in_bedprimerpairs(bedfile_path)

    # Get the primerpairs for this new MSA
    primerpairs_in_msa = [
        pp for pp in all_primerpairs if pp.chrom_name == msa_obj._chrom_name
    ]

    if len(primerpairs_in_msa) == 0:
        logger.critical(
            f"No primerpairs found for {msa_obj._chrom_name} in {bedfile_path}",
        )
        raise UsageError(
            f"No primerpairs found for {msa_obj._chrom_name} in {bedfile_path}"
        )

    # Get all the seqs in each pool
    seqs_in_pools = [[] for _ in range(config.n_pools)]
    for pp in primerpairs_in_msa:
        seqs_in_pools[pp.pool].extend([*pp.fprimer.seqs, *pp.rprimer.seqs])

    # Find the indexes in the MSA that the primerbed refer to
    assert msa_obj._mapping_array is not None
    mapping_list = list(msa_obj._mapping_array)

    # For primerpair in the bedfile, check if new seqs need to be added by digestion the MSA
    for pp in primerpairs_in_msa:
        logger.info(
            f"Checking {pp.amplicon_prefix}_{pp.amplicon_number}_LEFT",
        )
        msa_fkmer_end = msa_obj._ref_to_msa.get(pp.fprimer.end)

        if msa_fkmer_end is None:
            continue

        _end_col, fseq_counts = f_digest_to_count(
            msa_obj.array, config, msa_fkmer_end, config.min_base_freq
        )

        # Check for early return conditions
        if detect_early_return(fseq_counts):
            logger.warning(
                f"Early return for {pp.amplicon_prefix}_{pp.amplicon_number}_LEFT. Skipping",
            )
            continue

        # Thermo check each sequence
        seqstatuss: list[SeqStatus] = []
        for seq, count in fseq_counts.items():
            if isinstance(seq, DIGESTION_ERROR):
                thermo_status = seq
                seq = None
            else:
                thermo_status = thermo_check(seq, config=config)
            seqstatuss.append(SeqStatus(seq, count, thermo_status))
        seqstatuss = sorted(seqstatuss, key=lambda x: x.count, reverse=True)

        # Decide if the new seqs should be added
        for seqstatus in seqstatuss:
            if not report_check(
                seqstatus=seqstatus,
                current_primer_seqs=pp.fprimer.seqs,
                seqs_in_pools=seqs_in_pools,
                pool=pp.pool,
                dimerscore=config.dimer_score,
                logger=logger,
                config=config,
            ):
                continue

            # Add the new seq
            seqs_in_pools[pp.pool].append(seqstatus.seq)

        # Handle the right primer
        logger.info(
            f"Checking {pp.amplicon_prefix}_{pp.amplicon_number}_RIGHT",
        )
        msa_rkmer_start = mapping_list.index(pp.rprimer.start)
        _start_col, rseq_counts = r_digest_to_count(
            msa_obj.array, config, msa_rkmer_start, config.min_base_freq
        )
        # Check for early return conditions
        if detect_early_return(rseq_counts):
            logger.warning(
                "Early return for {pp.amplicon_prefix}_{pp.amplicon_number}_RIGHT",
            )
            continue

        # Valid seqs
        valid_rseqs = {
            reverse_complement(seq): count
            for seq, count in rseq_counts.items()
            if isinstance(seq, str)
        }

        # Thermo check each sequence
        rseqstatuss: list[SeqStatus] = []
        for seq, count in valid_rseqs.items():
            if isinstance(seq, DIGESTION_ERROR):
                thermo_status = seq
                seq = None
            else:
                thermo_status = thermo_check(seq, config=config)
            rseqstatuss.append(SeqStatus(seq, count, thermo_status))
        rseqstatuss = sorted(rseqstatuss, key=lambda x: x.count, reverse=True)

        # Decide if the new seqs should be added
        for rseqstatus in rseqstatuss:
            if not report_check(
                seqstatus=rseqstatus,
                current_primer_seqs=pp.rprimer.seqs,
                seqs_in_pools=seqs_in_pools,
                pool=pp.pool,
                dimerscore=config.dimer_score,
                logger=logger,
                config=config,
            ):
                continue

            # Add the new seq
            seqs_in_pools[pp.pool].append(rseqstatus.seq)

    # Write out the new bedfile
    with open(OUTPUT_DIR / "primer.bed", "w") as f:
        for pp in primerpairs_in_msa:
            pp.amplicon_prefix = msa_obj._uuid
            f.write(pp.to_bed() + "\n")

    # Amplicon and primertrimmed files should not have changed. Can be copied from the input dir
    # Not sure how to handle the amplicon names, as the primerstem has changed?
    ## Keep original names for now

    # Write the config dict to file
    with open(OUTPUT_DIR / "config.json", "w") as outfile:
        outfile.write(json.dumps(base_cfg, sort_keys=True))
