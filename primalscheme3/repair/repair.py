import hashlib
import json
import pathlib
import shutil
from enum import Enum

from click import UsageError
from primalschemers._core import do_pool_interact  # type: ignore

# Core imports
from primalscheme3.core.bedfiles import read_bedlines_to_bedprimerpairs
from primalscheme3.core.config import Config
from primalscheme3.core.digestion import (
    DIGESTION_ERROR,
    DIGESTION_RESULT,
    f_digest_to_result,
    r_digest_to_result,
)
from primalscheme3.core.logger import setup_rich_logger
from primalscheme3.core.msa import MSA
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.core.seq_functions import reverse_complement
from primalscheme3.core.thermo import THERMO_RESULT


class NewPrimerStatus(Enum):
    VALID = "valid"
    PRESENT = "present"
    FAILED = "failed"


class SeqStatus:
    seq: str | None
    count: int
    thermo_status: THERMO_RESULT | DIGESTION_ERROR

    def __init__(
        self,
        seq: str | None,
        count: int,
        thermo_status: THERMO_RESULT | DIGESTION_ERROR,
    ):
        self.seq = seq
        self.count = count
        self.thermo_status = thermo_status

    def __str__(self) -> str:
        return f"{self.seq}\t{self.count}\t{self.thermo_status}"


def detect_early_return(seq_counts: list[DIGESTION_RESULT]) -> bool:
    """
    Checks for an early return condition, will return True condition is met
    """
    # Check for early return conditions
    for dr in seq_counts:
        if dr.count == -1:
            return True
    return False


def report_check(
    seqstatus: DIGESTION_RESULT,
    current_primer_seqs: set[str],
    seqs_bytes_in_pools: list[list[str]],
    pool: int,
    dimerscore: float,
    logger,
    config: Config,
) -> bool:
    """
    Will carry out the checks and report the results via the logger. Will return False if the seq should not be added
    """

    report_seq = seqstatus.seq if isinstance(seqstatus.seq, str) else "DIGESTION_ERROR"
    report_seq = report_seq.rjust(config.primer_size_max + 5, " ")

    # Check it passed thermo
    if (
        seqstatus.thermo_check(config=config) != THERMO_RESULT.PASS
        or seqstatus.seq is None
    ):
        logger.warning(
            f"{report_seq}\t{round(seqstatus.count, 4)}\t[red]{NewPrimerStatus.FAILED.value}[/red]: {seqstatus.thermo_check(config=config)}",
        )
        return False

    # Check it is a new seq
    if seqstatus.seq in current_primer_seqs:
        logger.info(
            f"{report_seq}\t{round(seqstatus.count, 4)}\t[blue]{NewPrimerStatus.PRESENT.value}[/blue]: In scheme",
        )
        return False

    # Check for minor allele
    if seqstatus.count < config.min_base_freq:
        logger.warning(
            f"{report_seq}\t{round(seqstatus.count, 4)}\t[red]{NewPrimerStatus.FAILED.value}[/red]: Minor allele",
        )
        return False

    # Check for dimer with pool
    if do_pool_interact(
        [seqstatus.seq.encode()],  # type: ignore
        seqs_bytes_in_pools[pool],
        dimerscore,
    ):
        logger.warning(
            f"{report_seq}\t{round(seqstatus.count, 4)}\t[red]{NewPrimerStatus.FAILED.value}[/red]: Interaction with pool",
        )
        return False

    # Log the seq
    logger.info(
        f"{report_seq}\t{round(seqstatus.count, 4)}\t[green]{NewPrimerStatus.VALID.value}[/green]: Can be added",
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
        config=config,
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
    all_primerpairs, _header = read_bedlines_to_bedprimerpairs(bedfile_path)

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
    seqs_bytes_in_pools = [[] for _ in range(config.n_pools)]
    for pp in primerpairs_in_msa:
        seqs_bytes_in_pools[pp.pool].extend(
            [*pp.fprimer.seq_bytes(), *pp.rprimer.seq_bytes()]
        )

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

        _end_col, fseq_counts = f_digest_to_result(
            msa_obj.array, config, msa_fkmer_end, config.min_base_freq
        )

        # Change count to freq
        fseq_total_count = sum([dr.count for dr in fseq_counts])
        for dr in fseq_counts:
            dr.count = dr.count / fseq_total_count

        # Check for early return conditions
        if detect_early_return(fseq_counts):
            logger.warning(
                f"Early return for {pp.amplicon_prefix}_{pp.amplicon_number}_LEFT. Skipping",
            )
            continue

        seqstatuss = sorted(fseq_counts, key=lambda x: x.count, reverse=True)

        # Decide if the new seqs should be added
        for seqstatus in seqstatuss:
            if not report_check(
                seqstatus=seqstatus,
                current_primer_seqs=pp.fprimer.seqs_bytes(),
                seqs_bytes_in_pools=seqs_bytes_in_pools,
                pool=pp.pool,
                dimerscore=config.dimer_score,
                logger=logger,
                config=config,
            ):
                continue

            # Add the new seq
            seqs_bytes_in_pools[pp.pool].append(seqstatus.seq.encode())  # type: ignore

        # Handle the right primer
        logger.info(
            f"Checking {pp.amplicon_prefix}_{pp.amplicon_number}_RIGHT",
        )
        msa_rkmer_start = mapping_list.index(pp.rprimer.start)
        _start_col, rseq_counts = r_digest_to_result(
            msa_obj.array, config, msa_rkmer_start, config.min_base_freq
        )
        # Check for early return conditions
        if detect_early_return(rseq_counts):
            logger.warning(
                "Early return for {pp.amplicon_prefix}_{pp.amplicon_number}_RIGHT",
            )
            continue
        # Valid seqs

        rseq_total_count = sum([dr.count for dr in rseq_counts])
        for dr in rseq_counts:
            dr.count = dr.count / rseq_total_count

        for dr in rseq_counts:
            if dr.seq is not None and not isinstance(dr.seq, DIGESTION_ERROR):
                dr.seq = reverse_complement(dr.seq)

        rseq_counts = sorted(rseq_counts, key=lambda x: x.count, reverse=True)

        # Decide if the new seqs should be added
        for rseqstatus in rseq_counts:
            if not report_check(
                seqstatus=rseqstatus,
                current_primer_seqs=pp.rprimer.seq_bytes(),
                seqs_bytes_in_pools=seqs_bytes_in_pools,
                pool=pp.pool,
                dimerscore=config.dimer_score,
                logger=logger,
                config=config,
            ):
                continue

            # Add the new seq
            seqs_bytes_in_pools[pp.pool].append(rseqstatus.seq.encode())  # type: ignore

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
