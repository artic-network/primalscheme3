import hashlib
import json
import pathlib
from enum import Enum

from click import UsageError
from primalschemers import do_pool_interact  # type: ignore

# Core imports
from primalscheme3.core.bedfiles import read_bedlines_to_bedprimerpairs
from primalscheme3.core.config import Config
from primalscheme3.core.digestion import (
    DIGESTION_ERROR,
    DIGESTION_RESULT,
)
from primalscheme3.core.logger import setup_rich_logger
from primalscheme3.core.msa import MSA
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.core.thermo import THERMO_RESULT, thermo_check


def primers_equal(s1, s2) -> bool:
    for b1, b2 in zip(s1, s2, strict=False):
        if b2 != b1:
            return False
    return True


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
    seq: str,
    count: float,
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

    report_seq = seq if isinstance(seq, str) else "DIGESTION_ERROR"
    report_seq = report_seq.rjust(config.primer_size_max + 5, " ")

    thermo_status = thermo_check(seq, config)

    # Check it passed thermo
    if thermo_status != THERMO_RESULT.PASS or seq is None:
        logger.warning(
            f"{report_seq}\t{round(count, 4)}\t[red]{NewPrimerStatus.FAILED.value}[/red]: {thermo_status.name}",
        )
        return False

    # Check it is a new seq
    if any([primers_equal(seq, x) for x in current_primer_seqs]):
        logger.info(
            f"{report_seq}\t{round(count, 4)}\t[blue]{NewPrimerStatus.PRESENT.value}[/blue]: In scheme",
        )
        return False

    # Check for minor allele
    if count < config.min_base_freq:
        logger.warning(
            f"{report_seq}\t{round(count, 4)}\t[red]{NewPrimerStatus.FAILED.value}[/red]: Minor allele",
        )
        return False

    # Check for dimer with pool
    if do_pool_interact(
        [seq.encode()],  # type: ignore
        seqs_bytes_in_pools[pool],
        dimerscore,
    ):
        logger.warning(
            f"{report_seq}\t{round(count, 4)}\t[red]{NewPrimerStatus.FAILED.value}[/red]: Interaction with pool",
        )
        return False

    # Log the seq
    logger.info(
        f"{report_seq}\t{round(count, 4)}\t[green]{NewPrimerStatus.VALID.value}[/green]: Can be added",
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

    config.min_base_freq = 0.1

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

    # Read in the bedfile, and find chrom names
    all_primerpairs, _header = read_bedlines_to_bedprimerpairs(bedfile_path)
    bedfile_chroms: set[str] = {pp.chrom_name for pp in all_primerpairs}  # type: ignore

    if None in bedfile_chroms:
        logger.warning("Invalid bedfile. Chrom name has been parsed to None")
        exit()

    # Read in the MSA file
    msa_obj = MSA(
        name=msa_path.stem,
        path=msa_path,
        msa_index=0,
        logger=logger,
        progress_manager=pm,
        config=config,
    )

    # Search for the required chroms in msa._seq_dict
    msa_seq_names_to_index = {id: i for i, id in enumerate(msa_obj._seq_dict.keys())}
    names_in_both = bedfile_chroms.intersection(msa_seq_names_to_index.keys())
    if len(names_in_both) == 0:
        logger.warning(
            f"BedFile Chrom names ({', '.join(bedfile_chroms)}) not found in MSA seq IDs"
        )
        exit()
    elif len(names_in_both) > 1:
        logger.warning(
            f"Multiple BedFile Chrom names ({', '.join(names_in_both)}) found in MSA seq IDs"
        )
        exit()

    # Change the MSA chrom to the correct ref
    msa_obj.set_reference_genome(msa_seq_names_to_index[names_in_both.pop()])

    logger.info(
        f"Read in MSA: [blue]{msa_path.name}[/blue] ({msa_obj._chrom_name})\t"
        f"seqs:[green]{msa_obj.array.shape[0]}[/green]\t"
        f"cols:[green]{msa_obj.array.shape[1]}[/green]"
    )

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
    msa_obj.write_msa_to_file(local_msa_path)

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
            [*pp.fprimer.seqs_bytes(), *pp.rprimer.seqs_bytes()]
        )

    # Update the MSA Digester to use the MSA index
    msa_obj.create_digester(local_msa_path, config.ncores, remap=False)

    fp_indexes = [msa_obj._ref_to_msa[pp.fprimer.end] for pp in primerpairs_in_msa]
    rp_indexes = [msa_obj._ref_to_msa[pp.rprimer.start] for pp in primerpairs_in_msa]

    # Digest the positions, returning all seqs
    msa_obj.digest_rs(
        config=config,
        indexes=(fp_indexes, rp_indexes),
        rs_thermo=False,
        py_hairpin=False,
    )

    # For primerpair in the bedfile, check if new seqs need to be added by digestion the MSA
    for pp in primerpairs_in_msa:
        logger.info(
            f"Checking {pp.amplicon_prefix}_{pp.amplicon_number}_LEFT",
        )
        msa_fkmer_end = msa_obj._ref_to_msa[pp.fprimer.end]

        if msa_fkmer_end is None:
            continue

        # Find the FKmer object
        new_fkmers = [fk for fk in msa_obj.fkmers if fk.end == msa_fkmer_end]
        if len(new_fkmers) != 1:
            logger.critical(f"Digestion failed for FKmer:{msa_fkmer_end}")
            continue
        new_fkmer = new_fkmers[0]

        fseq_count = {
            seq: count
            for seq, count in zip(new_fkmer.seqs(), new_fkmer.counts(), strict=True)
        }

        # Decide if the new seqs should be added
        for seq in fseq_count:
            if not report_check(
                seq=seq,
                count=fseq_count[seq],
                current_primer_seqs=pp.fprimer.seqs(),
                seqs_bytes_in_pools=seqs_bytes_in_pools,
                pool=pp.pool,
                dimerscore=config.dimer_score,
                logger=logger,
                config=config,
            ):
                continue

            # Add the new seq
            seqs_bytes_in_pools[pp.pool].append(seq.encode())  # type: ignore

        # Handle the right primer
        logger.info(
            f"Checking {pp.amplicon_prefix}_{pp.amplicon_number}_RIGHT",
        )
        msa_rkmer_end = msa_obj._ref_to_msa[pp.rprimer.start]

        if msa_rkmer_end is None:
            continue

        # Find the FKmer object
        new_rkmers = [rk for rk in msa_obj.rkmers if rk.start == msa_rkmer_end]
        if len(new_rkmers) != 1:
            logger.critical(f"Digestion failed for RKmer:{msa_rkmer_end}")
            continue
        new_rkmer = new_rkmers[0]

        rseq_count = {
            seq: count
            for seq, count in zip(new_rkmer.seqs(), new_rkmer.counts(), strict=True)
        }

        # Decide if the new seqs should be added
        for seq in rseq_count:
            if not report_check(
                seq=seq,
                count=rseq_count[seq],
                current_primer_seqs=pp.rprimer.seqs(),
                seqs_bytes_in_pools=seqs_bytes_in_pools,
                pool=pp.pool,
                dimerscore=config.dimer_score,
                logger=logger,
                config=config,
            ):
                continue

            # Add the new seq
            seqs_bytes_in_pools[pp.pool].append(seq.encode())  # type: ignore

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
        outfile.write(json.dumps(config.to_dict(), sort_keys=True))
