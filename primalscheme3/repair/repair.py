import hashlib
import json
import pathlib
from enum import Enum

from primalschemers import FKmer, RKmer, do_pool_interact  # type: ignore

# Core imports
from primalscheme3.core.bedfiles import read_bedlines_to_bedprimerpairs
from primalscheme3.core.config import Config
from primalscheme3.core.errors import UsageError
from primalscheme3.core.logger import setup_rich_logger
from primalscheme3.core.msa import MSA, digest_kwargs
from primalscheme3.core.progress_tracker import ProgressManager

# Tags on a return_all digest's statuses and errors that are not failures
LOW_FREQ = "LowFreq"
IGNORED = "Ignored"


def required_n_pools(primerpairs, min_pools: int) -> int:
    """
    Number of pool buckets needed to index every primerpair by its (0-indexed)
    ``pool`` attribute.

    Sized by the highest pool index present, not the count of distinct pools, so
    sparse / non-contiguous pool numbering (e.g. pools {0, 2}) does not under-size
    the list and cause an IndexError when indexing ``pools[pp.pool]``.
    """
    max_pool = max((pp.pool for pp in primerpairs), default=-1)
    return max(min_pools, max_pool + 1)


class NewPrimerStatus(Enum):
    VALID = "valid"
    PRESENT = "present"
    FAILED = "failed"


def covered_by(seq: str, primer_seqs) -> bool:
    """Whether an existing primer has the same 3' end as ``seq``: equal, or
    one is a 3'-suffix of the other. Primers are 5'->3' with the 3' end last,
    for LEFT and RIGHT alike."""
    return any(seq.endswith(p) or p.endswith(seq) for p in primer_seqs)


def log_seq(
    level,
    report_seq: str,
    count: float,
    status: NewPrimerStatus,
    reason: str,
    colour: str,
):
    level(
        f"{report_seq}\t{round(count, 4)}\t[{colour}]{status.value}[/{colour}]: {reason}"
    )


def report_check(
    seq: str,
    count: float,
    statuses: list[str],
    current_primer_seqs: set[str],
    seqs_bytes_in_pools: list[list[bytes]],
    pool: int,
    dimerscore: float,
    logger,
    config: Config,
) -> bool:
    """
    Report whether ``seq`` can be added to the primer, and why not. Returns
    True if it should be added. ``statuses`` are the digest's verdicts on the
    sequence (failed checks, plus the LowFreq tag).
    """
    report_seq = seq.rjust(config.primer_size_max + 5, " ")
    failures = [s for s in statuses if s != LOW_FREQ]

    # Already covered: an existing primer shares this 3' end (one is a
    # 3'-suffix of the other), i.e. the same site walked to another length
    if covered_by(seq, current_primer_seqs):
        log_seq(
            logger.info, report_seq, count, NewPrimerStatus.PRESENT, "In scheme", "blue"
        )
        return False

    # Every check the sequence fails
    if failures:
        log_seq(
            logger.warning,
            report_seq,
            count,
            NewPrimerStatus.FAILED,
            ", ".join(failures),
            "red",
        )
        return False

    # Below min_base_freq
    if LOW_FREQ in statuses:
        log_seq(
            logger.warning,
            report_seq,
            count,
            NewPrimerStatus.FAILED,
            "Minor allele",
            "red",
        )
        return False

    # Dimer with the primers in its pool
    if do_pool_interact([seq.encode()], seqs_bytes_in_pools[pool], dimerscore):
        log_seq(
            logger.warning,
            report_seq,
            count,
            NewPrimerStatus.FAILED,
            "Interaction with pool",
            "red",
        )
        return False

    log_seq(
        logger.info, report_seq, count, NewPrimerStatus.VALID, "Can be added", "green"
    )
    return True


def check_side(
    pp,
    side: str,
    msa_pos: int | None,
    kmer_at,
    existing,
    kmer_cls,
    seqs_bytes_in_pools: list[list[bytes]],
    config: Config,
    logger,
):
    """Report every sequence at one primer site; return the primer with
    the valid new sequences added."""
    name = f"{pp.amplicon_prefix}_{pp.amplicon_number}_{side}"
    logger.info(f"Checking {name}")
    kmer = kmer_at.get(msa_pos) if msa_pos is not None else None
    if kmer is None:
        logger.critical(f"{name}: no sequences at this site in the MSA")
        return existing

    # Genomes with no usable sequence here, and why
    for tags, count in kmer.errors():
        logger.warning(f"{name}: {', '.join(tags)} in {count:g} genomes")
    for status in kmer.kmer_statuses():
        logger.warning(f"{name}: {status} between the sequences found here")

    current = set(existing.seqs())
    new_seqs, new_counts = [], []
    statuses = kmer.statuses() or [[] for _ in kmer.seqs()]
    for seq, count, tags in zip(kmer.seqs(), kmer.counts(), statuses, strict=True):
        if report_check(
            seq=seq,
            count=count,
            statuses=tags,
            current_primer_seqs=current,
            seqs_bytes_in_pools=seqs_bytes_in_pools,
            pool=pp.pool,
            dimerscore=config.dimer_score,
            logger=logger,
            config=config,
        ):
            new_seqs.append(seq)
            new_counts.append(count)
            seqs_bytes_in_pools[pp.pool].append(seq.encode())
    if not new_seqs:
        return existing

    # Add the new sequences to the primer
    old_counts = existing.counts()
    counts = None if old_counts is None else [*old_counts, *new_counts]
    pos = existing.end if kmer_cls is FKmer else existing.start
    return kmer_cls(
        [*existing.seqs_bytes(), *(s.encode() for s in new_seqs)], pos, counts
    )


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
    n_pools = required_n_pools(all_primerpairs, config.n_pools)
    seqs_bytes_in_pools = [[] for _ in range(n_pools)]
    for pp in primerpairs_in_msa:
        seqs_bytes_in_pools[pp.pool].extend(
            [*pp.fprimer.seqs_bytes(), *pp.rprimer.seqs_bytes()]
        )

    # Digest at the existing primer sites (MSA coordinates), returning every
    # sequence with the checks it fails, rather than only passing primers
    msa_obj.create_digester(local_msa_path, config.ncores, remap=False)

    ref_to_msa = msa_obj._ref_to_msa
    fp_sites = {
        pp.fprimer.end: ref_to_msa.get(pp.fprimer.end) for pp in primerpairs_in_msa
    }
    rp_sites = {
        pp.rprimer.start: ref_to_msa.get(pp.rprimer.start) for pp in primerpairs_in_msa
    }
    fkmers, rkmers, logs = msa_obj._digester.digest(
        findexes=sorted({i for i in fp_sites.values() if i is not None}),
        rindexes=sorted({i for i in rp_sites.values() if i is not None}),
        thermo_check=True,
        return_all=True,
        **digest_kwargs(config),
    )
    for line in logs:
        logger.debug(line)
    fkmer_at = {fk.end: fk for fk in fkmers}
    rkmer_at = {rk.start: rk for rk in rkmers}

    for pp in primerpairs_in_msa:
        pp.fprimer = check_side(
            pp,
            "LEFT",
            fp_sites[pp.fprimer.end],
            fkmer_at,
            pp.fprimer,
            FKmer,
            seqs_bytes_in_pools,
            config,
            logger,
        )
        pp.rprimer = check_side(
            pp,
            "RIGHT",
            rp_sites[pp.rprimer.start],
            rkmer_at,
            pp.rprimer,
            RKmer,
            seqs_bytes_in_pools,
            config,
            logger,
        )

    # Write out the new bedfile
    with open(OUTPUT_DIR / "primer.bed", "w") as f:
        for pp in primerpairs_in_msa:
            pp.amplicon_prefix = msa_obj._uuid
            f.write(pp.to_bed() + "\n")

    # Amplicon and primertrimmed files should not have changed. Can be copied from the input dir
    # Not sure how to handle the amplicon names, as the primerstem has changed?
    ## Keep original names for now

    # Write the config, with the new MSA's entry
    with open(OUTPUT_DIR / "config.json", "w") as outfile:
        outfile.write(json.dumps(base_cfg, sort_keys=True))
