import hashlib
import json
import pathlib
from collections import Counter
from enum import Enum

from click import UsageError

# Interaction checker
from primalscheme3.core.bedfiles import (
    read_bedlines_to_bedprimerpairs,
)
from primalscheme3.core.config import Config
from primalscheme3.core.logger import setup_rich_logger
from primalscheme3.core.mismatches import MatchDB
from primalscheme3.core.msa import MSA
from primalscheme3.core.multiplex import Multiplex, PrimerPairCheck
from primalscheme3.core.progress_tracker import ProgressManager


class ReplaceErrors(Enum):
    SamePrimer = "SamePrimer"
    ClashWithPool = "ClashWithPool"
    InteractWithPool = "InteractWithPool"


def replace(
    config_path: pathlib.Path,
    ampliconsizemax: int,
    ampliconsizemin: int,
    primerbed: pathlib.Path,
    primername: str,
    msapath: pathlib.Path,
    pm: ProgressManager | None,
    output: pathlib.Path,
    force,
):
    """
    List all replacements primers
    """
    # Read in the config file
    with open(config_path) as file:
        _cfg: dict = json.load(file)
    config = Config()
    config.assign_kwargs(**_cfg)
    config.min_overlap = 0
    config.in_memory_db = True  # Set up the db

    # See if the output dir already exists.
    if output.is_dir() and not force:
        raise UsageError(f"{output} already exists, please use --force to override")

    # Create the output dir and a work subdir
    pathlib.Path.mkdir(output, exist_ok=True)
    pathlib.Path.mkdir(output / "work", exist_ok=True)

    if pm is None:
        pm = ProgressManager()

    # Set up the logger without a file
    logger = setup_rich_logger(logfile=str((output / "work/file.log").absolute()))

    # Update the amplicon size if it is provided
    if ampliconsizemax:
        logger.info(
            f"Updated min/max amplicon size to {ampliconsizemin}/{ampliconsizemax}"
        )

        config.amplicon_size_max = ampliconsizemax
        config.amplicon_size_min = ampliconsizemin

    # If more than two pools are given throw error
    if config.n_pools > 2:
        raise UsageError("ERROR: repair is only supported with two pools")

    # Create a mapping of chromname/reference to msa_index
    msa_chrom_to_index: dict[str, int] = {
        msa_data["msa_chromname"]: msa_index
        for msa_index, msa_data in _cfg["msa_data"].items()
    }

    # Read in the bedfile
    bedprimerpairs, headers = read_bedlines_to_bedprimerpairs(primerbed)

    # Map each primer to an MSA index
    for primerpair in bedprimerpairs:
        msa_index = msa_chrom_to_index.get(str(primerpair.chrom_name), None)
        if msa_index is not None:
            primerpair.msa_index = msa_index
        elif _cfg["bedfile"]:
            # This case can happen when primers are added to the scheme via --bedfile.
            # Set the msa index to -1
            primerpair.msa_index = -1
        else:
            raise UsageError(f"ERROR: {primerpair.chrom_name} not found in MSA data")

    bedprimerpairs.sort(key=lambda x: (x.chrom_name, x.amplicon_number))

    # Extract the stem from the primername
    try:
        prefix, ampliconnumber = primername.split("_")[:2]
        primerstem = f"{ampliconnumber}_{prefix}"
    except ValueError:
        raise UsageError(
            f"ERROR: {primername} cannot be parsed using _ as delim"
        ) from None

    # Find primernumber from bedfile
    wanted_pp = None
    for pp in bedprimerpairs:
        if pp.match_primer_stem(primerstem):
            wanted_pp = pp
    if wanted_pp is None:
        raise UsageError(f"ERROR: {primername} not found in bedfile")
    else:
        logger.info("Found amplicon to replace:")
        logger.info(wanted_pp.to_bed())

    # Read in the MSAs from config["msa_data"]
    msa_data = _cfg["msa_data"].get(wanted_pp.msa_index)
    if msa_data is None:
        if wanted_pp.msa_index == -1:
            raise UsageError(f"ERROR: The Primer {primername} was added via --bedfile")
        else:
            raise UsageError(f"ERROR: MSA index {wanted_pp.msa_index} not found")

    msa = MSA(
        name=msa_data["msa_name"],
        path=msapath,
        msa_index=wanted_pp.msa_index,
        logger=None,
        progress_manager=pm,
        config=config,
    )
    # Check the hashes match
    with open(msa.path, "rb") as f:
        msa_checksum = hashlib.file_digest(f, "md5").hexdigest()
    if msa_checksum != msa_data["msa_checksum"]:
        raise UsageError("ERROR: MSA checksums do not match.")
    else:
        logger.info("MSA checksums match")

    # Create the multiplex object.
    msa_dict = {wanted_pp.msa_index: msa}
    match_db = MatchDB("", [], config)
    multiplex = Multiplex(config, match_db, msa_dict)

    # Add all primers into the multiplex
    for bpp in bedprimerpairs:
        multiplex.add_primer_pair_to_pool(bpp, bpp.pool, bpp.msa_index)

    # Remove the wanted_pp
    multiplex.remove_primerpair(wanted_pp)

    # Find any ols pp
    fp_ol = set(multiplex._lookup[wanted_pp.msa_index][:, wanted_pp.fprimer.end])
    rp_ol = set(multiplex._lookup[wanted_pp.msa_index][:, wanted_pp.rprimer.start - 1])

    # Parse the overlaps to find the coords for the primers
    left_ol_ref_index = None
    for overlapping_fp in fp_ol:
        # Skip None
        if overlapping_fp is None:
            continue
        if (
            left_ol_ref_index is None
            or overlapping_fp.rprimer.start - 1 > left_ol_ref_index
        ):
            left_ol_ref_index = wanted_pp.rprimer.start - 1

    # If no left ol set the index to something reasonable
    if left_ol_ref_index is None:
        left_ol_ref_index = max(0, wanted_pp.rprimer.start - config.amplicon_size_min)

    findexes = [
        *range(
            msa._ref_to_msa[max(0, left_ol_ref_index - (config.amplicon_size_max * 2))],
            msa._ref_to_msa[left_ol_ref_index],
        )
    ]

    right_ol_ref_index = None
    for overlapping_rpp in rp_ol:
        # Skip None
        if overlapping_rpp is None:
            continue
        if (
            right_ol_ref_index is None
            or overlapping_rpp.fprimer.end < right_ol_ref_index
        ):
            right_ol_ref_index = overlapping_rpp.fprimer.end

    # If no left ol set the index to something reasonable
    if right_ol_ref_index is None:
        right_ol_ref_index = max(0, wanted_pp.fprimer.end + config.amplicon_size_min)

    rindexes = [
        *range(
            msa._ref_to_msa[right_ol_ref_index],
            min(
                msa.array.shape[1],
                msa._ref_to_msa[right_ol_ref_index + config.amplicon_size_max * 2],
            ),
        )
    ]

    # Targeted digestion leads to a mismatch of the indexes.
    # Digest the MSA into FKmers and RKmers
    msa.digest_rs(config, (findexes, rindexes))  ## Primer are remapped at this point.
    logger.info(f"Digested into {len(msa.fkmers)} FKmers and {len(msa.rkmers)} RKmers")

    # Generate all primerpairs then interaction check
    msa.generate_primerpairs(
        amplicon_size_max=config.amplicon_size_max,
        amplicon_size_min=config.amplicon_size_min,
        dimerscore=config.dimer_score,
    )

    logger.info(f"Generated {len(msa.primerpairs)} possible amplicons")
    if len(msa.primerpairs) == 0:
        logger.critical("Failed to generate amplicons. Please increase amplicon size.")
        return

    # Throw all the possible pp at the multiplex and see what sticks
    valid_pp = []

    # Keep track of the statuses
    results_counter = Counter()
    for pp in msa.primerpairs:
        result = multiplex.check_primerpair_can_be_added(pp, wanted_pp.pool)
        results_counter.update([result])

        if result == PrimerPairCheck.OK:
            pp.amplicon_number = wanted_pp.amplicon_number
            pp.pool = wanted_pp.pool
            pp.msa_index = wanted_pp.msa_index
            valid_pp.append(pp)

    # Report the status
    for k, v in results_counter.most_common():
        logger.info(f"{k}: {v}")

    if len(valid_pp) == 0:
        logger.critical("No valid primerpairs can be added.")

    # Sort the pp depending on a score
    valid_pp.sort(key=lambda pp: len(pp.all_seqs()))

    multiplex.add_primer_pair_to_pool(
        valid_pp[0], valid_pp[0].pool, valid_pp[0].msa_index
    )

    logger.info(f"Added following primerpair to scheme. \n{valid_pp[0].to_bed()}")
    # Write primer bed file
    with open(output / "primer.bed", "w") as outfile:
        primer_bed_str = multiplex.to_bed()
        outfile.write(primer_bed_str)
    exit()
