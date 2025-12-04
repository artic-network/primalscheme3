import pathlib
from collections import Counter
from enum import Enum

from click import UsageError

# Interaction checker
from primalscheme3.core.bedfiles import (
    read_bedlines_to_bedprimerpairs,
)
from primalscheme3.core.config import Config
from primalscheme3.core.create_report_data import (
    generate_all_plotdata,
)
from primalscheme3.core.create_reports import generate_all_plots_html
from primalscheme3.core.logger import setup_rich_logger
from primalscheme3.core.mismatches import MatchDB
from primalscheme3.core.msa import MSA
from primalscheme3.core.multiplex import Multiplex, PrimerPairCheck
from primalscheme3.core.primer_visual import (
    primer_mismatch_heatmap,
)
from primalscheme3.core.progress_tracker import ProgressManager


class ReplaceRunModes(Enum):
    AddBest = "addbest"
    ListAll = "listall"


def replace(
    config: Config,
    primerbed: pathlib.Path,
    primername: str,
    msapath: pathlib.Path,
    pm: ProgressManager | None,
    output: pathlib.Path,
    force: bool,
    mode: ReplaceRunModes,
    mask_old_sites: bool = False,
):
    """
    List all replacements primers
    """
    offline_plots = False

    # See if the output dir already exists.
    if output.is_dir() and not force:
        raise UsageError(f"{output} already exists, please use --force to override")

    # Create the output dir and a work subdir
    pathlib.Path.mkdir(output, exist_ok=True)
    pathlib.Path.mkdir(output / "work", exist_ok=True)

    if pm is None:
        pm = ProgressManager()

    logger = setup_rich_logger(logfile=str((output / "work/file.log").absolute()))

    # Update the amplicon size
    logger.info(
        f"Updated min/max amplicon size to {config.amplicon_size_min}/{config.amplicon_size_max}"
    )

    # Read in the bedfile
    bedprimerpairs, headers = read_bedlines_to_bedprimerpairs(primerbed)
    chroms = {bpp.chrom_name for bpp in bedprimerpairs}

    # Read in the MSA
    MSA_INDEX = 0
    msa = MSA(
        name=msapath.stem,
        path=msapath,
        msa_index=MSA_INDEX,
        logger=logger,
        progress_manager=pm,
        config=config,
    )

    # Assign all primers with matching chroms to the msa_index
    correct_msa_bpp = 0
    for bpp in bedprimerpairs:
        if bpp.chrom_name == msa._chrom_name:
            correct_msa_bpp += 1
            bpp.msa_index = MSA_INDEX
        # else default is -1

    if correct_msa_bpp == 0:
        logger.critical(
            f"Chroms found in bedfile ({', '.join([c for c in chroms if c is not None])}) do not match MSA chromname ({msa._chrom_name})."
        )
        exit()

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
        left_ol_ref_index = max(0, wanted_pp.rprimer.start)

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
        right_ol_ref_index = max(0, wanted_pp.fprimer.end)

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

    # See if to hard mask old primer sites
    unwanted_f_ends = {wanted_pp.fprimer.end}
    unwanted_r_starts = {wanted_pp.rprimer.start}

    if mask_old_sites:
        unwanted_f_ends = unwanted_f_ends.union(
            range(min(wanted_pp.fprimer.starts()), wanted_pp.fprimer.end)
        )
        unwanted_r_starts = unwanted_f_ends.union(
            range(wanted_pp.rprimer.start, max(wanted_pp.rprimer.ends()))
        )

    # Remove the primers from old scheme
    msa.fkmers = [fk for fk in msa.fkmers if fk.end not in unwanted_f_ends]
    msa.rkmers = [rk for rk in msa.rkmers if rk.start not in unwanted_r_starts]

    logger.info(f"Digested into {len(msa.fkmers)} FKmers and {len(msa.rkmers)} RKmers")

    # Generate all primerpairs then interaction check
    msa.generate_primerpairs(
        amplicon_size_max=config.amplicon_size_max,
        amplicon_size_min=config.amplicon_size_min,
        dimerscore=config.dimer_score,
    )

    ## TODO ENSURE THE PP SPAN REQUIRED REGIONS IN REGION MODE

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
        return

    # Sort the pp depending on a score
    valid_pp.sort(key=lambda pp: len(pp.all_seqs()))

    # Write all valid primers to file
    if mode == ReplaceRunModes.ListAll:
        with open(output / "primer.bed", "w") as outfile:
            for i, pp in enumerate(valid_pp, 1):
                outfile.write(f"# PrimerPair Option: {i}\n")
                outfile.write(pp.to_bed())
        exit()

    multiplex.add_primer_pair_to_pool(
        valid_pp[0], valid_pp[0].pool, valid_pp[0].msa_index
    )

    logger.info(f"Added following primerpair to scheme. \n{valid_pp[0].to_bed()}")
    # Write primer bed file
    with open(output / "primer.bed", "w") as outfile:
        primer_bed_str = multiplex.to_bed()
        outfile.write(primer_bed_str)

    # Writing plot data
    plot_data = generate_all_plotdata(
        list(msa_dict.values()),
        output / "work",
        last_pp_added=multiplex._last_pp_added,
    )

    # Write the plot
    with open(output / "plot.html", "w") as outfile:
        outfile.write(
            generate_all_plots_html(plot_data, output, offline_plots=offline_plots)
        )
    with open(output / "primer.html", "w") as outfile:
        for i, msa_obj in enumerate(msa_dict.values()):
            try:
                outfile.write(
                    primer_mismatch_heatmap(
                        array=msa_obj.array,
                        seqdict=msa_obj._seq_dict,
                        bedfile=output / "primer.bed",
                        offline_plots=True if offline_plots and i == 0 else False,
                        mapping=config.mapping,
                    )
                )
            except UsageError:
                logger.warning(
                    f"No Primers found for {msa_obj._chrom_name}. Skipping Plot!"
                )
