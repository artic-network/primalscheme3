import hashlib
import json
import pathlib
from collections import Counter
from enum import Enum
from time import sleep

import dnaio
from click import UsageError

# Interaction checker
from primalschemers._core import do_pool_interact  # type: ignore

from primalscheme3.core.bedfiles import (
    read_bedlines_to_bedprimerpairs,
    read_in_extra_primers,
)
from primalscheme3.core.config import Config, MappingType
from primalscheme3.core.create_report_data import (
    generate_all_plotdata,
)
from primalscheme3.core.create_reports import generate_all_plots_html
from primalscheme3.core.errors import DigestionFailNoPrimerPairs
from primalscheme3.core.logger import setup_rich_logger
from primalscheme3.core.mapping import (
    generate_consensus,
    generate_reference,
)
from primalscheme3.core.mismatches import MatchDB
from primalscheme3.core.msa import MSA
from primalscheme3.core.primer_visual import (
    plot_primer_thermo_profile_html,
    primer_mismatch_heatmap,
)
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.scheme.classes import Scheme, SchemeReturn


class ReplaceErrors(Enum):
    SamePrimer = "SamePrimer"
    ClashWithPool = "ClashWithPool"
    InteractWithPool = "InteractWithPool"


def schemereplace(
    config_path: pathlib.Path,
    ampliconsizemax: int,
    ampliconsizemin: int,
    primerbed: pathlib.Path,
    primername: str,
    msapath: pathlib.Path,
    pm: ProgressManager | None,
    ncores: int = 1,
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

    if pm is None:
        pm = ProgressManager()

    # Set up the logger without a file
    logger = setup_rich_logger(logfile=None)

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
        mapping=config.mapping.value,
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

    # Target digest the MSA
    msa_fp_end = msa._ref_to_msa[wanted_pp.fprimer.region()[1]]
    msa_rp_start = msa._ref_to_msa[wanted_pp.fprimer.region()[0]]

    findexes = sorted(
        [
            x
            for x in range(
                msa_fp_end - config.amplicon_size_max,
                msa_fp_end + config.amplicon_size_max,
            )
            if x >= 0 and x < msa.array.shape[1]
        ]
    )
    rindexes = sorted(
        [
            x
            for x in range(
                msa_rp_start - config.amplicon_size_max,
                msa_rp_start + config.amplicon_size_max,
            )
            if x >= 0 and x < msa.array.shape[1]
        ]
    )

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

    # Find the primers on either side of the wanted primer
    if wanted_pp.amplicon_number == 0:
        left_pp = None
        cov_start = wanted_pp.rprimer.region()[0]
    else:
        left_pp = [
            x
            for x in bedprimerpairs
            if x.amplicon_number == wanted_pp.amplicon_number - 1
        ][0]
        cov_start = left_pp.rprimer.region()[0]

    if wanted_pp.amplicon_number == max([x.amplicon_number for x in msa.primerpairs]):
        right_pp = None
        cov_end = wanted_pp.fprimer.region()[1]
    else:
        right_pp = [
            x
            for x in bedprimerpairs
            if x.amplicon_number == wanted_pp.amplicon_number + 1
        ][0]
        cov_end = right_pp.fprimer.region()[1]

    # Find primerpairs that span the gap
    spanning_primerpairs = [
        x
        for x in msa.primerpairs
        if x.fprimer.region()[1] < cov_start - config.min_overlap
        and x.rprimer.region()[0] > cov_end + config.min_overlap
    ]

    if len(spanning_primerpairs) > 0:
        logger.info(f"Spanning amplicons found: {len(spanning_primerpairs)}")
    else:
        raise UsageError(
            "ERROR: No spanning primers found. Please increase amplicon size"
        )

    # Sort for number of primer pairs
    spanning_primerpairs.sort(key=lambda x: len(x.all_seqs()))

    # Find all primerpairs that in the same pool as the wanted primer
    same_pool_primerpairs = [x for x in bedprimerpairs if x.pool == wanted_pp.pool]
    seqs_bytes_in_same_pool = [
        seq for seq in (x.all_seq_bytes() for x in same_pool_primerpairs) for seq in seq
    ]

    # Find all primerpairs that in the same pool and same msa as the wanted primer
    spanning_pool_primerpairs = [
        x
        for x in msa.primerpairs
        if x.pool == wanted_pp.pool and x.msa_index == wanted_pp.msa_index
    ]

    accepted_primerpairs = []

    # Keep a counter to see what goes wrong
    result_counter = Counter()
    for pos_primerpair in spanning_primerpairs:
        # Make sure the new primerpair doesn't contain the primers we want to replace
        if (
            pos_primerpair.fprimer == wanted_pp.fprimer
            or pos_primerpair.rprimer == wanted_pp.rprimer
        ):
            result_counter.update([ReplaceErrors.SamePrimer])
            continue

        # Check for all clashes
        clash = False
        # If they overlap
        for current_pp in spanning_pool_primerpairs:
            # If clash with any skip
            if range(
                max(pos_primerpair.fprimer.region()[0], current_pp.fprimer.region()[0]),
                min(pos_primerpair.rprimer.region()[1], current_pp.rprimer.region()[1])
                + 1,
            ):
                clash = True
                break
        if clash:
            result_counter.update([ReplaceErrors.ClashWithPool])
            continue

        # If they dont overlap
        # Check for interactions
        if do_pool_interact(
            pos_primerpair.all_seq_bytes(), seqs_bytes_in_same_pool, config.dimer_score
        ):
            result_counter.update([ReplaceErrors.InteractWithPool])
            continue

        pos_primerpair.pool = wanted_pp.pool
        pos_primerpair.amplicon_number = wanted_pp.amplicon_number

        accepted_primerpairs.append(pos_primerpair)

    accepted_primerpairs.sort(key=lambda a: a.fprimer.end)

    logger.info(f"Found {len(accepted_primerpairs)} valid replacement amplicons")
    for pp_number, pp in enumerate(accepted_primerpairs, 1):
        print(f"Amplicon {pp_number}: {len(pp.all_seqs())} Primers")
        print(pp.to_bed())

    # If No valid primerpairs
    if len(accepted_primerpairs) == 0:
        logger.info("The following errors prevented primerpairs being added:")
        for k, v in result_counter.items():
            logger.info(f"\t- {k.value}:\t{v}")


def schemecreate(
    msa: list[pathlib.Path],
    output_dir: pathlib.Path,
    config: Config,
    pm: ProgressManager | None,
    force: bool = False,
    input_bedfile: pathlib.Path | None = None,
    offline_plots: bool = True,
):
    """
    Creates a scheme based on multiple sequence alignments (MSA).

    Args:
        msa (list[pathlib.Path]): A list of paths to MSA files.
        output_dir (pathlib.Path): The directory where the scheme will be saved.
        config (Config): Configuration settings for scheme creation.
        pm (ProgressManager | None): Optional. A progress manager instance for UI feedback.
        force (bool): If True, existing output directories will be overwritten.
        input_bedfile (pathlib.Path | None): Optional. A path to an input BED file for incorporating specific primer pairs.
        offline_plots (bool): If True, plots will be generated for offline use.

    Raises:
        SystemExit: If the output directory already exists and the force flag is not set.
    """
    ARG_MSA = msa
    OUTPUT_DIR = pathlib.Path(output_dir).absolute()  # Keep absolute path

    # Create the Config_dict
    cfg_dict = config.to_dict()

    # See if the output dir already exists.
    if OUTPUT_DIR.is_dir() and not force:
        raise UsageError(f"{OUTPUT_DIR} already exists, please use --force to override")

    # Create the output dir and a work subdir
    pathlib.Path.mkdir(OUTPUT_DIR, exist_ok=True)
    pathlib.Path.mkdir(OUTPUT_DIR / "work", exist_ok=True)

    # Set up the logger
    logger = setup_rich_logger(str(OUTPUT_DIR / "work/file.log"))

    if pm is None:
        pm = ProgressManager()

    # Create the mismatch db
    logger.info(
        "Creating the Mismatch Database",
    )
    mismatch_db = MatchDB(
        OUTPUT_DIR / "work/mismatch",
        [str(x) for x in ARG_MSA] if config.use_matchdb else [],
        config,
    )
    logger.info(
        f"[green]Created[/green]: "
        f"{OUTPUT_DIR.relative_to(OUTPUT_DIR.parent)}/work/mismatch.db",
    )

    # If the bedfile flag is given add the primers into the scheme
    if input_bedfile is not None:
        bedprimerpairs = read_in_extra_primers(input_bedfile, config, logger)

    # Create a dict full of msa data
    msa_data = {}
    msa_dict: dict[int, MSA] = {}

    msa_has_primerpairs_bool = {msa_index: False for msa_index, _ in enumerate(ARG_MSA)}

    for msa_index, msa_path in enumerate(ARG_MSA):
        msa_data[msa_index] = {}

        # Read in the MSA
        msa_obj = MSA(
            name=msa_path.stem,
            path=msa_path,
            msa_index=msa_index,
            mapping=config.mapping.value,
            logger=logger,
            progress_manager=pm,
            config=config,
        )

        # copy the msa into the output / work dir
        local_msa_path = OUTPUT_DIR / "work" / msa_path.name
        msa_obj.write_msa_to_file(local_msa_path)

        # Create MSA checksum
        with open(local_msa_path, "rb") as f:
            msa_data[msa_index]["msa_checksum"] = hashlib.file_digest(
                f, "md5"
            ).hexdigest()

        logger.info(
            f"Read in MSA: [blue]{msa_obj._chrom_name}[/blue]\t"
            f"seqs:[green]{msa_obj.array.shape[0]}[/green]\t"
            f"cols:[green]{msa_obj.array.shape[1]}[/green]"
        )

        # Add some msa data to the dict
        msa_data[msa_index]["msa_name"] = msa_obj.name
        msa_data[msa_index]["msa_path"] = str(
            "work/" + msa_path.name
        )  # Write local path
        msa_data[msa_index]["msa_chromname"] = msa_obj._chrom_name
        msa_data[msa_index]["msa_uuid"] = msa_obj._uuid

        # Add the msa to the dict
        msa_dict[msa_index] = msa_obj

    # Check for collisions in the MSA._chrom_names names
    if len({msa_obj._chrom_name for msa_obj in msa_dict.values()}) != len(msa_dict):
        logger.critical("Duplicate chrom names found in MSA data. Exiting.")
        raise UsageError("Duplicate chrom names found in MSA data")

    # Read in all MSAs before digestion
    for msa_index, msa_obj in msa_dict.items():
        # Digest the MSA into FKmers and RKmers
        msa_obj.digest_rs(config, None)  # Default to one core
        logger.info(
            f"[blue]{msa_obj._chrom_name}[/blue]: digested to "
            f"[green]{len(msa_obj.fkmers)}[/green] FKmers and "
            f"[green]{len(msa_obj.rkmers)}[/green] RKmers"
        )

        # Add the msa to the scheme
        msa_dict[msa_index] = msa_obj

        if len(msa_obj.fkmers) == 0 or len(msa_obj.rkmers) == 0:
            logger.critical(
                f"No valid FKmers or RKmers found for [blue]{msa_obj._chrom_name}[/blue]"
            )
            continue

        # Generate all primerpairs then interaction check
        msa_obj.generate_primerpairs(
            amplicon_size_max=config.amplicon_size_max,
            amplicon_size_min=config.amplicon_size_min,
            dimerscore=config.dimer_score,
        )
        logger.info(
            f"[blue]{msa_obj._chrom_name}[/blue]: Generated "
            f"[green]{len(msa_obj.primerpairs)}[/green] possible amplicons"
        )

        if len(msa_obj.primerpairs) == 0:
            logger.critical(
                f"No valid primers found for [blue]{msa_obj._chrom_name}[/blue]"
            )
        else:
            msa_has_primerpairs_bool[msa_index] = True

    # Add MSA data into cfg
    cfg_dict["msa_data"] = msa_data

    # If all MSAs have no primers, exit
    if not any(msa_has_primerpairs_bool.values()):
        logger.critical("No valid primers found in any MSA")
        raise DigestionFailNoPrimerPairs("No valid primerpairs found in any MSA")

    # Create the scheme object early
    scheme = Scheme(config=config, matchDB=mismatch_db, msa_dict=msa_dict)

    msa_chrom_to_index: dict[str, int] = {
        msa._chrom_name: msa_index for msa_index, msa in msa_dict.items()
    }
    # Add the bedprimerpairs into the scheme
    if input_bedfile is not None and bedprimerpairs:  # type: ignore
        # if input_bedfile != None then bedprimerpairs is assigned
        for bedpp in bedprimerpairs:
            # Map the primerpair to the msa via chromname
            bedpp.msa_index = msa_chrom_to_index.get(bedpp.chrom_name, -1)  # type: ignore
            scheme.add_primer_pair_to_pool(bedpp, bedpp.pool, bedpp.msa_index)
            logger.debug(
                f"Added {bedpp.amplicon_prefix} from [blue]{input_bedfile.name}[/blue]",
            )

    # Start the Scheme generation
    for msa_index, msa_obj in msa_dict.items():
        # Set up the pm for the MSA
        scheme_pt = pm.create_sub_progress(
            iter=None, chrom=msa_obj.name, process="Creating Scheme", leave=False
        )
        scheme_pt.manual_update(n=0, total=msa_obj.array.shape[1])

        while True:
            # Provide the coverage to the progress tracker
            coverage = scheme.get_coverage_percent(msa_index)
            scheme_pt.manual_update(count=coverage if coverage is not None else 0)
            # Update the progress tracker to the current state of the walk
            if scheme._last_pp_added:
                scheme_pt.manual_update(n=scheme._last_pp_added[-1].rprimer.region()[1])

            match scheme.try_ol_primerpairs(msa_obj.primerpairs, msa_index):
                case SchemeReturn.ADDED_OL_PRIMERPAIR:
                    last_pp_added = scheme._last_pp_added[-1]
                    logger.info(
                        "Added [green]overlapping[/green] amplicon for "
                        f"[blue]{msa_obj._chrom_name}[/blue]: {last_pp_added.fprimer.region()[0]}\t"
                        f"{last_pp_added.rprimer.region()[1]}\t{last_pp_added.pool + 1}"
                    )
                    continue
                case SchemeReturn.ADDED_FIRST_PRIMERPAIR:
                    last_pp_added = scheme._last_pp_added[-1]
                    logger.info(
                        "Added [green]first[/green] amplicon for "
                        f"[blue]{msa_obj._chrom_name}[/blue]: {last_pp_added.fprimer.region()[0]}\t"
                        f"{last_pp_added.rprimer.region()[1]}\t{last_pp_added.pool + 1}",
                    )
                    continue
                case SchemeReturn.NO_OL_PRIMERPAIR:
                    pass  # Do nothing move on to next step
                case SchemeReturn.NO_FIRST_PRIMERPAIR:
                    logger.warning(
                        f"No valid primerpairs found for [blue]{msa_obj._chrom_name}[/blue]",
                    )
                    break

            # Try to backtrack
            if config.backtrack:
                logger.info("Backtracking...")
                match scheme.try_backtrack(msa_obj.primerpairs, msa_index):
                    case SchemeReturn.ADDED_BACKTRACKED:
                        last_pp_added = scheme._last_pp_added[-1]
                        logger.info(
                            f"Backtracking allowed [green]overlapping[/green] amplicon for [blue]{msa_obj._chrom_name}[/blue]: "
                            f"{last_pp_added.fprimer.region()[0]}\t{last_pp_added.rprimer.region()[1]}\t{last_pp_added.pool + 1}",
                        )
                    case SchemeReturn.NO_BACKTRACK:
                        logger.info(
                            f"Could not backtrack for [blue]{msa_obj._chrom_name}[/blue]",
                        )
                        pass  # Do nothing move on to next step

            # Try and add a walking primer
            match scheme.try_walk_primerpair(msa_obj.primerpairs, msa_index):
                case SchemeReturn.ADDED_WALK_PRIMERPAIR:
                    last_pp_added = scheme._last_pp_added[-1]
                    logger.info(
                        f"Added [yellow]walking[/yellow] amplicon for [blue]{msa_obj._chrom_name}[/blue]: "
                        f"{last_pp_added.fprimer.region()[0]}\t{last_pp_added.rprimer.region()[1]}\t{last_pp_added.pool + 1}",
                    )
                case _:
                    break

        if config.circular:
            match scheme.try_circular(msa_obj):
                case SchemeReturn.ADDED_CIRCULAR:
                    last_pp_added = scheme._last_pp_added[-1]
                    logger.info(
                        f"Added [green]circular[/green] amplicon for [blue]{msa_obj._chrom_name}[/blue]: "
                        f"{last_pp_added.fprimer.region()[0]}\t{last_pp_added.rprimer.region()[1]}\t{last_pp_added.pool + 1}",
                    )
                case SchemeReturn.NO_CIRCULAR:
                    logger.info(
                        f"No [red]circular[/red] amplicon for [blue]{msa_obj._chrom_name}[/blue]"
                    )
        # Close the progress tracker
        scheme_pt.manual_update(n=scheme_pt.total)
        scheme_pt.close()

    # Create the progress tracker for the final steps
    logger.info("Writing output files")
    upload_steps = 8
    upload_pt = pm.create_sub_progress(
        iter=None,
        chrom="file creation",
        process="Creating output files",
        leave=False,
        bar_format="{l_bar}{bar}",
        unit="%",
    )
    upload_pt.manual_update(n=0, total=upload_steps)

    # Write primer bed file
    with open(OUTPUT_DIR / "primer.bed", "w") as outfile:
        primer_bed_str = scheme.to_bed()
        outfile.write(primer_bed_str)
    upload_pt.manual_update(n=1, update=True)

    # Write amplicon bed file
    with open(OUTPUT_DIR / "amplicon.bed", "w") as outfile:
        amp_bed_str = scheme.to_amplicons(trim_primers=False)
        outfile.write(amp_bed_str)
    with open(OUTPUT_DIR / "primertrim.amplicon.bed", "w") as outfile:
        outfile.write(scheme.to_amplicons(trim_primers=True))
    upload_pt.manual_update(n=2, update=True)

    # Write all the consensus sequences to a single file
    with dnaio.FastaWriter(
        OUTPUT_DIR / "reference.fasta", line_length=60
    ) as reference_outfile:
        for msa_obj in msa_dict.values():
            if config.mapping == MappingType.FIRST:
                seq_str = generate_reference(msa_obj.array)
            elif config.mapping == MappingType.CONSENSUS:
                seq_str = generate_consensus(msa_obj.array)
            else:
                raise ValueError("Mapping must be 'first' or 'consensus'")

            reference_outfile.write(
                dnaio.SequenceRecord(name=msa_obj._chrom_name, sequence=seq_str)
            )

    upload_pt.manual_update(n=3, update=True)

    # Create all hashes
    ## Generate the bedfile hash, and add it into the config
    primer_md5 = hashlib.md5("\n".join(primer_bed_str).encode()).hexdigest()
    cfg_dict["primer.bed.md5"] = primer_md5

    ## Generate the amplicon hash, and add it into the config
    amp_md5 = hashlib.md5(amp_bed_str.encode()).hexdigest()
    cfg_dict["amplicon.bed.md5"] = amp_md5

    ## Read in the reference file and generate the hash
    with open(OUTPUT_DIR / "reference.fasta") as reference_outfile:
        ref_md5 = hashlib.md5(reference_outfile.read().encode()).hexdigest()
    cfg_dict["reference.fasta.md5"] = ref_md5

    # Write the config dict to file
    with open(OUTPUT_DIR / "config.json", "w") as outfile:
        outfile.write(json.dumps(cfg_dict, sort_keys=True))

    upload_pt.manual_update(n=4, update=True)

    # Create qc coverage data
    qc_data = {}
    for msa_index, msa_obj in msa_dict.items():
        msa_data = {"coverage": scheme.get_coverage_percent(msa_index)}
        msa_data["n_amplicons"] = len(
            [x for x in scheme._last_pp_added if x.msa_index == msa_index]
        )
        msa_data["n_gaps"] = len(scheme.get_coverage_gaps(msa_index))
        qc_data[msa_obj.name] = msa_data

    with open(OUTPUT_DIR / "work" / "qc.json", "w") as outfile:
        outfile.write(json.dumps(qc_data, sort_keys=True))

    ## DO THIS LAST AS THIS CAN TAKE A LONG TIME

    # Create primer thermo profiles
    with open(OUTPUT_DIR / "work" / "primer_thermo.html", "w") as outfile:
        outfile.write(
            plot_primer_thermo_profile_html(
                OUTPUT_DIR / "primer.bed",
                config,
                offline_plots=offline_plots,
            )
        )

    # Writing plot data
    plot_data = generate_all_plotdata(
        list(msa_dict.values()),
        OUTPUT_DIR / "work",
        last_pp_added=scheme._last_pp_added,
    )
    upload_pt.manual_update(n=5, update=True)

    # Write the plot
    with open(OUTPUT_DIR / "plot.html", "w") as outfile:
        outfile.write(
            generate_all_plots_html(plot_data, OUTPUT_DIR, offline_plots=offline_plots)
        )
    upload_pt.manual_update(n=6, update=True)
    with open(OUTPUT_DIR / "primer.html", "w") as outfile:
        for i, msa_obj in enumerate(msa_dict.values()):
            try:
                outfile.write(
                    primer_mismatch_heatmap(
                        array=msa_obj.array,
                        seqdict=msa_obj._seq_dict,
                        bedfile=OUTPUT_DIR / "primer.bed",
                        offline_plots=True if offline_plots and i == 0 else False,
                        mapping=config.mapping,
                    )
                )
            except UsageError:
                logger.warning(
                    f"No Primers found for {msa_obj._chrom_name}. Skipping Plot!"
                )

    upload_pt.manual_update(n=7, update=True)

    # Close the progress tracker
    upload_pt.manual_update(n=upload_steps, update=True)
    sleep(0.5)  # Sleep to allow the progress bar to finish
    upload_pt.close()

    logger.info("Completed Successfully")
