# Core imports
import hashlib
import json

# General import
import pathlib
import shutil
import sys

import dnaio
from click import UsageError

# version import
from primalscheme3.core.bedfiles import read_in_extra_primers
from primalscheme3.core.config import Config, MappingType
from primalscheme3.core.create_report_data import generate_all_plotdata
from primalscheme3.core.create_reports import generate_all_plots_html
from primalscheme3.core.logger import setup_rich_logger
from primalscheme3.core.mapping import generate_consensus, generate_reference
from primalscheme3.core.mismatches import MatchDB
from primalscheme3.core.primer_visual import (
    plot_primer_thermo_profile_html,
    primer_mismatch_heatmap,
)
from primalscheme3.core.progress_tracker import ProgressManager

# Module imports
from primalscheme3.panel.panel_classes import (
    Panel,
    PanelMSA,
    PanelReturn,
    PanelRunModes,
    Region,
    RegionParser,
)


def mean_gc_diff(seqs: list[str] | set[str], target_gc=0.5) -> float:
    gc_diff = []
    for seq in seqs:
        gc_diff.append(abs(target_gc - ((seq.count("G") + seq.count("C")) / len(seq))))
    return sum(gc_diff) / len(seqs)


def read_region_bedfile(path) -> list[list[str]]:
    """
    Bedfiles need to be in the format:
    chrom start end name score
    """
    ## If bedfile given, parse it:
    bed_lines = []
    with open(path) as bedfile:
        for line in bedfile.readlines():
            line = line.strip()
            if not line or line.startswith("#"):  # Skip empty lines and header lines
                continue
            bed_lines.append(line.split("\t"))
    return bed_lines


def panelcreate(
    msa: list[pathlib.Path],
    output_dir: pathlib.Path,
    config: Config,
    pm: ProgressManager | None,
    force: bool = False,
    input_bedfile: pathlib.Path | None = None,
    region_bedfile: pathlib.Path | None = None,
    mode: PanelRunModes = PanelRunModes.ENTROPY,
    max_amplicons: int | None = None,
    max_amplicons_msa: int | None = None,
    max_amplicons_region_group: int | None = None,
    offline_plots: bool = True,
):
    ARG_MSA = msa
    OUTPUT_DIR = pathlib.Path(output_dir).absolute()

    # Config Dicts
    config_dict = config.to_dict()
    config_dict["max_amplicons"] = max_amplicons
    config_dict["mode"] = mode.value
    config_dict["max_amplicons_msa"] = max_amplicons_msa
    config_dict["max_amplicons_region_group"] = max_amplicons_region_group

    # Enforce region only has a region bedfile
    if mode == PanelRunModes.REGION_ONLY and region_bedfile is None:
        raise UsageError(
            "region-only mode requires a region bedfile. Provide with --region-bedfile"
        )

    # See if the output dir already exists
    if OUTPUT_DIR.is_dir() and not force:
        raise UsageError(f"{OUTPUT_DIR} already exists, please use --force to override")

    # Create the output dir and a work subdir
    pathlib.Path.mkdir(OUTPUT_DIR, exist_ok=True)
    pathlib.Path.mkdir(OUTPUT_DIR / "work", exist_ok=True)

    ## Set up the logger
    logger = setup_rich_logger(str(OUTPUT_DIR / "work" / "file.log"))

    # Check
    if mode != PanelRunModes.REGION_ONLY and max_amplicons_region_group is not None:
        logger.warning(
            "max-amplicons-group is set but mode is not region-only. Ignoring"
        )
        max_amplicons_region_group = None

    ## Set up the progress manager
    if pm is None:
        pm = ProgressManager()

    # Create the mismatch db
    logger.info(
        "Creating the Mismatch Database",
    )
    mismatch_db = MatchDB(
        OUTPUT_DIR / "work/mismatch",
        [str(x) for x in ARG_MSA] if config.use_matchdb else [],
        config.mismatch_kmersize,
    )
    logger.info(
        f"[green]Created[/green]: "
        f"{OUTPUT_DIR.relative_to(OUTPUT_DIR.parent)}/work/mismatch.db",
    )

    regions_mapping: dict[Region, str | None] | None = None
    # Read in the region_bedfile if given
    if region_bedfile is not None:
        region_bed_lines = read_region_bedfile(region_bedfile)

        regions_mapping = {
            RegionParser.from_list(bed_list): None for bed_list in region_bed_lines
        }
        try:
            shutil.copy(region_bedfile, OUTPUT_DIR / region_bedfile.name)  # type: ignore # Copy the bedfile
        except shutil.SameFileError:
            pass

    # If the bedfile flag is given add the primers into the scheme
    if input_bedfile is not None:
        bedprimerpairs = read_in_extra_primers(input_bedfile, config, logger)

    ## Read in the MSAs
    msa_dict: dict[int, PanelMSA] = {}
    msa_data: dict = {}
    for msa_index, msa_path in enumerate(ARG_MSA):
        msa_data[msa_index] = {}

        # Read in the MSA
        msa_obj = PanelMSA(
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
            f"Read in MSA: [blue]{local_msa_path.name}[/blue] ({msa_obj._chrom_name})\t"
            f"seqs:[green]{msa_obj.array.shape[0]}[/green]\t"
            f"cols:[green]{msa_obj.array.shape[1]}[/green]"
        )
        # Add the MSA to the dict
        msa_dict[msa_index] = msa_obj

        # Add some msa data to the dict
        msa_data[msa_index]["msa_name"] = msa_obj.name
        msa_data[msa_index]["msa_path"] = str(
            "work/" + msa_path.name
        )  # Write local path
        msa_data[msa_index]["msa_chromname"] = msa_obj._chrom_name
        msa_data[msa_index]["msa_uuid"] = msa_obj._uuid

        # Add the regions
        msa_regions = None
        if regions_mapping is not None:
            msa_regions = []
            for region in regions_mapping.keys():
                if region.chromname == msa_obj._chrom_name:
                    msa_regions.append(region)
                    regions_mapping[region] = msa_obj._chrom_name

            # Print Number mapped
            logger.info(
                f"[blue]{msa_obj._chrom_name}[/blue]: "
                f"[green]{len(msa_regions)}[/green] regions mapped",
            )

        # is mode is all msa_regions is None
        # Creates the score array
        msa_obj.add_regions(msa_regions)
        msa_obj.create_score_array(
            msa_regions, mode=mode
        )  # if mode is not REGION_ONLY, msa_regions is None

    # Print mapping stats
    if regions_mapping is not None:
        mapped_regions = sum(1 for x in regions_mapping.values() if x is not None)
        unmapped_regions = len(regions_mapping) - mapped_regions
        logger.info(
            f"Regions mapped: ([green]{mapped_regions}[/green]). Regions not mapped: ({unmapped_regions})",
        )
        if unmapped_regions > 0:
            logger.warning(
                f"Regions not mapped: {set([x.chromname for x, y in regions_mapping.items() if y is None])}"
            )
        if mapped_regions == 0:
            logger.critical(
                f"No regions mapped. Input genomes use chrom names ({', '.join(x._chrom_name for x in msa_dict.values())}). Regions use chrom names ({', '.join(set(x.chromname for x in regions_mapping.keys()))})"
            )
            sys.exit(1)

        # Check all genomes have regions
        for msa_obj in msa_dict.values():
            if msa_obj._chrom_name not in regions_mapping.values():  # type: ignore
                logger.critical(
                    f"Provided genome ({msa_obj._chrom_name}) has no regions. Please remove genome or add regions"
                )
                sys.exit(1)

    # Start the digestion loop
    for _msa_index, msa_obj in msa_dict.items():
        if msa_obj.regions is not None:
            indexes = set()
            for region in msa_obj.regions:
                indexes.update(
                    range(
                        msa_obj._ref_to_msa[region.start],
                        msa_obj._ref_to_msa[region.stop],
                    )
                )
            # Clean up
            findexes = sorted(
                {
                    fi
                    for fi in (
                        range(i - config.amplicon_size_max, i + 1) for i in indexes
                    )
                    for fi in fi
                    if fi >= 0 and fi < msa_obj.array.shape[1]
                }
            )
            rindexes = sorted(
                {
                    ri
                    for ri in (range(i, i + config.amplicon_size_max) for i in indexes)
                    for ri in ri
                    if ri >= 0 and ri < msa_obj.array.shape[1]
                }
            )
        # Split the logic for the different modes
        logger.info(
            f"Digesting [blue]{msa_obj._chrom_name}[/blue]",
        )
        match mode:
            case PanelRunModes.REGION_ONLY:
                msa_obj.digest_rs(
                    config=config,
                    indexes=(findexes, rindexes),  # type: ignore
                )
            case _:
                msa_obj.digest_rs(
                    config=config,
                    indexes=None,
                )

        # Log the digestion
        logger.info(
            f"[blue]{msa_obj._chrom_name}[/blue]: digested to "
            f"[green]{len(msa_obj.fkmers)}[/green] FKmers and "
            f"[green]{len(msa_obj.rkmers)}[/green] RKmers"
        )
        if len(msa_obj.fkmers) == 0 or len(msa_obj.rkmers) == 0:
            logger.critical(
                f"No valid FKmers or RKmers found for [blue]{msa_obj._chrom_name}[/blue]"
            )
            continue

        # Generate all primerpairs
        msa_obj.generate_primerpairs(
            amplicon_size_max=config.amplicon_size_max,
            amplicon_size_min=config.amplicon_size_max,
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

        match mode:
            case PanelRunModes.REGION_ONLY:
                # Filter the primerpairs to only include the ones with scores (cover regions)
                msa_obj.primerpairs = [
                    x for x in msa_obj.primerpairs if msa_obj.get_pp_score(x) > 0
                ]
            case _:
                continue

    # Add all the msa_data to the cfg
    config_dict["msa_data"] = msa_data

    ## Digestion finished, now create the panel

    # Create the regions for the regions_mapping
    chrom_name_regions_mapping: dict[str, dict[str, list[Region]]] = {}
    if regions_mapping is not None:
        for region, chrom_name in regions_mapping.items():
            # Add only mapped regions
            if chrom_name is None:
                continue
            if chrom_name not in chrom_name_regions_mapping:
                chrom_name_regions_mapping[chrom_name] = {}

            if region.name not in chrom_name_regions_mapping[chrom_name]:
                chrom_name_regions_mapping[chrom_name][region.name] = []
            chrom_name_regions_mapping[chrom_name][region.name].append(region)

    # Create a lookup dict for the msa index to name
    msa_index_to_name = {k: v._chrom_name for k, v in msa_dict.items()}

    # Create a dict to store how many amplicons have been added to each msa
    msa_index_to_amplicon_count = {k: 0 for k in msa_data.keys()}

    # Create the panel object
    panel: Panel = Panel(msa_dict, config=config, matchdb=mismatch_db, logger=logger)

    # MSA_INDEX_TO_CHROMNAME =
    msa_chromname_to_index = {
        msa._chrom_name: msa.msa_index for msa in msa_dict.values()
    }

    # Read in the input_bedfile if given
    if input_bedfile is not None:
        for bedpp in bedprimerpairs:  # type: ignore
            bedpp.msa_index = msa_chromname_to_index.get(
                bedpp.chrom_name,  # type: ignore
                -1,
            )
            panel._add_primerpair(
                bedpp,
                pool=bedpp.pool,
                msa_index=bedpp.msa_index,  # type: ignore
            )
            logger.debug(
                f"Added {bedpp.amplicon_prefix} from [blue]{input_bedfile.name}[/blue]",
            )

    # Add the first primerpair
    counter = 0

    # Randomise the order of primerpairs, sorting them by fprimer hash
    for msa_obj in msa_dict.values():
        msa_obj.primerpairs = sorted(
            msa_obj.primerpairs, key=lambda x: x.fprimer.__hash__()
        )

    while max_amplicons is None or counter < max_amplicons:
        # Check if the panel is complete
        if all(panel._is_msa_index_finished.values()):
            break

        match panel.add_next_primerpair(max_amplicons_group=max_amplicons_region_group):
            case PanelReturn.ADDED_PRIMERPAIR:
                added_pp = panel._last_pp_added[-1]
                # Update the amplicon count
                msa_index_to_amplicon_count[added_pp.msa_index] += 1
                logger.info(
                    f"Added amplicon ([green]{msa_index_to_amplicon_count.get(added_pp.msa_index)}[/green]) "
                    f"for [blue]{msa_index_to_name.get(added_pp.msa_index)}[/blue]: "
                    f"{added_pp.fprimer.region()[0]}\t{added_pp.rprimer.region()[1]}\t{added_pp.pool + 1}",
                )
                counter += 1

                # If number of amplicons per msa is reached. Set as finished
                if max_amplicons_msa is not None and (
                    msa_index_to_amplicon_count[added_pp.msa_index] >= max_amplicons_msa
                ):
                    panel._is_msa_index_finished[added_pp.msa_index] = True

                chrom_name = msa_index_to_name[added_pp.msa_index]
                continue
            case PanelReturn.NO_PRIMERPAIRS_IN_MSA:
                logger.debug("Skipping MSA")
                continue
            case _:  # pragma: no cover
                logger.error("Unknown return from add_next_primerpair")

    # Log that the panel is finished
    logger.info(
        f"Finished creating the panel. [green]{len(panel._last_pp_added)}[/green] amplicons total",
    )
    # Print the amplicons count for each msa
    for msa_index, amplicon_count in msa_index_to_amplicon_count.items():
        logger.info(
            f"[blue]{msa_index_to_name.get(msa_index)}[/blue]: [green]{amplicon_count}[/green] amplicons",
        )

    region_to_coverage = {}
    # If region bedfile given, check that all regions have been covered
    if region_bedfile is not None:
        for msa_obj in panel._msa_dict.values():
            assert msa_obj.regions is not None
            for region in msa_obj.regions:
                # Skip if region is in a group
                if region.group is not None:
                    continue

                region_coverage = panel._coverage[msa_obj.msa_index][
                    region.start : region.stop
                ]
                region_mean_coverage = region_coverage.mean()
                percent = round(region_mean_coverage * 100, 2)

                log_percent_col = "yellow" if percent < 100 else "green"
                logger.info(
                    f"[blue]{msa_obj._chrom_name}[/blue]:{region.name} "
                    f"{region.start}:{region.stop} "
                    f"[{log_percent_col}]{percent}%[/{log_percent_col}] covered",
                )
                region_to_coverage[region] = region_mean_coverage

    logger.info(
        "Writing outputs...",
    )
    # Write primer bed file
    with open(OUTPUT_DIR / "primer.bed", "w") as outfile:
        primer_bed_str = panel.to_bed()
        outfile.write(primer_bed_str)

    # Write amplicon bed file
    with open(OUTPUT_DIR / "amplicon.bed", "w") as outfile:
        amp_bed_str = panel.to_amplicons(trim_primers=False)
        outfile.write(amp_bed_str)
    with open(OUTPUT_DIR / "primertrim.amplicon.bed", "w") as outfile:
        outfile.write(panel.to_amplicons(trim_primers=True))

    # Write all the consensus sequences to a single file
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

    # Generate all the hashes
    ## Generate the bedfile hash, and add it into the config
    primer_md5 = hashlib.md5("\n".join(primer_bed_str).encode()).hexdigest()
    config_dict["primer.bed.md5"] = primer_md5

    ## Generate the amplicon hash, and add it into the config
    amp_md5 = hashlib.md5(amp_bed_str.encode()).hexdigest()
    config_dict["amplicon.bed.md5"] = amp_md5

    ## Read in the reference file and generate the hash
    with open(OUTPUT_DIR / "reference.fasta") as reference_outfile:
        ref_md5 = hashlib.md5(reference_outfile.read().encode()).hexdigest()
    config_dict["reference.fasta.md5"] = ref_md5

    # Write the config dict to file
    # Add the bedfile to the cfg
    config_dict["region_bedfile"] = str(region_bedfile)
    config_dict["input_bedfile"] = str(input_bedfile)
    with open(OUTPUT_DIR / "config.json", "w") as outfile:
        outfile.write(json.dumps(config_dict, sort_keys=True))

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
        last_pp_added=panel._last_pp_added,
    )

    # Write the plot
    with open(OUTPUT_DIR / "plot.html", "w") as outfile:
        outfile.write(
            generate_all_plots_html(plot_data, OUTPUT_DIR, offline_plots=offline_plots)
        )

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

    logger.info("Completed Successfully")
