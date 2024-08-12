# Core imports
import hashlib
import json

# General import
import pathlib
import shutil
import sys
from enum import Enum

from Bio import Seq, SeqIO, SeqRecord
from loguru import logger

# version import
from primalscheme3.core.bedfiles import read_in_bedprimerpairs
from primalscheme3.core.config import Config, MappingType
from primalscheme3.core.create_report_data import generate_all_plotdata
from primalscheme3.core.create_reports import generate_all_plots
from primalscheme3.core.logger import setup_logger
from primalscheme3.core.mapping import generate_consensus, generate_reference
from primalscheme3.core.mismatches import MatchDB
from primalscheme3.core.primer_visual import primer_mismatch_heatmap
from primalscheme3.core.progress_tracker import ProgressManager

# Module imports
from primalscheme3.panel.panel_classes import (
    Panel,
    PanelMSA,
    PanelReturn,
    Region,
)

logger = logger.opt(colors=True)


def mean_gc_diff(seqs: list[str] | set[str], target_gc=0.5) -> float:
    gc_diff = []
    for seq in seqs:
        gc_diff.append(abs(target_gc - ((seq.count("G") + seq.count("C")) / len(seq))))
    return sum(gc_diff) / len(seqs)


class PanelRunModes(Enum):
    ALL = "all"
    REGION_ONLY = "region-only"
    REGION_ALL = "region-all"


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
    inputbedfile: pathlib.Path | None = None,
    regionbedfile: pathlib.Path | None = None,
    mode: PanelRunModes = PanelRunModes.ALL,
    max_amplicons: int | None = None,
    offline_plots: bool = True,
    plot_mismatches: bool = True,
):
    ARG_MSA = msa
    OUTPUT_DIR = pathlib.Path(output_dir).absolute()

    # Config Dicts
    config_dict = config.to_json()
    config_dict["max_amplicons"] = max_amplicons
    config_dict["mode"] = mode.value

    # Enforce mapping
    if config.mapping != MappingType.FIRST:
        sys.exit("ERROR: mapping must be 'first'")

    # Enforce region only has a region bedfile
    if mode == PanelRunModes.REGION_ONLY and regionbedfile is None:
        sys.exit("ERROR: region-only mode requires a region bedfile")

    # See if the output dir already exists
    if OUTPUT_DIR.is_dir() and not force:
        sys.exit(f"ERROR: {OUTPUT_DIR} already exists, please use --force to override")

    # Create the output dir and a work subdir
    pathlib.Path.mkdir(OUTPUT_DIR, exist_ok=True)
    pathlib.Path.mkdir(OUTPUT_DIR / "work", exist_ok=True)

    ## Set up the logger
    logger = setup_logger(OUTPUT_DIR)

    ## Set up the progress manager
    if pm is None:
        pm = ProgressManager()

    # Create the mismatch db
    logger.info(
        "Creating the Mismatch Database",
    )
    mismatch_db = MatchDB(
        OUTPUT_DIR / "work/mismatch",
        [str(x) for x in ARG_MSA],
        config.mismatch_kmersize,
    )
    logger.info(
        "<green>Created:</> {path}",
        path=f"{OUTPUT_DIR.relative_to(OUTPUT_DIR.parent)}/work/mismatch.db",  # Write the path relative to the parent dir
    )

    regions_mapping: dict[Region, str | None] | None = None
    # Read in the regionbedfile if given
    if regionbedfile is not None:
        bed_lines = read_region_bedfile(regionbedfile)
        regions_mapping = {
            Region(
                bedline[0],
                int(bedline[1]),
                int(bedline[2]),
                bedline[3],
                int(bedline[4]),
            ): None
            for bedline in bed_lines
        }
        try:
            shutil.copy(regionbedfile, OUTPUT_DIR / regionbedfile.name)  # type: ignore # Copy the bedfile
        except shutil.SameFileError:
            pass

    ## Read in the MSAs
    msa_dict: dict[int, PanelMSA] = {}
    msa_data: dict = {}
    for msa_index, msa_path in enumerate(ARG_MSA):
        msa_data[msa_index] = {}
        # copy the msa into the output / work dir
        local_msa_path = OUTPUT_DIR / "work" / msa_path.name
        shutil.copy(msa_path, local_msa_path)

        # Create MSA checksum
        with open(local_msa_path, "rb") as f:
            msa_data[msa_index]["msa_checksum"] = hashlib.file_digest(
                f, "md5"
            ).hexdigest()

        # Read in the MSA
        msa_obj = PanelMSA(
            name=local_msa_path.stem,
            path=local_msa_path,
            msa_index=msa_index,
            mapping=config.mapping.value,
            logger=logger,
            progress_manager=pm,
        )
        logger.info(
            "Read in MSA: <blue>{msa_path} -> '{chromname}'</>\tseqs:<green>{msa_rows}</>\tcols:<green>{msa_cols}</>",
            msa_path=msa_path.name,
            chromname=msa_obj._chrom_name,
            msa_rows=msa_obj.array.shape[0],
            msa_cols=msa_obj.array.shape[1],
        )
        if "/" in msa_obj._chrom_name:
            new_chromname = msa_obj._chrom_name.split("/")[0]
            logger.warning(
                "<red>WARNING</>: Having a '/' in the chromname {msachromname} will cause issues with figure generation bedfile output. Parsing chromname <yellow>{msachromname}</> -> <green>{new_chromname}</>",
                msachromname=msa_obj._chrom_name,
                new_chromname=new_chromname,
            )
            msa_obj._chrom_name = new_chromname

        # Add the regions
        if regions_mapping is not None:
            msa_regions = []
            for region in regions_mapping.keys():
                if region.chromname == msa_obj._chrom_name:
                    msa_regions.append(region)
                    regions_mapping[region] = msa_obj._chrom_name
            msa_obj.add_regions(msa_regions)

            # Print Number mapped
            logger.info(
                "<blue>{msa_name}</>: <green>{mapped}</> regions mapped",
                msa_name=msa_obj._chrom_name,
                mapped=len(msa_regions),
            )

            # Create indexes from the regions
            indexes = set()
            for region in msa_regions:
                if max(region.start, region.stop) > len(msa_obj._mapping_array):
                    logger.error(
                        "Region {regionname} is out of bounds for {msa_name}",
                        regionname=region.name,
                        msa_name=msa_obj._chrom_name,
                    )
                    sys.exit(1)
                indexes.update(
                    [
                        msa_obj._ref_to_msa[x]
                        for x in range(region.start, region.stop)
                        if msa_obj._ref_to_msa[x] is not None
                    ]
                )

            findexes = list(
                {
                    fi
                    for fi in (
                        range(i - config.amplicon_size_max, i + 1) for i in indexes
                    )
                    for fi in fi
                    if fi >= 0 and fi < msa_obj.array.shape[1]
                }
            )
            findexes.sort()
            rindexes = list(
                {
                    ri
                    for ri in (range(i, i + config.amplicon_size_max) for i in indexes)
                    for ri in ri
                    if ri >= 0 and ri < msa_obj.array.shape[1]
                }
            )
            rindexes.sort()

        # Add some msa data to the dict
        msa_data[msa_index]["msa_name"] = msa_obj.name
        msa_data[msa_index]["msa_path"] = str(local_msa_path.absolute())
        msa_data[msa_index]["msa_chromname"] = msa_obj._chrom_name
        msa_data[msa_index]["msa_uuid"] = msa_obj._uuid

        logger.info(
            "Read in MSA: <blue>{msa_path}</>\tseqs:<green>{msa_rows}</>\tcols:<green>{msa_cols}</>",
            msa_path=msa_obj.name,
            msa_rows=msa_obj.array.shape[0],
            msa_cols=msa_obj.array.shape[1],
        )

        # Split the logic for the different modes
        match mode:
            case PanelRunModes.REGION_ONLY:
                msa_obj.digest(config=config, indexes=(findexes, rindexes))  # type: ignore
            case _:
                msa_obj.digest(config=config, indexes=None)

        # Log the digestion
        logger.info(
            "<blue>{msa_path}</>: <green>regions</> digested to <green>{num_fkmers}</> FKmers and <green>{num_rkmers}</> RKmers",
            msa_path=msa_obj._chrom_name,
            num_fkmers=len(msa_obj.fkmers),
            num_rkmers=len(msa_obj.rkmers),
        )

        # Generate all primerpairs
        msa_obj.generate_primerpairs(
            amplicon_size_max=config.amplicon_size_max,
            amplicon_size_min=config.amplicon_size_max,
            dimerscore=config.dimer_score,
        )
        logger.info(
            "<blue>{msa_path}</>: Generated <green>{num_pp}</> possible amplicons",
            msa_path=msa_obj._chrom_name,
            num_pp=len(msa_obj.primerpairs),
        )

        match mode:
            case PanelRunModes.REGION_ONLY:
                # Filter the primerpairs to only include the ones with scores (cover regions)
                msa_obj.primerpairs = [
                    x for x in msa_obj.primerpairs if msa_obj.get_pp_score(x) > 0
                ]
                msa_obj.primerpairs.sort(key=lambda x: x.fprimer.end)
                print("Filtered primerpairs", len(msa_obj.primerpairs))
            case _:
                continue

        # Add the MSA to the dict
        msa_dict[msa_index] = msa_obj

    # Add all the msa_data to the cfg
    config_dict["msa_data"] = msa_data

    ## Digestion finished, now create the panel

    # Create a lookup dict for the msa index to name
    msa_index_to_name = {k: v._chrom_name for k, v in msa_dict.items()}

    # Create a dict to store how many amplicons have been added to each msa
    msa_index_to_amplicon_count = {k: 0 for k in msa_data.keys()}

    # Create the panel object
    panel: Panel = Panel(msa_dict, config=config, matchdb=mismatch_db)

    # MSA_INDEX_TO_CHROMNAME =
    msa_chromname_to_index = {
        msa._chrom_name: msa.msa_index for msa in msa_dict.values()
    }

    # Read in the inputbedfile if given
    if inputbedfile is not None:
        bedprimerpairs, _headers = read_in_bedprimerpairs(inputbedfile)
        for bedprimerpair in bedprimerpairs:
            bedprimerpair.msa_index = msa_chromname_to_index.get(
                bedprimerpair.chrom_name,  # type: ignore
                -1,
            )
            panel._add_primerpair(
                bedprimerpair,
                pool=bedprimerpair.pool,
                msa_index=bedprimerpair.msa_index,  # type: ignore
            )
            logger.info(
                "Added primerpair from inputbedfile: {bedprimerpair}",
                bedprimerpair=bedprimerpair.amplicon_prefix,
            )
    # if pool > 1 then should be an ol scheme
    if config.n_pools > 1:
        # Sort primerpairs start position
        for msa_obj in panel._msa_dict.values():
            msa_obj.primerpairs.sort(key=lambda x: x.fprimer.end)

    # Add the first primerpair

    counter = 0
    while max_amplicons is None or counter < max_amplicons:
        match panel.add_next_primerpair():
            case PanelReturn.ADDED_PRIMERPAIR:
                # Update the amplicon count
                msa_index_to_amplicon_count[panel._last_pp_added[-1].msa_index] += 1
                logger.info(
                    "Added <blue>amplicon</> (<green>{msaampliconnumber}</>) for <blue>{msa_name}</>: {primer_start}\t{primer_end}",
                    primer_start=panel._last_pp_added[-1].fprimer.region()[0],
                    primer_end=panel._last_pp_added[-1].rprimer.region()[1],
                    msa_name=msa_index_to_name.get(panel._last_pp_added[-1].msa_index),
                    msaampliconnumber=msa_index_to_amplicon_count.get(
                        panel._last_pp_added[-1].msa_index
                    ),
                )
                counter += 1
                continue
            case PanelReturn.NO_PRIMERPAIRS:
                logger.info(
                    "No more valid <red>amplicons</> for <blue>{msa_name}</>",
                    msa_name=msa_index_to_name.get(panel._current_msa_index),
                )
                break

    # Try adding all remaining primerpairs
    while max_amplicons is None or counter < max_amplicons:
        match panel.keep_adding():
            case PanelReturn.ADDED_PRIMERPAIR:
                # Update the amplicon count
                msa_index_to_amplicon_count[panel._last_pp_added[-1].msa_index] += 1
                logger.info(
                    "Added <blue>amplicon</> (<green>{msaampliconnumber}</>) for <blue>{msa_name}</>: {primer_start}\t{primer_end}",
                    primer_start=panel._last_pp_added[-1].fprimer.region()[0],
                    primer_end=panel._last_pp_added[-1].rprimer.region()[1],
                    msa_name=msa_index_to_name.get(panel._last_pp_added[-1].msa_index),
                    msaampliconnumber=msa_index_to_amplicon_count.get(
                        panel._last_pp_added[-1].msa_index
                    ),
                )
                counter += 1
                continue
            case PanelReturn.NO_PRIMERPAIRS:
                logger.info("All primerpairs added")
                break
            case PanelReturn.NO_MORE_PRIMERPAIRS_IN_MSA:
                logger.info(
                    "No more valid <red>amplicons</> for <blue>{msa_name}</>",
                    msa_name=msa_index_to_name.get(panel._current_msa_index),
                )
            case PanelReturn.MOVING_TO_NEW_MSA:
                logger.debug(
                    "Skipping <blue>{msa_name}</>",
                    msa_name=msa_index_to_name.get(panel._current_msa_index),
                )  # Catch case but no nothing
            case _:
                print("ERROR?")

    # Log that the panel is finished
    logger.info(
        "Finished creating the panel. <green>{num_amplicons}</> amplicons created",
        num_amplicons=len(panel._last_pp_added),
    )
    # Print the amplicons count for each msa
    for msa_index, amplicon_count in msa_index_to_amplicon_count.items():
        logger.info(
            "<blue>{msa_name}</>: <green>{amplicon_count}</> amplicons",
            msa_name=msa_index_to_name.get(msa_index),
            amplicon_count=amplicon_count,
        )

    region_to_coverage = {}
    # If region bedfile given, check that all regions have been covered
    if regionbedfile is not None:
        for msa_obj in panel._msa_dict.values():
            assert msa_obj.regions is not None
            for region in msa_obj.regions:
                region_coverage = panel._coverage[msa_obj.msa_index][
                    region.start : region.stop
                ]
                region_mean_coverage = region_coverage.mean()

                logger.info(
                    "<blue>{msa_name}</>:{regionname} <yellow>{regionstart}:{regionstop}</> {percent} covered",
                    msa_name=msa_obj._chrom_name,
                    regionstart=region.start,
                    regionstop=region.stop,
                    percent=f"{round(region_mean_coverage * 100, 2)}%",
                    regionname=region.name,
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
    with open(OUTPUT_DIR / "reference.fasta", "w") as reference_outfile:
        reference_records = []
        for msa_obj in msa_dict.values():
            if config.mapping == MappingType.FIRST:
                seq_str = generate_reference(msa_obj.array)
            elif config.mapping == MappingType.CONSENSUS:
                seq_str = generate_consensus(msa_obj.array)
            else:
                raise ValueError("Mapping must be 'first' or 'consensus'")

            reference_records.append(
                SeqRecord.SeqRecord(
                    seq=Seq.Seq(seq_str),
                    id=msa_obj._chrom_name,
                )
            )

        SeqIO.write(reference_records, reference_outfile, "fasta")

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
    config_dict["regionbedfile"] = str(regionbedfile)
    config_dict["inputbedfile"] = str(inputbedfile)
    with open(OUTPUT_DIR / "config.json", "w") as outfile:
        outfile.write(json.dumps(config_dict, sort_keys=True))

    ## DO THIS LAST AS THIS CAN TAKE A LONG TIME
    # Writing plot data
    plot_data = generate_all_plotdata(
        list(msa_dict.values()),
        OUTPUT_DIR / "work",
        last_pp_added=panel._last_pp_added,
    )
    generate_all_plots(plot_data, OUTPUT_DIR, offline_plots=offline_plots)

    with open(OUTPUT_DIR / "primer.html", "w") as outfile:
        for i, msa_obj in enumerate(msa_dict.values()):
            outfile.write(
                primer_mismatch_heatmap(
                    array=msa_obj.array,
                    seqdict=msa_obj._seq_dict,
                    bedfile=OUTPUT_DIR / "primer.bed",
                    offline_plots=True if offline_plots and i == 0 else False,
                )
            )
