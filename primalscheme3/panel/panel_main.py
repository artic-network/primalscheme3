# Core imports
import hashlib
import json

# General import
import pathlib
import shutil
import sys
from enum import Enum
from math import sqrt

from Bio import Seq, SeqIO, SeqRecord
from loguru import logger

# version import
from primalscheme3 import __version__
from primalscheme3.core.bedfiles import read_in_bedprimerpairs
from primalscheme3.core.config import config_dict
from primalscheme3.core.create_report_data import generate_all_plotdata
from primalscheme3.core.create_reports import generate_all_plots
from primalscheme3.core.logger import setup_loger
from primalscheme3.core.mapping import generate_consensus, generate_reference
from primalscheme3.core.mismatches import MatchDB

# Module imports
from primalscheme3.panel.minimal_scheme_classes import (
    Panel,
    PanelMSA,
    PanelReturn,
    Region,
)

logger = logger.opt(colors=True)


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
    argmsa: list[pathlib.Path],
    outputdir: pathlib.Path,
    primer_gc_max: float,
    primer_gc_min: float,
    primer_tm_max: float,
    primer_tm_min: float,
    dimerscore: float,
    cores: int,
    ampliconsizemax: int,
    ampliconsizemin: int,
    force: bool,
    npools: int,
    reducekmers: bool,
    minbasefreq: float,
    regionbedfile: pathlib.Path | None,
    inputbedfile: pathlib.Path,
    mapping: str,
    maxamplicons: int,
    mode: PanelRunModes,
    ignore_n: bool,
):
    ARG_MSA = argmsa
    OUTPUT_DIR = pathlib.Path(outputdir).absolute()

    cfg = config_dict
    # Primer Digestion settings
    cfg["primer_gc_min"] = primer_gc_min
    cfg["primer_gc_max"] = primer_gc_max
    cfg["primer_tm_min"] = primer_tm_min
    cfg["primer_tm_max"] = primer_tm_max
    cfg["dimerscore"] = dimerscore
    cfg["n_cores"] = cores
    cfg["output_dir"] = str(OUTPUT_DIR)
    cfg["amplicon_size_max"] = ampliconsizemax
    cfg["amplicon_size_min"] = ampliconsizemin
    cfg["force"] = force
    cfg["npools"] = npools
    cfg["reducekmers"] = reducekmers
    cfg["minbasefreq"] = minbasefreq
    cfg["mode"] = mode.value

    # Add the mismatch params to the cfg
    cfg["mismatch_fuzzy"] = True
    cfg["mismatch_kmersize"] = 14
    cfg["mismatch_product_size"] = ampliconsizemax

    # Add plots to the cfg
    cfg["plot"] = True
    cfg["disable_progress_bar"] = False

    # Add the bedfile to the cfg
    cfg["regionbedfile"] = str(regionbedfile) if regionbedfile is not None else None
    cfg["inputbedfile"] = str(inputbedfile)

    # Set the mapping
    cfg["mapping"] = mapping

    # Set the max amount of amplicons
    cfg["maxamplicons"] = maxamplicons

    # Add the version
    cfg["algorithmversion"] = f"primalscheme3:{__version__}"
    cfg["primerclass"] = "primerpanels"

    cfg["primer_size_min"] = 14

    # Add ignore_n to the cfg
    cfg["ignore_n"] = ignore_n

    # Enforce only one pool
    if npools != 1:
        sys.exit("ERROR: primalpanel only supports one pool")

    # Enforce region only has a region bedfile
    if mode == PanelRunModes.REGION_ONLY and regionbedfile is None:
        sys.exit("ERROR: region-only mode requires a region bedfile")

    # See if the output dir already exsits
    if OUTPUT_DIR.is_dir() and not force:
        sys.exit(f"ERROR: {OUTPUT_DIR} already exists, please use --force to override")

    # Create the output dir and a work subdir
    pathlib.Path.mkdir(OUTPUT_DIR, exist_ok=True)
    pathlib.Path.mkdir(OUTPUT_DIR / "work", exist_ok=True)

    ## Set up the logger
    logger = setup_loger(OUTPUT_DIR)

    # Create the mismatch db
    logger.info(
        "Creating the Mismatch Database",
    )
    mismatch_db = MatchDB(
        OUTPUT_DIR / "work/mismatch", [str(x) for x in ARG_MSA], cfg["primer_size_min"]
    )
    logger.info(
        "<green>Created:</> {path}",
        path=f"{OUTPUT_DIR.relative_to(OUTPUT_DIR.parent)}/work/mismatch.db",  # Write the path relative to the parent dir
    )

    regions_mapping: dict[Region, str | None] | None = None
    # Read in the regionbedfile if given
    if cfg["regionbedfile"] is not None:
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
        msa = PanelMSA(
            name=local_msa_path.stem,
            path=local_msa_path,
            msa_index=msa_index,
            mapping=cfg["mapping"],
            logger=logger,
        )
        if "/" in msa._chrom_name:
            new_chromname = msa._chrom_name.split("/")[0]
            logger.warning(
                "<red>WARNING</>: Having a '/' in the chromname {msachromname} will cause issues with figure generation bedfile output. Parsing chromname <yellow>{msachromname}</> -> <green>{new_chromname}</>",
                msachromname=msa._chrom_name,
                new_chromname=new_chromname,
            )
            msa._chrom_name = new_chromname

        # Add the regions
        if regions_mapping is not None:
            msa_regions = []
            for region in regions_mapping.keys():
                if region.chromname == msa._chrom_name:
                    msa_regions.append(region)
                    regions_mapping[region] = msa._chrom_name
            msa.add_regions(msa_regions)

            # Print Number mapped
            logger.info(
                "<blue>{msa_name}</>: <green>{mapped}</> regions mapped",
                msa_name=msa._chrom_name,
                mapped=len(msa_regions),
            )

            # Create indexes from the regions
            indexes = set()
            for region in msa_regions:
                indexes.update(range(region.start, region.stop))

            findexes = {
                fi
                for fi in (range(i - cfg["amplicon_size_max"], i + 1) for i in indexes)
                for fi in fi
                if fi >= 0 and fi < msa.array.shape[1]
            }
            rindexes = {
                ri
                for ri in (range(i, i + cfg["amplicon_size_max"]) for i in indexes)
                for ri in ri
                if ri >= 0 and ri < msa.array.shape[1]
            }

        # Add some msa data to the dict
        msa_data[msa_index]["msa_name"] = msa.name
        msa_data[msa_index]["msa_path"] = str(local_msa_path.absolute())
        msa_data[msa_index]["msa_chromname"] = msa._chrom_name
        msa_data[msa_index]["msa_uuid"] = msa._uuid

        logger.info(
            "Read in MSA: <blue>{msa_path}</>\tseqs:<green>{msa_rows}</>\tcols:<green>{msa_cols}</>",
            msa_path=msa_path.name,
            msa_rows=msa.array.shape[0],
            msa_cols=msa.array.shape[1],
        )

        # Split the logic for the different modes
        match mode:
            case PanelRunModes.REGION_ONLY:
                msa.digest(cfg, indexes=(findexes, rindexes))  # type: ignore
                msa.remove_kmers_that_clash_with_regions()
            case _:
                msa.digest(cfg, indexes=None)

        # Log the digestion
        logger.info(
            "<blue>{msa_path}</>: <green>regions</> digested to <green>{num_fkmers}</> FKmers and <green>{num_rkmers}</> RKmers",
            msa_path=msa_path.stem,
            num_fkmers=len(msa.fkmers),
            num_rkmers=len(msa.rkmers),
        )

        # Generate all primerpairs
        msa.generate_primerpairs(
            amplicon_size_max=cfg["amplicon_size_max"],
            amplicon_size_min=cfg["amplicon_size_min"],
            dimerscore=cfg["dimerscore"],
        )
        logger.info(
            "<blue>{msa_path}</>: Generated <green>{num_pp}</> possible amplicons",
            msa_path=msa_path.stem,
            num_pp=len(msa.primerpairs),
        )

        if mode == PanelRunModes.ALL:
            msa.primerpairs.sort(
                key=lambda pp: msa.get_pp_entropy(pp) ** 2 / sqrt(len(pp.all_seqs())),  # type: ignore
                reverse=True,
            )
        if mode == PanelRunModes.REGION_ONLY:
            # Sort the primerpairs by the number of SNPs in the amplicon
            msa.primerpairs.sort(
                key=lambda pp: (
                    msa.get_pp_score(pp),
                    -sqrt(len(pp.all_seqs())),  # type: ignore
                    pp.fprimer.__hash__(),
                ),  # type: ignore # Use a HASH Prevent sequential primerpairs being added
                reverse=True,
            )
            # Remove all primerpairs with a score of 0
            msa.primerpairs = [pp for pp in msa.primerpairs if msa.get_pp_score(pp) > 0]

        # Add the MSA to the dict
        msa_dict[msa_index] = msa

    # Add all the msa_data to the cfg
    cfg["msa_data"] = msa_data

    ## Digestion finished, now create the panel

    # Create a lookup dict for the msa index to name
    msa_index_to_name = {k: v for k, v in enumerate([msa.name for msa in ARG_MSA])}

    # Create a dict to store how many amplicons have been addded to each msa
    msa_index_to_amplicon_count = {k: 0 for k in msa_index_to_name.keys()}

    # Create the panel object
    panel: Panel = Panel([x for x in msa_dict.values()], cfg, mismatch_db)

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
            logger.debug(
                "Added primerpair from inputbedfile: {bedprimerpair}",
                bedprimerpair=bedprimerpair,
            )

    counter = 0
    while counter < cfg["maxamplicons"]:
        match panel.try_add_primerpair():
            case PanelReturn.ADDED_PRIMERPAIR:
                # Update the amplicon count
                msa_index_to_amplicon_count[panel._last_pp_added[-1].msa_index] += 1
                logger.info(
                    "Added <blue>amplicon</> (<green>{msaampliconnumber}</>) for <blue>{msa_name}</>: {primer_start}\t{primer_end}",
                    primer_start=panel._last_pp_added[-1].start,
                    primer_end=panel._last_pp_added[-1].end,
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
    while counter < cfg["maxamplicons"]:
        match panel.keep_adding():
            case PanelReturn.ADDED_PRIMERPAIR:
                # Update the amplicon count
                msa_index_to_amplicon_count[panel._last_pp_added[-1].msa_index] += 1
                logger.info(
                    "Added <blue>amplicon</> (<green>{msaampliconnumber}</>) for <blue>{msa_name}</>: {primer_start}\t{primer_end}",
                    primer_start=panel._last_pp_added[-1].start,
                    primer_end=panel._last_pp_added[-1].end,
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

    # If region bedfile given, check that all regions have been covered
    if cfg["regionbedfile"] is not None:
        for msa in panel.msas:
            regions_covered = []
            regions_not_covered = []
            # Find all primerpairs in this msa
            primer_pairs = [
                x for x in panel._last_pp_added if x.msa_index == msa.msa_index
            ]
            covered_positions = {
                x
                for x in (
                    range(pp.fprimer.end, pp.rprimer.start) for pp in primer_pairs
                )
                for x in x
            }
            if msa.regions is None:
                continue
            # See which regions are completely covered by the covered_positions
            for region in msa.regions:
                all_regions = {x for x in range(region.start, region.stop)}

                if all_regions.issubset(covered_positions):
                    regions_covered.append(region)
                else:
                    regions_not_covered.append(region)

            logger.info(
                "<blue>{msa_name}</>: <green>{covered}</> covered regions, <red>{not_covered}</> not covered regions",
                msa_name=msa.name,
                covered=len(regions_covered),
                not_covered=len(regions_not_covered),
            )

            # Write which regions are not covered to logger
            for region in regions_not_covered:
                logger.info(
                    "<blue>{msa_name}</>: <red>{regionstart}:{regionstop}</> not covered",
                    msa_name=msa.name,
                    regionstart=region.start,
                    regionstop=region.stop,
                )

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
        if cfg["mapping"] == "first":
            for msa in msa_dict.values():
                reference_records.append(
                    SeqRecord.SeqRecord(
                        seq=Seq.Seq(generate_reference(msa.array)),
                        id=msa._chrom_name,
                    )
                )
        elif cfg["mapping"] == "consensus":
            for msa in msa_dict.values():
                reference_records.append(
                    SeqRecord.SeqRecord(
                        seq=Seq.Seq(generate_consensus(msa.array)),
                        id=msa._chrom_name,
                    )
                )
        SeqIO.write(reference_records, reference_outfile, "fasta")

    # Generate all the hashes
    ## Generate the bedfile hash, and add it into the config
    primer_md5 = hashlib.md5("\n".join(primer_bed_str).encode()).hexdigest()
    cfg["primer.bed.md5"] = primer_md5

    ## Generate the amplicon hash, and add it into the config
    amp_md5 = hashlib.md5(amp_bed_str.encode()).hexdigest()
    cfg["amplicon.bed.md5"] = amp_md5

    ## Read in the reference file and generate the hash
    with open(OUTPUT_DIR / "reference.fasta") as reference_outfile:
        ref_md5 = hashlib.md5(reference_outfile.read().encode()).hexdigest()
    cfg["reference.fasta.md5"] = ref_md5

    # Write the config dict to file
    # Add the bedfile to the cfg
    cfg["regionbedfile"] = str(regionbedfile)
    cfg["inputbedfile"] = str(inputbedfile)
    with open(OUTPUT_DIR / "config.json", "w") as outfile:
        outfile.write(json.dumps(cfg, sort_keys=True))

    ## DO THIS LAST AS THIS CAN TAKE A LONG TIME
    # Writing plot data
    plot_data = generate_all_plotdata(
        list(msa_dict.values()),
        OUTPUT_DIR / "work",
        last_pp_added=panel._last_pp_added,
    )
    generate_all_plots(plot_data, OUTPUT_DIR)
