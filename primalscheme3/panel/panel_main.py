# Core imports
from primalscheme3.core.config import config_dict
from primalscheme3.core.mismatches import MatchDB
from primalscheme3.core.create_reports import generate_all_plots
from primalscheme3.core.create_report_data import generate_all_plotdata
from primalscheme3.core.mapping import generate_consensus, generate_reference
from primalscheme3.core.bedfiles import read_in_bedprimerpairs
from primalscheme3.core.logger import setup_loger

# version import
from primalscheme3 import __version__

# Module imports
from primalscheme3.panel.minimal_scheme_classes import (
    PanelMSA,
    Panel,
    PanelReturn,
    Region,
)

# General import
import pathlib
import sys
from loguru import logger
import json
from math import sqrt
from Bio import SeqIO, SeqRecord, Seq
import shutil
import hashlib


logger = logger.opt(colors=True)


def read_region_bedfile(path) -> list[list[str]]:
    ## If bedfile given, parse it:
    bed_lines = []
    try:
        with open(path, "r") as bedfile:
            for line in bedfile.readlines():
                line = line.strip()
                if line:
                    data = line.split()[:3]

                    # Validate bedfile
                    if int(data[1]) > int(data[2]):
                        raise ValueError(
                            f"Start position is greater than end position: {line}"
                        )
                    bed_lines.append(line.split()[:3])
    except:
        raise ValueError(f"Cannot parse bedfile: {path}")
    return bed_lines


def panelcreate(args):
    ARG_MSA = args.msa
    OUTPUT_DIR = pathlib.Path(args.output).absolute()

    cfg = config_dict
    # Primer Digestion settings
    cfg["primer_gc_min"] = args.primer_gc_min
    cfg["primer_gc_max"] = args.primer_gc_max
    cfg["primer_tm_min"] = args.primer_tm_min
    cfg["primer_tm_max"] = args.primer_tm_max
    cfg["dimerscore"] = args.dimerscore
    cfg["n_cores"] = args.cores
    cfg["output_dir"] = str(OUTPUT_DIR)
    cfg["amplicon_size_max"] = args.ampliconsizemax
    cfg["amplicon_size_min"] = args.ampliconsizemin
    cfg["force"] = args.force
    cfg["npools"] = args.npools
    cfg["reducekmers"] = args.reducekmers
    cfg["minbasefreq"] = args.minbasefreq

    # Add the mismatch params to the cfg
    cfg["mismatch_fuzzy"] = True
    cfg["mismatch_kmersize"] = 20
    cfg["mismatch_product_size"] = args.ampliconsizemax

    # Add plots to the cfg
    cfg["plot"] = True
    cfg["disable_progress_bar"] = False

    # Add the bedfile to the cfg
    cfg["regionbedfile"] = (
        str(args.regionbedfile) if args.regionbedfile is not None else None
    )
    cfg["inputbedfile"] = str(args.inputbedfile)

    # Set the mapping
    cfg["mapping"] = args.mapping

    # Set the max amount of amplicons
    cfg["maxamplicons"] = args.maxamplicons

    # Add the version
    cfg["algorithmversion"] = f"primalscheme3:{__version__}"
    cfg["primerclass"] = "primerpanels"

    # Enforce only one pool
    if args.npools != 1:
        sys.exit("ERROR: primalpanel only supports one pool")

    # Enforce region only has a region bedfile
    if args.mode == "region-only" and args.regionbedfile is None:
        sys.exit("ERROR: region-only mode requires a region bedfile")

    # See if the output dir already exsits
    if OUTPUT_DIR.is_dir() and not args.force:
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
    mismatch_db = MatchDB(OUTPUT_DIR / "work/mismatch", ARG_MSA, cfg["primer_size_min"])
    logger.info(
        "<green>Created:</> {path}",
        path=f"{OUTPUT_DIR.relative_to(OUTPUT_DIR.parent)}/work/mismatch.db",  # Write the path relative to the parent dir
    )

    # Read in the regionbedfile if given
    if cfg["regionbedfile"] is not None:
        msa_name_to_regions: dict[str, list[Region]] = {
            msa_path.stem: [] for msa_path in ARG_MSA
        }
        bed_lines = read_region_bedfile(args.regionbedfile)

        # Map each bedfile line to an MSA regions
        msa_name_to_indexes = dict()
        failed_region_count = 0
        mapped_region_count = 0

        for bedline in bed_lines:
            # If the region can be mapped to an MSA using file name
            if bedline[0] in msa_name_to_regions.keys():
                msa_name_to_regions[bedline[0]].append(
                    Region(bedline[0], int(bedline[1]), int(bedline[2]))
                )
                mapped_region_count += 1
            else:
                logger.info(
                    "Cannot map <red>{bedline_name}</>:<blue>{bedline_start}</>:<blue>{bedline_end}</>",
                    bedline_name=bedline[0],
                    bedline_start=bedline[1],
                    bedline_end=bedline[2],
                    msa=bedline[0],
                )
                failed_region_count += 1

        logger.info(
            "Mapped (<green>{mapped_region_count}</>) regions. Failed to map (<red>{failed_region_count}</>) regions",
            mapped_region_count=mapped_region_count,
            failed_region_count=failed_region_count,
        )

        # Turn the regions into a range of indexes
        for msa_name, regions in msa_name_to_regions.items():
            findexes = {
                fi
                for fi in (
                    range(r.start - cfg["amplicon_size_max"], r.start - 1)
                    for r in regions
                )
                for fi in fi
            }
            rindexes = {
                ri
                for ri in (
                    range(r.stop + 1, r.stop + cfg["amplicon_size_max"])
                    for r in regions
                )
                for ri in ri
            }
            msa_name_to_indexes[msa_name] = (findexes, rindexes)
    else:
        msa_name_to_regions = {msa_path.stem: [] for msa_path in ARG_MSA}

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
            regions=msa_name_to_regions.get(msa_path.stem, []),
        )
        if "/" in msa._chrom_name:
            new_chromname = msa._chrom_name.split("/")[0]
            logger.warning(
                "<red>WARNING</>: Having a '/' in the chromname {msachromname} will cause issues with figure generation bedfile output. Parsing chromname <yellow>{msachromname}</> -> <green>{new_chromname}</>",
                msachromname=msa._chrom_name,
                new_chromname=new_chromname,
            )
            msa._chrom_name = new_chromname

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
        if args.mode in ["all", "region-all"]:
            # Digest all of each MSAs
            msa.digest(cfg, indexes=None)
        elif args.mode == "region-only":
            # Digest only the wanted regions
            msa.digest(cfg, indexes=msa_name_to_indexes.get(msa_path.stem))  # type: ignore
            msa.remove_kmers_that_clash_with_regions()  # Remove kmers that clash with regions
        else:
            sys.exit(f"ERROR: {args.mode} is not a valid mode")

        # Log the digestion
        logger.info(
            "<blue>{msa_path}</>: <green>regions</> digested to <green>{num_fkmers}</> FKmers and <green>{num_rkmers}</> RKmers",
            msa_path=msa_path.stem,
            num_fkmers=len(msa.fkmers),
            num_rkmers=len(msa.rkmers),
        )

        # Generate all primerpairs
        msa.generate_primerpairs(cfg=cfg)
        logger.info(
            "<blue>{msa_path}</>: Generated <green>{num_pp}</> possible amplicons",
            msa_path=msa_path.stem,
            num_pp=len(msa._primerpairs),
        )

        if args.mode == "all":
            msa._primerpairs.sort(
                key=lambda x: msa.get_pp_entropy(x) ** 2 / sqrt(len(x.all_seqs())),  # type: ignore
                reverse=True,
            )
        if args.mode == "region-only":
            # Sort the primerpairs by the number of SNPs in the amplicon
            msa._primerpairs.sort(
                key=lambda pp: msa.get_pp_snp_score(pp),  # type: ignore
                reverse=True,
            )

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

    # Read in the inputbedfile if given
    if args.inputbedfile is not None:
        for bedprimerpair in read_in_bedprimerpairs(args.inputbedfile):
            # Add the primerpair into the panel, with a fake msaindex
            panel._add_primerpair(bedprimerpair, pool=bedprimerpair.pool, msa_index=-1)

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
    with open(OUTPUT_DIR / f"primer.bed", "w") as outfile:
        primer_bed_str = []
        for pp in panel._last_pp_added:
            # If there is an corrasponding msa
            ## Primers parsed in via bed do not have an msa_index
            if msa := msa_dict.get(pp.msa_index):
                chrom_name = msa._chrom_name
                primer_prefix = msa._uuid

                st = pp.__str__(
                    chrom_name,
                    primer_prefix,
                )
            else:
                # This Chrom name is not used in the bedfile
                # As Bedprimerpairs have there own name/prefix
                chrom_name = "scheme"
                primer_prefix = "scheme"
                st = pp.__str__(chrom_name, primer_prefix)

            primer_bed_str.append(st.strip())
        outfile.write("\n".join(primer_bed_str))

    # Write out the amplcion bed file
    with open(OUTPUT_DIR / "amplicon.bed", "w") as outfile:
        amp_bed_str = []
        for pp in panel._last_pp_added:
            if msa := msa_dict.get(pp.msa_index):
                chrom_name = msa._chrom_name
                primer_prefix = msa._uuid
            else:
                chrom_name = "scheme"
                primer_prefix = "scheme"

            region_list = []

            for regions in msa_name_to_regions.get(
                msa_index_to_name.get(pp.msa_index, -1), []
            ):
                if regions.start > pp.start and regions.stop < pp.end:
                    region_list.append(f"{regions.start}:{regions.stop}")

            regions_str = "\t".join(region_list)

            amp_bed_str.append(
                f"{chrom_name}\t{pp.start}\t{pp.end}\t{primer_prefix}_{pp.amplicon_number}\t{pp.pool + 1}\t{regions_str}"
            )
        outfile.write("\n".join(amp_bed_str))

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
    amp_md5 = hashlib.md5("\n".join(amp_bed_str).encode()).hexdigest()
    cfg["amplicon.bed.md5"] = amp_md5

    ## Read in the reference file and generate the hash
    with open(OUTPUT_DIR / "reference.fasta", "r") as reference_outfile:
        ref_md5 = hashlib.md5(reference_outfile.read().encode()).hexdigest()
    cfg["reference.fasta.md5"] = ref_md5

    # Write the config dict to file
    with open(OUTPUT_DIR / f"config.json", "w") as outfile:
        outfile.write(json.dumps(cfg, sort_keys=True))

    # Write the config dict to file
    # Add the bedfile to the cfg
    cfg["regionbedfile"] = str(args.regionbedfile)
    cfg["inputbedfile"] = str(args.inputbedfile)
    with open(OUTPUT_DIR / f"config.json", "w") as outfile:
        outfile.write(json.dumps(cfg, sort_keys=True))

    ## DO THIS LAST AS THIS CAN TAKE A LONG TIME
    # Writing plot data
    plot_data = generate_all_plotdata(
        list(msa_dict.values()),
        OUTPUT_DIR / "work",
        last_pp_added=panel._last_pp_added,
    )
    generate_all_plots(plot_data, OUTPUT_DIR)
