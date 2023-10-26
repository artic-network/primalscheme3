# Interaction checker
from primaldimer_py import do_pools_interact_py

# Module imports
from primal_digest.thermo import *
from primal_digest.config import config_dict
from primal_digest.classes import Scheme, SchemeReturn
from primal_digest.bedfiles import (
    parse_bedfile,
    calc_median_bed_tm,
    BedPrimerPair,
    read_in_bedfile,
)
from primal_digest.mismatches import MatchDB
from primal_digest import __version__
from primal_digest.create_reports import generate_plot
from primal_digest.msa import MSA
from primal_digest.mapping import generate_consensus, generate_reference
from primal_digest.create_report_data import generate_data


# Added
import hashlib
import sys
import pathlib
from loguru import logger
import json
from Bio import SeqIO, SeqRecord, Seq
import shutil

logger = logger.opt(colors=True)

"""
This is a test of a new dynamic digestion algo
"""


def replace(args):
    """
    List all replacements primers
    """
    # Read in the config file
    with open(args.config, "r") as file:
        cfg: dict = json.load(file)

    # If more than two pools are given throw error
    if cfg["npools"] > 2:
        raise ValueError("ERROR: repair is only surported with two pools")

    # Create a mapping of chromname/referance to msa_index
    msa_chrom_to_index: dict[int:str] = {
        msa_data["msa_chromname"]: msa_index
        for msa_index, msa_data in cfg["msa_data"].items()
    }

    # Read in the bedfile
    bedprimerpairs: list[BedPrimerPair] = read_in_bedfile(args.primerbed)
    # Map each primer to an MSA index
    for primerpair in bedprimerpairs:
        primerpair.msa_index = msa_chrom_to_index[primerpair.chromname]

    bedprimerpairs.sort(key=lambda x: (x.chromname, x.amplicon_number))

    # Extract the stem from the primername
    try:
        prefix, ampliconnumber = args.primername.split("_")[:2]
        primerstem = f"{ampliconnumber}_{prefix}"
    except ValueError:
        raise ValueError(f"ERROR: {args.primername} cannot be parsed using _ as delim")

    # Find primernumber from bedfile
    wanted_pp = None
    for pp in bedprimerpairs:
        if pp.match_primer_stem(primerstem):
            wanted_pp = pp
    if wanted_pp is None:
        raise ValueError(f"ERROR: {args.primername} not found in bedfile")
    else:
        print(wanted_pp.__str__())

    # Read in the MSAs from config["msa_data"]
    # Only digest the wanted msa TODO
    msa_data = cfg["msa_data"][wanted_pp.msa_index]
    msa = MSA(
        name=msa_data["msa_name"],
        path=msa_data["msa_path"],
        msa_index=wanted_pp.msa_index,
        mapping=cfg["mapping"],
    )
    # Check the hashes match
    with open(msa.path, "rb") as f:
        msa_checksum = hashlib.file_digest(f, "md5").hexdigest()
    if msa_checksum != msa_data["msa_checksum"]:
        raise ValueError(f"ERROR: MSA checksums do not match.")

    # Digest the MSA into FKmers and RKmers
    msa.digest(cfg)

    # Generate all primerpairs then interaction check
    msa.generate_primerpairs(cfg)

    # Find the primers on either side of the wanted primer
    if wanted_pp.amplicon_number == 0:
        left_pp = None
        cov_start = wanted_pp.rprimer.start
    else:
        left_pp = [
            x
            for x in bedprimerpairs
            if x.amplicon_number == wanted_pp.amplicon_number - 1
        ][0]
        cov_start = left_pp.rprimer.start

    if wanted_pp.amplicon_number == max([x.amplicon_number for x in msa.primerpairs]):
        right_pp = None
        cov_end = wanted_pp.fprimer.end
    else:
        right_pp = [
            x
            for x in bedprimerpairs
            if x.amplicon_number == wanted_pp.amplicon_number + 1
        ][0]
        cov_end = right_pp.fprimer.end

    # Find primerpairs that span the gap
    spanning_primerpairs = [
        x
        for x in msa.primerpairs
        if x.fprimer.end < cov_start - cfg["min_overlap"]
        and x.rprimer.start > cov_end + cfg["min_overlap"]
    ]

    if not spanning_primerpairs:
        raise ValueError(f"ERROR: No spanning primers found")

    print(len(spanning_primerpairs))

    # Sort for number of primer pairs
    spanning_primerpairs.sort(key=lambda x: len(x.all_seqs()))

    # Find all primerpairs that in the same pool as the wanted primer
    same_pool_primerpairs = [x for x in bedprimerpairs if x.pool == wanted_pp.pool]
    seqs_in_same_pool = [
        seq for seq in (x.all_seqs() for x in same_pool_primerpairs) for seq in seq
    ]

    # Find all primerpairs that in the same pool and same msa as the wanted primer
    spanning_pool_primerpairs = [
        x
        for x in msa.primerpairs
        if x.pool == wanted_pp.pool and x.msa_index == wanted_pp.msa_index
    ]

    for pos_primerpair in spanning_primerpairs:
        # Make sure the new primerpair doesnt contain the primers we want to replace
        if (
            pos_primerpair.fprimer == wanted_pp.fprimer
            or pos_primerpair.rprimer == wanted_pp.rprimer
        ):
            continue

        # Check for all clashes
        clash = False
        # If they overlap
        for current_pp in spanning_pool_primerpairs:
            # If clash with any skip
            if range(
                max(pos_primerpair.start, current_pp.start),
                min(pos_primerpair.end, current_pp.end) + 1,
            ):
                clash = True
                break
        if clash:
            continue

        # If they dont overlap
        # Check for interactions
        if do_pools_interact_py(
            list(pos_primerpair.all_seqs()), seqs_in_same_pool, cfg["dimerscore"]
        ):
            continue

        # primer passes!
        print(
            f"fprimer: {pos_primerpair.fprimer.end}\trprimer: {pos_primerpair.rprimer.start}"
        )
        print(pos_primerpair.__str__(msa._chrom_name, msa._uuid))

    ## TODO check for clashes between the last primer in the same pool and the first primer in the spanning primer


def create(args):
    ARG_MSA = args.msa
    OUTPUT_DIR = pathlib.Path(args.output).absolute()

    cfg = config_dict

    # Add version to config
    cfg["primaldigest_version"] = f"primaldigest:{__version__}"
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
    cfg["min_overlap"] = args.minoverlap
    cfg["force"] = args.force
    cfg["npools"] = args.npools
    cfg["reducekmers"] = args.reducekmers
    cfg["minbasefreq"] = args.minbasefreq

    # Add the mismatch params to the cfg
    cfg["mismatch_fuzzy"] = True
    cfg["mismatch_kmersize"] = 20
    cfg["mismatch_product_size"] = args.ampliconsizemax

    # Add plots to the cfg
    cfg["plot"] = args.plot
    cfg["disable_progress_bar"] = False

    # Add the mapping to the cfg
    cfg["mapping"] = args.mapping

    # Add circular to the cfg
    cfg["circular"] = args.circular

    # Add the backtrack to the cfg
    cfg["backtrack"] = args.backtrack

    # Add the bedfile path if given
    if args.bedfile:
        cfg["bedfile"] = str(args.bedfile)
    else:
        cfg["bedfile"] = False

    # See if the output dir already exsits
    if OUTPUT_DIR.is_dir() and not args.force:
        sys.exit(f"ERROR: {OUTPUT_DIR} already exists, please use --force to override")

    # Create the output dir and a work subdir
    pathlib.Path.mkdir(OUTPUT_DIR, exist_ok=True)
    pathlib.Path.mkdir(OUTPUT_DIR / "work", exist_ok=True)

    ## Set up the loggers
    logger.remove()  # Remove default stderr logger
    # Add the deep log file
    logger.add(
        OUTPUT_DIR / "work/file.log",
        colorize=False,
        format="{time:YYYY-MM-DD at HH:mm:ss} | {level} | {message}",
        enqueue=True,
    )
    # Add the nice stdout progress
    logger.add(sys.stdout, colorize=True, format="{message}", level="INFO")

    # Create the mismatch db
    logger.info(
        "Creating the Mismatch Database",
    )
    mismatch_db = MatchDB(OUTPUT_DIR / "work/mismatch", ARG_MSA, cfg["primer_size_min"])
    logger.info(
        "<green>Created:</> {path}",
        path=f"{OUTPUT_DIR}/work/mismatch.db",
    )

    # Create the scheme object early
    scheme = Scheme(cfg=cfg, matchDB=mismatch_db)

    # If the bedfile flag is given add the primers into the scheme
    if args.bedfile:
        current_pools = parse_bedfile(args.bedfile, cfg["npools"])
        # Assign the bedfile generated pool to the scheme in a hacky way
        scheme._pools = current_pools
        primer_tms = calc_median_bed_tm(current_pools, cfg)
        logger.info(
            "Read in bedfile: <blue>{msa_path}</>: <green>{num_pp}</> Primers with median Tm of <green>{tm}</>",
            msa_path=args.bedfile.name,
            num_pp=len(
                [
                    primer
                    for primer in (pool for pool in current_pools)
                    for primer in primer
                ]
            ),
            tm=round(primer_tms[1], 2),
        )

        # If the bedfile tm is dif to the config Throw error
        if (
            primer_tms[0] < cfg["primer_tm_min"] - 5
            or primer_tms[2] > cfg["primer_tm_max"] + 5
        ):
            logger.warning(
                "Tm in bedfile <green>{tm_min}</>:<green>{tm_median}</>:<green>{tm_max}</> (min, median, max) fall outside range of <green>{conf_tm_min}</> to <green>{conf_tm_max}</>)",
                tm_min=round(primer_tms[0], 1),
                tm_median=round(primer_tms[1], 1),
                tm_max=round(primer_tms[2], 1),
                conf_tm_min=cfg["primer_tm_min"],
                conf_tm_max=cfg["primer_tm_max"],
            )

    # Create a dict full of msa data
    msa_data = {}

    # Read in the MSAs
    msa_dict: dict[int:MSA] = {}
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
        msa = MSA(
            name=local_msa_path.stem,
            path=msa_path,
            msa_index=msa_index,
            mapping=cfg["mapping"],
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
        msa_data[msa_index]["msa_path"] = str(
            "work/" + msa_path.name
        )  # Write localpath
        msa_data[msa_index]["msa_chromname"] = msa._chrom_name
        msa_data[msa_index]["msa_uuid"] = msa._uuid

        logger.info(
            "Read in MSA: <blue>{msa_path}</>\tseqs:<green>{msa_rows}</>\tcols:<green>{msa_cols}</>",
            msa_path=msa.name,
            msa_rows=msa.array.shape[0],
            msa_cols=msa.array.shape[1],
        )

        # Digest the MSA into FKmers and RKmers
        msa.digest(cfg)
        logger.info(
            "<blue>{msa_path}</>: digested to <green>{num_fkmers}</> FKmers and <green>{num_rkmers}</> RKmers",
            msa_path=msa.name,
            num_fkmers=len(msa.fkmers),
            num_rkmers=len(msa.rkmers),
        )

        # Generate all primerpairs then interaction check
        msa.generate_primerpairs(cfg)
        logger.info(
            "<blue>{msa_path}</>: Generated <green>{num_pp}</> possible amplicons",
            msa_path=msa.name,
            num_pp=len(msa.primerpairs),
        )

        # Add the msa to the scheme
        msa_dict[msa_index] = msa

    # Add MSA data into cfg
    cfg["msa_data"] = msa_data

    # Start the Scheme generation
    for msa_index, msa in msa_dict.items():
        # Add the first primer, and if no primers can be added move to next msa
        if (
            scheme.add_first_primer_pair(msa.primerpairs, msa_index)
            == SchemeReturn.ADDED_FIRST_PRIMERPAIR
        ):
            logger.info(
                "Added <green>first</> amplicon for <blue>{msa_name}</>: {primer_start}\t{primer_end}\t{primer_pool}",
                primer_start=scheme._last_pp_added[-1].start,
                primer_end=scheme._last_pp_added[-1].end,
                primer_pool=scheme._last_pp_added[-1].pool + 1,
                msa_name=msa.name,
            )
        else:
            logger.warning(
                "No valid primers found for <blue>{msa_name}</>",
                msa_name=msa.name,
            )
            continue

        while True:
            # Try and add an overlapping primer
            ol_result = scheme.try_ol_primerpairs(msa.primerpairs, msa_index)
            if ol_result == SchemeReturn.ADDED_OL_PRIMERPAIR:
                logger.info(
                    "Added <green>overlapping</> amplicon for <blue>{msa_name}</>: {primer_start}\t{primer_end}\t{primer_pool}",
                    primer_start=scheme._last_pp_added[-1].start,
                    primer_end=scheme._last_pp_added[-1].end,
                    primer_pool=scheme._last_pp_added[-1].pool + 1,
                    msa_name=msa.name,
                )
                continue
            elif ol_result == SchemeReturn.NO_OL_PRIMERPAIR:
                pass  # Do nothing move on to next step

            # Try to backtrack
            if cfg["backtrack"]:
                backtrack_result = scheme.try_backtrack(msa.primerpairs, msa_index)
                # If successful log and continue
                if backtrack_result == SchemeReturn.ADDED_BACKTRACKED:
                    logger.info(
                        "Backtracking allowed <green>overlapping</> amplicon to be added for <blue>{msa_name}</>: {primer_start}\t{primer_end}\t{primer_pool}",
                        primer_start=scheme._last_pp_added[-1].start,
                        primer_end=scheme._last_pp_added[-1].end,
                        primer_pool=scheme._last_pp_added[-1].pool + 1,
                        msa_name=msa.name,
                    )
                    continue
                # If cannot backtrack and backtrack is enabled, log then move on
                elif backtrack_result == SchemeReturn.NO_BACKTRACK:
                    logger.info(
                        "Could not backtrack for <blue>{msa_name}</>",
                        msa_name=msa.name,
                    )

            # Try and add a walking primer
            if (
                scheme.try_walk_primerpair(msa.primerpairs, msa_index)
                == SchemeReturn.ADDED_WALK_PRIMERPAIR
            ):
                logger.info(
                    "Added <yellow>walking</> amplicon for <blue>{msa_name}</>: {primer_start}\t{primer_end}\t{primer_pool}",
                    primer_start=scheme._last_pp_added[-1].start,
                    primer_end=scheme._last_pp_added[-1].end,
                    primer_pool=scheme._last_pp_added[-1].pool + 1,
                    msa_name=msa.name,
                )
            else:
                break

        if cfg["circular"] and scheme.try_circular(msa) == SchemeReturn.ADDED_CIRULAR:
            logger.info(
                "Added <green>circular</> amplicon for <blue>{msa_name}</>: {primer_start}\t{primer_end}\t{primer_pool}",
                primer_start=scheme._last_pp_added[-1].start,
                primer_end=scheme._last_pp_added[-1].end,
                primer_pool=scheme._last_pp_added[-1].pool + 1,
                msa_name=msa.name,
            )

    logger.info("Writting output files")

    # Write primer bed file
    with open(OUTPUT_DIR / "primer.bed", "w") as outfile:
        primer_bed_str = []
        for pp in scheme.all_primers():
            # If there is an corrasponding msa
            ## Primers parsed in via bed do not have an msa_index
            if msa := msa_dict.get(pp.msa_index):
                chrom_name = msa._chrom_name
                primer_prefix = msa._uuid
            else:
                # This Chrom name is not used in the bedfile
                # As BedLines have there own name/prefix
                chrom_name = "scheme"
                primer_prefix = "scheme"

            st = pp.__str__(
                chrom_name,
                primer_prefix,
            )
            primer_bed_str.append(st.strip())
        outfile.write("\n".join(primer_bed_str))

    # Write amplicon bed file
    with open(OUTPUT_DIR / "amplicon.bed", "w") as outfile:
        amp_bed_str = []
        for pp in scheme.all_primers():
            if msa := msa_dict.get(pp.msa_index):
                chrom_name = msa._chrom_name
                primer_prefix = msa._uuid
            else:
                chrom_name = "scheme"
                primer_prefix = "scheme"

            amp_bed_str.append(
                f"{chrom_name}\t{pp.start}\t{pp.end}\t{primer_prefix}_{pp.amplicon_number}\t{pp.pool + 1}"
            )
        outfile.write("\n".join(amp_bed_str))

    # Create the fancy plots
    if cfg["plot"]:
        for msa in msa_dict.values():
            generate_plot(msa, scheme._pools, OUTPUT_DIR)

    # Writing plot data
    for msa in msa_dict.values():
        generate_data(msa, OUTPUT_DIR / "work")

    # Write all the consensus sequences to a single file
    with open(OUTPUT_DIR / "reference.fasta", "w") as reference_outfile:
        reference_records = []
        if cfg["mapping"] == "first":
            for msa in msa_dict.values():
                reference_records.append(
                    SeqRecord.SeqRecord(
                        seq=Seq.Seq(generate_reference(msa.array)),
                        id=msa._chrom_name,
                        description="",
                    )
                )
        elif cfg["mapping"] == "consensus":
            for msa in msa_dict.values():
                reference_records.append(
                    SeqRecord.SeqRecord(
                        seq=Seq.Seq(generate_consensus(msa.array)),
                        id=msa._chrom_name,
                        description="",
                    )
                )
        SeqIO.write(reference_records, reference_outfile, "fasta")

    # Create all hashes
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
