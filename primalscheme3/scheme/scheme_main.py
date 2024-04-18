import hashlib
import json
import pathlib
import shutil
import sys

from Bio import Seq, SeqIO, SeqRecord

# Interaction checker
from primaldimer_py import do_pools_interact_py  # type: ignore


from primalscheme3.__init__ import __version__
from primalscheme3.core.bedfiles import (
    read_in_bedprimerpairs,
)
from primalscheme3.core.config import config_dict
from primalscheme3.core.create_report_data import (
    generate_all_plotdata,
)
from primalscheme3.core.create_reports import generate_all_plots
from primalscheme3.core.logger import setup_loger
from primalscheme3.core.mapping import generate_consensus, generate_reference
from primalscheme3.core.mismatches import MatchDB
from primalscheme3.core.msa import MSA
from primalscheme3.scheme.classes import Scheme, SchemeReturn


def schemereplace(
    config: pathlib.Path,
    ampliconsizemax: int,
    ampliconsizemin: int,
    primerbed: pathlib.Path,
    primername: str,
    msapath: pathlib.Path,
):
    """
    List all replacements primers
    """
    # Read in the config file
    with open(config) as file:
        cfg: dict = json.load(file)

    # Update the amplicon size if it is provided
    if ampliconsizemax:
        print(
            f"Updating amplicon size max to {ampliconsizemax} and min to {ampliconsizemin}"
        )

        cfg["amplicon_size_max"] = ampliconsizemax
        cfg["amplicon_size_min"] = ampliconsizemin

    # Add the newer primer_max_walk to the config if not already there
    if "primer_max_walk" not in cfg:
        cfg["primer_max_walk"] = 100

    # If more than two pools are given throw error
    if cfg["npools"] > 2:
        raise ValueError("ERROR: repair is only surported with two pools")

    # Create a mapping of chromname/referance to msa_index
    msa_chrom_to_index: dict[str, int] = {
        msa_data["msa_chromname"]: msa_index
        for msa_index, msa_data in cfg["msa_data"].items()
    }

    # Read in the bedfile
    bedprimerpairs, headers = read_in_bedprimerpairs(primerbed)
    # Map each primer to an MSA index
    for primerpair in bedprimerpairs:
        msa_index = msa_chrom_to_index.get(str(primerpair.chrom_name), None)
        if msa_index is not None:
            primerpair.msa_index = msa_index
        elif cfg["bedfile"]:
            # This case can happen when primers are added to the scheme via --bedfile.
            # Set the msa index to -1
            primerpair.msa_index = -1
        else:
            raise ValueError(f"ERROR: {primerpair.chrom_name} not found in MSA data")

    bedprimerpairs.sort(key=lambda x: (x.chrom_name, x.amplicon_number))

    # Extract the stem from the primername
    try:
        prefix, ampliconnumber = primername.split("_")[:2]
        primerstem = f"{ampliconnumber}_{prefix}"
    except ValueError:
        raise ValueError(f"ERROR: {primername} cannot be parsed using _ as delim")

    # Find primernumber from bedfile
    wanted_pp = None
    for pp in bedprimerpairs:
        if pp.match_primer_stem(primerstem):
            wanted_pp = pp
    if wanted_pp is None:
        raise ValueError(f"ERROR: {primername} not found in bedfile")
    else:
        print(wanted_pp.__str__())

    # Read in the MSAs from config["msa_data"]
    msa_data = cfg["msa_data"].get(wanted_pp.msa_index)
    if msa_data is None:
        if wanted_pp.msa_index == -1:
            raise ValueError(f"ERROR: The Primer {primername} was added via --bedfile")
        else:
            raise ValueError(f"ERROR: MSA index {wanted_pp.msa_index} not found")

    msa = MSA(
        name=msa_data["msa_name"],
        path=msapath,
        msa_index=wanted_pp.msa_index,
        mapping=cfg["mapping"],
    )
    # Check the hashes match
    with open(msa.path, "rb") as f:
        msa_checksum = hashlib.file_digest(f, "md5").hexdigest()
    if msa_checksum != msa_data["msa_checksum"]:
        raise ValueError("ERROR: MSA checksums do not match.")
    else:
        print("MSA checksums match")

    # Targeted digestion leads to a mismatch of the indexes.
    # Digest the MSA into FKmers and RKmers
    msa.digest(cfg)  ## Primer are remapped at this point.
    print(f"Found {len(msa.fkmers)} FKmers and {len(msa.rkmers)} RKmers")

    # Generate all primerpairs then interaction check
    msa.generate_primerpairs(
        amplicon_size_max=cfg["amplicon_size_max"],
        amplicon_size_min=cfg["amplicon_size_min"],
        dimerscore=cfg["dimerscore"],
    )

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

    if len(spanning_primerpairs) > 0:
        print(f"Spanning amplicons found: {len(spanning_primerpairs)}")
    else:
        raise ValueError("ERROR: No spanning primers found")

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

    accepted_primerpairs = []
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

        pos_primerpair.pool = wanted_pp.pool
        pos_primerpair.amplicon_number = wanted_pp.amplicon_number

        accepted_primerpairs.append(pos_primerpair)

    accepted_primerpairs.sort(key=lambda a: a.fprimer.end)

    print(f"Found {len(accepted_primerpairs)} valid replacement amplicons")
    for pp_number, pp in enumerate(accepted_primerpairs, 1):
        print(f"Amplicon {pp_number}: {len(pp.all_seqs())} Primers")
        print(pp.to_bed())


def schemecreate(
    argmsa: list[pathlib.Path],
    primer_gc_min: float,
    primer_gc_max: float,
    primer_tm_min: float,
    primer_tm_max: float,
    dimerscore: float,
    ncores: int,
    output_dir: pathlib.Path,
    ampliconsizemax: int,
    ampliconsizemin: int,
    minoverlap: int,
    npools: int,
    reducekmers: bool,
    minbasefreq: int,
    circular: bool,
    backtrack: bool,
    ignore_n: bool,
    bedfile: pathlib.Path | None = None,
    force: bool = False,
    mapping: str = "first",
    no_plot: bool = False,
):
    ARG_MSA = argmsa
    OUTPUT_DIR = pathlib.Path(output_dir).absolute()  # Keep absolute path

    cfg = config_dict

    # Add version to config
    cfg["algorithmversion"] = f"primalscheme3:{__version__}"
    cfg["primerclass"] = "primerschemes"
    # Primer Digestion settings
    cfg["primer_gc_min"] = primer_gc_min
    cfg["primer_gc_max"] = primer_gc_max
    cfg["primer_tm_min"] = primer_tm_min
    cfg["primer_tm_max"] = primer_tm_max
    cfg["dimerscore"] = dimerscore
    cfg["n_cores"] = ncores
    cfg["output_dir"] = str(OUTPUT_DIR)  # Write localpath
    cfg["amplicon_size_max"] = ampliconsizemax
    cfg["amplicon_size_min"] = ampliconsizemin
    cfg["min_overlap"] = minoverlap
    cfg["force"] = force
    cfg["npools"] = npools
    cfg["reducekmers"] = reducekmers
    cfg["minbasefreq"] = minbasefreq

    # Add the mismatch params to the cfg
    cfg["mismatch_fuzzy"] = True
    cfg["mismatch_kmersize"] = cfg["primer_size_min"]
    cfg["mismatch_product_size"] = ampliconsizemax

    # Add plots to the cfg
    cfg["plot"] = no_plot
    cfg["disable_progress_bar"] = False

    # Add the mapping to the cfg
    cfg["mapping"] = mapping

    # Add circular to the cfg
    cfg["circular"] = circular

    # Add the backtrack to the cfg
    cfg["backtrack"] = backtrack

    # Add the bedfile path if given
    if bedfile is not None:
        cfg["bedfile"] = str(bedfile)
    else:
        cfg["bedfile"] = False

    # Add ignore_n to the cfg
    cfg["ignore_n"] = ignore_n

    # See if the output dir already exsits
    if OUTPUT_DIR.is_dir() and not force:
        sys.exit(f"ERROR: {OUTPUT_DIR} already exists, please use --force to override")

    # Create the output dir and a work subdir
    pathlib.Path.mkdir(OUTPUT_DIR, exist_ok=True)
    pathlib.Path.mkdir(OUTPUT_DIR / "work", exist_ok=True)

    # Set up the logger
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
        path=f"{OUTPUT_DIR.relative_to(OUTPUT_DIR.parent)}/work/mismatch.db",
    )

    # Create the scheme object early
    scheme = Scheme(cfg=cfg, matchDB=mismatch_db)

    # If the bedfile flag is given add the primers into the scheme
    if bedfile is not None:
        bedprimerpairs, _headers = read_in_bedprimerpairs(bedfile)
        # Check the number of pools in the given bedfile, is less or equal to npools arg
        pools_in_bed = {primer.pool for primer in bedprimerpairs}
        if max(pools_in_bed) > cfg["npools"]:
            sys.exit(
                f"ERROR: The number of pools in the bedfile is greater than --npools: {max(pools_in_bed)} > {cfg['npools']}"
            )

        # Calculate the bedfile tm
        primer_tms = [
            tm for tm in (pp.calc_tm(cfg) for pp in bedprimerpairs) for tm in tm
        ]
        logger.info(
            "Read in bedfile: <blue>{msa_path}</>: <green>{num_pp}</> PrimersPairs with mean Tm of <green>{tm}</>",
            msa_path=bedfile.name,
            num_pp=len(bedprimerpairs),
            tm=round(sum(primer_tms) / len(primer_tms), 2),
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
    msa_dict: dict[int, MSA] = {}
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

        if len(msa.fkmers) == 0 or len(msa.rkmers) == 0:
            logger.critical(
                "No valid FKmers or RKmers found for <blue>{msa_name}</>",
                msa_name=msa.name,
            )
            sys.exit(1)

        # Generate all primerpairs then interaction check
        msa.generate_primerpairs(
            amplicon_size_max=cfg["amplicon_size_max"],
            amplicon_size_min=cfg["amplicon_size_min"],
            dimerscore=cfg["dimerscore"],
        )
        logger.info(
            "<blue>{msa_path}</>: Generated <green>{num_pp}</> possible amplicons",
            msa_path=msa.name,
            num_pp=len(msa.primerpairs),
        )

        if len(msa.primerpairs) == 0:
            logger.critical(
                "No valid primers found for <blue>{msa_name}</>",
                msa_name=msa.name,
            )
            sys.exit(1)

        # Add the msa to the scheme
        msa_dict[msa_index] = msa

    # Add MSA data into cfg
    cfg["msa_data"] = msa_data

    msa_chrom_to_index: dict[str, int] = {
        msa._chrom_name: msa_index for msa_index, msa in msa_dict.items()
    }
    # Add the bedprimerpairs into the scheme
    if bedfile is not None and bedprimerpairs:
        for pp in bedprimerpairs:
            # Map the primerpair to the msa via chromname
            pp.msa_index = msa_chrom_to_index.get(pp.chrom_name, -1)  # type: ignore
            scheme.add_primer_pair_to_pool(pp, pp.pool, pp.msa_index)

    # Start the Scheme generation
    for msa_index, msa in msa_dict.items():
        while True:
            match scheme.try_ol_primerpairs(msa.primerpairs, msa_index):
                case SchemeReturn.ADDED_OL_PRIMERPAIR:
                    logger.info(
                        "Added <green>overlapping</> amplicon for <blue>{msa_name}</>: {primer_start}\t{primer_end}\t{primer_pool}",
                        primer_start=scheme._last_pp_added[-1].start,
                        primer_end=scheme._last_pp_added[-1].end,
                        primer_pool=scheme._last_pp_added[-1].pool + 1,
                        msa_name=msa.name,
                    )
                    continue
                case SchemeReturn.ADDED_FIRST_PRIMERPAIR:
                    logger.info(
                        "Added <green>first</> amplicon for <blue>{msa_name}</>: {primer_start}\t{primer_end}\t{primer_pool}",
                        primer_start=scheme._last_pp_added[-1].start,
                        primer_end=scheme._last_pp_added[-1].end,
                        primer_pool=scheme._last_pp_added[-1].pool + 1,
                        msa_name=msa.name,
                    )
                    continue
                case SchemeReturn.NO_OL_PRIMERPAIR:
                    pass  # Do nothing move on to next step
                case SchemeReturn.NO_FIRST_PRIMERPAIR:
                    logger.warning(
                        "No valid primers found for <blue>{msa_name}</>",
                        msa_name=msa.name,
                    )
                    break

            # Try to backtrack
            if cfg["backtrack"]:
                logger.info("Backtracking...")
                match scheme.try_backtrack(msa.primerpairs, msa_index):
                    case SchemeReturn.ADDED_BACKTRACKED:
                        logger.info(
                            "Backtracking allowed <green>overlapping</> amplicon to be added for <blue>{msa_name}</>: {primer_start}\t{primer_end}\t{primer_pool}",
                            primer_start=scheme._last_pp_added[-1].start,
                            primer_end=scheme._last_pp_added[-1].end,
                            primer_pool=scheme._last_pp_added[-1].pool + 1,
                            msa_name=msa.name,
                        )
                    case SchemeReturn.NO_BACKTRACK:
                        logger.info(
                            "Could not backtrack for <blue>{msa_name}</>",
                            msa_name=msa.name,
                        )
                        pass  # Do nothing move on to next step

            # Try and add a walking primer
            match scheme.try_walk_primerpair(msa.primerpairs, msa_index):
                case SchemeReturn.ADDED_WALK_PRIMERPAIR:
                    logger.info(
                        "Added <yellow>walking</> amplicon for <blue>{msa_name}</>: {primer_start}\t{primer_end}\t{primer_pool}",
                        primer_start=scheme._last_pp_added[-1].start,
                        primer_end=scheme._last_pp_added[-1].end,
                        primer_pool=scheme._last_pp_added[-1].pool + 1,
                        msa_name=msa.name,
                    )
                case _:
                    break

        if cfg["circular"]:
            match scheme.try_circular(msa):
                case SchemeReturn.ADDED_CIRULAR:
                    logger.info(
                        "Added <green>circular</> amplicon for <blue>{msa_name}</>: {primer_start}\t{primer_end}\t{primer_pool}",
                        primer_start=scheme._last_pp_added[-1].start,
                        primer_end=scheme._last_pp_added[-1].end,
                        primer_pool=scheme._last_pp_added[-1].pool + 1,
                        msa_name=msa.name,
                    )
                case SchemeReturn.NO_CIRCULAR:
                    logger.info(
                        "No <red>circular</> amplicon for <blue>{msa_name}</>",
                        msa_name=msa.name,
                    )

    logger.info("Writting output files")

    # Write primer bed file
    with open(OUTPUT_DIR / "primer.bed", "w") as outfile:
        primer_bed_str = scheme.to_bed()
        outfile.write(primer_bed_str)

    # Write amplicon bed file
    with open(OUTPUT_DIR / "amplicon.bed", "w") as outfile:
        amp_bed_str = scheme.to_amplicons(trim_primers=False)
        outfile.write(amp_bed_str)
    with open(OUTPUT_DIR / "primertrim.amplicon.bed", "w") as outfile:
        outfile.write(scheme.to_amplicons(trim_primers=True))

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
    amp_md5 = hashlib.md5(amp_bed_str.encode()).hexdigest()
    cfg["amplicon.bed.md5"] = amp_md5

    ## Read in the reference file and generate the hash
    with open(OUTPUT_DIR / "reference.fasta") as reference_outfile:
        ref_md5 = hashlib.md5(reference_outfile.read().encode()).hexdigest()
    cfg["reference.fasta.md5"] = ref_md5

    # Write the config dict to file
    with open(OUTPUT_DIR / "config.json", "w") as outfile:
        outfile.write(json.dumps(cfg, sort_keys=True))

    ## DO THIS LAST AS THIS CAN TAKE A LONG TIME
    # Writing plot data
    plot_data = generate_all_plotdata(
        list(msa_dict.values()),
        OUTPUT_DIR / "work",
        last_pp_added=scheme._last_pp_added,
    )
    generate_all_plots(plot_data, OUTPUT_DIR)

    logger.info("Completed Successfully")
