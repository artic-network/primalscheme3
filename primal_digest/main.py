# Module imports
from primal_digest.thermo import *
from primal_digest.cli import cli
from primal_digest.config import config_dict
from primal_digest.classes import Scheme
from primal_digest.bedfiles import parse_bedfile, calc_median_bed_tm
from primal_digest.mismatches import MatchDB
from primal_digest import __version__
from primal_digest.create_reports import generate_plot
from primal_digest.msa import MSA
from primal_digest.mapping import generate_consensus


# Added
import hashlib
import sys
import pathlib
from loguru import logger
import json
from Bio import SeqIO, SeqRecord, Seq

logger = logger.opt(colors=True)

"""
This is a test of a new dynamic digestion algo
"""


def main():
    args = cli()
    ARG_MSA = args.msa
    OUTPUT_DIR = pathlib.Path(args.output).absolute()

    cfg = config_dict

    # Add version to config
    cfg["primal_digest_version"] = __version__
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

        # Create MSA checksum
        with open(msa_path, "rb") as f:
            msa_data[msa_index]["msa_checksum"] = hashlib.file_digest(
                f, "md5"
            ).hexdigest()

        # Read in the MSA
        msa = MSA(
            name=msa_path.stem,
            path=msa_path,
            msa_index=msa_index,
            mapping=cfg["mapping"],
        )

        # Add some msa data to the dict
        msa_data[msa_index]["msa_name"] = msa.name
        msa_data[msa_index]["msa_path"] = str(msa_path.absolute())
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
        if scheme.add_first_primer_pair(msa.primerpairs, msa_index):
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
            if scheme.try_ol_primerpairs(msa.primerpairs, msa_index):
                logger.info(
                    "Added <green>overlapping</> amplicon for <blue>{msa_name}</>: {primer_start}\t{primer_end}\t{primer_pool}",
                    primer_start=scheme._last_pp_added[-1].start,
                    primer_end=scheme._last_pp_added[-1].end,
                    primer_pool=scheme._last_pp_added[-1].pool + 1,
                    msa_name=msa.name,
                )
                continue
            # Try to backtrack
            # elif scheme.try_backtrack(msa.primerpairs, msa_index):
            #     logger.info(
            #         "Added <yellow>backtracking</> amplicon for <blue>{msa_name}</>: {primer_start}\t{primer_end}\t{primer_pool}",
            #         primer_start=scheme._last_pp_added[-1].start,
            #         primer_end=scheme._last_pp_added[-1].end,
            #         primer_pool=scheme._last_pp_added[-1].pool + 1,
            #         msa_name=msa.name,
            #     )
            #     continue
            # Try and add a walking primer
            elif scheme.try_walk_primerpair(msa.primerpairs, msa_index):
                logger.info(
                    "Added <yellow>walking</> amplicon for <blue>{msa_name}</>: {primer_start}\t{primer_end}\t{primer_pool}",
                    primer_start=scheme._last_pp_added[-1].start,
                    primer_end=scheme._last_pp_added[-1].end,
                    primer_pool=scheme._last_pp_added[-1].pool + 1,
                    msa_name=msa.name,
                )
            else:
                break

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
                # As BedPrimers have there own name/prefix
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
            if pp_name := msa_dict.get(pp.msa_index):
                pp_name = pp_name.name
            else:
                pp_name = pp.ref
            amp_bed_str.append(
                f"{pp_name}\t{pp.start}\t{pp.end}\tAMP_{pp.amplicon_number}\t{pp.pool + 1}"
            )
        outfile.write("\n".join(amp_bed_str))

    # Create the fancy plots
    if cfg["plot"]:
        for msa in msa_dict.values():
            generate_plot(msa, scheme._pools, OUTPUT_DIR)

    # Write all the consensus sequences to a single file
    with open(OUTPUT_DIR / "referance.fasta", "w") as referance_outfile:
        referance_records = []
        for msa in msa_dict.values():
            referance_records.append(
                SeqRecord.SeqRecord(
                    seq=Seq.Seq(generate_consensus(msa.array)),
                    id=msa._chrom_name,
                )
            )
        SeqIO.write(referance_records, referance_outfile, "fasta")

    # Create all hashes
    ## Generate the bedfile hash, and add it into the config
    primer_md5 = hashlib.md5("\n".join(primer_bed_str).encode()).hexdigest()
    cfg["primer_checksum"] = primer_md5

    ## Generate the amplicon hash, and add it into the config
    amp_md5 = hashlib.md5("\n".join(amp_bed_str).encode()).hexdigest()
    cfg["amplicon_checksum"] = amp_md5

    ## Read in the referance file and generate the hash
    with open(OUTPUT_DIR / "referance.fasta", "r") as referance_outfile:
        ref_md5 = hashlib.md5(referance_outfile.read().encode()).hexdigest()
    cfg["referance_checksum"] = ref_md5

    ## Write the info.json file
    info_dict = {
        "algorithm_version": f"primaldigest:{__version__}",
        "scheme_version": "v0.0.0",
        "ampliconsize": cfg["amplicon_size_max"],
        "primer_checksum": cfg["primer_checksum"],
        "referance_checksum": cfg["referance_checksum"],
        "amplicon_checksum": cfg["amplicon_checksum"],
    }

    # Write the config dict to file
    with open(OUTPUT_DIR / f"config.json", "w") as outfile:
        outfile.write(json.dumps(cfg, sort_keys=True))


if __name__ == "__main__":
    main()
