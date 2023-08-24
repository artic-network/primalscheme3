# Module imports
from primal_digest.thermo import *
from primal_digest.cli import cli
from primal_digest.config import config_dict
from primal_digest.classes import FKmer, RKmer, Scheme
from primal_digest.digestion import digest, generate_valid_primerpairs
from primal_digest.bedfiles import parse_bedfile, calc_median_bed_tm
from primal_digest.seq_functions import remove_end_insertion
from primal_digest.mismatches import MatchDB

import numpy as np
from Bio import SeqIO
import json

# Added
import hashlib
import sys
import pathlib
from loguru import logger

logger = logger.opt(colors=True)

"""
This is a test of a new dynamic digestion algo
"""


def main():
    args = cli()
    ARG_MSA = args.msa
    OUTPUT_DIR = pathlib.Path(args.output).absolute()

    cfg = config_dict
    # Primer Digestion settings
    cfg["primer_gc_min"] = args.primer_gc_min
    cfg["primer_gc_max"] = args.primer_gc_max
    cfg["primer_tm_min"] = args.primer_tm_min
    cfg["primer_tm_max"] = args.primer_tm_max
    cfg["dimerscore"] = args.dimerscore
    cfg["refname"] = args.refnames
    cfg["n_cores"] = args.cores
    cfg["output_prefix"] = args.prefix
    cfg["output_dir"] = str(OUTPUT_DIR)
    cfg["msa_paths"] = [str(x) for x in ARG_MSA]
    cfg["mismatches_self"] = args.mismatches_self
    cfg["mismatches_alt"] = args.mismatches_alt
    cfg["amplicon_size_max"] = args.ampliconsizemax
    cfg["amplicon_size_min"] = args.ampliconsizemin
    cfg["min_overlap"] = args.minoverlap
    cfg["force"] = args.force
    cfg["npools"] = args.npools
    cfg["reducekmers"] = args.reducekmers
    cfg["minbasefreq"] = args.minbasefreq
    cfg["msa_index_to_ref_name"] = {
        index: msa_name for index, msa_name in enumerate(cfg["refname"])
    }
    # Add the mismatch params to the cfg
    cfg["mismatch_fuzzy"] = True
    cfg["mismatch_kmersize"] = 20
    cfg["mismatch_product_size"] = args.ampliconsizemax

    # Add the bedfile path if given
    if args.bedfile:
        cfg["bedfile"] = str(args.bedfile)
    else:
        cfg["bedfile"] = False

    # See if the output dir already exsits
    if OUTPUT_DIR.is_dir():
        if not args.force:
            sys.exit(
                f"ERROR: {OUTPUT_DIR} already exists, please use --force to override"
            )
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

    # Read in the MSAs
    msa_list: list[np.ndarray] = []
    for msa_path in ARG_MSA:
        records = SeqIO.parse(msa_path, "fasta")
        align_array = np.array([record.seq.upper() for record in records], dtype=str)
        align_array = remove_end_insertion(align_array)
        msa_list.append(align_array)

        logger.info(
            "Read in MSA: <blue>{msa_path}</>\tseqs:<green>{msa_rows}</>\tcols:<green>{msa_cols}</>",
            msa_path=msa_path.name,
            msa_rows=align_array.shape[0],
            msa_cols=align_array.shape[1],
        )

    raw_msa_list: list[np.ndarray] = []
    for msa_path in ARG_MSA:
        records = SeqIO.parse(msa_path, "fasta")
        align_array = np.array([record.seq.upper() for record in records], dtype=str)
        raw_msa_list.append(align_array)

    # Generate the Kmers for each array in msa_list
    unique_f_r_msa: list[list[list[FKmer], list[RKmer]]] = []

    for msa_index, msa_array in enumerate(msa_list):
        mp_thermo_pass_fkmers, mp_thermo_pass_rkmers = digest(msa_array, cfg)

        # Append them to the msa
        unique_f_r_msa.append((mp_thermo_pass_fkmers, mp_thermo_pass_rkmers))
        logger.info(
            "<blue>{msa_path}</>: digested to <green>{num_fkmers}</> FKmers and <green>{num_rkmers}</> RKmers",
            msa_path=msa_path.name,
            num_fkmers=len(mp_thermo_pass_fkmers),
            num_rkmers=len(mp_thermo_pass_rkmers),
        )

    ## TODO use the matchdb to find mispriming

    msa_primerpairs = []
    ## Generate all valid primerpairs for each msa
    for msa_index, unique_fr_kmers in enumerate(unique_f_r_msa):
        # Generate all primerpairs then interaction check
        msa_primerpairs.append(
            generate_valid_primerpairs(
                unique_fr_kmers[0], unique_fr_kmers[1], cfg, msa_index=msa_index
            )
        )
        # Log some stats
        logger.info(
            "<blue>{msa_path}</>: Generated <green>{num_pp}</> possible amplicons",
            msa_path=msa_path.name,
            num_pp=len(msa_primerpairs[msa_index]),
        )

    for msa_index, primerpairs_in_msa in enumerate(msa_primerpairs):
        # Add the first primer, and if no primers can be added move to next msa
        if not scheme.add_first_primer_pair(primerpairs_in_msa, msa_index):
            continue

        while True:
            # Try and add an overlapping primer
            if scheme.try_ol_primerpairs(primerpairs_in_msa, cfg, msa_index):
                continue
            # Try and add a walking primer
            elif scheme.try_walk_primerpair(primerpairs_in_msa, cfg, msa_index):
                continue
            else:
                break

    logger.info("Writting output files")

    # Write primer bed file
    with open(OUTPUT_DIR / f"{cfg['output_prefix']}.primer.bed", "w") as outfile:
        primer_bed_str = []
        for pp in scheme.all_primers():
            st = pp.__str__(
                cfg["msa_index_to_ref_name"].get(pp.msa_index, "NA"),
                cfg["msa_index_to_ref_name"].get(pp.msa_index, "NA"),
            )
            primer_bed_str.append(st.strip())
        outfile.write("\n".join(primer_bed_str))

    # Write amplicon bed file
    with open(OUTPUT_DIR / f"{cfg['output_prefix']}.amplicon.bed", "w") as outfile:
        amp_bed_str = []
        for pp in scheme.all_primers():
            amp_bed_str.append(
                f"{cfg['msa_index_to_ref_name'].get(pp.msa_index, 'NA')}\t{pp.start}\t{pp.end}\tAMP_{pp.amplicon_number}\t{pp.pool + 1}"
            )
        outfile.write("\n".join(amp_bed_str))

    # Generate the bedfile hash, and add it into the config
    bed_md5 = hashlib.md5("\n".join(primer_bed_str).encode())
    cfg["md5_bed"] = bed_md5.hexdigest()

    # Generate the amplicon hash, and add it into the config
    bed_md5 = hashlib.md5("\n".join(amp_bed_str).encode())
    cfg["md5_amp"] = bed_md5.hexdigest()

    # Write the config dict to file
    with open(OUTPUT_DIR / f"config.json", "w") as outfile:
        outfile.write(json.dumps(cfg, sort_keys=True))


if __name__ == "__main__":
    main()
