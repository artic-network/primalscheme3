# Module imports
from primal_panel.classes2 import Region, Scheme
from primal_panel.cli import cli

# Submodule imports
from primal_digest.digestion import digest
from primal_digest.seq_functions import remove_end_insertion
from primal_digest.classes import PrimerPair
from primal_digest.config import config_dict, AMBIGUOUS_DNA
from primal_digest.primaldimer_py import do_seqs_interact_py, do_pools_interact_py

# Use the new interaction checker
import random


# General import
import pathlib
import sys
from Bio import SeqIO
import numpy as np
from loguru import logger
from uuid import uuid4
import json


logger = logger.opt(colors=True)


def network_opto(scheme):
    iteration = 1
    while scheme.total_interactions() > 0:
        logger.info(
            "Iteration: {i}\tScore: <blue>{score}</>\r",
            score=scheme.total_interactions(),
            i=iteration,
        )
        scheme.replace_worst_primerpair()
        iteration += 1
    scheme.view()
    print(scheme.to_bed_format().strip())


def does_overlap(pp: PrimerPair, covered: tuple[int, int]) -> bool:
    """
    Will check for overlaps between all given primers
    False means no overlap (good)
    """
    pp_start = pp.start()
    pp_end = pp.end()

    for s, e in covered:
        if s <= pp_start <= e or s <= pp_end <= e:
            # Check if new primers fall into covered region
            return True
        elif pp_start <= s <= pp_end:
            # Check of the new primers overlaps a covered region
            return True

    return False


def main():
    # Create a run UUID
    run_uuid = uuid4()

    logger.remove()  # Remove default stderr logger
    logger.add(sys.stderr, colorize=True, format="{message}", level="INFO")

    logger.info("Starting run <blue>{ref}</>", ref=run_uuid)

    args = cli()
    ARG_MSA: list[pathlib.Path] = [pathlib.Path(x) for x in args.msa]
    OUTPUT_DIR = pathlib.Path(args.output).absolute()

    cfg = config_dict

    # Run Settings
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
    cfg["dimerscore"] = args.dimerscore
    cfg["ampliconnumber"] = args.ampliconnumber
    cfg["use_cache"] = args.use_cache

    # See if the output dir already exsits
    if OUTPUT_DIR.is_dir() and not args.force:
        sys.exit(f"ERROR: {OUTPUT_DIR} already exists, please use --force to override")

    ## If bedfile given, parse it:
    if args.bedfile:
        bed_lines = []
        try:
            with open(args.bedfile, "r") as bedfile:
                for line in bedfile.readlines():
                    line = line.strip()
                    if line:
                        bed_lines.append(line.split()[:3])

            print(bed_lines)
        except:
            sys.exit("Could not parse bedfile")
    else:
        print("no bed")

    ## Generate the regions
    regions = []
    for bedline in bed_lines:
        region = Region(bedline)
        regions.append(region)
        logger.debug(
            "Found region: <blue>{region_refname}</>\t{region_start}\t{region_end}",
            region_refname=region.ref_name,
            region_start=region.start,
            region_end=region.end,
        )
    ## Ensure the regions for the same genome do not overlap

    # Read in the MSAs
    msa_list: list[np.ndarray] = []
    for msa_path in ARG_MSA:
        records = SeqIO.parse(msa_path, "fasta")
        align_array = np.array([record.seq.upper() for record in records], dtype=str)
        align_array = remove_end_insertion(align_array)
        msa_list.append(align_array)

    # Get all msa file_names:
    msa_file_names = [x.name for x in ARG_MSA]

    ## Map each region to an MSA
    for region in regions:
        logger.info(
            "For region: <blue>{region_refname}</>\t{region_start}\t{region_end}",
            region_refname=region.ref_name,
            region_start=region.start,
            region_end=region.end,
        )
        # If region can be mapped
        if region.ref_name in msa_file_names:
            region.add_msa(
                msa_list[msa_file_names.index(region.ref_name)],
                cfg["amplicon_size_max"],
            )
            logger.info(
                "\tMapped to <green>{msa_path}</>",
                msa_path=ARG_MSA[msa_file_names.index(region.ref_name)],
            )
        # If region cannot be mapped
        else:
            logger.error(
                "\tCannot map to <red>MSA</>",
            )
            continue

        # Digest each region
        region.digest(cfg=cfg, thermo_cfg=thermo_cfg)
        logger.info(
            "\tDigested region into <green>{nfkmer} fkmers</> and <green>{nrkmer} rkmers</>",
            region_refname=region.ref_name,
            region_start=region.start,
            region_end=region.end,
            nfkmer=len(region.fkmers),
            nrkmer=len(region.rkmers),
        )

        # Generate all PrimerPairs
        region.generate_all_primerpairs(cfg=cfg)

        # Find all amplcions which span the region
        region.find_single_pp_coverage(cfg)

        logger.info(
            "\tOf <green>{nt_pp}</> PrimerPairs, <green>{ns_pp}</> fully span region",
            nt_pp=len(region.primer_pairs),
            ns_pp=len(region._single_primer_pairs),
        )

        if len(region._single_primer_pairs) == 0:
            logger.info("\t<red>Removed</> region due to no PrimerPairs")

    # Remove regions with no good amplicons
    # scheme = Scheme([x for x in regions if len(x._single_primer_pairs) > 0])

    regions = [x for x in regions if len(x._single_primer_pairs) > 0]

    if cfg["npools"] != 1:
        sys.exit("TODO: Just single pools for now :)")

    # The trys to find a global solution :)
    # network_opto(scheme)

    # Quick and dirty greedy solution
    regions.sort(key=lambda r: len(r._single_primer_pairs))

    for i in range(0, 100):
        print(f"Attempt: {i}")
        pool: list[PrimerPair] = []
        n_errors = 0

        for region in regions:
            all_seqs_in_pool = [
                x for sublist in (x.all_seqs() for x in pool) for x in sublist
            ]
            added = False
            for primerpair in region._single_primer_pairs:
                if not all_seqs_in_pool or not do_pools_interact_py(
                    all_seqs_in_pool,
                    list(primerpair.all_seqs()),
                    cfg["dimer_threshold"],
                ):
                    # If they do not interact
                    primerpair.refname = region.ref_name
                    added = True
                    pool.append(primerpair)
                    break

            if added == False:
                print(f"NO AMPS FOR {region.__str__()}")
                n_errors += 1

        # If amplicons are added for all regions
        if n_errors == 0:
            break
        else:
            random.shuffle(regions)

    print("".join([pp.__str__(pp.refname, "panel") for pp in pool]))

    exit()

    # Add the Run UUID into the config
    cfg["run_uuid"] = run_uuid
    # Write the config file, combining the cfg and thermo_cfg
    with open(OUTPUT_DIR / f"config.json", "w") as outfile:
        comb = thermo_cfg | cfg
        outfile.write(json.dumps(comb, sort_keys=True))


if __name__ == "__main__":
    main()
