import json
import pathlib
import sys

# Core imports
from primalscheme3.core.msa import MSA
from primalscheme3.core.bedfiles import read_in_bedprimerpairs
from primalscheme3.core.logger import setup_loger

from primaldimer_py import do_pools_interact_py  # type: ignore


def repair(args):
    OUTPUT_DIR = pathlib.Path(args.output).absolute()  # Keep absolute path

    # Read in the config file
    with open(args.config) as f:
        base_cfg = json.load(f)

    # Overwrite the config file
    base_cfg["n_cores"] = args.cores  # Overwrite the number of cores
    base_cfg["output_dir"] = str(OUTPUT_DIR)
    base_cfg["primer_max_walk"] = 50

    # See if the output dir already exsits
    if OUTPUT_DIR.is_dir() and not args.force:
        sys.exit(f"ERROR: {OUTPUT_DIR} already exists, please use --force to override")

    # Create the output dir and a work subdir
    pathlib.Path.mkdir(OUTPUT_DIR, exist_ok=True)
    pathlib.Path.mkdir(OUTPUT_DIR / "work", exist_ok=True)

    # Set up the logger
    logger = setup_loger(OUTPUT_DIR)

    # Read in the MSA file
    msa = MSA(
        name=args.msa.stem,
        path=args.msa,
        msa_index=0,
        mapping=base_cfg["mapping"],
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

    logger.info(
        "Read in MSA: <blue>{msa_path}</>\tseqs:<green>{msa_rows}</>\tcols:<green>{msa_cols}</>",
        msa_path=msa.name,
        msa_rows=msa.array.shape[0],
        msa_cols=msa.array.shape[1],
    )

    # Read in the bedfile
    all_primerpairs = read_in_bedprimerpairs(args.bedfile)

    # Get the primerpairs for this new MSA
    primerpairs_in_msa = [
        pp for pp in all_primerpairs if pp.chromname == msa._chrom_name
    ]
    if len(primerpairs_in_msa) == 0:
        raise ValueError(
            f"No primerpairs found for {msa._chrom_name} in {args.bedfile}"
        )

    # Digest the MSA
    # Go lower to ensure digestion products are always returned even if error
    msa.digest(
        base_cfg,
        indexes=(
            [x.fprimer.end for x in primerpairs_in_msa],
            [x.rprimer.start for x in primerpairs_in_msa],
        ),
    )

    # Test
    from primalscheme3.core.digestion import (
        mp_f_digest,
        mp_r_digest,
        walk_left,
        wrap_walk,
    )
    from primalscheme3.core.thermo import calc_tm, thermo_check_kmers, forms_hairpin

    for x in primerpairs_in_msa:
        total_col_seqs = set()
        for row_index in range(0, msa.array.shape[0]):
            # Get the start and end col
            start_array = msa.array[
                row_index, x.fprimer.end - base_cfg["primer_size_min"] : x.fprimer.end
            ]
            start_seq = "".join(start_array).replace("-", "")

            results = wrap_walk(
                walk_left,
                msa.array,
                x.fprimer.end,
                x.fprimer.end - base_cfg["primer_size_min"],
                row_index=row_index,
                seq_str=start_seq,
                cfg=base_cfg,
            )
            total_col_seqs.update(results)
        # Thermo check each sequence
        print(x.fprimer.end, total_col_seqs)

        seqs_results = {}
        for seq in total_col_seqs:
            if type(seq) != str:
                continue

            seqs_results[seq] = "PASS"

            # Check for hairpins
            if forms_hairpin([seq], base_cfg):
                seqs_results[seq] = "FAIL_HAIRPIN"
                continue
            # Check for THERMO
            if not thermo_check_kmers([seq], base_cfg):
                seqs_results[seq] = "FAIL_THERMO"
                continue
            # Check for dimer
            if do_pools_interact_py([seq], [seq], base_cfg["dimerscore"]):
                seqs_results[seq] = "FAIL_DIMER"
                continue

        print(seqs_results)

    logger.info(
        "<blue>{msa_path}</>: digested to <green>{num_fkmers}</> FKmers and <green>{num_rkmers}</> RKmers",
        msa_path=msa.name,
        num_fkmers=len(msa.fkmers),
        num_rkmers=len(msa.rkmers),
    )

    # Create a dodgy pool object to run the interaction checker on
    seqs_in_pools = [[] for _ in range(base_cfg["npools"])]
    for pp in all_primerpairs:
        seqs_in_pools[pp.pool].extend([*pp.fprimer.seqs, *pp.rprimer.seqs])

    # Keep track of what has been added
    addedseqs = set()
    # Repair the fprimer class in the primerpairs
    for pp in primerpairs_in_msa:
        # Get the fprimer from the msa for the same index
        newfprimer = [fkmer for fkmer in msa.fkmers if fkmer.end == pp.fprimer.end]
        if len(newfprimer) != 1:
            print(f"no new fprimer for {pp.fprimer.end}")
            continue

        # Find what new seqs need to be added
        new_fprimer_seqs = newfprimer[0].seqs.difference(pp.fprimer.seqs)
        # If new seqs run the interaction checker
        if new_fprimer_seqs:
            # guard for interaction
            if do_pools_interact_py(
                [*new_fprimer_seqs], seqs_in_pools[pp.pool], base_cfg["dimerscore"]
            ):
                logger.warning(
                    "<red>WARNING</>: Interaction detected for new fprimer {pp_chromname}:{pp_ampliconprefix}:{pp_amplicon_number}",
                    pp_chromname=pp.chromname,
                    pp_ampliconprefix=pp.ampliconprefix,
                    pp_amplicon_number=pp.amplicon_number,
                )
                raise ValueError()

            # Update the fprimer
            pp.fprimer.seqs.update(new_fprimer_seqs)
            addedseqs.update(new_fprimer_seqs)
            # Update the pool seqs
            seqs_in_pools[pp.pool].extend(new_fprimer_seqs)
            logger.info(
                "Added new seqs {seqs} to <green>fprimer</> {pp_chromname}:{pp_ampliconprefix}:{pp_amplicon_number}",
                pp_chromname=pp.chromname,
                pp_ampliconprefix=pp.ampliconprefix,
                pp_amplicon_number=pp.amplicon_number,
                seqs=",".join(new_fprimer_seqs),
            )

        # Find what rprimer seqs need to be added
        new_rprimer = [rkmer for rkmer in msa.rkmers if rkmer.start == pp.rprimer.start]
        if len(new_rprimer) != 1:
            raise ValueError(
                f"Could not find a new rprimer for {pp.chromname}:{pp.ampliconprefix}:{pp.amplicon_number}:RIGHT"
            )
        new_rprimer_seqs = new_rprimer[0].seqs.difference(pp.rprimer.seqs)
        # If new seqs run the interaction checker
        if new_rprimer_seqs:
            # guard for interaction
            if do_pools_interact_py(
                [*new_rprimer_seqs], seqs_in_pools[pp.pool], base_cfg["dimerscore"]
            ):
                logger.warning(
                    "<red>WARNING</>: Interaction detected for new rprimer {pp_chromname}:{pp_ampliconprefix}:{pp_amplicon_number}",
                    pp_chromname=pp.chromname,
                    pp_ampliconprefix=pp.ampliconprefix,
                    pp_amplicon_number=pp.amplicon_number,
                )
                raise ValueError()

            # Update the rprimer
            pp.rprimer.seqs.update(new_rprimer_seqs)
            addedseqs.update(new_rprimer_seqs)
            # Update the pool seqs
            seqs_in_pools[pp.pool].extend(new_rprimer_seqs)
            logger.info(
                "Added new seqs {seqs} to <green>rprimer</> {pp_chromname}:{pp_ampliconprefix}:{pp_amplicon_number}",
                pp_chromname=pp.chromname,
                pp_ampliconprefix=pp.ampliconprefix,
                pp_amplicon_number=pp.amplicon_number,
                seqs=",".join(new_rprimer_seqs),
            )

    logger.info("Added {n} new seqs", n=len(addedseqs))
