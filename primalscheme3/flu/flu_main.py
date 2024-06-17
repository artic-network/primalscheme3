import hashlib
import json
import pathlib
import shutil
import sys

from Bio import Seq, SeqIO, SeqRecord

from primalscheme3.core.bedfiles import read_in_bedlines
from primalscheme3.core.classes import FKmer, RKmer
from primalscheme3.core.config import config_dict
from primalscheme3.core.create_report_data import (
    generate_all_plotdata,
)
from primalscheme3.core.create_reports import generate_all_plots
from primalscheme3.core.digestion import generate_valid_primerpairs
from primalscheme3.core.logger import setup_loger
from primalscheme3.core.mapping import generate_consensus, generate_reference
from primalscheme3.core.mismatches import MatchDB
from primalscheme3.core.msa import MSA
from primalscheme3.core.progress_tracker import ProgressManager

## Submodules
from primalscheme3.scheme.classes import Scheme, SchemeReturn


def create_flu(args):
    # Read in the bedfile
    bedlines, _header = read_in_bedlines(args.bedfile)

    # Create the output dir
    OUTPUT_DIR = pathlib.Path(args.output).absolute()  # Keep absolute path
    ARG_MSA = args.msa
    # See if the output dir already exsits
    if OUTPUT_DIR.is_dir() and not args.force:
        sys.exit(f"ERROR: {OUTPUT_DIR} already exists, please use --force to override")
    # Create the output dir and a work subdir
    pathlib.Path.mkdir(OUTPUT_DIR, exist_ok=True)
    pathlib.Path.mkdir(OUTPUT_DIR / "work", exist_ok=True)

    # Set up the logger
    logger = setup_loger(OUTPUT_DIR)
    # Set up the progress manager
    pm = ProgressManager()

    ## Hard code some params for now
    cfg = config_dict
    cfg["amplicon_size_max"] = 1200
    cfg["amplicon_size_min"] = 900
    cfg["npools"] = 2
    cfg["mapping"] = "first"
    cfg["minbasefreq"] = 0.01
    cfg["n_cores"] = 8
    cfg["mismatch_product_size"] = 1100
    cfg["backtrack"] = True

    # Create an empty Mismatch
    mismatch_db = MatchDB(OUTPUT_DIR / "work/mismatch", [], cfg["primer_size_min"])
    logger.info(
        "<green>Created:</> {path}",
        path=f"{OUTPUT_DIR.relative_to(OUTPUT_DIR.parent)}/work/mismatch.db",
    )

    # Create a dict full of msa data
    msa_data = {}

    # Create the scheme
    scheme = Scheme(cfg=cfg, matchDB=mismatch_db)

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
            progress_manager=pm,
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
        # Add the msa to the scheme
        msa_dict[msa_index] = msa

        # Digest the MSA into FKmers and RKmers
        msa.digest(cfg)
        logger.info(
            "<blue>{msa_path}</>: digested to <green>{num_fkmers}</> FKmers and <green>{num_rkmers}</> RKmers",
            msa_path=msa.name,
            num_fkmers=len(msa.fkmers),
            num_rkmers=len(msa.rkmers),
        )

        # Generate the left_primerpairs
        left_primer_pairs = generate_valid_primerpairs(
            fkmers=[
                FKmer(
                    seqs={"GTTACGCGCCAGCAAAAGCAGG", "GTTACGCGCCAGCGAAAGCAGG"}, end=22
                )  ## Correct start index for all genes
            ],
            rkmers=msa.rkmers,
            dimerscore=cfg["dimerscore"],
            amplicon_size_max=cfg["amplicon_size_max"],
            amplicon_size_min=cfg["amplicon_size_min"],
            msa_index=msa.msa_index,
            progress_manager=pm,
        )
        left_primer_pairs.sort(
            key=lambda x: len(x.all_seqs())
            / (x.rprimer.start - cfg["amplicon_size_min"] + 100)
        )
        for lpp in left_primer_pairs:
            lpp.chrom_name = msa._chrom_name
            lpp.amplicon_prefix = msa._uuid

        # Generate the right_primerpairs
        msa_length = len(generate_reference(msa.array))
        right_primer_site = msa_length - 22
        right_primer_pairs = generate_valid_primerpairs(
            fkmers=msa.fkmers,
            rkmers=[RKmer(seqs={"GTTACGCGCCAGTAGAAACAAGG"}, start=right_primer_site)],
            dimerscore=cfg["dimerscore"],
            amplicon_size_max=cfg["amplicon_size_max"],
            amplicon_size_min=cfg["amplicon_size_min"],
            msa_index=msa.msa_index,
            progress_manager=pm,
        )
        for rpp in right_primer_pairs:
            rpp.chrom_name = msa._chrom_name
            rpp.amplicon_prefix = msa._uuid

        all_primerpairs = [*left_primer_pairs, *right_primer_pairs]

        while True:
            # Try and add an overlapping primer
            match scheme.try_ol_primerpairs(all_primerpairs, msa_index):
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
                    print("NO_OL_PRIMERPAIR")
                    break

            # Try and backtrack
            if cfg["backtrack"]:
                match scheme.try_backtrack(all_primerpairs, msa_index):
                    case SchemeReturn.ADDED_BACKTRACKED:
                        logger.info(
                            "Backtracking allowed <green>overlapping</> amplicon to be added for <blue>{msa_name}</>: {primer_start}\t{primer_end}\t{primer_pool}",
                            primer_start=scheme._last_pp_added[-1].start,
                            primer_end=scheme._last_pp_added[-1].end,
                            primer_pool=scheme._last_pp_added[-1].pool + 1,
                            msa_name=msa.name,
                        )
                        continue
                    case SchemeReturn.NO_BACKTRACK:
                        logger.info(
                            "Could not backtrack for <blue>{msa_name}</>",
                            msa_name=msa.name,
                        )
                        pass  # Do nothing move on to next step

                # Try and add a walking primer
            if (
                scheme.try_walk_primerpair(right_primer_pairs, msa_index)
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

    logger.info("Writting output files")
    # Write primer bed file
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
    amp_md5 = hashlib.md5("\n".join(amp_bed_str).encode()).hexdigest()
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
