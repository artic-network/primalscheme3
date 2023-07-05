from primal_digest.thermo import *
from primal_digest.cli import cli
from primal_digest.config import (
    config_dict,
    thermo_config,
)
from primaldimer_py import do_pools_interact_py
from primal_digest.classes import FKmer, RKmer, PrimerPair, Scheme, BedPrimer, BedRecord
from primal_digest.digestion import *
from primal_digest.get_window import *
from primal_digest.bedfiles import parse_bedfile

import numpy as np
from Bio import SeqIO
from multiprocessing import Pool
import kmertools
import json

# Added
import hashlib
import sys
import pathlib


"""
This is a test of a new dynamic digestion algo
"""


def mp_pp_inter_free(data: tuple[PrimerPair, dict]) -> bool:
    """
    True means interaction
    """
    pp = data[0]
    cfg = data[1]
    return do_pools_interact_py(
        list(pp.fprimer.seqs), list(pp.rprimer.seqs), cfg["dimerscore"]
    )


def main():
    args = cli()
    ARG_MSA = args.msa
    OUTPUT_DIR = pathlib.Path(args.output).absolute()

    cfg = config_dict
    thermo_cfg = thermo_config

    # Primer Digestion settings
    thermo_cfg["kmer_size_max"] = 36
    thermo_cfg["kmer_size_min"] = 14
    thermo_cfg["primer_gc_min"] = args.primer_gc_min
    thermo_cfg["primer_gc_max"] = args.primer_gc_max
    thermo_cfg["primer_tm_min"] = args.primer_tm_min
    thermo_cfg["primer_tm_max"] = args.primer_tm_max
    thermo_cfg["dimerscore"] = args.dimerscore

    # Run Settings
    cfg["refname"] = args.refnames
    cfg["n_cores"] = args.cores
    cfg["output_prefix"] = args.prefix
    cfg["output_dir"] = str(OUTPUT_DIR)
    cfg["msa_paths"] = ARG_MSA
    cfg["mismatches_self"] = args.mismatches_self
    cfg["mismatches_alt"] = args.mismatches_alt
    cfg["amplicon_size_max"] = args.ampliconsizemax
    cfg["amplicon_size_min"] = args.ampliconsizemin
    cfg["min_overlap"] = args.minoverlap
    cfg["force"] = args.force
    cfg["npools"] = args.npools
    cfg["dimerscore"] = args.dimerscore

    msa_index_to_ref_name = {
        index: msa_name for index, msa_name in enumerate(cfg["refname"])
    }

    cfg["msa_index_to_ref_name"] = msa_index_to_ref_name

    # See if the output dir already exsits
    if OUTPUT_DIR.is_dir():
        if not args.force:
            sys.exit(
                f"ERROR: {OUTPUT_DIR} already exists, please use --force to override"
            )

    # Read in the MSAs
    msa_list: list[np.ndarray] = []
    for msa_path in ARG_MSA:
        records = SeqIO.parse(msa_path, "fasta")
        align_array = np.array([record.seq.upper() for record in records], dtype=str)
        align_array = remove_end_insertion(align_array)
        msa_list.append(align_array)

    raw_msa_list: list[np.ndarray] = []
    for msa_path in ARG_MSA:
        records = SeqIO.parse(msa_path, "fasta")
        align_array = np.array([record.seq.upper() for record in records], dtype=str)
        raw_msa_list.append(align_array)

    # Generate the Kmers for each array in msa_list
    unique_f_r_msa: list[list[list[FKmer], list[RKmer]]] = []

    for msa_index, msa_array in enumerate(msa_list):
        mp_thermo_pass_fkmers, mp_thermo_pass_rkmers = digest(
            msa_array, cfg, thermo_cfg
        )
        # Use the custom Rust unique checker

        f_kmer_bools = []
        r_kmer_bools = []

        for index, raw_msa_array in enumerate(raw_msa_list):
            referance_seqs = ["".join(x) for x in raw_msa_array]
            if msa_index == index:
                r_results = kmertools.rkmer_is_unique(
                    rkmers=[(x.start, list(x.seqs)) for x in mp_thermo_pass_rkmers],
                    referance_seqs=referance_seqs,
                    n_cores=cfg["n_cores"],
                    mismatches=cfg["mismatches_self"],
                    detect_expected=False,
                )  # Should be false
                f_results = kmertools.fkmer_is_unique(
                    fkmers=[(x.end, list(x.seqs)) for x in mp_thermo_pass_fkmers],
                    referance_seqs=referance_seqs,
                    n_cores=cfg["n_cores"],
                    mismatches=cfg["mismatches_self"],
                    detect_expected=False,
                )
                f_kmer_bools.append(f_results)
                r_kmer_bools.append(r_results)

            else:
                r_results = kmertools.rkmer_is_unique(
                    rkmers=[(x.start, list(x.seqs)) for x in mp_thermo_pass_rkmers],
                    referance_seqs=referance_seqs,
                    n_cores=cfg["n_cores"],
                    mismatches=cfg["mismatches_alt"],
                    detect_expected=True,
                )
                f_results = kmertools.fkmer_is_unique(
                    fkmers=[(x.end, list(x.seqs)) for x in mp_thermo_pass_fkmers],
                    referance_seqs=referance_seqs,
                    n_cores=cfg["n_cores"],
                    mismatches=cfg["mismatches_alt"],
                    detect_expected=True,
                )
                f_kmer_bools.append(f_results)
                r_kmer_bools.append(r_results)

        ## Combine the data for each kmer for each msa
        unique_f_bools = [all(x) for x in zip(*f_kmer_bools)]
        unique_r_bools = [all(x) for x in zip(*r_kmer_bools)]

        ## Get the unique kmers
        unique_fkmers = [
            fkmer
            for (bool, fkmer) in zip(unique_f_bools, mp_thermo_pass_fkmers)
            if bool
        ]
        unique_rkmers = [
            rkmer
            for (bool, rkmer) in zip(unique_r_bools, mp_thermo_pass_rkmers)
            if bool
        ]

        # Append them to the msa
        unique_f_r_msa.append((unique_fkmers, unique_rkmers))

    msa_primerpairs = []
    # Generate all valid primerpairs for each msa
    for msa_index, unique_fr_kmers in enumerate(unique_f_r_msa):
        wanted_fkmers = unique_fr_kmers[0]
        wanted_rkmers = unique_fr_kmers[1]

        # Generate all primer pairs
        primer_pairs: list[PrimerPair] = []
        for f in wanted_fkmers:
            pos_r = get_r_window_FAST2(
                kmers=wanted_rkmers,
                start=min(f.starts()) + cfg["amplicon_size_min"],
                end=min(f.starts()) + cfg["amplicon_size_max"],
            )
            for r in pos_r:
                primer_pairs.append(PrimerPair(f, r))

        # Filter out primerpairs if f primer seqs interact with r primer seqs
        with Pool(cfg["n_cores"]) as p:
            mp_pp_bool = p.map(
                mp_pp_inter_free, [(pp, thermo_cfg) for pp in primer_pairs]
            )

        iter_free_primer_pairs: list[PrimerPair] = [
            pp for (bool, pp) in zip(mp_pp_bool, primer_pairs) if not bool
        ]
        iter_free_primer_pairs.sort(key=lambda pp: (pp.fprimer.end, -pp.rprimer.start))

        msa_primerpairs.append(iter_free_primer_pairs)

    scheme = Scheme(cfg=cfg)

    # If the bedfile flag is given add the primers into the scheme
    if args.bedfile:
        current_pools = parse_bedfile(args.bedfile, cfg["npools"])
        # Assign the bedfile generated pool to the scheme in a hacky way
        scheme._pools = current_pools

    for msa_index, primerpairs_in_msa in enumerate(msa_primerpairs):
        # Add the first primer, and if no primers can be added move to next msa
        if not scheme.add_first_primer_pair(primerpairs_in_msa, msa_index):
            continue

        while True:
            # Try and add an overlapping primer
            if scheme.try_ol_primerpairs(primerpairs_in_msa, thermo_cfg, msa_index):
                continue
            # Try and add a walking primer
            elif scheme.try_walk_primerpair(primerpairs_in_msa, thermo_cfg, msa_index):
                continue
            else:
                break

    # Create the output dir
    pathlib.Path.mkdir(OUTPUT_DIR, exist_ok=True)

    # Write primer bed file
    with open(OUTPUT_DIR / f"{cfg['output_prefix']}.primer.bed", "w") as outfile:
        primer_bed_str = []
        for pp in scheme.all_primers():
            st = pp.__str__(
                msa_index_to_ref_name.get(pp.msa_index, "NA"),
                msa_index_to_ref_name.get(pp.msa_index, "NA"),
            )
            primer_bed_str.append(st.strip())
        outfile.write("\n".join(primer_bed_str))

    # Write amplicon bed file
    with open(OUTPUT_DIR / f"{cfg['output_prefix']}.amplicon.bed", "w") as outfile:
        amp_bed_str = []
        for pp in scheme.all_primers():
            amp_bed_str.append(
                f"{msa_index_to_ref_name.get(pp.msa_index, 'NA')}\t{pp.start}\t{pp.end}\tAMP_{pp.amplicon_number}\t{pp.pool + 1}"
            )
        outfile.write("\n".join(amp_bed_str))

    # Generate the bedfile hash, and add it into the config
    bed_md5 = hashlib.md5("\n".join(primer_bed_str).encode())
    cfg["md5_bed"] = bed_md5.hexdigest()

    # Generate the amplicon hash, and add it into the config
    bed_md5 = hashlib.md5("\n".join(amp_bed_str).encode())
    cfg["md5_amp"] = bed_md5.hexdigest()

    # Write the config file, combining the cfg and thermo_cfg
    with open(OUTPUT_DIR / f"config.json", "w") as outfile:
        comb = thermo_cfg | cfg
        outfile.write(json.dumps(comb, sort_keys=True))


if __name__ == "__main__":
    main()
