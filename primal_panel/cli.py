#!/usr/bin/python3

import argparse
import sys
import pathlib


def check_valid_freq(value):
    fvalue = float(value)
    if 0 <= fvalue <= 1:
        return fvalue
    else:
        raise argparse.ArgumentTypeError(
            "--minbasefreq must be 0 <= value <= 1. value: %s" % value
        )


def cli():
    description = "Generates a primerscheme panel from an MSA"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-m",
        "--msa",
        help="Paths to the MSA files",
        type=pathlib.Path,
        default=[
            "/Users/kentcg/Downloads/GCF_000005845.2_ASM584v2_genomic.fna",
            "/Users/kentcg/quick_lab_fork/SARS-CoV-2/400/v5.0.0_400/MN908947.3.fasta",
            "/Users/kentcg/Downloads/138954_final.fasta",
        ],  # required=True,
        nargs="+",
    )
    parser.add_argument(
        "-c",
        "--cores",
        help="The number of cores to use in Kmer digestion and thermo checking",
        type=int,
        default=8,
    )
    parser.add_argument(
        "--ampliconsizemax",
        help="The max size of an amplicon [100<=x<=2000]",
        type=int,
        default=1000,
    )
    parser.add_argument(
        "--ampliconsizemin",
        help="The min size of an amplicon [100<=x<=2000] Default = 0 (0.9 * ampliconsizemax)",
        type=int,
        default=800,
    )
    parser.add_argument(
        "--ampliconnumber",
        help="The number of amplicons for each msa",
        type=int,
        default=2,
    )
    parser.add_argument(
        "-p",
        "--prefix",
        help="The prefix used in the bedfile name",
        type=str,
        default="output",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="The output directory of the primer.bed file",
        type=str,
        default="/Users/kentcg/primal-nf/data/dengue_virus_type_2--11060/11060_scheme",
    )
    parser.add_argument(
        "--force",
        help="Over ride the output directory",
        action="store_true",
    )
    parser.add_argument(
        "--minoverlap",
        help="min amount of overlap between primers",
        type=int,
        default=20,
    )
    parser.add_argument(
        "--primer_gc_min",
        type=float,
        default=30,
    )
    parser.add_argument(
        "--primer_gc_max",
        type=float,
        default=55,
    )
    parser.add_argument(
        "--primer_tm_min",
        type=float,
        default=59.5,
    )
    parser.add_argument(
        "--primer_tm_max",
        type=float,
        default=62.5,
    )
    parser.add_argument(
        "--use_cache",
        help="Save / Load data from the caches",
        action="store_true",
    )
    parser.add_argument(
        "--dev",
        help="Enable dev options",
        action="store_true",
    )
    parser.add_argument(
        "-b",
        "--bedfile",
        help="A bedfile containing regions of intrest. Please ensure the region tab is the same as the filename of the msa it corrasponds to",
        type=pathlib.Path,
        default=False,
    )
    parser.add_argument("--npools", help="Number of pools to use", default=2, type=int)
    parser.add_argument(
        "--dimerscore", help="Threshold for dimer interaction", default=-26, type=float
    )
    parser.add_argument(
        "--reducekmers",
        help="An existing bedfile to add primers to",
        type=bool,
        default=False,
    )
    parser.add_argument(
        "--minbasefreq",
        help="Min freq to be included,[0<=x<=1]",
        type=check_valid_freq,
        default=0.0,
    )

    args = parser.parse_args()
    # Parse number of pools
    if args.npools < 1:
        sys.exit(f"ERROR: Needs at least one pool, {args.npools} requested")

    # If there is an invalid number of cores
    if int(args.cores) <= 0:
        sys.exit(f"ERROR: {int(args.cores)} is not a valid core count")

    if 100 <= int(args.ampliconsizemax) <= 2000:
        pass
    else:
        sys.exit(
            f"ERROR: {int(args.ampliconsizemax)} is outside of the range for --amplicon-size-max [100<=x<=2000]"
        )

    if args.minoverlap < 0:
        sys.exit(
            f"ERROR: {int(args.minoverlap)} is outside of the range for --minoverlap [0<=x]"
        )

    if args.ampliconsizemin == 0:
        args.ampliconsizemin = 0.9 * args.ampliconsizemax

    if 100 <= int(args.ampliconsizemin) <= int(args.ampliconsizemax):
        pass
    else:
        sys.exit(
            f"ERROR: {int(args.ampliconsizemin)} is outside of the range for --amplicon-size-min [100<=x<={int(args.ampliconsizemax)}]"
        )

    # Check gc min is less than GC max
    if args.primer_gc_max <= args.primer_gc_min:
        sys.exit(f"ERROR: --primer_gc_max cannot be smaller than --primer_gc_min")
    # Check Tms
    if args.primer_tm_max <= args.primer_tm_min:
        sys.exit(f"ERROR: --primer_tm_max cannot be smaller than --primer_tm_min")

    return args
