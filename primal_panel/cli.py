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
        required=True,
        nargs="+",
    )
    parser.add_argument(
        "--regionbedfile",
        help="Path to the bedfile containing the wanted regions",
        type=pathlib.Path,
        default=None,
    )
    parser.add_argument(
        "--inputbedfile",
        help="Path to a primer.bedfile containing the precalculated primers",
        type=pathlib.Path,
        default=None,
    )
    parser.add_argument(
        "--mode",
        help="Select what mode for selecting regions in --regionbedfile. \n'region-only': only create primers for the regions \n'region-all': select regions first then keep adding amplicons \n'all': add amplicons based on entropy",
        choices=["region-only", "region-all", "all"],
        type=str,
        default="region-only",
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
        default="output",
    )
    parser.add_argument(
        "--force",
        help="Override the output directory",
        action="store_true",
    )
    parser.add_argument(
        "--minoverlap",
        help="Min amount of coverage overlap between amplicons",
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
    parser.add_argument(
        "--mapping",
        choices=["first", "consensus"],
        default="first",
        type=str,
    )
    parser.add_argument(
        "--maxamplicons",
        help="Max number of amplicons to create",
        default=100,
        type=lambda x: int(x) if int(x) > 0 else sys.exit("ERROR: Must be > 0"),
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

    # Deal with conflict in params
    if args.regionbedfile is not None and args.mode == "all":
        sys.exit(
            f"ERROR: Cannot use --regionbedfile and --mode {args.mode} at the same time, please select a region-* mode or not use --regionbedfile"
        )
    if args.regionbedfile is None and args.mode in ["region-only", "region-all"]:
        sys.exit(f"ERROR: Cannot use --mode '{args.mode}' without --regionbedfile")

    return args
