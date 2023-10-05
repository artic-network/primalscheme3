#!/usr/bin/python3

import argparse
import sys
import pathlib

# Module imports
from primal_digest.__init__ import __version__


def check_valid_freq(value):
    fvalue = float(value)
    if 0 <= fvalue <= 1:
        return fvalue
    else:
        raise argparse.ArgumentTypeError(
            f"--minbasefreq must be 0 <= value <= 1. value: {value}"
        )


def cli():
    description = "Generates a primerscheme from an MSA"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-v", "--version", action="version", version=__version__)
    parser.add_argument(
        "-m",
        "--msa",
        help="Paths to the MSA files",
        type=pathlib.Path,
        required=True,
        nargs="+",
    )
    parser.add_argument(
        "-c",
        "--cores",
        help="The number of cores to use in Kmer digestion and thermo checking",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--ampliconsize",
        help="The size of an amplicon. Use single value for ± 10 percent, or two values to set min, max [100<=x<=2000]",
        type=int,
        default=400,
        nargs="+",
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
    parser.add_argument("--npools", help="Number of pools to use", default=2, type=int)
    parser.add_argument(
        "--dimerscore", help="Threshold for dimer interaction", default=-26, type=float
    )
    parser.add_argument(
        "--bedfile", help="An existing bedfile to add primers to", type=pathlib.Path
    )
    parser.add_argument(
        "--reducekmers",
        help="Should number of sequences in each Kmer be reduced",
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
        "--plot",
        type=bool,
        default=True,
        help="Should HTML plots be generated",
    )
    parser.add_argument(
        "--mapping",
        help="How should the primers in the bedfile be mapped",
        choices=["consensus", "first"],
        default="consensus",
    )

    args = parser.parse_args()

    # Check the bedfile exsists if given
    if args.bedfile and not args.bedfile.is_file():
        sys.exit(f"ERROR: No file found at: '{str(args.bedfile.absolute())}'")

    # If there is an invalid number of cores
    if int(args.cores) <= 0:
        sys.exit(f"ERROR: {int(args.cores)} is not a valid core count")

    if args.minoverlap < 0:
        sys.exit(
            f"ERROR: {int(args.minoverlap)} is outside of the range for --minoverlap [0<=x]"
        )

    # Validate amplicon size
    if len(args.ampliconsize) == 1:
        args.ampliconsizemin = int(args.ampliconsize[0] * 0.9)
        args.ampliconsizemax = int(args.ampliconsize[0] * 1.1)
    elif len(args.ampliconsize) == 2:
        args.ampliconsizemin = int(args.ampliconsize[0])
        args.ampliconsizemax = int(args.ampliconsize[1])
    else:
        sys.exit(
            f"ERROR: --ampliconsize must be a single value or two values [100<=x<=2000]"
        )

    # Check amplicon size
    if args.ampliconsizemin >= args.ampliconsizemax:
        sys.exit(f"ERROR: --ampliconsize min cannot be greater than max [100<=x<=2000]")

    # Check gc min is less than GC max
    if args.primer_gc_max <= args.primer_gc_min:
        sys.exit(f"ERROR: --primer_gc_max cannot be smaller than --primer_gc_min")
    # Check Tms
    if args.primer_tm_max <= args.primer_tm_min:
        sys.exit(f"ERROR: --primer_tm_max cannot be smaller than --primer_tm_min")
    # Check npools
    if args.npools < 1:
        sys.exit(f"ERROR: --npools cannot be less than 1")

    return args
