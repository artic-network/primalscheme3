#!/usr/bin/python3

import argparse
import sys
import pathlib


def cli():
    description = "Generates a primerscheme from an MSA"
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
        "-c",
        "--cores",
        help="The number of cores to use in Kmer digestion and thermo checking",
        type=int,
        default=1,
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
        default=900,
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
        "--refnames",
        help="The name of the referance used in the outputed bed file",
        type=str,
        default=["MSA"],
        nargs="+",
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
        "--mismatches_self",
        help="How many mismatches between the primer sequence and its own genome should still count as a hit",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--mismatches_alt",
        help="How many mismatches between the primer sequence and its other genomes should still count as a hit",
        type=int,
        default=2,
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
        help="An existing bedfile to add primers to",
        type=bool,
        default=False,
    )
    args = parser.parse_args()

    # Check the bedfile exsists if given
    if args.bedfile and not args.bedfile.is_file():
        sys.exit(f"ERROR: No file found at: '{str(args.bedfile.absolute())}'")

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

    ## Check the number of msa and the number of msa names match
    if len(args.refnames) != len(args.msa) and args.refnames != ["msa"]:
        alt_names = [f"msa_{x}" for x in range(0, len(args.msa))]
        print(
            f"ISSUE: The number of args in --msa ({len(args.msa)}) != the number of args in --refnames ({len(args.refnames)}).\tUsing refnames of {alt_names}"
        )
        args.refnames = alt_names

    # Check gc min is less than GC max
    if args.primer_gc_max <= args.primer_gc_min:
        sys.exit(f"ERROR: --primer_gc_max cannot be smaller than --primer_gc_min")
    # Check Tms
    if args.primer_tm_max <= args.primer_tm_min:
        sys.exit(f"ERROR: --primer_tm_max cannot be smaller than --primer_tm_min")

    # Check the mismatches are both posative values
    if args.mismatches_self < 0:
        sys.exit(
            f"ERROR: --mismatches_self must be > 0. Got value {args.mismatches_self}"
        )
    if args.mismatches_alt < 0:
        sys.exit(
            f"ERROR: --mismatches_alt must be > 0. Got value {args.mismatches_alt}"
        )

    return args
