#!/usr/bin/python3
import argparse
import pathlib
import sys

# Module imports
from primalscheme3.__init__ import __version__
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.flu.flu_main import create_flu
from primalscheme3.interaction.interaction import visulise_interactions
from primalscheme3.panel.panel_main import PanelRunModes, panelcreate
from primalscheme3.remap.remap import remap
from primalscheme3.repair.repair import repair

# Import main functions
from primalscheme3.scheme.scheme_main import schemecreate, schemereplace

## Commands are in the format of
# {pclass}-{mode}
# pclass = pannel or scheme

# Example to create a scheme
# scheme-create

# To repair a scheme
# scheme-repair

# To create a pannel
# pannel-create


def validate_panel_create_args(args):
    # Validate amplicon size
    if len(args.ampliconsize) == 1:
        args.ampliconsizemin = int(args.ampliconsize[0] * 0.9)
        args.ampliconsizemax = int(args.ampliconsize[0] * 1.1)
    elif len(args.ampliconsize) == 2:
        args.ampliconsizemin = int(args.ampliconsize[0])
        args.ampliconsizemax = int(args.ampliconsize[1])
    else:
        raise ValueError(
            "ERROR: --ampliconsize must be a single value or two values [100<=x<=2000]"
        )


def validate_scheme_create_args(args) -> None:
    # Validate backtrack is not used with more than 2 pools as it is not supported
    if args.backtrack and args.npools > 2:
        raise ValueError("backtrack is not supported with >2 pools")

    # Validate amplicon size
    if len(args.ampliconsize) == 1:
        args.ampliconsizemin = int(args.ampliconsize[0] * 0.9)
        args.ampliconsizemax = int(args.ampliconsize[0] * 1.1)
    elif len(args.ampliconsize) == 2:
        args.ampliconsizemin = int(args.ampliconsize[0])
        args.ampliconsizemax = int(args.ampliconsize[1])
    else:
        raise ValueError(
            "ERROR: --ampliconsize must be a single value or two values [100<=x<=2000]"
        )

    # Check amplicon size
    if args.ampliconsizemin >= args.ampliconsizemax:
        raise ValueError(
            "ERROR: --ampliconsize min cannot be greater than max [100<=x<=2000]"
        )
    # Check npools
    if args.npools < 1:
        raise ValueError("ERROR: --npools cannot be less than 1")

    # Check the bedfile exsists if given
    if args.bedfile and not args.bedfile.is_file():
        raise ValueError(f"ERROR: No file found at: '{str(args.bedfile.absolute())}'")
    # If there is an invalid number of cores
    if int(args.cores) <= 0:
        raise ValueError(f"ERROR: {int(args.cores)} is not a valid core count")

    if args.minoverlap < 0:
        raise ValueError(
            f"ERROR: {int(args.minoverlap)} is outside of the range for --minoverlap [0<=x]"
        )


def validate_scheme_replace_args(args):
    # Validate amplicon size
    if args.ampliconsize is None:  # Optional argument as it can be read from the config
        pass
    elif len(args.ampliconsize) == 1:
        args.ampliconsizemin = int(args.ampliconsize[0] * 0.9)
        args.ampliconsizemax = int(args.ampliconsize[0] * 1.1)
    elif len(args.ampliconsize) == 2:
        args.ampliconsizemin = int(args.ampliconsize[0])
        args.ampliconsizemax = int(args.ampliconsize[1])
    else:
        raise ValueError(
            "ERROR: --ampliconsize must be a single value or two values [100<=x<=2000]"
        )

        # Check amplicon size
    if args.ampliconsizemin >= args.ampliconsizemax:
        raise ValueError(
            "ERROR: --ampliconsize min cannot be greater than max [100<=x<=2000]"
        )


def check_valid_freq(value):
    fvalue = float(value)
    if 0 <= fvalue <= 1:
        return fvalue
    else:
        raise argparse.ArgumentTypeError(
            f"--minbasefreq must be 0 <= value <= 1. value: {value}"
        )


def check_path_is_file(value: str | pathlib.Path) -> pathlib.Path:
    if isinstance(value, str):
        value = pathlib.Path(value)
    if not value.is_file():
        raise argparse.ArgumentTypeError(f"No file found at: '{str(value.absolute())}'")
    return value


def cli():
    description = "PrimerScheme3: A tool for designing highly multiplexed amplicon sequencing primers"
    global_parser = argparse.ArgumentParser(
        description=description,
    )
    global_parser.add_argument("-v", "--version", action="version", version=__version__)
    # Add some global args
    global_parser.add_argument(
        "--primer_gc_min",
        type=float,
        default=30,
    )
    global_parser.add_argument(
        "--primer_gc_max",
        type=float,
        default=55,
    )
    global_parser.add_argument(
        "--primer_tm_min",
        type=float,
        default=59.5,
    )
    global_parser.add_argument(
        "--primer_tm_max",
        type=float,
        default=62.5,
    )
    global_parser.add_argument(
        "-o",
        "--output",
        help="The output directory of the primer.bed file",
        type=str,
        default="output",
    )
    global_parser.add_argument(
        "--force",
        help="Override the output directory",
        action="store_true",
    )
    global_parser.add_argument(
        "--ignore-n", help="Ignore Ns in the MSA", action="store_true"
    )
    subparsers = global_parser.add_subparsers(
        title="subcommands", help="scheme types", required=True
    )

    # Add the panel create subparser
    panel_create_parser = subparsers.add_parser(
        "panel-create",
        help="create a primer panel",
    )
    panel_create_parser.add_argument(
        "-m",
        "--msa",
        help="Paths to the MSA files",
        type=pathlib.Path,
        required=True,
        nargs="+",
    )
    panel_create_parser.add_argument(
        "--regionbedfile",
        help="Path to the bedfile containing the wanted regions",
        type=pathlib.Path,
        default=None,
    )
    panel_create_parser.add_argument(
        "--inputbedfile",
        help="Path to a primer.bedfile containing the precalculated primers",
        type=pathlib.Path,
        default=None,
    )
    panel_create_parser.add_argument(
        "--mode",
        help="Select what mode for selecting regions in --regionbedfile. \n'region-only': only create primers for the regions \n'region-all': select regions first then keep adding amplicons \n'all': add amplicons based on entropy",
        type=PanelRunModes,
        default="region-only",
    )
    panel_create_parser.add_argument(
        "-c",
        "--cores",
        help="The number of cores to use in Kmer digestion and thermo checking",
        type=int,
        default=8,
    )
    panel_create_parser.add_argument(
        "--ampliconsize",
        help="The size of an amplicon. Use single value for ± 10 percent, or two values to set min, max [100<=x<=2000]",
        type=int,
        default=[400],
        nargs="+",
    )
    panel_create_parser.add_argument(
        "--ampliconnumber",
        help="The number of amplicons for each msa",
        type=int,
        default=2,
    )
    panel_create_parser.add_argument(
        "--npools", help="Number of pools to use", default=1, type=int
    )
    panel_create_parser.add_argument(
        "--dimerscore", help="Threshold for dimer interaction", default=-26, type=float
    )
    panel_create_parser.add_argument(
        "--reducekmers",
        help="An existing bedfile to add primers to",
        action="store_true",
    )
    panel_create_parser.add_argument(
        "--minbasefreq",
        help="Min freq to be included,[0<=x<=1]",
        type=check_valid_freq,
        default=0.0,
    )
    panel_create_parser.add_argument(
        "--mapping",
        choices=["first", "consensus"],
        default="first",
        type=str,
    )
    panel_create_parser.add_argument(
        "--maxamplicons",
        help="Max number of amplicons to create",
        default=100,
        type=lambda x: int(x) if int(x) > 0 else sys.exit("ERROR: Must be > 0"),
    )
    panel_create_parser.set_defaults(func=panelcreate)

    # Parser to create a primer schemes
    scheme_create_parser = subparsers.add_parser(
        "scheme-create", help="generate an overlapping primerscheme"
    )

    scheme_create_parser.add_argument(
        "-m",
        "--msa",
        help="Paths to the MSA files",
        type=check_path_is_file,
        required=True,
        nargs="+",
    )
    scheme_create_parser.add_argument(
        "-c",
        "--cores",
        help="The number of cores to use in Kmer digestion and thermo checking",
        type=int,
        default=1,
    )
    scheme_create_parser.add_argument(
        "--ampliconsize",
        help="The size of an amplicon. Use single value for ± 10 percent, or two values to set min, max [100<=x<=2000]",
        type=int,
        default=[400],
        nargs="+",
    )
    scheme_create_parser.add_argument(
        "--minoverlap",
        help="min amount of overlap between primers",
        type=int,
        default=20,
    )
    scheme_create_parser.add_argument(
        "--npools", help="Number of pools to use", default=2, type=int
    )
    scheme_create_parser.add_argument(
        "--dimerscore", help="Threshold for dimer interaction", default=-26, type=float
    )
    scheme_create_parser.add_argument(
        "--bedfile", help="An existing bedfile to add primers to", type=pathlib.Path
    )
    scheme_create_parser.add_argument(
        "--reducekmers",
        help="Should number of sequences in each Kmer be reduced",
        action="store_true",
    )
    scheme_create_parser.add_argument(
        "--minbasefreq",
        help="Min freq to be included,[0<=x<=1]",
        type=check_valid_freq,
        default=0.0,
    )
    scheme_create_parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Prevents the creation of a plots of the primer scheme",
    )
    scheme_create_parser.add_argument(
        "--mapping",
        help="How should the primers in the bedfile be mapped",
        choices=["consensus", "first"],
        default="consensus",
    )
    scheme_create_parser.add_argument(
        "--circular",
        help="Should a circular amplicon be added (vv experimental)",
        action="store_true",
    )
    scheme_create_parser.add_argument(
        "--backtrack",
        help="Should the algorythm backtrack (vv experimental)",
        action="store_true",
    )
    scheme_create_parser.set_defaults(func=schemecreate)

    # Add the replace subparser
    scheme_replace_parser = subparsers.add_parser(
        "scheme-replace",
        help="replace a primer in a bedfile",
    )
    scheme_replace_parser.add_argument(
        "--primername",
        help="The name of the primer to replace, as found in the bedfile",
        type=str,
        required=True,
    )
    scheme_replace_parser.add_argument(
        "--primerbed",
        help="The bedfile containing the primer to replace",
        type=check_path_is_file,
        required=True,
    )
    scheme_replace_parser.add_argument(
        "--msa",
        help="The msa used to create the original primer scheme",
        type=check_path_is_file,
        required=True,
    )
    scheme_replace_parser.add_argument(
        "--config",
        help="The config.json used to create the original primer scheme",
        type=check_path_is_file,
        required=True,
    )
    scheme_replace_parser.add_argument(
        "--ampliconsize",
        help="The size of an amplicon. Use single value for ± 10 percent, or two values to set min, max [100<=x<=2000]",
        type=int,
        default=[400],
        nargs="+",
    )
    scheme_replace_parser.set_defaults(func=schemereplace)

    interactions_parser = subparsers.add_parser(
        "interactions", help="Shows all interactions within a bedfile"
    )
    interactions_parser.add_argument(
        "--bedfile", help="Path to the bedfile", type=check_path_is_file, required=True
    )
    interactions_parser.add_argument(
        "--threshold",
        help="Only show interactions more severe (Lower score) than this value",
        type=float,
        default=-26,
    )
    interactions_parser.set_defaults(func=visulise_interactions)

    ##############################
    #   Parser for repair mode   #
    ##############################
    repair_parser = subparsers.add_parser(
        "repair",
        help="repair a primer scheme via adding more primers to account for new mutations",
    )
    repair_parser.add_argument(
        "--bedfile", help="Path to the bedfile", type=pathlib.Path
    )
    repair_parser.add_argument(
        "--msa",
        help="An MSA, with the reference.fasta, aligned to any new genomes with mutations",
        type=pathlib.Path,
        required=True,
    )
    repair_parser.add_argument(
        "--config", help="Path to the config.json", type=pathlib.Path, required=True
    )
    repair_parser.add_argument(
        "-c",
        "--cores",
        help="The number of cores to use in Kmer digestion and thermo checking",
        type=int,
        default=1,
    )
    repair_parser.set_defaults(func=repair)

    ############################
    #   Parser for flu mode   #
    ############################

    flu_parser = subparsers.add_parser(
        "flu",
        help="create a flu panel",
    )
    flu_parser.add_argument("--bedfile", help="Path to the bedfile", type=pathlib.Path)
    flu_parser.add_argument(
        "--msa",
        help="An MSA, with the reference.fasta, aligned to any new genomes with mutations",
        type=pathlib.Path,
        required=True,
        nargs="+",
    )
    flu_parser.set_defaults(func=create_flu)

    ############################
    #   Parser for remap mode  #
    ############################

    remap_parser = subparsers.add_parser(
        "remap", help="Remap a primer scheme to a new reference genome"
    )
    remap_parser.add_argument(
        "-b", "--bedfile", help="Path to the bedfile", type=pathlib.Path, required=True
    )
    remap_parser.add_argument(
        "-i",
        "--id-to-remap-to",
        help="The ID of the reference genome to remap to",
        type=str,
        required=True,
    )
    remap_parser.add_argument(
        "-m",
        "--msa",
        help="Path to the MSA file",
        type=pathlib.Path,
        required=True,
    )
    remap_parser.set_defaults(func=remap)

    args = global_parser.parse_args()

    # Validate some global args
    if args.primer_gc_max <= args.primer_gc_min:
        raise ValueError(
            f"ERROR: --primer_gc_max ({args.primer_gc_max}) cannot be smaller than --primer_gc_min ({args.primer_gc_min})",
        )

    # Check Tms
    if args.primer_tm_max <= args.primer_tm_min:
        raise ValueError(
            f"ERROR: --primer_tm_max ({args.primer_tm_max}) cannot be smaller than --primer_tm_min ({args.primer_tm_min})"
        )

    # if output directory exists and force is not set
    if pathlib.Path(args.output).is_dir() and not args.force:
        raise ValueError(
            f"ERROR: Output directory '{args.output}' already exists. Use --force to override"
        )

    pm = ProgressManager()

    if args.func == schemecreate:
        validate_scheme_create_args(args)
        schemecreate(
            argmsa=args.msa,
            primer_gc_min=args.primer_gc_min,
            primer_gc_max=args.primer_gc_max,
            primer_tm_min=args.primer_tm_min,
            primer_tm_max=args.primer_tm_max,
            force=args.force,
            ignore_n=args.ignore_n,
            ampliconsizemin=args.ampliconsizemin,
            ampliconsizemax=args.ampliconsizemax,
            minoverlap=args.minoverlap,
            npools=args.npools,
            dimerscore=args.dimerscore,
            bedfile=args.bedfile,
            reducekmers=args.reducekmers,
            minbasefreq=args.minbasefreq,
            mapping=args.mapping,
            no_plot=args.no_plot,
            ncores=args.cores,
            output_dir=args.output,
            backtrack=args.backtrack,
            circular=args.circular,
            progress_manager=pm,
        )
    elif args.func == schemereplace:
        validate_scheme_replace_args(args)
        schemereplace(
            config=args.config,
            primername=args.primername,
            ampliconsizemin=args.ampliconsizemin,
            ampliconsizemax=args.ampliconsizemax,
            primerbed=args.primerbed,
            msapath=args.msa,
        )
    elif args.func == panelcreate:
        validate_panel_create_args(args)
        panelcreate(
            argmsa=args.msa,
            outputdir=args.output,
            primer_gc_max=args.primer_gc_max,
            primer_gc_min=args.primer_gc_min,
            primer_tm_max=args.primer_tm_max,
            primer_tm_min=args.primer_tm_min,
            force=args.force,
            ignore_n=args.ignore_n,
            regionbedfile=args.regionbedfile,
            inputbedfile=args.inputbedfile,
            mode=args.mode,
            cores=args.cores,
            ampliconsizemin=args.ampliconsizemin,
            ampliconsizemax=args.ampliconsizemax,
            npools=args.npools,
            dimerscore=args.dimerscore,
            reducekmers=args.reducekmers,
            minbasefreq=args.minbasefreq,
            mapping=args.mapping,
            maxamplicons=args.maxamplicons,
        )
    elif args.func == visulise_interactions:
        visulise_interactions(bedpath=args.bedfile, threshold=args.threshold)
    elif args.func == repair:
        repair(
            config_path=args.config,
            msa_path=args.msa,
            cores=args.cores,
            bedfile_path=args.bedfile,
            output_dir=args.output,
            force=args.force,
        )
    elif args.func == create_flu:
        create_flu(args)
    elif args.func == remap:
        remap(
            bedfile_path=args.bedfile,
            id_to_remap_to=args.id_to_remap_to,
            msa_path=args.msa,
            output_dir=args.output,
        )
    else:
        raise ValueError(f"ERROR: No function found for {args.func}")


if __name__ == "__main__":
    cli()
