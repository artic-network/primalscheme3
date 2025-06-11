#!/usr/bin/python3
import argparse
import pathlib
from typing import Annotated, Optional

import typer

# Module imports
from primalscheme3.__init__ import __version__
from primalscheme3.core.config import Config, MappingType
from primalscheme3.core.msa import parse_msa
from primalscheme3.core.primer_visual import bedfile_plot_html, primer_mismatch_heatmap
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.interaction.interaction import (
    visualise_interactions,
)
from primalscheme3.panel.panel_main import PanelRunModes, panelcreate
from primalscheme3.repair.repair import repair

# Import main functions
from primalscheme3.scheme.scheme_main import schemecreate, schemereplace

## Commands are in the format of
# {pclass}-{mode}
# pclass = panel or scheme

# Example to create a scheme
# scheme-create

# To repair a scheme
# scheme-repair

# To create a panel
# panel-create


def check_path_is_file(value: str | pathlib.Path) -> pathlib.Path:
    if isinstance(value, str):
        value = pathlib.Path(value)
    if not value.is_file():
        raise argparse.ArgumentTypeError(f"No file found at: '{str(value.absolute())}'")
    return value


# Create the main app
app = typer.Typer(name="primalscheme3", no_args_is_help=True)


def check_output_dir(output: pathlib.Path, force: bool):
    if output.exists() and not force:
        raise typer.BadParameter(
            f"--output '{output}' directory already exists. Use --force to overwrite"
        )


def typer_callback_version(value: bool):
    if value:
        version = typer.style(__version__, fg=typer.colors.GREEN, bold=True)
        typer.echo("PrimalScheme3 version: " + version)
        raise typer.Exit()


@app.callback()
def primalscheme3(
    value: Annotated[bool, typer.Option] = typer.Option(
        False, "--version", callback=typer_callback_version
    ),
):
    pass


@app.command(no_args_is_help=True)
def scheme_create(
    msa: Annotated[
        list[pathlib.Path],
        typer.Option(
            help="The MSA to design against. To use multiple MSAs, use multiple --msa flags. (--msa 1.fasta --msa 2.fasta)",
            exists=True,
            readable=True,
            resolve_path=True,
        ),
    ],
    output: Annotated[
        pathlib.Path,
        typer.Option(
            help="The output directory",
            resolve_path=True,
        ),
    ],
    amplicon_size: Annotated[
        int,
        typer.Option(
            help="The size of an amplicon. Min / max size are ± 10 percent [100<=x<=2000]",
            min=100,
            max=2000,
        ),
    ] = Config.amplicon_size,
    bedfile: Annotated[
        Optional[pathlib.Path],
        typer.Option(
            help="An existing bedfile to add primers to",
            exists=True,
            readable=True,
            resolve_path=True,
        ),
    ] = None,
    min_overlap: Annotated[
        int,
        typer.Option(help="min amount of overlap between primers", min=0),
    ] = Config.min_overlap,
    n_pools: Annotated[
        int, typer.Option(help="Number of pools to use", min=1)
    ] = Config.n_pools,
    dimer_score: Annotated[
        float,
        typer.Option(
            help="Threshold for dimer interaction",
        ),
    ] = Config.dimer_score,
    min_base_freq: Annotated[
        float,
        typer.Option(help="Min freq to be included,[0<=x<=1]", min=0.0, max=1.0),
    ] = Config.min_base_freq,
    mapping: Annotated[
        MappingType,
        typer.Option(
            help="How should the primers in the bedfile be mapped",
        ),
    ] = Config.mapping.value,  # type: ignore
    circular: Annotated[
        bool, typer.Option(help="Should a circular amplicon be added")
    ] = Config.circular,
    backtrack: Annotated[
        bool, typer.Option(help="Should the algorithm backtrack")
    ] = Config.backtrack,
    ignore_n: Annotated[
        bool,
        typer.Option(help="Should N in the input genomes be ignored"),
    ] = Config.ignore_n,
    force: Annotated[
        bool, typer.Option(help="Override the output directory")
    ] = Config.force,
    input_bedfile: Annotated[
        Optional[pathlib.Path],
        typer.Option(
            help="Path to a primer.bedfile containing the pre-calculated primers"
        ),
    ] = Config.input_bedfile,
    high_gc: Annotated[bool, typer.Option(help="Use high GC primers")] = Config.high_gc,
    offline_plots: Annotated[
        bool,
        typer.Option(
            "--offline-plots/--online-plots",
            help="Offline plots includes 3Mb of dependencies, so they can be viewed offline",
        ),
    ] = True,
    use_matchdb: Annotated[
        bool,
        typer.Option(
            help="Create and use a mispriming database",
        ),
    ] = Config.use_matchdb,
    ncores: Annotated[
        int,
        typer.Option(help="Number of CPU cores to use during digestion", min=1),
    ] = Config.ncores,
    use_annealing: Annotated[
        bool,
        typer.Option(
            "--use-annealing/--use-tm",
            help="Using annealing proportion rather than Tm to calculate primers",
        ),
    ] = Config.use_annealing,
):
    """
    Creates a tiling overlap scheme for each MSA file
    """
    # Update the config with CLI params
    config = Config(**locals())

    # Check the output directory
    check_output_dir(output, force)

    # Set up the progress manager
    pm = ProgressManager()
    schemecreate(
        msa=msa,
        output_dir=output,
        config=config,
        pm=pm,
        force=force,
        input_bedfile=input_bedfile,
        offline_plots=offline_plots,
    )


@app.command(no_args_is_help=True)
def scheme_replace(
    primername: Annotated[
        str, typer.Argument(help="The name of the primer to replace")
    ],
    primerbed: Annotated[
        pathlib.Path,
        typer.Argument(
            help="The bedfile containing the primer to replace",
            exists=True,
            readable=True,
            resolve_path=True,
        ),
    ],
    msa: Annotated[
        pathlib.Path,
        typer.Argument(
            help="The msa used to create the original primer scheme",
            exists=True,
            readable=True,
            resolve_path=True,
        ),
    ],
    amplicon_size_max: Annotated[
        int,
        typer.Option(
            help="The max size of an amplicon [100<=x<=2000]", min=100, max=2000
        ),
    ],
    amplicon_size_min: Annotated[
        int,
        typer.Option(
            help="The min size of an amplicon. [100<=x<=2000]", min=100, max=2000
        ),
    ],
    config: Annotated[
        pathlib.Path,
        typer.Option(
            help="The config.json used to create the original primer scheme",
            exists=True,
            readable=True,
            resolve_path=True,
        ),
    ],
):
    """
    Replaces a primerpair in a bedfile
    """

    # Set up the progress manager
    pm = ProgressManager()

    schemereplace(
        config_path=config,
        primername=primername,
        ampliconsizemax=amplicon_size_max,
        ampliconsizemin=amplicon_size_min,
        primerbed=primerbed,
        msapath=msa,
        pm=pm,
    )


@app.command(no_args_is_help=True)
def panel_create(
    msa: Annotated[
        list[pathlib.Path],
        typer.Option(
            help="Paths to the MSA files", exists=True, readable=True, resolve_path=True
        ),
    ],
    output: Annotated[
        pathlib.Path,
        typer.Option(
            help="The output directory",
            resolve_path=True,
        ),
    ],
    region_bedfile: Annotated[
        Optional[pathlib.Path],
        typer.Option(
            help="Path to the bedfile containing the wanted regions",
            readable=True,
            dir_okay=False,
            file_okay=True,
            exists=True,
        ),
    ] = None,
    input_bedfile: Annotated[
        Optional[pathlib.Path],
        typer.Option(
            help="Path to a primer.bedfile containing the pre-calculated primers",
            readable=True,
            dir_okay=False,
            file_okay=True,
            exists=True,
        ),
    ] = None,
    mode: Annotated[
        PanelRunModes,
        typer.Option(
            help="Select what run mode",
        ),
    ] = PanelRunModes.REGION_ONLY.value,  # type: ignore
    amplicon_size: Annotated[
        int, typer.Option(help="The size of an amplicon")
    ] = Config.amplicon_size,
    n_pools: Annotated[
        int, typer.Option(help="Number of pools to use", min=1)
    ] = Config.n_pools,
    dimer_score: Annotated[
        float, typer.Option(help="Threshold for dimer interaction")
    ] = Config.dimer_score,
    min_base_freq: Annotated[
        float,
        typer.Option(help="Min freq to be included,[0<=x<=1]", min=0.0, max=1.0),
    ] = Config.min_base_freq,
    mapping: Annotated[
        MappingType,
        typer.Option(
            help="How should the primers in the bedfile be mapped",
        ),
    ] = Config.mapping.value,  # type: ignore
    max_amplicons: Annotated[
        Optional[int], typer.Option(help="Max number of amplicons to create", min=1)
    ] = None,
    max_amplicons_msa: Annotated[
        Optional[int], typer.Option(help="Max number of amplicons for each MSA", min=1)
    ] = None,
    max_amplicons_region_group: Annotated[
        Optional[int],
        typer.Option(help="Max number of amplicons for each region", min=1),
    ] = None,
    force: Annotated[bool, typer.Option(help="Override the output directory")] = False,
    high_gc: Annotated[bool, typer.Option(help="Use high GC primers")] = Config.high_gc,
    offline_plots: Annotated[
        bool,
        typer.Option(
            "--offline-plots/--online-plots",
            help="Offline plots includes 3Mb of dependencies, so they can be viewed offline",
        ),
    ] = True,
    use_matchdb: Annotated[
        bool,
        typer.Option(
            help="Create and use a mispriming database",
        ),
    ] = Config.use_matchdb,
    ncores: Annotated[
        int,
        typer.Option(help="Number of CPU cores to use during digestion", min=1),
    ] = Config.ncores,
    use_annealing: Annotated[
        bool,
        typer.Option(
            "--use-annealing/--use-tm",
            help="Using annealing proportion rather than Tm to calculate primers",
        ),
    ] = Config.use_annealing,
):
    """
    Creates a primer panel
    """
    # Update the config with CLI params
    config = Config(**locals())

    # Check the output directory
    check_output_dir(output, force)

    # Set up the progress manager
    pm = ProgressManager()

    panelcreate(
        msa=msa,
        output_dir=output,
        region_bedfile=region_bedfile,
        input_bedfile=input_bedfile,
        mode=mode,
        config=config,
        pm=pm,
        max_amplicons=max_amplicons,
        max_amplicons_msa=max_amplicons_msa,
        max_amplicons_region_group=max_amplicons_region_group,
        force=force,
        offline_plots=offline_plots,
    )


@app.command(no_args_is_help=True)
def interactions(
    bedfile: Annotated[
        pathlib.Path,
        typer.Argument(
            help="Path to the bedfile",
            exists=True,
            readable=True,
            resolve_path=True,
        ),
    ],
    threshold: Annotated[
        float,
        typer.Option(
            help="Only show interactions more severe (Lower score) than this value",
        ),
    ] = -26.0,
):
    """
    Shows all the primer-primer interactions within a bedfile
    """
    visualise_interactions(bedfile, threshold)


@app.command(no_args_is_help=True)
def repair_mode(
    bedfile: Annotated[
        pathlib.Path,
        typer.Option(
            help="Path to the bedfile",
            exists=True,
            readable=True,
            resolve_path=True,
        ),
    ],
    msa: Annotated[
        pathlib.Path,
        typer.Option(
            help="An MSA, with the reference.fasta, aligned to any new genomes with mutations",
            exists=True,
            readable=True,
            resolve_path=True,
        ),
    ],
    config: Annotated[
        pathlib.Path,
        typer.Option(
            help="Path to the config.json",
            exists=True,
            readable=True,
            resolve_path=True,
        ),
    ],
    output: Annotated[
        pathlib.Path, typer.Option(help="The output directory", dir_okay=True)
    ],
    force: Annotated[bool, typer.Option(help="Override the output directory")] = False,
):
    """
    Repairs a primer scheme via adding more primers to account for new mutations
    """
    # Set up the progress manager
    pm = ProgressManager()

    repair(
        config_path=config,
        bedfile_path=bedfile,
        force=force,
        pm=pm,
        output_dir=output,
        msa_path=msa,
    )


@app.command(no_args_is_help=True)
def visualise_primer_mismatches(
    msa: Annotated[
        pathlib.Path,
        typer.Argument(
            help="The MSA used to design the scheme",
            readable=True,
            exists=True,
            dir_okay=False,
            resolve_path=True,
        ),
    ],
    bedfile: Annotated[
        pathlib.Path,
        typer.Argument(
            help="The bedfile containing the primers",
            readable=True,
            exists=True,
            dir_okay=False,
            resolve_path=True,
        ),
    ],
    output: Annotated[
        pathlib.Path,
        typer.Option(help="Output location of the plot", dir_okay=False, writable=True),
    ] = pathlib.Path("primer.html"),
    include_seqs: Annotated[
        bool,
        typer.Option(help="Reduces plot filesize, by excluding primer sequences"),
    ] = True,
    offline_plots: Annotated[
        bool,
        typer.Option(
            help="Includes 3Mb of dependencies into the plots, so they can be viewed offline"
        ),
    ] = True,
):
    """
    Visualise mismatches between primers and the input genomes
    """

    array, seqdict = parse_msa(msa)

    with open(output, "w") as outfile:
        outfile.write(
            primer_mismatch_heatmap(
                array=array,
                seqdict=seqdict,
                bedfile=bedfile,
                offline_plots=offline_plots,
                include_seqs=include_seqs,
            )
        )


@app.command(no_args_is_help=True)
def visualise_bedfile(
    bedfile: Annotated[
        pathlib.Path,
        typer.Argument(
            help="The bedfile containing the primers",
            readable=True,
            exists=True,
            dir_okay=False,
            resolve_path=True,
        ),
    ],
    ref_id: Annotated[str, typer.Option(help="The reference genome ID")],
    ref_path: Annotated[
        pathlib.Path,
        typer.Argument(
            help="The bedfile containing the primers",
            readable=True,
            exists=True,
            dir_okay=False,
            resolve_path=True,
        ),
    ],
    output: Annotated[
        pathlib.Path,
        typer.Option(help="Output location of the plot", dir_okay=False, writable=True),
    ] = pathlib.Path("bedfile.html"),
):
    """
    Visualise the bedfile
    """
    import dnaio

    refs = []
    ref_genome: str | None = None

    # Find the wanted file
    with dnaio.open(ref_path, mode="r") as ref_file:
        for record in ref_file:
            refs.append(record.id)

            # See if wanted genome
            if record.id == ref_id:
                ref_genome = record.sequence

    if ref_genome is None:
        raise typer.BadParameter(
            f"Reference genome ID '{ref_id}' not found in '{ref_path}'. Options: {', '.join(refs)}"
        )

    with open(output, "w") as outfile:
        outfile.write(
            bedfile_plot_html(bedfile=bedfile, ref_name=ref_id, ref_seq=ref_genome)
        )


if __name__ == "__main__":
    app()
