import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.offline.offline import plot

import numpy as np
from collections import Counter
from itertools import groupby
from operator import itemgetter
import pathlib

# Module imports
from primalscheme3.core.classes import PrimerPair
from primalscheme3.core.msa import MSA
from primalscheme3.core.seq_functions import entropy_score_array


def reduce_data(results: list[tuple[int, float]]) -> list[tuple[int, float]]:
    """
    Reduce the size of data by merging consecutive points
    """
    reduced_results = []
    for iindex, (index, oc) in enumerate(results):
        # Add first point
        if iindex == 0:
            reduced_results.append((index, oc))
            continue
        # Add the last point
        if iindex == len(results) - 1:
            reduced_results.append((index, oc))
            continue

        # If the previous point is the same, and the next point is the same
        if results[iindex - 1][1] == oc and results[iindex + 1][1] == oc:
            continue
        else:
            reduced_results.append((index, oc))
    return reduced_results


def calc_base_consensus(align_array) -> list[float]:
    results = []
    # Calculate the base proportions
    for index, column in enumerate(align_array.T):
        counts = Counter(column)
        results.append(
            (
                index,
                counts.most_common()[0][0],
                counts.most_common()[0][1] / len(column),
            )
        )
    return results


def calc_occupancy(align_array) -> list[tuple[int, float]]:
    results = []
    # Calculate the base proportions
    for index, column in enumerate(align_array.T):
        gaps = np.count_nonzero(column == "-")
        gaps += np.count_nonzero(column == "")
        results.append((index, 1 - (gaps / len(column))))
    return reduce_data(results)


def calc_gc(align_array, kmer_size=30) -> list[tuple[int, float]]:
    results = []
    # Calculate the base proportions
    for col_index in range(0, align_array.shape[1] - kmer_size, 15):
        slice = align_array[:, col_index : col_index + kmer_size]
        ng = np.count_nonzero(slice == "G")
        nc = np.count_nonzero(slice == "C")

        n_invalid = np.count_nonzero(slice == "-")
        gc_prop = round((ng + nc) / ((len(slice) * kmer_size) - n_invalid), 2)

        results.append((col_index, gc_prop))
    return reduce_data(results)


def calc_variance(align_array, kmer_size=30) -> list[float]:
    results = []
    # Calculate the base proportions
    for col_index in range(0, align_array.shape[1] - kmer_size, 5):
        slice = align_array[:, col_index : col_index + kmer_size]
        seqs = {"".join(x) for x in slice}
        results.append((col_index, len(seqs)))
    return results


def generate_plot(msa: MSA, scheme_pools: list[list[PrimerPair]], outdir: pathlib.Path):
    """Generate a plot for a single MSA"""
    npools = len(scheme_pools)
    length = msa.array.shape[1]

    # Filter out non MSA-Primerpairs, into a flat list
    msa_primers: list[PrimerPair] = []
    for pool in scheme_pools:
        for primerpair in pool:
            if primerpair.msa_index == msa.msa_index:
                msa_primers.append(primerpair)

    # Filter primers that are circular
    circular_pp = []
    included_primers = []
    for pp in msa_primers:
        if pp.start > pp.end:
            circular_pp.append(pp)
        else:
            included_primers.append(pp)

    # Remap the included primers to the MSA if they have been mapped to an genome
    if msa._mapping_array is not None:
        mapping_list = list(msa._mapping_array)
        for fkmer in msa.fkmers:
            fkmer.end = mapping_list.index(fkmer.end)
            fkmer._starts = {fkmer.end - len(x) for x in fkmer.seqs}
        for rkmer in msa.rkmers:
            rkmer.start = mapping_list.index(rkmer.start)
            rkmer._ends = {rkmer.start + len(x) for x in rkmer.seqs}

    # Generate the uncovered_regions
    # Add regions of no coverage
    uncovered_indexes = {x for x in range(0, length)}
    for primerpair in included_primers:
        uncovered_indexes -= set(
            range(primerpair.fprimer.end, primerpair.rprimer.start)
        )
    for cpp in circular_pp:
        uncovered_indexes -= set(range(cpp.fprimer.end, length))
        uncovered_indexes -= set(range(0, cpp.rprimer.start))

    # Plot the uncovered regions
    uncovered_indexes_list = sorted(uncovered_indexes)
    # Generate continous regions
    uncovered_regions = []
    for k, g in groupby(enumerate(uncovered_indexes_list), lambda ix: ix[0] - ix[1]):
        uncovered_regions.append(list(map(itemgetter(1), g)))
    uncovered_regions = [(min(x), max(x)) for x in uncovered_regions]

    # Create an empty figre with the Fprimer hovers
    fig = make_subplots(
        cols=1,
        rows=4,
        shared_xaxes=True,
        vertical_spacing=0.02,
        row_heights=[0.2, 0.2, 0.2, 0.1],
        specs=[
            [{"secondary_y": False}],
            [{"secondary_y": True}],
            [{"secondary_y": False}],
            [{"secondary_y": False}],
        ],
    )
    fig.add_trace(
        go.Scatter(
            x=[x.start for x in included_primers],
            y=[x.pool + 1 for x in included_primers],
            opacity=0,
            name="FPrimers",
            hovertext=[f"AMP_{x.amplicon_number}_LEFT" for x in included_primers],
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=[x.end for x in included_primers],
            y=[x.pool + 1 for x in included_primers],
            opacity=0,
            name="RPrimers",
            hovertext=[f"AMP_{x.amplicon_number}_RIGHT" for x in included_primers],
        ),
        row=1,
        col=1,
    )

    # Plot the amplicons lines
    for amplicon in included_primers:
        fig.add_shape(
            type="line",
            y0=amplicon.pool + 1,
            y1=amplicon.pool + 1,
            x0=amplicon.fprimer.end,
            x1=amplicon.rprimer.start,
            line=dict(color="LightSeaGreen", width=5),
            row=1,
            col=1,
        )
        fig.add_shape(
            type="rect",
            y0=amplicon.pool + 1 - 0.05,
            y1=amplicon.pool + 1 + 0.05,
            x0=min(amplicon.fprimer.starts()),
            x1=amplicon.fprimer.end,
            fillcolor="LightSalmon",
            line=dict(color="LightSalmon", width=2),
            row=1,
            col=1,
        )
        fig.add_shape(
            type="rect",
            y0=amplicon.pool + 1 - 0.05,
            y1=amplicon.pool + 1 + 0.05,
            x0=amplicon.rprimer.start,
            x1=max(amplicon.rprimer.ends()),
            fillcolor="LightSalmon",
            line=dict(color="LightSalmon", width=2),
            row=1,
            col=1,
        )
    # Plot the circular primers
    for pp in circular_pp:
        # Add the left side line
        fig.add_shape(
            type="line",
            y0=pp.pool + 1,
            y1=pp.pool + 1,
            x1=pp.rprimer.start,
            x0=0,
            line=dict(color="LightSeaGreen", width=5),
            row=1,
            col=1,
        )
        # Add the right side line
        fig.add_shape(
            type="line",
            y0=pp.pool + 1,
            y1=pp.pool + 1,
            x0=pp.fprimer.end,
            x1=length,
            line=dict(color="LightSeaGreen", width=5),
            row=1,
            col=1,
        )
        fig.add_shape(
            type="rect",
            y0=pp.pool + 1 - 0.05,
            y1=pp.pool + 1 + 0.05,
            x0=min(pp.fprimer.starts()),
            x1=pp.fprimer.end,
            fillcolor="LightSalmon",
            line=dict(color="LightSalmon", width=2),
            row=1,
            col=1,
        )
        fig.add_shape(
            type="rect",
            y0=pp.pool + 1 - 0.05,
            y1=pp.pool + 1 + 0.05,
            x0=pp.rprimer.start,
            x1=max(pp.rprimer.ends()),
            fillcolor="LightSalmon",
            line=dict(color="LightSalmon", width=2),
            row=1,
            col=1,
        )

    # Add the uncovered regions
    for region in uncovered_regions:
        fig.add_vrect(
            x0=region[0],
            x1=region[1],
            fillcolor="#F0605D",
            line=dict(width=0),
            opacity=0.5,
            row=1,  # type: ignore
            col=1,  # type: ignore
        )

    # Add the base occupancy
    occupancy = calc_occupancy(msa.array)
    fig.add_trace(
        go.Scatter(
            x=[x[0] for x in occupancy],
            y=[x[1] for x in occupancy],
            mode="lines",
            name="Base Occupancy",
            line=dict(color="#F0605D", width=2),
            fill="tozeroy",
            opacity=0.5,
        ),
        row=2,
        col=1,
    )
    gc_prop = calc_gc(msa.array)
    fig.add_trace(
        go.Scatter(
            x=[x[0] for x in gc_prop],
            y=[x[1] for x in gc_prop],
            mode="lines",
            name="GC Prop",
            line=dict(color="#005c68", width=2),
        ),
        row=2,
        col=1,
        secondary_y=True,
    )

    # Add the entropy plot
    entropy = entropy_score_array(msa.array)
    fig.add_trace(
        go.Scatter(
            x=[x for x in range(len(entropy))],
            y=[x for x in entropy],
            opacity=1,
            name="Sequence Entropy",
            mode="lines",
        ),
        row=3,
        col=1,
    )

    # Add all posible Fkmers
    fig.add_trace(
        go.Scatter(
            x=[fkmer.end for fkmer in msa.fkmers],
            y=[1 for _ in msa.rkmers],
            hovertext=[f"Number Seqs: {len(fkmer.seqs)}" for fkmer in msa.fkmers],
            marker=dict(symbol="triangle-right", size=10),
            mode="markers",
        ),
        row=4,
        col=1,
    )
    # Add all posible Rkmers
    fig.add_trace(
        go.Scatter(
            x=[rkmer.start for rkmer in msa.rkmers],
            y=[0.5 for _ in msa.rkmers],
            hovertext=[f"Number Seqs: {len(rkmer.seqs)}" for rkmer in msa.rkmers],
            marker=dict(symbol="triangle-left", size=10),
            mode="markers",
        ),
        row=4,
        col=1,
    )
    # Add the base plot settings
    fig.update_xaxes(
        showline=True,
        mirror=True,
        ticks="outside",
        linewidth=2,
        linecolor="black",
        tickformat=",d",
        title_font=dict(size=18, family="Arial", color="Black"),
        range=[0, length],
        title="Position",
    )
    fig.update_yaxes(
        showline=True,
        mirror=True,
        ticks="outside",
        linewidth=2,
        linecolor="black",
        fixedrange=True,
        title_font=dict(size=18, family="Arial", color="Black"),
    )

    # Update the top plot
    fig.update_yaxes(
        range=[0.5, npools + 0.5],
        title="pool",
        tickmode="array",
        tickvals=sorted({x.pool + 1 for x in included_primers}),
        row=1,
        col=1,
    )
    # Update the second plot
    fig.update_yaxes(
        range=[-0.1, 1.1],
        title="Base Occupancy",
        tickmode="array",
        tickvals=[0, 0.25, 0.5, 0.75, 1],
        row=2,
        col=1,
        secondary_y=False,
    )
    fig.update_yaxes(
        range=[-0.1, 1.1],
        title="GC%",
        tickmode="array",
        tickvals=[0, 0.25, 0.5, 0.75, 1],
        ticktext=[0, 25, 50, 75, 100],
        row=2,
        col=1,
        secondary_y=True,
        side="right",
    )
    # Update the third plot
    fig.update_yaxes(
        title="Entropy",
        tickmode="array",
        row=3,
        col=1,
    )
    # Update the fourth plot
    fig.update_yaxes(
        title="Thermo-passing Primers",
        range=[0.5 - 0.1, 1 + 0.1],
        tickmode="array",
        row=4,
        col=1,
        secondary_y=False,
    )

    # fig.update_layout(paper_bgcolor="#000000")
    fig.update_layout(height=900, title_text=msa._chrom_name, showlegend=False)
    # plot_bgcolor="rgba(246, 237, 202, 0.5)",

    plot(
        fig,
        filename=str(outdir.absolute() / (msa._chrom_name + ".html")),
        auto_open=False,
    )
    fig.write_image(
        str(outdir.absolute() / (msa._chrom_name + ".png")),
        format="png",
        height=900,
        width=1600,
    )
