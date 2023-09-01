import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.offline.offline import plot

import numpy as np
from collections import Counter
from Bio import SeqIO
from itertools import groupby
from operator import itemgetter
import pathlib


class Primer:
    chrom: str
    start: int
    end: int
    name: str
    pool: int
    strand: str
    sequence: str
    _amplicon_number: int
    _primer_number: int

    def __init__(self, bed_line: list[str]):
        self.chrom = bed_line[0]
        self.start = int(bed_line[1])
        self.end = int(bed_line[2])
        self.name = bed_line[3]
        self.pool = int(bed_line[4])
        self.strand = bed_line[5]
        self.sequence = bed_line[6]

        # Auto generate the amplicon number
        self._amplicon_number = int(self.name.split("_")[1])
        self._primer_number = int(self.name.split("_")[3])

    def to_bed(self) -> str:
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.name}\t{self.pool}\t{self.strand}\t{self.sequence}"

    def __hash__(self) -> int:
        return hash(self.to_bed())


class Amplicon:
    fprimer: list[Primer]
    rprimer: list[Primer]

    def __init__(self, fprimer, rprimer) -> None:
        self.fprimer = fprimer
        self.rprimer = rprimer

    @property
    def start(self) -> int:
        return min(x.start for x in self.fprimer)

    @property
    def end(self) -> int:
        return max(x.end for x in self.rprimer)

    @property
    def coverage_start(self) -> int:
        return max(x.end for x in self.fprimer)

    @property
    def coverage_end(self) -> int:
        return min(x.start for x in self.rprimer)

    @property
    def number(self) -> int:
        return self.fprimer[0]._amplicon_number

    @property
    def chrom(self) -> int:
        return self.fprimer[0].chrom

    @property
    def pool(self) -> int:
        return self.fprimer[0].pool


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


def calc_occupancy(align_array) -> list[float]:
    results = []
    # Calculate the base proportions
    for index, column in enumerate(align_array.T):
        gaps = np.count_nonzero(column == "-")
        results.append((index, 1 - (gaps / len(column))))
    return results


def calc_gc(align_array, kmer_size=30) -> list[float]:
    results = []
    # Calculate the base proportions
    for col_index in range(0, align_array.shape[1] - kmer_size, 15):
        slice = align_array[:, col_index : col_index + kmer_size]
        ng = np.count_nonzero(slice == "G")
        nc = np.count_nonzero(slice == "C")

        n_invalid = np.count_nonzero(slice == "-")
        gc_prop = (ng + nc) / ((len(slice) * kmer_size) - n_invalid)

        results.append((col_index, gc_prop))
    return results


def calc_variance(align_array, kmer_size=30) -> list[float]:
    results = []
    # Calculate the base proportions
    for col_index in range(0, align_array.shape[1] - kmer_size, 5):
        slice = align_array[:, col_index : col_index + kmer_size]
        seqs = {"".join(x) for x in slice}
        results.append((col_index, len(seqs)))
    return results


def read_bedfile(bedfile: str) -> list[Primer]:
    primers = []
    with open(bedfile, "r") as f:
        for line in f:
            primers.append(Primer(line.strip().split("\t")))
    return primers


def generate_all_amplicons(
    primers: list[Primer],
) -> list[Amplicon]:
    """Will group primers into amplicons based on the primer name."""
    return_amplicons = []

    # Group the primers by chrom
    chroms = dict()
    for primer in primers:
        if primer.chrom not in chroms:
            chroms[primer.chrom] = []
        chroms[primer.chrom].append(primer)

    # Group the primers by amplicon
    amplicons = dict()
    for _chrom, primers in chroms.items():
        amplicon_numbers = dict()
        for primer in primers:
            if primer._amplicon_number not in amplicon_numbers:
                amplicon_numbers[primer._amplicon_number] = []
            amplicon_numbers[primer._amplicon_number].append(primer)

        # Create the amplicons
        for _amplicon_number, amplicon_primers in amplicon_numbers.items():
            fprimers = [x for x in amplicon_primers if x.strand == "+"]
            rprimers = [x for x in amplicon_primers if x.strand == "-"]
            return_amplicons.append(Amplicon(fprimers, rprimers))

    return_amplicons.sort(key=lambda x: (x.chrom, x.number))
    return return_amplicons


def create_plots(
    primer_bed_path: str, outdir: pathlib.Path, chrom_to_msapath: dict[str:str]
):
    # Read the bedfile
    primers = read_bedfile(primer_bed_path)
    # Generate the amplicons
    amplicons = generate_all_amplicons(primers)

    # Group amplicons by Chrom
    amplicons_by_chrom = dict()
    for amplicon in amplicons:
        if amplicon.chrom not in amplicons_by_chrom:
            amplicons_by_chrom[amplicon.chrom] = []
        amplicons_by_chrom[amplicon.chrom].append(amplicon)

    # Generate the plots
    for chrom, amplicons in amplicons_by_chrom.items():
        generate_plot(chrom, amplicons, chrom_to_msapath[chrom], outdir)


def generate_plot(
    msa_name: str, amplicons: list[Amplicon], msa_path: str, outdir: pathlib.Path
):
    # Read the MSA
    seq_index = SeqIO.index(msa_path, "fasta")
    align_array = np.array(
        [record.seq.upper() for id, record in seq_index.items()], dtype=str
    )
    npools = max({x.pool for x in amplicons})
    length = list({len(record.seq) for id, record in seq_index.items()})

    # Generate the uncovered_regions
    # Add regions of no coverage
    uncovered_indexes = {x for x in range(0, length[0])}
    for amplicon in amplicons:
        uncovered_indexes -= set(range(amplicon.coverage_start, amplicon.coverage_end))
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
        rows=3,
        shared_xaxes=True,
        vertical_spacing=0.02,
        row_heights=[0.2, 0.2, 0.2],
        specs=[
            [{"secondary_y": False}],
            [{"secondary_y": True}],
            [{"secondary_y": False}],
        ],
    )
    fig.add_trace(
        go.Scatter(
            x=[x.start for x in amplicons],
            y=[x.pool for x in amplicons],
            opacity=0,
            name="FPrimers",
            hovertext=[x.fprimer[0].name for x in amplicons],
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=[x.end for x in amplicons],
            y=[x.pool for x in amplicons],
            opacity=0,
            name="RPrimers",
            hovertext=[x.rprimer[0].name for x in amplicons],
        ),
        row=1,
        col=1,
    )

    # Plot the amplicons lines
    for amplicon in amplicons:
        fig.add_shape(
            type="line",
            y0=amplicon.pool,
            y1=amplicon.pool,
            x0=amplicon.coverage_start,
            x1=amplicon.coverage_end,
            line=dict(color="LightSeaGreen", width=5),
            row=1,
            col=1,
        )
        fig.add_shape(
            type="rect",
            y0=amplicon.pool - 0.05,
            y1=amplicon.pool + 0.05,
            x0=amplicon.start,
            x1=amplicon.coverage_start,
            fillcolor="LightSalmon",
            line=dict(color="LightSalmon", width=2),
            row=1,
            col=1,
        )
        fig.add_shape(
            type="rect",
            y0=amplicon.pool - 0.05,
            y1=amplicon.pool + 0.05,
            x0=amplicon.coverage_end,
            x1=amplicon.end,
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
            row=1,
            col=1,
        )

    # Add the base occupancy
    occupancy = calc_occupancy(align_array)
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
    gc_prop = calc_gc(align_array)
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

    # Add the bottom plot
    varience = calc_variance(align_array)
    fig.add_trace(
        go.Scatter(
            x=[x[0] for x in varience],
            y=[x[1] for x in varience],
            opacity=1,
            name="Sequence variance",
            mode="lines",
        ),
        row=3,
        col=1,
    )

    # Update the top plot
    fig.update_xaxes(
        title="position",
        range=[0, length[0]],
        showline=True,
        mirror=True,
        ticks="outside",
        linewidth=2,
        linecolor="black",
        row=1,
        col=1,
        title_font=dict(size=18, family="Arial", color="Black"),
    )
    fig.update_yaxes(
        range=[0.5, npools + 0.5],
        title="pool",
        tickmode="array",
        tickvals=sorted({x.pool for x in amplicons}),
        showline=True,
        mirror=True,
        ticks="outside",
        linewidth=2,
        linecolor="black",
        row=1,
        col=1,
        title_font=dict(size=18, family="Arial", color="Black"),
    )
    # Update the bottom plot
    fig.update_xaxes(
        title="position",
        range=[0, length[0]],
        showline=True,
        mirror=True,
        ticks="outside",
        linewidth=2,
        linecolor="black",
        row=2,
        col=1,
        title_font=dict(size=18, family="Arial", color="Black"),
    )
    fig.update_yaxes(
        range=[-0.1, 1.1],
        title="Base Occupancy",
        tickmode="array",
        tickvals=[0, 0.25, 0.5, 0.75, 1],
        showline=True,
        mirror=True,
        ticks="outside",
        linewidth=2,
        linecolor="black",
        row=2,
        col=1,
        secondary_y=False,
        title_font=dict(size=18, family="Arial", color="#F0605D"),
    )
    fig.update_yaxes(
        range=[-0.1, 1.1],
        title="GC%",
        tickmode="array",
        tickvals=[0, 0.25, 0.5, 0.75, 1],
        ticktext=[0, 25, 50, 75, 100],
        showline=True,
        linewidth=2,
        linecolor="black",
        row=2,
        col=1,
        secondary_y=True,
        side="right",
        title_font=dict(size=18, family="Arial", color="#005c68"),
    )
    # Update the bottom plot
    fig.update_xaxes(
        title="position",
        range=[0, length[0]],
        showline=True,
        mirror=True,
        ticks="outside",
        linewidth=2,
        linecolor="black",
        row=3,
        col=1,
        title_font=dict(size=18, family="Arial", color="Black"),
    )
    fig.update_yaxes(
        title="Base Variance",
        tickmode="array",
        showline=True,
        mirror=True,
        ticks="outside",
        linewidth=2,
        linecolor="black",
        row=3,
        col=1,
        secondary_y=False,
        title_font=dict(size=18, family="Arial", color="Black"),
    )
    # fig.update_layout(paper_bgcolor="#000000")
    fig.update_layout(height=900, title_text=msa_name, showlegend=False)
    # plot_bgcolor="rgba(246, 237, 202, 0.5)",

    plot(fig, filename=str(outdir.absolute() / (msa_name + ".html")), auto_open=False)
    fig.write_image(
        str(outdir.absolute() / (msa_name + ".png")),
        format="png",
        height=900,
        width=1600,
    )


if __name__ == "__main__":
    primer_bed_path = "/Users/kentcg/schemes/lassa-v2/lassa-v2-0.02/output.primer.bed"
    chrom_to_msapath = {
        "L": "/Users/kentcg/schemes/lassa-v2/L.align.mod.fasta",
        "S": "/Users/kentcg/schemes/lassa-v2/S.align.fasta",
    }
    outdir = pathlib.Path("/Users/kentcg/schemes/lassa-v2/lassa-v2-0.02/")
    create_plots(primer_bed_path, outdir, chrom_to_msapath)
