import pathlib

import numpy as np
import plotly.graph_objects as go
from click import UsageError
from plotly.subplots import make_subplots
from primalbedtools.bedfiles import BedLine, BedLineParser
from primalbedtools.primerpairs import create_primerpairs

# Create in the classes from primalscheme3
from primalscheme3.core.config import Config, MappingType
from primalscheme3.core.create_report_data import calc_gc
from primalscheme3.core.mapping import (
    check_for_end_on_gap,
    create_mapping,
    fix_end_on_gap,
    ref_index_to_msa,
)
from primalscheme3.core.seq_functions import extend_ambiguous_base, reverse_complement
from primalscheme3.core.thermo import calc_annealing_profile


class PlotlyText:
    """
    A class to hold the text for a plotly heatmap.
    """

    primer_name: str
    primer_seq: str
    genome_seq: str

    def __init__(
        self,
        primer_name: str,
        primer_seq: str,
        genome_seq: str,
    ):
        self.primer_name = primer_name
        self.primer_seq = primer_seq
        self.genome_seq = genome_seq

    def format_str(self) -> str:
        # parsedseqs
        cigar = []
        for p, g in zip(self.primer_seq[::-1], self.genome_seq[::-1]):
            if p == g:
                cigar.append("|")
            else:
                cigar.append(".")
        cigar = "".join(cigar)[::-1]
        return f"5'{self.primer_seq}: {self.primer_name}<br>5'{cigar}<br>5'{self.genome_seq[-len(self.primer_seq) :]}"


def get_primers_from_msa(
    array: np.ndarray, index: int, forward: bool = True, length: int = 20, row=0
) -> dict[int, str | None]:
    """
    Get a primer from an MSA array.
    """
    row_data = {}
    if forward:
        for row in range(array.shape[0]):
            gaps = 0
            row_data[row] = None
            while index - length - gaps >= 0:
                # Get slice
                initial_slice = array[row, index - length - gaps : index]
                # Check for gaps on set base
                if initial_slice[-1] == "-":
                    break
                sequence = "".join(initial_slice).replace("-", "")
                # Covered removed gaps
                if "" in initial_slice:
                    break
                # Check for gaps in the slice
                if len(sequence) == length:
                    row_data[row] = sequence
                    break
                # Walk left
                gaps += 1
    else:
        for row in range(array.shape[0]):
            gaps = 0
            row_data[row] = None
            while index + length + gaps <= array.shape[1]:
                # Get slice
                initial_slice = array[row, index : index + length + gaps]
                # Check for gaps on set base
                if initial_slice[0] == "-":
                    break
                sequence = "".join(initial_slice).replace("-", "")
                # Covered removed gaps
                if "" in initial_slice:
                    break
                # Check for gaps in the slice
                if len(sequence) == length:
                    row_data[row] = reverse_complement(sequence)
                    break
                # Walk right
                gaps += 1
    return row_data


def calc_primer_hamming(seq1, seq2) -> int:
    """
    Calculate the hamming distance between two sequences of equal length. Ignores N.
    :param seq1: The primer sequence in 5' to 3' orientation.
    :param seq2: The primer sequence in 5' to 3' orientation.
    :return: The number of mismatches between the two sequences.
    """
    dif = 0
    for seq1b, seq2b in zip(seq1[::-1], seq2[::-1]):
        seq1b_exp = set(extend_ambiguous_base(seq1b))
        seq2b_exp = set(extend_ambiguous_base(seq2b))

        if not seq1b_exp & seq2b_exp and (seq1b != "N" and seq2b != "N"):
            dif += 1

    return dif


def primer_mismatch_heatmap(
    array: np.ndarray,
    seqdict: dict,
    bedfile: pathlib.Path,
    include_seqs: bool = True,
    offline_plots: bool = True,
    mapping: MappingType = MappingType.FIRST,
) -> str:
    """
    Create a heatmap of primer mismatches in an MSA.
    :param array: The MSA array.
    :param seqdict: The sequence dictionary.
    :param bedfile: The bedfile of primers.
    :param include_seqs: Reduces plot size by removing hovertext.
    :raises: click.UsageError
    """
    # Read in the bedfile

    _header, bedlines = BedLineParser.from_file(bedfile)

    # Find the mapping genome
    bed_chrom_names = {bedline.chrom for bedline in bedlines}

    # Reference genome
    primary_ref = bed_chrom_names.intersection(seqdict.keys())

    if len(primary_ref) == 0:
        # Try to fix a common issue with Jalview
        parsed_seqdict = {"_".join(k.split("/")): v for k, v in seqdict.items()}
        primary_ref = bed_chrom_names.intersection(parsed_seqdict.keys())
        seqdict = parsed_seqdict

    # handle errors if mapping is set to first
    if len(primary_ref) == 0 and mapping == MappingType.FIRST:
        raise UsageError(
            f"Primer chrom names ({', '.join(bed_chrom_names)}) not found in MSA ({', '.join(seqdict.keys())})"
        )
    # If consensus mapping ensure only one chrom in bedfile
    elif mapping == MappingType.CONSENSUS:
        primary_ref = ["Consensus"]
        if len(bed_chrom_names) > 1:
            raise UsageError(
                f"Primer chrom names ({', '.join(bed_chrom_names)}) not found in MSA ({', '.join(seqdict.keys())})"
            )
    else:  # mapping == MappingType.FIRST & len(primary_ref) > 0
        # Filter the bedlines for only the reference genome
        bedlines = [bedline for bedline in bedlines if bedline.chrom in primary_ref]

    kmers_names = [bedline.primername for bedline in bedlines]

    # Create mapping array
    # Find index of primary ref
    if mapping == MappingType.FIRST:
        mapping_index = [
            i for i, (k, v) in enumerate(seqdict.items()) if k in primary_ref
        ][0]
        mapping_array, array = create_mapping(array, mapping_index)
    else:
        mapping_array = np.array([x for x in range(array.shape[1])])

    ref_index_to_msa_dict = ref_index_to_msa(mapping_array)

    # Group Primers by basename
    basename_to_line: dict[str, set[BedLine]] = {
        "_".join(name.split("_")[:-1]): set() for name in kmers_names
    }
    for bedline in bedlines:
        basename = "_".join(bedline.primername.split("_")[:-1])
        basename_to_line[basename].add(bedline)

    basename_to_index = {bn: i for i, bn in enumerate(basename_to_line.keys())}

    seq_to_primername = {line.sequence: line.primername for line in bedlines}

    # Create the scoremap
    scoremap = np.empty((array.shape[0], len(basename_to_line)))
    scoremap.fill(None)
    textmap = np.empty((array.shape[0], len(basename_to_line)), dtype="str")
    textmap.fill("None")
    textmap = textmap.tolist()

    # get FPrimer sequences for each basename
    for bn, lines in basename_to_line.items():
        # Get primer size
        primer_len_max = max(len(line.sequence) for line in lines)

        # Set the direction
        if "LEFT" in bn:
            forward = True
            primer_end = list(lines)[0].end
            # Check for the end on a gap edge case and fix it
            if check_for_end_on_gap(ref_index_to_msa_dict, primer_end):
                msa_index = fix_end_on_gap(ref_index_to_msa_dict, primer_end)
            else:
                msa_index = ref_index_to_msa_dict[list(lines)[0].end]
        else:
            forward = False
            msa_index = ref_index_to_msa_dict[list(lines)[0].start]

        # Get the primer sequences
        msa_data = get_primers_from_msa(array, msa_index, forward, primer_len_max)

        # Get the score for each genome
        primer_seqs = {line.sequence for line in lines}

        for genome_index, genome_seq in msa_data.items():
            # Caused by gaps in the msa
            if genome_seq is None:
                if forward:
                    slice = array[genome_index, msa_index - primer_len_max : msa_index]
                    slice[slice == ""] = "-"
                    genome_seq = "".join(slice)
                else:
                    slice = array[genome_index, msa_index : msa_index + primer_len_max]
                    slice[slice == ""] = "-"
                    genome_seq = reverse_complement("".join(slice))

                textmap[genome_index][basename_to_index[bn]] = PlotlyText(
                    primer_seq=[x for x in primer_seqs][0],
                    genome_seq="".join(slice),
                    primer_name=bn,
                ).format_str()
                continue
            # Quick check for exact match
            if genome_seq in primer_seqs:
                scoremap[genome_index, basename_to_index[bn]] = 0
                primer_seq = "".join(primer_seqs.intersection({genome_seq}))
                textmap[genome_index][basename_to_index[bn]] = PlotlyText(
                    primer_seq=primer_seq,
                    genome_seq=genome_seq,
                    primer_name=seq_to_primername.get(primer_seq, "Unknown"),
                ).format_str()
                continue
            # Calculate the hamming distance between all
            seq_to_scores: dict[str, int] = {}
            for primer_seq in primer_seqs:
                seq_to_scores[primer_seq] = calc_primer_hamming(primer_seq, genome_seq)
            scoremap[genome_index, basename_to_index[bn]] = min(
                seq_to_scores.values(),  # type: ignore
            )
            primer_seq = "".join(
                [
                    k
                    for k, v in seq_to_scores.items()
                    if v == scoremap[genome_index, basename_to_index[bn]]
                ][0]
            )
            textmap[genome_index][basename_to_index[bn]] = PlotlyText(
                genome_seq=genome_seq,
                primer_seq=primer_seq,
                primer_name=seq_to_primername.get(primer_seq, "Unknown"),
            ).format_str()

    # Hovertemplate string
    if include_seqs:
        hovertemplatestr = "%{text}<br>" + "<b>Mismatches: %{z}</b><br>"
    else:
        hovertemplatestr = ""

    # Create the heatmap
    fig = go.Figure(
        data=go.Heatmap(
            z=scoremap,
            x=list(basename_to_line.keys()),
            y=[x for x in seqdict.keys()],
            colorscale="Viridis",
            text=textmap if include_seqs else None,  # only show text if not minimal
            hovertemplate=hovertemplatestr,
            xgap=0.1,
            ygap=0.1,
            name="Primer Mismatches",
        )
    )
    fig.update_layout(
        font=dict(family="Courier New, monospace"),
        hoverlabel=dict(font_family="Courier New, monospace"),
        title_text=f"Primer Mismatches: {list(primary_ref)[0]}",
        coloraxis=dict(cmax=10, cmin=0),
    )
    fig.update_yaxes(autorange="reversed")

    # Remove unnecessary plot elements
    fig.update_layout(
        modebar_remove=[
            "select2d",
            "lasso2d",
            "select",
            "autoScale2d",
            "zoom",
            "toImage",
        ]
    )

    return fig.to_html(
        include_plotlyjs=True if offline_plots else "cdn", full_html=False
    )


def bedfile_plot_html(
    bedfile: pathlib.Path, ref_name: str, ref_seq: str, offline_plots: bool = False
) -> str:
    """
    Create a plotly heatmap from a bedfile.
    """
    # Read in the bedfile
    _header, bedlines = BedLineParser.from_file(bedfile)
    primerpairs = create_primerpairs(bedlines)

    # Filter primerpairs for the reference genome
    wanted_primerspairs = [pp for pp in primerpairs if pp.chrom == ref_name]
    if len(wanted_primerspairs) == 0:
        raise ValueError(f"No primers found for {ref_name}")

    # Calculate the GC data
    ref_array = np.array(
        [list(ref_seq.upper())],
        dtype="U1",
        ndmin=2,
    )
    gc = calc_gc(ref_array)

    # Create the plot
    fig = make_subplots(
        cols=1,
        rows=2,
        shared_xaxes=True,
        row_heights=[1, 0.5],
        specs=[
            [{"secondary_y": False}],
            [{"secondary_y": False}],
        ],
    )

    # Add the GC data
    fig.add_trace(
        go.Scattergl(
            x=[x[0] for x in gc],
            y=[x[1] for x in gc],
            mode="lines",
            name="GC Prop",
            line=dict(color="#005c68", width=2),
            fill="tozeroy",
        ),
        row=2,
        col=1,
    )

    # Add the primer lines
    for pp in wanted_primerspairs:
        print(f"Adding primer: {pp.amplicon_name}")
        fig.add_shape(
            type="rect",
            y0=pp.pool - 0.05,
            y1=pp.pool + 0.05,
            x0=pp.amplicon_start,
            x1=pp.coverage_start,
            fillcolor="LightSalmon",
            line=dict(color="darksalmon", width=2),
            row=1,
            col=1,
        )
        fig.add_shape(
            type="rect",
            y0=pp.pool - 0.05,
            y1=pp.pool + 0.05,
            x0=pp.coverage_end,
            x1=pp.amplicon_end,
            fillcolor="LightSalmon",
            line=dict(color="darksalmon", width=2),
            row=1,
            col=1,
        )
        # Handle circular genomes
        is_circular = pp.is_circular

        fig.add_shape(
            type="line",
            y0=pp.pool,
            y1=pp.pool,
            x0=pp.coverage_start,
            x1=pp.coverage_end if not is_circular else len(ref_seq),
            line=dict(color="LightSeaGreen", width=5),
            row=1,
            col=1,
        )
        if is_circular:
            fig.add_shape(
                type="line",
                y0=pp.pool,
                y1=pp.pool,
                x0=0,
                x1=pp.coverage_end,
                line=dict(color="LightSeaGreen", width=5),
                row=1,
                col=1,
            )

    fig.update_xaxes(
        showline=True,
        mirror=True,
        ticks="outside",
        linewidth=2,
        linecolor="black",
        tickformat=",d",
        title_font=dict(size=18, family="Arial", color="Black"),
        range=[0, len(ref_seq)],
        title="",  # Blank title for all x-axes
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
    pools = sorted({x.pool for x in wanted_primerspairs})
    fig.update_yaxes(
        range=[pools[0] - 0.5, pools[-1] + 0.5],
        title="pool",
        tickmode="array",
        tickvals=pools,
        row=1,
        col=1,
    )
    fig.update_yaxes(
        range=[-0.1, 1.1],
        title="GC%",
        tickmode="array",
        tickvals=[0, 0.25, 0.5, 0.75, 1],
        ticktext=[0, 25, 50, 75, 100],
        row=2,
        col=1,
    )

    # Remove unnecessary plot elements
    fig.update_layout(
        modebar_remove=[
            "select2d",
            "lasso2d",
            "select",
            "autoScale2d",
            "zoom",
            "toImage",
        ]
    )
    fig.update_layout(height=400, title_text=ref_name, showlegend=False)

    # Write a html version of the plot
    return fig.to_html(
        include_plotlyjs=True if offline_plots else "cdn", full_html=False
    )


def plot_primer_thermo_profile_html(
    bedfile: pathlib.Path, config: Config, offline_plots: bool = False
):
    """
    Plot primer annealing temperature profile
    """
    # Read in the bedfile
    _header, bedlines = BedLineParser.from_file(bedfile)
    fig = go.Figure()
    # Create a thermo profile for each primer
    for bedline in bedlines:
        profile = calc_annealing_profile(
            bedline.sequence,
            config.mv_conc,
            config.dv_conc,
            config.dntp_conc,
            config.dna_conc,
        )
        fig.add_trace(
            go.Scatter(
                y=list(profile.values()),
                x=list(profile.keys()),
                mode="lines",
                name="",
                line=dict(color="rgba(0, 0, 0, 0.3)", width=1),
                hovertemplate="%{y:.2f} @ %{x}°C <br>" + bedline.primername,
            ),
        )

    fig.update_xaxes(
        showline=True,
        mirror=True,
        ticks="outside",
        linewidth=2,
        linecolor="black",
        tickformat=",d",
        title_font=dict(size=18, family="Arial", color="Black"),
        title="Temperature (°C)",  # Blank title for all x-axes
    )
    fig.update_yaxes(
        showline=True,
        mirror=True,
        ticks="outside",
        linewidth=2,
        linecolor="black",
        fixedrange=True,
        title_font=dict(size=18, family="Arial", color="Black"),
        title="Percentage of Primer Annealed",
    )
    fig.update_layout(showlegend=False)
    return fig.to_html(
        include_plotlyjs=True if offline_plots else "cdn", full_html=False
    )
