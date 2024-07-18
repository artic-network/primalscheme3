import pathlib

import numpy as np
import plotly.graph_objects as go

# Create in the classes from primalscheme3
from primalscheme3.core.bedfiles import BedLine, read_in_bedlines
from primalscheme3.core.mapping import create_mapping
from primalscheme3.core.seq_functions import reverse_complement


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
        return f"5'{self.primer_seq}: {self.primer_name}<br>5'{cigar}<br>5'{self.genome_seq[-len(self.primer_seq):]}"


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
                inital_slice = array[row, index - length - gaps : index]
                # Check for gaps on set base
                if inital_slice[-1] == "-":
                    break
                sequence = "".join(inital_slice).replace("-", "")
                # Covered removed gaps
                if "" in inital_slice:
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
                inital_slice = array[row, index : index + length + gaps]
                # Check for gaps on set base
                if inital_slice[0] == "-":
                    break
                sequence = "".join(inital_slice).replace("-", "")
                # Covered removed gaps
                if "" in inital_slice:
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
        if seq1b != seq2b and (seq1b != "N" and seq2b != "N"):
            dif += 1
    return dif


def primer_mismatch_heatmap(
    array: np.ndarray,
    seqdict: dict,
    bedfile: pathlib.Path,
    outpath: pathlib.Path | None = None,
    show_plot: bool = False,
    include_seqs: bool = True,
    offline_plots: bool = True,
):
    """
    Create a heatmap of primer mismatches in an MSA.
    :param array: The MSA array.
    :param seqdict: The sequence dictionary.
    :param bedfile: The bedfile of primers.
    :param outpath: The output path for the plot.
    :param show_plot: Show the plot.
    :param minimal: Reduces plot size by removing hovertext.
    """
    # Read in the bedfile
    bedlines, _header = read_in_bedlines(bedfile)

    # Find the mapping genome
    bed_chrom_names = {bedline.chrom_name for bedline in bedlines}

    # Reference genome
    primary_ref = bed_chrom_names.intersection(seqdict.keys())
    if len(primary_ref) == 0:
        # Try to fix a common issue with Jalview
        parsed_seqdict = {k.split("/")[0]: v for k, v in seqdict.items()}
        primary_ref = bed_chrom_names.intersection(parsed_seqdict.keys())

        if len(primary_ref) == 0:
            raise ValueError(
                f"Chrom names ({', '.join(bed_chrom_names)}) not found in MSA"
            )
        else:  # Update seqdict
            seqdict = parsed_seqdict

    # Filter the bedlines for only the reference genome
    bedlines = [bedline for bedline in bedlines if bedline.chrom_name in primary_ref]

    kmers_names = [bedline.primername for bedline in bedlines]

    # Create mapping array
    # Find index of primary ref
    mapping_index = [i for i, (k, v) in enumerate(seqdict.items()) if k in primary_ref][
        0
    ]

    mapping_array, array = create_mapping(array, mapping_index)
    ref_index_to_msa = {
        x: i for i, x in enumerate(list(mapping_array)) if x is not None
    }

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
            msa_index = ref_index_to_msa[list(lines)[0].end]
        else:
            forward = False
            msa_index = ref_index_to_msa[list(lines)[0].start]

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
        )
    )
    fig.update_layout(
        font=dict(family="Courier New, monospace"),
        hoverlabel=dict(font_family="Courier New, monospace"),
    )
    fig.update_yaxes(autorange="reversed")

    # Output the plot
    if outpath:
        if outpath.suffix != ".html":
            outpath = outpath.with_suffix(".html")
        fig.write_html(str(outpath), include_plotlyjs=True if offline_plots else "cdn")

    if show_plot:
        fig.show()
