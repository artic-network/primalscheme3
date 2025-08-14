# Downsample to count
import itertools
import pathlib

import matplotlib.pyplot as plt
import networkx as nx
from primalschemers import FKmer, RKmer  # type: ignore

# type: ignore
from primalscheme3.core.config import Config
from primalscheme3.core.msa import MSA
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.core.thermo import THERMO_RESULT, calc_annealing_hetro, thermo_check

MIN_ANNEALING = 3
NAIVE_ANNEALING = 20
TARGET_ANNEALING = 20


VISUALISE_LABEL_LEN = 20


def visualise_graph(G, seq_counts, seq_thermo_dict, pos, added_seqs=None):
    # Create node sizes based on sequence counts
    # Scale the counts to reasonable node sizes (multiply by a scaling factor)
    node_sizes = [
        seq_counts[node] * 10 for node in G.nodes()
    ]  # Adjust scaling factor as needed

    # Create node colors based on thermodynamic check results and added status
    node_colors = []
    node_edge_colors = []
    node_edge_widths = []

    for node in G.nodes():
        thermo_result = seq_thermo_dict[node]
        is_added = added_seqs is not None and node in added_seqs

        if is_added:
            # Added sequences get a special color scheme
            if thermo_result == THERMO_RESULT.PASS:
                node_colors.append("darkgreen")  # Dark green for added + good thermo
            else:
                node_colors.append("darkred")  # Dark red for added + poor thermo
            node_edge_colors.append("gold")  # Gold border for added sequences
            node_edge_widths.append(3.0)  # Thicker border for added sequences
        else:
            # Non-added sequences use the original color scheme
            if thermo_result == THERMO_RESULT.PASS:
                node_colors.append("lightgreen")  # Light green for good thermo
            else:
                node_colors.append("lightcoral")  # Light coral for poor thermo
            node_edge_colors.append("black")  # Black border for non-added sequences
            node_edge_widths.append(0.5)  # Thinner border for non-added sequences

    # Create edge weights for visualization
    edge_weights = [G[u][v]["score"] for u, v in G.edges()]

    # Plot the graph
    plt.figure(figsize=(12, 8))

    # Draw nodes with sizes based on counts and colors based on thermo + added status
    nx.draw_networkx_nodes(
        G,
        pos,
        node_size=node_sizes,  # type: ignore
        node_color=node_colors,  # type: ignore
        alpha=0.8,
        edgecolors=node_edge_colors,  # type: ignore
        linewidths=node_edge_widths,  # type: ignore
    )

    # Draw edges with thickness based on annealing scores
    if edge_weights:  # Only draw edges if they exist
        nx.draw_networkx_edges(
            G,
            pos,
            width=[(w / MIN_ANNEALING) ** 2 for w in edge_weights],  # type: ignore
            alpha=0.5,
            edge_color="gray",
        )

    # Add labels (you might want to truncate long sequences)
    node_labels = {
        node: f"{node[:VISUALISE_LABEL_LEN]}..."
        if len(node) > VISUALISE_LABEL_LEN
        else node
        for node in G.nodes()
    }
    nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=8)

    # Add edge labels showing annealing scores
    if G.edges():  # Only add edge labels if edges exist
        edge_labels = {(u, v): f"{G[u][v]['score']:.1f}" for u, v in G.edges()}
        nx.draw_networkx_edge_labels(G, pos, edge_labels, font_size=6)

    plt.title(
        "Sequence Annealing Network\n(Node size = sequence count, Node color = thermo check, Gold border = added sequences)"
    )
    plt.axis("off")

    # Add legend for color scheme
    if added_seqs is not None:
        from matplotlib.patches import Patch

        legend_elements = [
            Patch(
                facecolor="darkgreen",
                edgecolor="gold",
                linewidth=2,
                label="Added + Thermo Pass",
            ),
            Patch(
                facecolor="darkred",
                edgecolor="gold",
                linewidth=2,
                label="Added + Thermo Fail",
            ),
            Patch(
                facecolor="lightgreen",
                edgecolor="black",
                linewidth=1,
                label="Not Added + Thermo Pass",
            ),
            Patch(
                facecolor="lightcoral",
                edgecolor="black",
                linewidth=1,
                label="Not Added + Thermo Fail",
            ),
        ]
        plt.legend(handles=legend_elements, loc="upper right", bbox_to_anchor=(1, 1))

    plt.tight_layout()
    plt.show()


def calc_scores(
    scores: dict[str, dict[str, float]], seq_counts, annealing_scores, failed_seqs
) -> dict[str, int]:
    new_scores = {seq: 0 for seq in seq_counts.keys() if seq not in failed_seqs}

    for seq, counts in seq_counts.items():
        if seq in failed_seqs:
            continue

        new_scores[seq] += (
            max((NAIVE_ANNEALING - annealing_scores.get(seq, 0)), 0) * counts
        )

        for s2, score in scores.get(seq, {}).items():
            new_scores[seq] += (
                max((score - annealing_scores.get(seq, 0)), 0) * seq_counts[s2]
            )

    # Return sorted dictionary by value descending
    return {
        k: v for k, v in sorted(new_scores.items(), key=lambda x: x[1], reverse=True)
    }


def add_annealing(
    g, seq, scores, annealing_scores: dict[str, float]
) -> dict[str, float]:
    annealing_scores[seq] += NAIVE_ANNEALING  # Add self annealing
    for n in g.neighbors(seq):
        annealing_scores[n] += scores[seq][n]  # Add neighbour annealing
    return annealing_scores


def update_annealing(
    seq, annealing_scores: dict[str, float], scores: dict[str, dict[str, float]]
):
    """
    Update the annealing scores based on the current graph state.
    This function is a placeholder for future logic to update annealing scores.
    """
    annealing_scores[seq] += NAIVE_ANNEALING  # Add self annealing
    for n in scores.get(seq, {}):
        annealing_scores[n] += scores[seq][n]  # Add neighbour annealing
    return annealing_scores


def downsample_kmer(
    kmer: FKmer | RKmer, config: Config, visualise: bool = False
) -> list[str] | None:
    """
    Takes a raw Kmer object from PrimalSchemers Kmer. Contains the un-thermo-checked sequence and
    the count of how often that sequence appeared in the MSA
    """
    # Sort seqs in desc
    seq_counts = {
        k: v
        for k, v in sorted(
            zip(kmer.seqs(), kmer.counts(), strict=True), key=lambda x: -x[1]
        )
    }

    # Thermo check each sequence.
    seq_thermo_dict = {seq: thermo_check(seq, config) for seq in seq_counts.keys()}

    # Create a graph
    G = nx.Graph()
    G.add_nodes_from(seq_counts.keys())

    # Calc annealing between all combinations of sequences
    scores: dict[str, dict[str, float]] = {}
    tups = itertools.combinations(seq_counts.keys(), 2)
    for s1, s2 in tups:
        an = calc_annealing_hetro(s1, s2, config)
        # Ignore arbitrary small annealing
        if an < MIN_ANNEALING:
            continue

        scores.setdefault(s1, {})[s2] = an
        scores.setdefault(s2, {})[s1] = an

        # Add edge to graph
        G.add_edge(s1, s2, score=an)

    if visualise:
        pos = nx.spring_layout(G, seed=42)
        visualise_graph(G, seq_counts, seq_thermo_dict, pos, added_seqs=set())

    # Start downsampling
    annealing_scores = {seq: 0.0 for seq in seq_counts.keys()}

    added_seqs = set()

    failed_seqs = set(
        seq for seq, status in seq_thermo_dict.items() if status != THERMO_RESULT.PASS
    )

    annealing_weight_sum = 0
    while (
        sum(
            {
                seq: min(TARGET_ANNEALING, score) * seq_counts[seq]
                for seq, score in annealing_scores.items()
            }.values()
        )
        < TARGET_ANNEALING * sum(seq_counts.values()) * 0.95
    ):
        # Calculate scores based on current annealing
        tmp_scores = calc_scores(scores, seq_counts, annealing_scores, failed_seqs)

        if len(tmp_scores) == 0:
            break

        for seq, _score in tmp_scores.items():
            # Add seq to fail seq so its score is not considered again
            if seq in added_seqs:
                failed_seqs.add(seq)

            added_seqs.add(seq)
            annealing_scores = update_annealing(
                seq,
                annealing_scores,
                scores,
            )
            annealing_weight_sum = sum(
                {
                    seq: min(TARGET_ANNEALING, score) * seq_counts[seq]
                    for seq, score in annealing_scores.items()
                }.values()
            )
            break

        # Show final visualization with added sequences highlighted
    if visualise:
        visualise_graph(G, seq_counts, seq_thermo_dict, pos, added_seqs)  # type: ignore

    # If we have not reached the target annealing, add more sequences
    if annealing_weight_sum / (TARGET_ANNEALING * sum(seq_counts.values())) >= 0.95:
        return sorted(added_seqs)
    else:
        return None


if __name__ == "__main__":
    # Example usage
    config = Config()
    config.primer_annealing_tempc = 65
    config.primer_annealing_prop = 10
    config.min_base_freq = 0
    config.ncores = 8

    pm = ProgressManager()
    msa_obj = MSA(
        "example",
        pathlib.Path(
            "/Users/kentcg/primerschemes/primerschemes/artic-measles/400/v1.0.0/work/all_genomes.align.ds.align.repaired.fasta"
        ),
        1,
        "first",
        pm,
        config=config,
    )

    msa_obj.digest_rs(config, ([*range(0, 1000)], []))

    for fk in msa_obj.fkmers:
        downsample_kmer(fk, config, True)
