# Downsample to count
import csv
import itertools
import math
import os

import matplotlib.pyplot as plt
import networkx as nx
from primalbedtools.bedfiles import group_amplicons
from primalbedtools.scheme import Scheme
from primalschemers import FKmer, RKmer

from primalscheme3.core.config import PRIMER_COUNT_ATTR_STRING, Config
from primalscheme3.core.primerscorer import PrimerScorer
from primalscheme3.core.seq_functions import primer_edit_dist
from primalscheme3.core.thermo import THERMO_RESULT, calc_annealing_hetro, thermo_check

MIN_ANNEALING = 3
NAIVE_ANNEALING = 20
TARGET_ANNEALING = 20


VISUALISE_LABEL_LEN = 20

# Edge-length multiplier for the spring layout. Higher = longer edges => connected
# nodes sit further apart, so the fixed-size markers stop overlapping. ~10 clears
# overlaps for typical clouds; lower it if clusters look too sparse.
SPRING_K = 14.0

# Node marker AREA is proportional to the raw read count of that sequence, so
# abundance is comparable across clouds and rare sequences read as tiny dots.
# area = NODE_MIN_AREA + NODE_AREA_PER_COUNT * count  (matplotlib size == area).
NODE_MIN_AREA = 8.0
NODE_AREA_PER_COUNT = 1.0


def _layout(G):
    """
    Organic force-directed (Fruchterman-Reingold) layout. A large SPRING_K makes
    edges long relative to the marker size, which is what actually stops nodes
    overlapping; disconnected clusters push apart on their own.
    """
    n = G.number_of_nodes()
    if n == 0:
        return {}
    return nx.spring_layout(
        G,
        seed=42,
        k=SPRING_K / math.sqrt(n),
        iterations=400,
    )


def visualise_graph(
    G,
    seq_counts,
    seq_thermo_dict,
    pos,
    added_seqs=None,
    accounted_seqs=None,
    save_path=None,
    score_cache=None,
    max_dct=None,
    covered_prop=None,
):
    # Colourblind-safe palette (Okabe-Ito)
    C_PASS = "#009E73"  # bluish green  -> thermo pass
    C_FAIL = "#D55E00"  # vermillion    -> thermo fail
    C_ADDED_EDGE = "#E69F00"  # orange   -> selected/added primer
    C_ACC_EDGE = "#0072B2"  # blue       -> accounted for
    C_EDGE = "#000000"

    # Node area == raw read count (see NODE_* constants), so a dominant sequence
    # is a big disc and rare ones are small dots.
    node_sizes = [
        NODE_MIN_AREA + NODE_AREA_PER_COUNT * seq_counts[node] for node in G.nodes()
    ]

    node_colors = []
    node_edge_colors = []
    node_edge_widths = []
    for node in G.nodes():
        passed = seq_thermo_dict[node] == THERMO_RESULT.PASS
        is_added = added_seqs is not None and node in added_seqs
        is_accounted = accounted_seqs is not None and node in accounted_seqs

        # Fill = thermo result; darker/opaque when added
        node_colors.append(C_PASS if passed else C_FAIL)
        if is_added:
            node_edge_colors.append(C_ADDED_EDGE)
            node_edge_widths.append(2.5)
        elif is_accounted:
            node_edge_colors.append(C_ACC_EDGE)
            node_edge_widths.append(2.0)
        else:
            node_edge_colors.append("white")
            node_edge_widths.append(0.6)

    fig, ax = plt.subplots(figsize=(8, 8), dpi=200)

    directional = score_cache is not None and max_dct is not None
    if directional:
        # Directed edges: an arrow primer -> template shows how well that primer
        # would amplify that (edit-distance-1) variant. Scores are asymmetric, so
        # each undirected edge becomes two arrows. Solid = amplifies within the
        # ΔCt cutoff (accounts for the target); dashed = below cutoff.
        DG = nx.DiGraph()
        DG.add_nodes_from(G.nodes())
        covered, uncovered = [], []
        for query in G.nodes():
            for ref, dct in score_cache.get(query, {}).items():
                if ref == query or ref not in G or dct is None:
                    continue
                (covered if dct >= max_dct else uncovered).append((query, ref))
        for elist, style, col, a in (
            (covered, "solid", "#333333", 0.6),
            (uncovered, "dashed", "#c00000", 0.45),
        ):
            if elist:
                nx.draw_networkx_edges(
                    DG,
                    pos,
                    ax=ax,
                    edgelist=elist,
                    edge_color=col,
                    style=style,
                    width=0.9,
                    alpha=a,
                    arrows=True,
                    arrowsize=8,
                    arrowstyle="-|>",
                    connectionstyle="arc3,rad=0.08",
                    node_size=node_sizes,  # type: ignore
                )
    else:
        edge_weights = [G[u][v].get("score", 5) for u, v in G.edges()]
        if edge_weights:
            nx.draw_networkx_edges(
                G,
                pos,
                ax=ax,
                width=[0.4 + (w / MIN_ANNEALING) * 0.6 for w in edge_weights],  # type: ignore
                alpha=0.35,
                edge_color=C_EDGE,
            )

    nx.draw_networkx_nodes(
        G,
        pos,
        ax=ax,
        node_size=node_sizes,  # type: ignore
        node_color=node_colors,  # type: ignore
        edgecolors=node_edge_colors,  # type: ignore
        linewidths=node_edge_widths,  # type: ignore
    )

    ax.set_aspect("equal")  # keep markers circular + spacing uniform in x/y
    ax.set_axis_off()
    ax.margins(0.08)

    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    def _dot(edge, lw, label):
        return Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#bbbbbb",
            markeredgecolor=edge,
            markeredgewidth=lw,
            markersize=11,
            label=label,
        )

    legend_elements = [
        Patch(facecolor=C_PASS, edgecolor="none", label="Thermo pass"),
        Patch(facecolor=C_FAIL, edgecolor="none", label="Thermo fail"),
        _dot(C_ADDED_EDGE, 2.5, "Selected (added)"),
        _dot(C_ACC_EDGE, 2.0, "Accounted for"),
    ]
    if directional:
        legend_elements += [
            Line2D(
                [0], [0], color="#333333", lw=1.2, label="primer→variant: amplifies"
            ),
            Line2D(
                [0],
                [0],
                color="#c00000",
                lw=1.2,
                ls="--",
                label="primer→variant: below cutoff",
            ),
        ]
    ax.legend(
        handles=legend_elements,
        loc="lower right",
        frameon=False,
        fontsize=9,
        handletextpad=0.6,
    )
    footer = "Node area ∝ read count"
    if covered_prop is not None:
        footer += f"   (total diversity covered = {covered_prop:.1%})"
    ax.text(
        0.0,
        0.0,
        footer,
        transform=ax.transAxes,
        fontsize=8,
        color="#666666",
    )

    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path, format="svg", bbox_inches="tight")
        plt.close()
    else:
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


def update_delta_ct(old: dict[str, None | float], new: dict[str, None | float]):
    """
    Merge delta-Ct values, ignoring None (unknown) scores.
    """
    for k, v in new.items():
        if v is None:
            continue
        if k not in old or old[k] is None:
            old[k] = v
        else:
            old[k] = max([old[k], v])  # type: ignore  TODO: creating the additive scoring model
    return old


def accounted_count(
    seq_counts: dict[str, float],
    seq_delta_ct: dict[str, None | float],
    max_dct: float,
) -> float:
    """
    Sum counts for sequences with known delta-Ct at or above the threshold.
    """
    return sum(
        seq_counts[s] for s, v in seq_delta_ct.items() if v is not None and v >= max_dct
    )


def coverage_prop(
    seq_counts: dict[str, float],
    seq_delta_ct: dict[str, None | float],
    max_dct: float,
) -> float:
    """
    Coverage proportion = accounted_count / total_count.
    """
    total = sum(seq_counts.values())
    if total == 0:
        return 0.0
    return accounted_count(seq_counts, seq_delta_ct, max_dct) / total


def delta_ct_downsample_seqs(
    seq_counts: dict[str, float],
    config: Config,
    ps: PrimerScorer,
    visualise: bool = False,
    outdir: str | None = None,
    name: str = "downsample",
):
    """
    Takes a raw Kmer object from PrimalSchemers Kmer. Contains the un-thermo-checked sequence and
    the count of how often that sequence appeared in the MSA
    """
    # Sort seq_counts so the most common sequences are first
    seq_counts = {
        k: v for k, v in sorted(seq_counts.items(), key=lambda x: x[1], reverse=True)
    }

    # Thermo check each sequence.
    seq_thermo_dict = {seq: thermo_check(seq, config) for seq in seq_counts.keys()}

    valid_edit_dist = []
    for n1, n2 in itertools.combinations(seq_counts.keys(), 2):
        if primer_edit_dist(n1, n2) == 1:
            valid_edit_dist.append((n1, n2))

    # Create a graph
    G = nx.Graph()
    G.add_nodes_from(seq_counts.keys())
    # Add edges for passing edit distance
    G.add_edges_from(valid_edit_dist)
    pos = None

    # Pre-calculate all mismatch scores for edges in the graph
    edge_list = list(G.edges())
    # For each edge (u, v), we need to score both ways: (u, v) and (v, u)
    score_pairs = []
    for u, v in edge_list:
        score_pairs.append((u, v))  # ref=u, query=v
        score_pairs.append((v, u))  # ref=v, query=u

    scores_flat = ps.score_seq_batch(score_pairs, config)

    # Store scores in a nested dictionary for fast lookup: score_cache[query][ref]
    score_cache: dict[str, dict[str, float | None]] = {
        s: {s: 0.0} for s in seq_counts.keys()
    }
    for (ref, query), score in zip(score_pairs, scores_flat, strict=True):
        if query not in score_cache:
            score_cache[query] = {query: 0.0}
        score_cache[query][ref] = score

    all_seqs_count = sum(seq_counts.values())

    always_add_count = config.downsample_always_add_prop * all_seqs_count

    # Start downsampling
    added_seqs = set()
    # Sequences that cannot be added to the reaction
    failed_seqs = set(
        seq for seq, status in seq_thermo_dict.items() if status != THERMO_RESULT.PASS
    )
    ordered_seqs = list(
        seq_counts.keys()
    )  # Preserve count-sorted order for deterministic ties
    order_index = {s: i for i, s in enumerate(ordered_seqs)}
    possible_seqs = set(seq_thermo_dict.keys()) - failed_seqs

    # Keep track of deltaCt
    seq_delta_ct: dict[str, None | float] = {s: None for s in seq_thermo_dict.keys()}

    # Add all valid sequences above threshold.
    for seq, count in seq_counts.items():
        # Skip failed or already added
        if seq not in failed_seqs and count >= always_add_count:
            # Calc dct for neighbors
            s = score_cache[seq]
            seq_delta_ct = update_delta_ct(seq_delta_ct, s)
            added_seqs.add(seq)
            possible_seqs.remove(seq)

    # Loop over other sequences to add the best for each iter.
    # While proportion of seqs covered above MAX_DCT is less than required
    while (
        coverage_prop(
            seq_counts,
            seq_delta_ct,
            config.downsample_max_dct,
        )
        < config.downsample_target
    ):
        # Sort the sequences to find the ones with highest count-weighted dCt
        tmp_scores = {}
        for s in ordered_seqs:
            if s not in possible_seqs:
                continue
            # Score the primer from cache
            scores = score_cache[s]

            # Filter scores to above threshold
            filtered_scores = {
                ref: dct
                for ref, dct in scores.items()
                if dct is not None and dct >= config.downsample_max_dct
            }

            # Calculate the score based on unaccounted count
            score = 0
            for i in filtered_scores:
                if (
                    seq_delta_ct[i] is None
                    or seq_delta_ct[i] < config.downsample_max_dct
                ):  # type: ignore
                    score += seq_counts[i]

            if score > 0:
                tmp_scores[s] = score

        if not tmp_scores:
            break

        # Deterministic tie-break: prefer earlier sequences by count order
        tmp_scores_list = [
            s
            for s, _ in sorted(
                tmp_scores.items(), key=lambda kv: (-kv[1], order_index[kv[0]])
            )
        ]

        added_seq = False
        for tmp in tmp_scores_list:
            if tmp_scores[tmp] > 0:
                s = score_cache.get(tmp, {tmp: 0.0})
                seq_delta_ct = update_delta_ct(seq_delta_ct, s)
                added_seqs.add(tmp)
                possible_seqs.remove(tmp)
                added_seq = True
                break

        if added_seq:
            continue
        else:
            break

    accounted_count_val = accounted_count(
        seq_counts,
        seq_delta_ct,
        config.downsample_max_dct,
    )
    tmp_target_prop = accounted_count_val / all_seqs_count

    # Quick + dirty per-cloud summary CSV (pre -> post + coverage prop)
    if outdir:
        summary_path = os.path.join(outdir, "downsample_summary.csv")
        write_header = not os.path.exists(summary_path)
        with open(summary_path, "a", newline="") as fh:
            w = csv.writer(fh)
            if write_header:
                w.writerow(
                    [
                        "name",
                        "n_seqs_pre",
                        "n_seqs_post",
                        "total_count",
                        "coverage_prop",
                        "passed",
                    ]
                )
            w.writerow(
                [
                    name,
                    len(seq_counts),
                    len(added_seqs),
                    all_seqs_count,
                    round(tmp_target_prop, 4),
                    tmp_target_prop >= config.downsample_target,
                ]
            )

    if visualise:
        account_seqs = {
            s
            for s, v in seq_delta_ct.items()
            if v is not None and v > config.downsample_max_dct
        }
        save_path = os.path.join(outdir, f"{name}.svg") if outdir else None
        if pos is None:
            pos = _layout(G)
        visualise_graph(
            G,
            seq_counts,
            seq_thermo_dict,
            pos,
            added_seqs,
            account_seqs,
            save_path=save_path,  # type: ignore
            score_cache=score_cache,
            max_dct=config.downsample_max_dct,
            covered_prop=tmp_target_prop,
        )
    if tmp_target_prop >= config.downsample_target:
        return sorted(added_seqs)
    else:
        return None


def downsample_seqs(
    seq_counts: dict[str, float],
    config: Config,
    visualise: bool = False,
    outdir: str | None = None,
    name: str = "downsample",
):
    """
    Takes a raw Kmer object from PrimalSchemers Kmer. Contains the un-thermo-checked sequence and
    the count of how often that sequence appeared in the MSA
    """

    # Sort the seq_counts in descending order by value
    seq_counts = {
        k: v for k, v in sorted(seq_counts.items(), key=lambda x: x[1], reverse=True)
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
        an = calc_annealing_hetro(
            s1, s2, config
        )  # Add NN to prevent dangling end bonuses which might inhibit PCR
        # Ignore arbitrary small annealing
        if an < MIN_ANNEALING:
            continue

        scores.setdefault(s1, {})[s2] = an
        scores.setdefault(s2, {})[s1] = an

        # Add edge to graph
        G.add_edge(s1, s2, score=an)

    if visualise:
        pos = _layout(G)

    # Start downsampling
    annealing_scores = {seq: 0.0 for seq in seq_counts.keys()}
    added_seqs = set()

    failed_seqs = set(
        seq for seq, status in seq_thermo_dict.items() if status != THERMO_RESULT.PASS
    )

    annealing_weight_sum = 0
    # Add all primers above a set prop
    total_count = sum(seq_counts.values())
    count_threshold = total_count * config.downsample_always_add_prop

    for seq, count in seq_counts.items():
        if count >= count_threshold and seq not in failed_seqs:
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

    while (
        sum(
            {
                seq: min(TARGET_ANNEALING, score) * seq_counts[seq]
                for seq, score in annealing_scores.items()
            }.values()
        )
        < TARGET_ANNEALING * total_count * config.downsample_target
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
        save_path = os.path.join(outdir, f"{name}.svg") if outdir else None
        visualise_graph(
            G,
            seq_counts,
            seq_thermo_dict,
            pos,  # type: ignore
            added_seqs,
            save_path=save_path,
            covered_prop=annealing_weight_sum
            / (TARGET_ANNEALING * sum(seq_counts.values())),
        )

    # If we have not reached the target annealing, add more sequences
    if (
        annealing_weight_sum / (TARGET_ANNEALING * sum(seq_counts.values()))
        >= config.downsample_target
    ):
        return sorted(added_seqs)
    else:
        return None
    pass


def visualise_edit_dist(seq_counts: dict[str, float], seq_thermo_dict, edit_dists):
    G = nx.Graph()
    # Add all nodes
    G.add_nodes_from(seq_counts.keys())

    # Add edit dists
    for n1, n2 in edit_dists:
        G.add_edge(n1, n2)

    node_sizes = [
        seq_counts[node] for node in G.nodes()
    ]  # Adjust scaling factor as needed

    # Create node colors based on thermodynamic check results and added status
    node_colors = []

    for node in G.nodes():
        thermo_result = seq_thermo_dict[node]

        # Non-added sequences use the original color scheme
        if thermo_result == THERMO_RESULT.PASS:
            node_colors.append("lightgreen")  # Light green for good thermo
        else:
            node_colors.append("lightcoral")  # Light coral for poor thermo

    pos = nx.spring_layout(G, seed=42)
    # Plot the graph
    plt.figure(figsize=(12, 8))

    # Draw nodes with sizes based on counts and colors based on thermo + added status
    nx.draw_networkx_nodes(
        G,
        pos,
        node_size=node_sizes,  # type: ignore
        node_color=node_colors,  # type: ignore
        alpha=0.8,
    )
    nx.draw_networkx_edges(G, pos)
    plt.show()


def downsample_kmer(
    kmer: FKmer | RKmer,
    config: Config,
    ps: PrimerScorer,
    visualise: bool = False,
    outdir: str | None = None,
    name: str = "downsample",
) -> list[str] | None:
    """
    Takes a raw Kmer object from PrimalSchemers Kmer. Contains the un-thermo-checked sequence and
    the count of how often that sequence appeared in the MSA
    """
    # Sort seqs in desc
    seq_counts = {k: v for k, v in zip(kmer.seqs(), kmer.counts(), strict=True)}
    return delta_ct_downsample_seqs(
        seq_counts=seq_counts,
        config=config,
        visualise=visualise,
        ps=ps,
        outdir=outdir,
        name=name,
    )


_WORKER_PS = None
_WORKER_CFG = None


def _ds_pool_init(model_path, config):
    global _WORKER_PS, _WORKER_CFG
    _WORKER_PS = PrimerScorer(model_path)
    _WORKER_CFG = config


def _ds_pool_worker(task):
    idx, seqs, counts = task
    seq_counts = {s: c for s, c in zip(seqs, counts, strict=True)}
    kept = delta_ct_downsample_seqs(
        seq_counts=seq_counts,
        config=_WORKER_CFG,
        ps=_WORKER_PS,
        visualise=False,
        outdir=None,
    )
    return idx, kept


def downsample_scheme(path, config, visualise: bool = False):
    scheme = Scheme.from_file(path)

    # Save the per-amplicon graphs as SVG into the run's output directory
    outdir = str(config.output) if visualise else None
    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)

    new_scheme = Scheme([], [])

    for i, amps in enumerate(group_amplicons(scheme.bedlines)):
        f_bls = amps["LEFT"]

        # Build a unique, human-readable name per amplicon for the SVG files
        amp_name = f"amp{i:03d}"
        if f_bls and getattr(f_bls[0], "chrom", None):
            amp_name = f"{f_bls[0].chrom}_amp{i:03d}"

        fseq_counts = {
            bl.sequence: float(bl.attributes.get(PRIMER_COUNT_ATTR_STRING, 1.0))
            for bl in f_bls
        }
        kept_fseqs = downsample_seqs(
            seq_counts=fseq_counts,
            config=config,
            visualise=visualise,
            outdir=outdir,
            name=f"{amp_name}_LEFT",
        )

        # Find which bedlines are being kept
        for fbl in f_bls:
            if kept_fseqs is not None and fbl.sequence in kept_fseqs:
                new_scheme.bedlines.append(fbl)

        r_bls = amps["RIGHT"]
        rseq_counts = {
            bl.sequence: float(bl.attributes.get(PRIMER_COUNT_ATTR_STRING, 1.0))
            for bl in r_bls
        }
        kept_rseqs = downsample_seqs(
            seq_counts=rseq_counts,
            config=config,
            visualise=visualise,
            outdir=outdir,
            name=f"{amp_name}_RIGHT",
        )

        # Find which bedlines are being kept
        for rbl in r_bls:
            if kept_rseqs is not None and rbl.sequence in kept_rseqs:
                new_scheme.bedlines.append(rbl)

    print(new_scheme.to_str())
