import json
import pathlib
import gzip
from itertools import groupby
from operator import itemgetter

# Module imports
from primalscheme3.core.msa import MSA
from primalscheme3.core.classes import PrimerPair
from primalscheme3.core.seq_functions import entropy_score_array
from primalscheme3.core.create_reports import calc_gc, calc_occupancy, reduce_data

# Panel imports
from primalscheme3.panel.minimal_scheme_classes import PanelMSA

# Plot format
# plot1: scheme coverage plot. Can be parsed from the bedfile
# plot2 Base occupancy + genome gc
# plot3: Entropy plot
# plot4: Thermo pass Fkmer and Rkmer plot


def generate_uncovered_data(length, primerpairs: list[PrimerPair]) -> dict[int, int]:
    uncovered_indexes = {x for x in range(0, length)}
    for primerpair in primerpairs:
        uncovered_indexes -= set(
            range(primerpair.fprimer.end, primerpair.rprimer.start)
        )

    # Plot the uncovered regions
    uncovered_indexes_list = sorted(uncovered_indexes)
    # Generate continous regions
    uncovered_regions = []
    for k, g in groupby(enumerate(uncovered_indexes_list), lambda ix: ix[0] - ix[1]):
        uncovered_regions.append(list(map(itemgetter(1), g)))

    data: dict[int, int] = dict()
    uncovered_regions = [(min(x), max(x)) for x in uncovered_regions]

    for start, end in uncovered_regions:
        data[start] = end
    return data


def generate_genome_gc_data(msa: MSA | PanelMSA, kmersize: int) -> dict[int, float]:
    """Creates a dict of the genome GC% with key as position and value as GC%"""
    gc_data = dict()
    for index, gc in calc_gc(msa.array, kmersize):
        gc_data[index] = gc
    return gc_data


def generate_genome_occupancy_data(msa: MSA | PanelMSA) -> dict[int, float]:
    """Creates a dict of the genome occupancy with key as position and value as occupancy"""
    occupancy_data = dict()
    for x, y in calc_occupancy(msa.array):
        occupancy_data[x] = y
    return occupancy_data


def generate_genome_entropy_data(msa: MSA | PanelMSA) -> dict[int, float]:
    """Creates a dict of the genome entropy with key as position and value as entropy"""
    results = []
    entropy_data = dict()

    try:
        entropy_array = list(msa._entropy_array)  # type: ignore
    except AttributeError:
        entropy_array = entropy_score_array(msa.array)

    # Calculate the entropy score for each position
    for x, y in enumerate(entropy_array):
        results.append((x, y))
    # Reduce the data
    reduced_data = reduce_data(results)
    for x, y in reduced_data:
        entropy_data[x] = y
    return entropy_data


def generate_thermo_pass_primer_data(msa: MSA | PanelMSA) -> dict[int, str]:
    primer_data = dict()

    fprimer_data = dict()
    for fkmer in msa.fkmers:
        fprimer_data[fkmer.end] = len(fkmer.seqs)
    primer_data["F"] = fprimer_data
    rprimer_data = dict()
    for rkmer in msa.rkmers:
        rprimer_data[rkmer.start] = len(rkmer.seqs)
    primer_data["R"] = rprimer_data
    return primer_data


def generate_amplicon_data(
    primerpairs: list[PrimerPair], msa: MSA | PanelMSA
) -> dict[str, dict[str, int | str]]:
    """
    Creates the amplicon plot data
    :param primerpairs: list of PrimerPair objects
    :return: dict of amplicon data
    """
    amplicon_data = dict()

    for primerpair in primerpairs:
        amplicon_data[primerpair.amplicon_number] = {
            "s": min(primerpair.fprimer.starts()),
            "cs": primerpair.fprimer.end,
            "ce": primerpair.rprimer.start,
            "e": max(primerpair.rprimer.ends()),
            "p": primerpair.pool + 1,
            "n": f"{msa._uuid}_{primerpair.amplicon_number}",
        }

    return amplicon_data


def generate_data(msa: MSA | PanelMSA, last_pp_added: list[PrimerPair]) -> dict:
    """
    Generate all the plot data for a single MSA
    :param msa: MSA object
    :param pools: The pools object
    :return: dict of all the plot data
    """
    # Filter the last primerpair added to the multiplex
    msa_pp: list[PrimerPair] = [
        x for x in last_pp_added if x.msa_index == msa.msa_index
    ]

    # Remap the included primers to the MSA if they have been mapped to an genome
    if msa._mapping_array is not None:
        mapping_list = list(msa._mapping_array)
        for fkmer in msa.fkmers:
            fkmer.end = mapping_list.index(fkmer.end)
            fkmer._starts = {fkmer.end - len(x) for x in fkmer.seqs}
        for rkmer in msa.rkmers:
            rkmer.start = mapping_list.index(rkmer.start)
            rkmer._ends = {rkmer.start + len(x) for x in rkmer.seqs}

    # Write all data to a single json file
    data = dict()
    data["gc"] = generate_genome_gc_data(msa, 30)
    data["entropy"] = generate_genome_entropy_data(msa)
    data["occupancy"] = generate_genome_occupancy_data(msa)
    data["thermo_pass"] = generate_thermo_pass_primer_data(msa)
    data["amplicons"] = generate_amplicon_data(msa_pp, msa)
    data["dims"] = [x for x in msa.array.shape]
    data["uncovered"] = generate_uncovered_data(msa.array.shape[1], msa_pp)

    return data


def generate_all_plotdata(
    msas: list[MSA] | list[PanelMSA],
    output_path: pathlib.Path,
    last_pp_added: list[PrimerPair],
) -> dict:
    """
    Generate all the plot data for all MSAs to plotdata.json.gz
    :param msa: list of MSA objects
    :param last_pp_added: list of PrimerPair objects added to the multiplex
    :param output_path: pathlib.Path to write the plotdata.json to
    :return: None
    """
    # Write all data to a single json file
    data = dict()
    for msa in msas:
        data[msa.name] = generate_data(msa, last_pp_added)

    # Write the data to a json file
    json_bytes = json.dumps(data, sort_keys=True).encode("utf-8")
    with gzip.open(output_path / "plotdata.json.gz", "wb") as fout:
        fout.write(json_bytes)

    # Return the data for use in plotting
    return data