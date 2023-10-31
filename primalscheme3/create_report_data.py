import json
import pathlib

# Module imports
from primal_digest.msa import MSA
from primal_digest.seq_functions import entropy_score_array
from primal_digest.create_reports import calc_gc, calc_occupancy

# Plot format
# plot1: scheme coverage plot. Can be parsed from the bedfile
# plot2 Base occupancy + genome gc
# plot3: Entropy plot
# plot4: Thermo pass Fkmer and Rkmer plot


def generate_genome_gc_data(msa: MSA, kmersize) -> dict[int, float]:
    """Creates a dict of the genome GC% with key as position and value as GC%"""
    gc_data = dict()
    for index, gc in calc_gc(msa.array, kmersize):
        gc_data[index] = gc
    return gc_data


def generate_genome_occupancy_data(msa: MSA) -> dict[int, float]:
    """Creates a dict of the genome occupancy with key as position and value as occupancy"""
    occupancy_data = dict()
    for x, y in enumerate(calc_occupancy(msa.array)):
        occupancy_data[x] = y
    return occupancy_data


def generate_genome_entropy_data(msa: MSA) -> dict[int, float]:
    """Creates a dict of the genome entropy with key as position and value as entropy"""
    entropy_data = dict()
    for x, y in enumerate(entropy_score_array(msa.array)):
        entropy_data[x] = y
    return entropy_data


def generate_thermo_pass_primer_data(msa: MSA) -> dict[int, str]:
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


def generate_data(msa: MSA, outdir: pathlib.Path):
    """Generate a plot for a single MSA"""

    # Write all data to a single json file
    data = dict()
    data["genome_gc"] = generate_genome_gc_data(msa, 30)
    data["genome_entropy"] = generate_genome_entropy_data(msa)
    data["genome_occupancy"] = generate_genome_occupancy_data(msa)
    data["thermo_pass"] = generate_thermo_pass_primer_data(msa)

    # Write all data to individual json files
    with open(outdir / f"{msa.name}_gc.json", "w") as f:
        json.dump(data["genome_gc"], f, sort_keys=True, indent=4)
    with open(outdir / f"{msa.name}_entropy.json", "w") as f:
        json.dump(data["genome_entropy"], f, sort_keys=True, indent=4)
    with open(outdir / f"{msa.name}_occupancy.json", "w") as f:
        json.dump(data["genome_occupancy"], f, sort_keys=True, indent=4)
    with open(outdir / f"{msa.name}_thermo_pass.json", "w") as f:
        json.dump(data["thermo_pass"], f, sort_keys=True, indent=4)
