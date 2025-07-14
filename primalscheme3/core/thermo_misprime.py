from math import exp, sqrt

from Bio import SeqIO
from primalbedtools.bedfiles import BedLineParser

from primalscheme3.core.config import Config
from primalscheme3.core.create_report_data import reduce_data
from primalscheme3.core.seq_functions import reverse_complement
from primalscheme3.core.thermo import R, calc_thermo_raw

annealing_temp_k = 65 + 273


def find_mispriming(
    bl,
    ref,
    config: Config,
    threshold: float = -5_000,
    region: None | tuple[int, int] = None,
    kmer_size: int = 120,
):
    matches = []

    step = kmer_size - ((len(bl.sequence) / 2) + 1)
    if region is None:
        start = 0
        end = len(ref) - kmer_size
    else:
        start = region[0] if region[0] >= 0 else 0
        end = region[1] if region[1] <= len(ref) - kmer_size else len(ref) - kmer_size

    f_data = []
    r_data = []

    # Do thermo align
    for i in range(start, end, int(step)):
        # Primer3 rcs the sequence
        kmer = ref[i : i + kmer_size]
        tr = calc_thermo_raw(bl.sequence, kmer, config=config, with_struct=True)
        # add all annealing
        k = exp((tr.ds / R) - (tr.dh / (R * annealing_temp_k)))

        r_data.append(
            (
                i,
                (1 / (1 + sqrt(1 / ((config.dna_conc / 4e9) * k)))) * 100,
            )
        )

        kmer_rc = reverse_complement(kmer)
        tr = calc_thermo_raw(bl.sequence, kmer_rc, config=config, with_struct=True)
        k = exp((tr.ds / R) - (tr.dh / (R * annealing_temp_k)))

        f_data.append(
            (
                i,
                (1 / (1 + sqrt(1 / ((config.dna_conc / 4e9) * k)))) * 100,
            )
        )

    reduce_f_data = reduce_data(f_data, 1)
    reduce_r_data = reduce_data(r_data, 1)

    for fd in reduce_f_data:
        print(bl.primername, fd[0], fd[1], "+")

    for rd in reduce_r_data:
        print(bl.primername, rd[0], rd[1], "-")

    return matches


if __name__ == "__main__":
    ref = SeqIO.read(
        "/Users/kentcg/schemes/tb/v4/v4-400-28-all-2/reference.fasta",
        "fasta",
    )
    _headers, primers = BedLineParser.from_file(
        "/Users/kentcg/schemes/tb/v4/v4-400-28-all-2/primer.bed"
    )
    config = Config()

    for bl in primers[:10]:
        mps = find_mispriming(bl, str(ref.seq), config, kmer_size=50, region=None)
