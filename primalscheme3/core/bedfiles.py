import logging
import pathlib
import sys

import primalbedtools.bedfiles as bf
from primalschemers._core import FKmer, RKmer  # type: ignore

from primalscheme3.core.classes import PrimerPair
from primalscheme3.core.config import Config


class BedPrimerPair(PrimerPair):
    """Class to contain a single primercloud from a bedfile, which contains the extra info parsed from the bedfile"""

    amplicon_prefix: str | None
    # Calc values
    _primername: str

    def __init__(
        self,
        fprimer: FKmer,
        rprimer: RKmer,
        msa_index: int,
        chrom_name: str,
        amplicon_prefix: str,
        amplicon_number: int,
        pool: int,
    ) -> None:
        self.fprimer = fprimer
        self.rprimer = rprimer
        self.chrom_name = chrom_name
        self.amplicon_prefix = amplicon_prefix
        self.msa_index = msa_index
        self.amplicon_number = amplicon_number
        self.pool = pool

        #
        self._primername = f"{self.amplicon_number}_{self.amplicon_prefix}"

    def match_primer_stem(self, primernamestem: str) -> bool:
        return self._primername == primernamestem


def read_bedlines_to_bedprimerpairs(
    path: pathlib.Path,
) -> tuple[list[BedPrimerPair], list[str]]:
    """
    uses primalbedtools to read in the bedlines. Parses all bedlines into primerclouds.
    """
    _headers, bedlines = bf.BedLineParser.from_file(path)
    grouped_bedlines = bf.group_primer_pairs(bedlines)

    primerpairs = []

    for fks, rks in grouped_bedlines:
        fk_end = {fk.end for fk in fks}
        if len(fk_end) != 1:
            raise ValueError("Cannot combine into Single FP")

        rk_start = {rk.start for rk in rks}
        if len(rk_start) != 1:
            raise ValueError("Cannot combine into Single RP")

        primerpairs.append(
            BedPrimerPair(
                fprimer=FKmer([fk.sequence.encode() for fk in fks], fk_end.pop()),
                rprimer=RKmer([rk.sequence.encode() for rk in rks], rk_start.pop()),
                msa_index=None,  # This is set later # type: ignore
                chrom_name=fks[0].chrom,
                amplicon_number=fks[0].amplicon_number,
                amplicon_prefix=fks[0].amplicon_prefix,
                pool=fks[0].ipool,
            )
        )

    primerpairs.sort(key=lambda x: (x.chrom_name, x.amplicon_number))

    return (primerpairs, _headers)


def create_bedfile_str(
    headers: list[str] | None, primerpairs: list[PrimerPair | BedPrimerPair]
) -> str:
    """
    Returns the multiplex as a bed file
    :return: str
    """
    primer_bed_str: list[str] = []

    # Ensure headers are commented and valid
    if headers is not None:
        for headerline in headers:
            if not headerline.startswith("#"):
                headerline = "# " + headerline
            primer_bed_str.append(headerline.strip())

    # Add the primerpairs to the bed file
    for pp in primerpairs:
        primer_bed_str.append(pp.to_bed().strip())

    return "\n".join(primer_bed_str) + "\n"


def create_amplicon_str(
    primerpairs: list[PrimerPair | BedPrimerPair], trim_primers: bool = False
) -> str:
    amplicon_str: list[str] = []
    # Add the amplicons to the string
    for pp in primerpairs:
        if trim_primers:
            amplicon_str.append(
                f"{pp.chrom_name}\t{pp.fprimer.region()[1]}\t{pp.rprimer.region()[0] - 1}\t{pp.amplicon_prefix}_{pp.amplicon_number}\t{pp.pool + 1}"
            )
        else:
            amplicon_str.append(
                f"{pp.chrom_name}\t{pp.fprimer.region()[0]}\t{pp.rprimer.region()[1]}\t{pp.amplicon_prefix}_{pp.amplicon_number}\t{pp.pool + 1}"
            )
    return "\n".join(amplicon_str) + "\n"


def read_in_extra_primers(
    input_bedfile: pathlib.Path, config: Config, logger: logging.Logger
) -> list[BedPrimerPair]:
    """
    Reads in Primers from a bedfile, and QC checks them for Tm and Pools
    """
    bedprimerpairs, _headers = read_bedlines_to_bedprimerpairs(input_bedfile)

    logger.info(
        f"Read in bedfile: [blue]{input_bedfile.name}[/blue]: "
        f"[green]{len(bedprimerpairs)}[/green] PrimersPairs containing "
        f"{len([primer for primers in (bedprimerpair.all_seqs() for bedprimerpair in bedprimerpairs) for primer in primers])} primers",
    )

    # Check the primers for Tm
    primer_tms = [
        tm for tm in (pp.calc_tm(config) for pp in bedprimerpairs) for tm in tm
    ]
    if min(primer_tms) < config.primer_tm_min or max(primer_tms) > config.primer_tm_max:
        logger.warning(
            f"Primer Tm outside range: {round(min(primer_tms), 2)} : {round(max(primer_tms), 2)} (range: {config.primer_tm_min} : {config.primer_tm_max})"
        )

    else:
        logger.info(
            f"Primer Tm range: [green]{min(primer_tms)}[/green] : [green]{max(primer_tms)}[/green]"
        )

    # Check pools are within the range
    pools_in_bed = {primer.pool for primer in bedprimerpairs}
    if max(pools_in_bed) > config.n_pools:
        logger.critical(
            f"The number of pools in the bedfile is greater than --npools: "
            f"{max(pools_in_bed)} > {config.n_pools}"
        )
        sys.exit(1)

    return bedprimerpairs
