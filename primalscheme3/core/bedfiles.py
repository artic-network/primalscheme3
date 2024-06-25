import pathlib
import re

from primalscheme3.core.classes import FKmer, PrimerPair, RKmer

# Module imports
from primalscheme3.core.seq_functions import expand_ambs

REGEX_PATTERN_PRIMERNAME = re.compile("\\d+(_RIGHT|_LEFT|_R|_L)")


def re_primer_name(string) -> list[str] | None:
    """
    Will return (amplicon_number, R/L) or None
    """
    match = REGEX_PATTERN_PRIMERNAME.search(string)
    if match:
        return match.group().split("_")
    return None


class BedPrimerPair(PrimerPair):
    """Class to contain a single primercloud from a bedfile, which contains the extra info parsed from the bedfile"""

    amplicon_prefix: str
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


class BedLine:
    """
    Contains a single line from a bedfile
    self.pool is stored as a 0 based index
    """

    chrom_name: str
    _start: int
    _end: int
    primername: str
    pool: int
    direction: str
    sequence: str
    # Calc values
    amplicon_number: int

    def __init__(self, bedline: list[str]) -> None:
        self.chrom_name = bedline[0]
        self._start = int(bedline[1])
        self._end = int(bedline[2])
        self.primername = bedline[3]
        self.pool = int(bedline[4]) - 1
        self.direction = bedline[5]
        self.sequence = bedline[6]

        # Check primername is valid
        if len(self.primername.split("_")) != 4:
            raise ValueError(
                f"Invalid primername: {self.primername} in bedline: {bedline}"
            )

        # Calc values
        self.amplicon_number = int(self.primername.split("_")[1])

    def all_seqs(self) -> set[str] | None:
        "Expands ambs bases"
        return expand_ambs([self.sequence])

    @property
    def msa_index(self) -> str:
        return self.chrom_name

    @property
    def start(self) -> int:
        return self._start

    @property
    def end(self) -> int:
        return self._end

    def __str__(self, *kwargs) -> str:
        # I use *kwargs so that it can have the same behavor as PrimerPairs
        return f"{self.chrom_name}\t{self.start}\t{self.end}\t{self.primername}\t{self.pool + 1}\t{self.direction}\t{self.sequence}"


def read_in_bedlines(path: pathlib.Path) -> tuple[list[BedLine], list[str]]:
    """
    Read in bedlines from a file.

    :param path: The path to the bed file.
    :type path: pathlib.Path
    :return: A list of BedLine objects.
    :rtype: tuple(list[BedLine], list[str])
    """
    bed_primers: list[BedLine] = []
    bed_headers: list[str] = []
    with open(path) as bedfile:
        for line in bedfile.readlines():
            line = line.strip()
            if not line:  # Skip empty lines
                continue
            elif line.startswith("#"):  # Store header lines
                bed_headers.append(line.strip())
            else:  # Store primer lines
                line = line.strip().split()
                bed_primers.append(BedLine(line))
    return (bed_primers, bed_headers)


def read_in_bedprimerpairs(path: pathlib.Path) -> tuple[list[BedPrimerPair], list[str]]:
    """
    Read in a bedfile and return a list of BedPrimerPairs, MSA index is set to None as it is not known at this point

    :param path: The path to the bed file.
    :type path: pathlib.Path
    :return: A list of BedPrimerPair objects, and the header lines from the bedfile.
    :rtype: tuple(list[BedPrimerPair], list[str])
    """

    # Read in the bedfile
    primerpairs = []
    primerlines, _headers = read_in_bedlines(path)  # Ignore headers for now

    # Group primers by referance
    ref_to_bedlines: dict[str, list[BedLine]] = dict()
    for ref in {bedline.chrom_name for bedline in primerlines}:
        ref_to_bedlines[ref] = [x for x in primerlines if x.chrom_name == ref]

    for ref, ref_bed_lines in ref_to_bedlines.items():
        # Group the bedlines by amplicon number
        for ampliconnumber in {
            int(bedline.amplicon_number) for bedline in ref_bed_lines
        }:
            amplicon_prefix = ref_bed_lines[0].primername.split("_")[0]
            ampliconlines = [
                x for x in ref_bed_lines if x.amplicon_number == ampliconnumber
            ]
            pool = ampliconlines[0].pool
            # Group the ampliconlines by direction
            fkmer = FKmer(
                [x.end for x in ampliconlines if x.direction == "+"][0],
                [x.sequence for x in ampliconlines if x.direction == "+"],
            )
            rkmer = RKmer(
                [x.start for x in ampliconlines if x.direction == "-"][0],
                [x.sequence for x in ampliconlines if x.direction == "-"],
            )
            primerpairs.append(
                BedPrimerPair(
                    fprimer=fkmer,
                    rprimer=rkmer,
                    msa_index=None,  # This is set later # type: ignore
                    chrom_name=ref,
                    amplicon_number=int(ampliconnumber),
                    amplicon_prefix=amplicon_prefix,
                    pool=pool,
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
