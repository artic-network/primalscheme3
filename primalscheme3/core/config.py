import pathlib
from enum import Enum
from typing import Any

from primalscheme3 import __version__


# Writen by Andy Smith, modified by: Chris Kent
class MappingType(Enum):
    """
    Enum for the mapping type
    """

    FIRST = "first"
    CONSENSUS = "consensus"


class Config:
    """
    PrimalScheme3 configuration.
    Class properties are defaults, can be overriden
    on instantiation (and will shadow class defaults)
    """

    # Run Settings
    output: pathlib.Path = pathlib.Path("./output")
    force: bool = False
    high_gc: bool = False
    input_bedfile: pathlib.Path | None = None
    version: str = __version__
    # Scheme Settings
    n_pools: int = 2
    min_overlap: int = 10
    mapping: MappingType = MappingType.FIRST
    circular: bool = False
    backtrack: bool = False
    min_base_freq: float = 0.0
    ignore_n: bool = False
    # Amplicon Settings
    amplicon_size: int = 400
    amplicon_size_min: int = 0
    amplicon_size_max: int = 0
    # Primer Settings
    _primer_size_default_min: int = 19
    _primer_size_default_max: int = 34
    _primer_size_hgc_min: int = 17
    _primer_size_hgc_max: int = 30
    _primer_gc_default_min: int = 30
    _primer_gc_default_max: int = 55
    _primer_gc_hgc_min: int = 40
    _primer_gc_hgc_max: int = 65
    primer_tm_min: float = 59.5
    primer_tm_max: float = 62.5
    primer_hairpin_th_max: float = 47.0
    primer_homopolymer_max: int = 5
    primer_max_walk: int = 80
    # MatchDB Settings
    editdist_max: int = 1
    mismatch_fuzzy: bool = True
    mismatch_kmersize: int  # Same as primer_size_min
    mismatch_product_size: int = 0
    # Thermodynamic Parameters
    mv_conc: float = 100.0
    dv_conc: float = 2.0
    dntp_conc: float = 0.8
    dna_conc: float = 15.0
    dimer_score: float = -26.0

    def __init__(self, **kwargs: Any) -> None:
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
        # Set amplicon size
        if self.amplicon_size_min == 0:
            self.amplicon_size_min = int(self.amplicon_size * 0.9)
        if self.amplicon_size_max == 0:
            self.amplicon_size_max = int(self.amplicon_size * 1.1)
        if self.high_gc:
            self.primer_size_min = self._primer_size_hgc_min
            self.primer_size_max = self._primer_size_hgc_max
            self.primer_gc_min = self._primer_gc_hgc_min
            self.primer_gc_max = self._primer_gc_hgc_max
        else:
            self.primer_size_min = self._primer_size_default_min
            self.primer_size_max = self._primer_size_default_max
            self.primer_gc_min = self._primer_gc_default_min
            self.primer_gc_max = self._primer_gc_default_max
        # Set MisMatch Kmer Size
        self.mismatch_kmersize = self.primer_size_min

    def items(self) -> dict[str, Any]:
        """
        Return a dict (key, val) for non-private, non-callable members
        """
        items = {}
        for key, val in self.__dict__.items():
            if not (callable(getattr(self, key)) or key.startswith("_")):
                items[key] = val
        return items

    def to_json(self) -> dict[str, Any]:
        """
        Return a dict (key, val) for non-private, non-callable members
        """
        json = {}
        for key, val in self.items().items():
            if isinstance(val, Enum):
                json[key] = val.value
            elif isinstance(val, pathlib.Path):
                json[key] = str(
                    val.relative_to(pathlib.Path.cwd()) if val.is_absolute() else val
                )
            else:
                json[key] = val
        return json

    def __str__(self) -> str:
        return "\n".join(f"{key}: {val}" for key, val in self.items())


config_dict = {
    "n_cores": 1,
    "amplicon_size_min": 360,
    "amplicon_size_max": 400,
    "min_overlap": 10,
    "mv_conc": 100.0,
    "dv_conc": 2.0,
    "dntp_conc": 0.8,
    "dna_conc": 15.0,
    "primer_gc_min": 30,
    "primer_gc_max": 55,
    "primer_tm_min": 59.5,
    "primer_tm_max": 62.5,
    "primer_homopolymer_max": 5,
    "primer_hairpin_th_max": 47.0,  #
    "editdist_max": 1,
    "editdist_end3p": 6,
    "primer_size_max": 36,
    "primer_size_min": 18,  # primer_size_min
    "primer_max_walk": 80,
    "minbasefreq": 0,
    "reducekmers": False,
    "dimerscore": -27,
    "mismatch_fuzzy": True,
    "mismatch_kmersize": 18,  # primer_size_min
    "ignore_n": False,
}

SIMPLE_BASES = {"A", "C", "G", "T"}

AMBIGUOUS_DNA = {
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
}
ALL_DNA: dict[str, str] = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
}
ALL_BASES: set[str] = {
    "A",
    "C",
    "G",
    "T",
    "M",
    "R",
    "W",
    "S",
    "Y",
    "K",
    "V",
    "H",
    "D",
    "B",
}
AMB_BASES = {"Y", "W", "R", "B", "H", "V", "D", "K", "M", "S"}
AMBIGUOUS_DNA_COMPLEMENT = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "M": "K",
    "R": "Y",
    "W": "W",
    "S": "S",
    "Y": "R",
    "K": "M",
    "V": "B",
    "H": "D",
    "D": "H",
    "B": "V",
    "X": "X",
    "N": "N",
}
