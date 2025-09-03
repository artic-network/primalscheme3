import pathlib
from enum import Enum
from importlib.metadata import version
from typing import Any


# Written by Andy Smith, modified by: Chris Kent
class MappingType(Enum):
    """
    Enum for the mapping type
    """

    FIRST = "first"
    CONSENSUS = "consensus"


class Config:
    """
    PrimalScheme3 configuration.
    Class properties are defaults, can be overridden
    on instantiation (and will shadow class defaults)
    """

    # Run Settings
    output: pathlib.Path = pathlib.Path("./output")
    force: bool = False
    high_gc: bool = False
    input_bedfile: pathlib.Path | None = None
    version: str = version("primalscheme3")
    ncores = 1
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
    # Thermo
    primer_tm_min: float = 59.5
    primer_tm_max: float = 62.5
    primer_annealing_tempc = 65
    primer_annealing_prop: float | None = None
    _primer_annealing_prop_default = 10
    use_annealing: bool = False

    primer_hairpin_th_max: float = 47.0
    primer_homopolymer_max: int = 5
    primer_max_walk: int = 80
    # MatchDB Settings
    use_matchdb: bool = True
    editdist_max: int = 1
    mismatch_fuzzy: bool = True
    mismatch_kmersize: int  # Same as primer_size_min
    mismatch_product_size: int = 0
    mismatch_in_memory: bool = False  # Use an in memory dict rather than a dbm file
    # Thermodynamic Parameters
    mv_conc: float = 100.0
    dv_conc: float = 2.0
    dntp_conc: float = 0.8
    dna_conc: float = 15.0
    dimer_score: float = -26.0

    def __init__(self, **kwargs: Any) -> None:
        self.assign_kwargs(**kwargs)
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
        # Set annealing
        if self.use_annealing:
            self.primer_annealing_prop = self._primer_annealing_prop_default

    def items(self) -> dict[str, Any]:
        """
        Return a dict (key, val) for non-private, non-callable members
        """
        items = {}
        for key in [x for x in dir(self) if not x.startswith("_")]:
            if not callable(getattr(self, key)):  # prevent functions
                value = getattr(self, key)
                # Convert Enum and Path objects to their values
                if isinstance(value, Enum):
                    items[key] = value.value
                if isinstance(value, pathlib.Path):
                    items[key] = str(value)
                else:
                    items[key] = value

        return items

    def to_dict(self) -> dict[str, Any]:
        """
        Return a dict (key, val) for non-private, non-callable members
        """
        dict = {}
        for key, val in self.items().items():
            if isinstance(val, Enum):
                dict[key] = val.value
            elif isinstance(val, pathlib.Path):
                dict[key] = str(val)
            else:
                dict[key] = val
        return dict

    def __str__(self) -> str:
        return "\n".join(f"{key}: {val}" for key, val in self.items())

    def assign_kwargs(self, **kwargs: Any) -> None:
        for key, value in kwargs.items():
            # Check if key is valid
            if value is None:
                continue

            if key == "input_bedfile":
                setattr(self, key, pathlib.Path(value))
                continue

            if hasattr(self, key):
                # Convert to expected type
                if isinstance(getattr(self, key), MappingType):
                    setattr(self, key, MappingType(value))
                elif isinstance(getattr(self, key), pathlib.Path):
                    setattr(self, key, pathlib.Path(value))
                elif isinstance(
                    getattr(self, key), bool
                ):  # Need to check bool before int
                    parsed_bool = str(value).lower()
                    setattr(self, key, parsed_bool == "true")
                elif isinstance(getattr(self, key), int):
                    setattr(self, key, int(value))
                elif isinstance(getattr(self, key), float):
                    setattr(self, key, float(value))
                elif isinstance(getattr(self, key), str):
                    setattr(self, key, str(value))
                else:
                    print(f"Could not parse {key} with value {value} ({type(value)}")


# All bases allowed in the input MSA
IUPAC_ALL_ALLOWED_DNA = {
    "A",
    "G",
    "K",
    "Y",
    "B",
    "S",
    "N",
    "H",
    "C",
    "W",
    "D",
    "R",
    "M",
    "T",
    "V",
    "-",
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
ALL_DNA_WITH_N: dict[str, str] = {
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
    "N": "ACGT",
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
ALL_BASES_WITH_N: set[str] = {
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
    "N",
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
    "-": "-",
}
