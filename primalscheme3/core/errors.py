from click import UsageError


class CustomErrors(Exception):
    """Base class for custom errors"""

    def __eq__(self, other: object) -> bool:
        return type(self) == type(other)

    def __hash__(self) -> int:
        return hash(type(self))


class GapOnSetBase(CustomErrors):
    """
    Defining position contains a gap
    """

    pass


class ContainsInvalidBase(CustomErrors):
    """
    Contains an invalid base
    """

    pass


class WalksOut(CustomErrors):
    """Walks out of the index of the MSA"""

    pass


class WalksTooFar(CustomErrors):
    """Primer walks more than primer_max_walk"""

    pass


class CustomRecursionError(CustomErrors):
    """Walks out of the index of the MSA"""

    pass


ERROR_SET = {
    WalksOut(),
    CustomRecursionError(),
    ContainsInvalidBase(),
    CustomErrors(),
    GapOnSetBase(),
    WalksTooFar(),
}


# MSA file errors
class MSAFileInvalid(UsageError):
    """Error raised when the MSA file is invalid"""

    pass


class MSAFileInvalidBase(MSAFileInvalid):
    """Error raised when the MSA file contains invalid bases"""

    pass


class MSAFileInvalidLength(MSAFileInvalid):
    """Error raised when the MSA file contains sequences of different lengths"""

    pass


class MSAFileDuplicateID(MSAFileInvalid):
    """Error raised when the MSA file contains sequences of different lengths"""

    pass


# Bed file errors
class BEDFileInvalid(UsageError):
    """Error raised when the BED file is invalid"""

    pass


# Digestion Errors
class DigestionFail(UsageError):
    """Base class for digestion failures"""

    pass


class DigestionFailNoPrimerPairs(DigestionFail):
    """Error raised when no primer pairs are found"""

    pass
