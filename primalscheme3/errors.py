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


class CustomRecursionError(CustomErrors):
    """Walks out of the index of the MSA"""

    pass


ERROR_SET = {
    WalksOut(),
    CustomRecursionError(),
    ContainsInvalidBase(),
    CustomErrors(),
    GapOnSetBase(),
}
