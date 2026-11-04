"""
Why primer sites fail: every check each rejected site fails, from a
return_all digest, and a table of them for the log.

A site often fails several checks (e.g. HighGC and HighAnnealing on the same
sequence), so each check is counted at every site it fails ("Fails at"), and
separately where it is the only failure ("Only failure": the sites relaxing
that check alone could rescue). Neither depends on the order checks run in.
"""

from collections import Counter

from rich import box
from rich.bar import Bar
from rich.table import Table

# Tags on a return_all kmer's statuses/errors that are not failures
LOW_FREQ = "LowFreq"
IGNORED = "Ignored"

PRIMER_DIMER = "PrimerDimer"
# Every sequence at the site is below min_base_freq
ONLY_LOW_FREQ = "OnlyLowFreq"

# Table sections, in order; any other check goes in the last
THERMO_CHECKS = {
    "HighGC",
    "LowGC",
    "HighTm",
    "LowTm",
    "HighAnnealing",
    "LowAnnealing",
    "Homopolymer",
    "Hairpin",
    "Fail",
    "ToLong",
}
INTERACTION_CHECKS = {PRIMER_DIMER}


def untagged(tags: list[str]) -> bool:
    return LOW_FREQ not in tags and IGNORED not in tags


def site_failures(kmer) -> frozenset[str]:
    """
    Every check a return_all ``kmer`` fails at or above min_base_freq: its
    sequences' failed checks, alignment errors and dimers. Empty if the kmer
    passes.
    """
    if kmer.passed():
        return frozenset()
    checks = set()
    for tags in kmer.statuses() or []:
        if LOW_FREQ not in tags:
            checks.update(tags)
    for tags, _ in kmer.errors():
        if untagged(tags):
            checks.add(tags[0])
    if kmer.dimers():
        checks.add(PRIMER_DIMER)
    return frozenset(checks or {ONLY_LOW_FREQ})


def failure_table(
    title: str,
    rejected: list[frozenset[str]],
    n_sites: int,
    caption: str | None = None,
) -> Table:
    """
    One row per check, grouped (thermodynamic, interaction, alignment):
    the rejected sites it fails at, and those where it is the only failure.
    ``rejected`` holds the checks each rejected site fails.
    """
    n_rejected = len(rejected)
    n_usable = n_sites - n_rejected
    fails_at: Counter = Counter()
    only: Counter = Counter()
    for checks in rejected:
        fails_at.update(checks)
        if len(checks) == 1:
            only.update(checks)

    table = Table(
        title=f"{title}: [green]{n_usable:,}[/green] of {n_sites:,} primer sites "
        f"usable" + (f" ({n_usable / n_sites:.0%})" if n_sites else ""),
        title_justify="left",
        caption=caption,
        caption_justify="left",
        box=box.SIMPLE_HEAD,
        show_footer=True,
    )
    table.add_column("Check", footer="Rejected sites")
    table.add_column("Fails at", justify="right", footer=f"{n_rejected:,}")
    table.add_column("%", justify="right")
    table.add_column("", width=20)
    table.add_column("Only failure", justify="right", footer=f"{only.total():,}")

    sections: list[list[str]] = [[], [], []]
    for check, _ in sorted(fails_at.items(), key=lambda kv: (-kv[1], kv[0])):
        if check in THERMO_CHECKS:
            sections[0].append(check)
        elif check in INTERACTION_CHECKS:
            sections[1].append(check)
        else:
            sections[2].append(check)

    for checks in (s for s in sections if s):
        if table.row_count:
            table.add_section()
        for check in checks:
            n = fails_at[check]
            share = n / n_rejected
            table.add_row(
                check,
                f"{n:,}",
                f"{share:.0%}" if share >= 0.01 else "<1%",
                Bar(n_rejected, 0, n, color="red"),
                f"{only[check]:,}",
            )
    return table
