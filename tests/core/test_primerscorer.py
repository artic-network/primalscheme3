"""
Tests for primerscorer.find_mismatch.

`find_mismatch` receives the two sequences ALREADY padded with a leading and
trailing "." (see PrimerScorer._get_features, which calls it as
``find_mismatch("." + ref + ".", "." + query + ".")``). The padding is
deliberate: the trailing "." lets a 3'-terminal substitution be detected, and
the leading "." marks the 5' end.

Primer clouds contain 3'-anchored primers of *different lengths*. Two such
primers are joined by a graph edge when ``primer_edit_dist == 1`` (a
3'-anchored Hamming distance that ignores 5' length differences). For those
edges the model must be able to score the single substitution.

The bug: because the padded strings then have different lengths, the 3'-anchored
comparison lines a padding "." up against a real base at the 5' boundary and
counts it as an extra mismatch, so ``find_mismatch`` reports != 1 mismatch and
returns None (or, for a pure length difference, indexes out of range and
crashes). The result is that every length-variant interaction scores None and is
invisible to the downsample coverage logic.

These tests keep the padding and assert the intended behaviour.
"""

import unittest

from primalscheme3.core.primerscorer import find_mismatch
from primalscheme3.core.seq_functions import primer_edit_dist


def fm(ref: str, query: str):
    """Call find_mismatch exactly like _get_features does: padded inputs."""
    return find_mismatch("." + ref + ".", "." + query + ".")


# A fixed primer (query) and a template variant (ref) differing by ONE
# substitution in the middle of the 3' overlap.
QUERY = "ACGTACGTACGTACGTACGT"
REF_SAME = "ACGTACGTACGTATGTACGT"  # one substitution vs QUERY (interior)
# The substitution sits at the same 3'-anchored position no matter how far REF
# is extended at the 5' end, so the mismatch descriptor must be identical.
EXPECTED_INTERIOR = (7, "AC/TA", "CA/GC", "T>C")

# A 3'-terminal substitution (last base differs).
REF_TERMINAL = "ACGTACGTACGTACGTACGA"
EXPECTED_TERMINAL = (1, "GT/CT", ".T/.T", "A>T")


class TestFindMismatchSameLength(unittest.TestCase):
    """Baseline: equal-length pairs must keep working exactly as before."""

    def test_interior_substitution(self):
        self.assertEqual(fm(REF_SAME, QUERY), EXPECTED_INTERIOR)

    def test_terminal_substitution(self):
        self.assertEqual(fm(REF_TERMINAL, QUERY), EXPECTED_TERMINAL)

    def test_more_than_one_substitution_is_none(self):
        two_subs = "ACGTACGTACGTATGTATGT"  # two diffs vs QUERY
        self.assertEqual(primer_edit_dist(two_subs, QUERY), 2)
        self.assertIsNone(fm(two_subs, QUERY))


class TestFindMismatchLengthVariant(unittest.TestCase):
    """
    The failure case: a length-variant neighbour (edge by primer_edit_dist)
    must still resolve to its single 3' substitution.
    """

    def test_ref_extended_by_two_is_an_edge(self):
        ref_long = "GG" + REF_SAME
        self.assertEqual(primer_edit_dist(ref_long, QUERY), 1)  # it IS a graph edge

    def test_ref_extended_by_two_is_scored(self):
        ref_long = "GG" + REF_SAME
        # Currently returns None -> the interaction is dropped. It must be scored.
        self.assertIsNotNone(fm(ref_long, QUERY))

    def test_ref_extended_matches_equivalent_same_length(self):
        # The 5' extension does not touch the 3' context, so the mismatch
        # descriptor must be identical to the equal-length pair.
        self.assertEqual(fm("GG" + REF_SAME, QUERY), EXPECTED_INTERIOR)
        self.assertEqual(fm("G" + REF_SAME, QUERY), EXPECTED_INTERIOR)

    def test_query_extended_is_scored(self):
        # Symmetric: the primer (query) being the longer one must also work.
        query_long = "GG" + QUERY
        self.assertEqual(primer_edit_dist(REF_SAME, query_long), 1)
        self.assertIsNotNone(fm(REF_SAME, query_long))

    def test_terminal_substitution_with_length_variant(self):
        # Terminal substitution should still be detected as terminal (idx == 1)
        # even when ref is 5'-extended.
        result = fm("GG" + REF_TERMINAL, QUERY)
        self.assertIsNotNone(result)
        self.assertEqual(result[0], 1)  # still terminal


class TestFindMismatchPureLengthDifference(unittest.TestCase):
    """
    A pure 5' length difference with no substitution is edit_dist 0 (never a
    graph edge). It must return None cleanly rather than crash.
    """

    def test_pure_extension_is_not_an_edge(self):
        self.assertEqual(primer_edit_dist("GG" + QUERY, QUERY), 0)

    def test_pure_extension_returns_none_without_crashing(self):
        # Currently raises IndexError.
        self.assertIsNone(fm("GG" + QUERY, QUERY))


class TestEdgeInvariant(unittest.TestCase):
    """
    Core invariant tying the graph to the scorer: any pair that forms a graph
    edge (primer_edit_dist == 1) must be scorable (find_mismatch is not None),
    regardless of length difference.
    """

    def test_every_edit_dist_one_pair_is_scorable(self):
        pairs = [
            (REF_SAME, QUERY),  # same length interior
            (REF_TERMINAL, QUERY),  # same length terminal
            ("G" + REF_SAME, QUERY),  # ref longer by 1
            ("GG" + REF_SAME, QUERY),  # ref longer by 2
            (REF_SAME, "GG" + QUERY),  # query longer by 2
        ]
        for ref, query in pairs:
            with self.subTest(ref=ref, query=query):
                self.assertEqual(primer_edit_dist(ref, query), 1)
                self.assertIsNotNone(
                    fm(ref, query),
                    msg=f"edge {ref!r}->{query!r} is not scorable",
                )


if __name__ == "__main__":
    unittest.main()
