import unittest

from primalscheme3.repair.repair import required_n_pools


class _PP:
    """Minimal stand-in for a PrimerPair: only the 0-indexed ``pool`` matters."""

    def __init__(self, pool: int):
        self.pool = pool


class Test_RequiredNPools(unittest.TestCase):
    def test_contiguous_pools(self):
        # pools {0, 1} -> need 2 buckets
        pps = [_PP(0), _PP(1), _PP(0)]
        self.assertEqual(required_n_pools(pps, min_pools=2), 2)

    def test_sparse_pools_do_not_undersize(self):
        # Regression: pools {0, 2} (pool 1 empty) must yield 3 buckets so that
        # indexing pools[pp.pool] for pp.pool == 2 does not raise IndexError.
        pps = [_PP(0), _PP(2)]
        self.assertEqual(required_n_pools(pps, min_pools=2), 3)
        # The highest index is addressable in a list of this size.
        buckets: list[list] = [[] for _ in range(required_n_pools(pps, 2))]
        for pp in pps:
            buckets[pp.pool].append(pp)  # must not raise

    def test_min_pools_floor(self):
        # More configured pools than are used -> keep the configured floor.
        pps = [_PP(0)]
        self.assertEqual(required_n_pools(pps, min_pools=4), 4)

    def test_empty_primerpairs(self):
        self.assertEqual(required_n_pools([], min_pools=2), 2)


if __name__ == "__main__":
    unittest.main()
