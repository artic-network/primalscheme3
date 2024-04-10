import unittest

from primaldimer_py import do_pools_interact_py  # type: ignore


class TestDoPoolsInteract(unittest.TestCase):
    def test_do_pools_interact_py(self):
        """
        Does this version detect known interactions found between 18_LEFT and 76_RIGHT (SARs-CoV-2:v3)?
        18_LEFT: TGGAAATACCCACAAGTTAATGGTTTAAC
        76_RIGHT: ACACCTGTGCCTGTTAAACCAT
        """
        dimerscore = -26

        pool1 = ["TGGAAATACCCACAAGTTAATGGTTTAAC"]
        pool2 = ["ACACCTGTGCCTGTTAAACCAT"]
        self.assertTrue(do_pools_interact_py(pool1, pool2, dimerscore))

    def test_not_do_pools_interact_py(self):
        """
        Does this version not detect known noninteractions
        AGCGTGGTTATTGGATGGGTTTG	AGCAAATCTTTACTAAAAAAAATTTACCTT
        """
        dimerscore = -26

        pool1 = ["AGCGTGGTTATTGGATGGGTTTG"]
        pool2 = ["AGCAAATCTTTACTAAAAAAAATTTACCTT"]
        self.assertFalse(do_pools_interact_py(pool1, pool2, dimerscore))

    def test_pool_do_pools_interact_py(self):
        """
        Can this version detect interactions from a pool also containing non iteractions.
        AGCGTGGTTATTGGATGGGTTTG	AGCAAATCTTTACTAAAAAAAATTTACCTT
        """
        dimerscore = -26

        pool1 = ["AGCGTGGTTATTGGATGGGTTTG", "TGGAAATACCCACAAGTTAATGGTTTAAC"]
        pool2 = ["AGCAAATCTTTACTAAAAAAAATTTACCTT", "ACACCTGTGCCTGTTAAACCAT"]
        self.assertTrue(do_pools_interact_py(pool1, pool2, dimerscore))


if __name__ == "__main__":
    unittest.main()
