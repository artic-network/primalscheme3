import unittest

from primalschemers import do_pool_interact  # type: ignore


class TestDoPoolsInteract(unittest.TestCase):
    def test_do_pools_interact_py(self):
        """
        Does this version detect known interactions found between 18_LEFT and 76_RIGHT (SARs-CoV-2:v3)?
        18_LEFT: TGGAAATACCCACAAGTTAATGGTTTAAC
        76_RIGHT: ACACCTGTGCCTGTTAAACCAT
        """
        dimerscore = -26

        pool1 = [b"TGGAAATACCCACAAGTTAATGGTTTAAC"]
        pool2 = [b"ACACCTGTGCCTGTTAAACCAT"]
        self.assertTrue(do_pool_interact(pool1, pool2, dimerscore))

    def test_not_do_pools_interact_py(self):
        """
        Does this version not detect known noninteractions
        AGCGTGGTTATTGGATGGGTTTG	AGCAAATCTTTACTAAAAAAAATTTACCTT
        """
        dimerscore = -26

        pool1 = [b"AGCGTGGTTATTGGATGGGTTTG"]
        pool2 = [b"AGCAAATCTTTACTAAAAAAAATTTACCTT"]
        self.assertFalse(do_pool_interact(pool1, pool2, dimerscore))

    def test_pool_do_pools_interact_py(self):
        """
        Can this version detect interactions from a pool also containing non iteractions.
        AGCGTGGTTATTGGATGGGTTTG	AGCAAATCTTTACTAAAAAAAATTTACCTT
        """
        dimerscore = -26

        pool1 = [b"AGCGTGGTTATTGGATGGGTTTG", b"TGGAAATACCCACAAGTTAATGGTTTAAC"]
        pool2 = [b"AGCAAATCTTTACTAAAAAAAATTTACCTT", b"ACACCTGTGCCTGTTAAACCAT"]
        self.assertTrue(do_pool_interact(pool1, pool2, dimerscore))


if __name__ == "__main__":
    unittest.main()
