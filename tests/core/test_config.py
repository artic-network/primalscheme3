import json
import pathlib
import unittest

from primalscheme3.core.config import Config, MappingType


class TestConfig(unittest.TestCase):
    def test_blank_config_init(self):
        config = Config()

        amplicon_size = 400

        self.assertEqual(config.amplicon_size, amplicon_size)
        # Check min and max are set correctly
        self.assertEqual(config.amplicon_size_min, int(amplicon_size * 0.9))
        self.assertEqual(config.amplicon_size_max, int(amplicon_size * 1.1))

        # Check Primer size and GC are set correctly
        self.assertFalse(config.high_gc)
        self.assertEqual(config.primer_size_min, config._primer_size_default_min)
        self.assertEqual(config.primer_size_max, config._primer_size_default_max)
        self.assertEqual(config.primer_gc_min, config._primer_gc_default_min)
        self.assertEqual(config.primer_gc_max, config._primer_gc_default_max)

        # Check mismatch kmer size is set correctly
        self.assertEqual(config.mismatch_kmersize, config.primer_size_min)

        # Check defaults
        self.assertFalse(config.circular)
        self.assertFalse(config.backtrack)
        self.assertEqual(config.min_base_freq, 0.0)
        self.assertFalse(config.ignore_n)

    def test_config_high_gc(self):
        config = Config(high_gc=True)

        self.assertTrue(config.high_gc)
        self.assertEqual(config.primer_size_min, config._primer_size_hgc_min)
        self.assertEqual(config.primer_size_max, config._primer_size_hgc_max)
        self.assertEqual(config.primer_gc_min, config._primer_gc_hgc_min)
        self.assertEqual(config.primer_gc_max, config._primer_gc_hgc_max)

    def test_config_custom(self):
        amplicon_size = 1000
        config = Config(amplicon_size=amplicon_size, backtrack=True, circular=True)

        self.assertEqual(config.amplicon_size, amplicon_size)
        self.assertEqual(config.amplicon_size_min, int(amplicon_size * 0.9))
        self.assertEqual(config.amplicon_size_max, int(amplicon_size * 1.1))

        self.assertTrue(config.backtrack)
        self.assertTrue(config.circular)

    def test_to_json(self):
        """
        Check that the to_json method returns a dictionary
        """
        default_config = Config()

        default_json = default_config.to_dict()
        self.assertIsInstance(default_json, dict)

    def test_assign_kwargs(self):
        config = Config()
        config.assign_kwargs(amplicon_size=1000, circular=True, mapping="first")

        self.assertEqual(config.amplicon_size, 1000)
        self.assertTrue(config.circular)

    def test_assign_kwargs_invalid(self):
        config = Config()

        # Test str
        config.assign_kwargs(amplicon_size="1000")
        self.assertEqual(config.amplicon_size, 1000)

        # Test Enum
        config.assign_kwargs(mapping="first")
        self.assertEqual(config.mapping, MappingType.FIRST)
        with self.assertRaises(ValueError):
            config.assign_kwargs(mapping="invalid")

        # test bool
        config.assign_kwargs(circular="True")
        self.assertTrue(config.circular)

        # Test Path
        config.assign_kwargs(output="test")
        self.assertEqual(config.output, pathlib.Path("test"))

    def test_json_round_trip(self):
        config = Config(
            amplicon_size=1000,
            output="test",
            ncores=1,
            high_gc=True,
            n_pools=3,
            min_overlap=20,
            mapping="first",
            circular=True,
            backtrack=True,
            min_base_freq=0.1,
            ignore_n=True,
        )

        json_str = json.dumps(config.to_dict())
        config_dict = json.loads(json_str)
        new_config = Config()
        new_config.assign_kwargs(**config_dict)

        for new, old in zip(
            new_config.items().items(), config.items().items(), strict=False
        ):
            self.assertEqual(new, old)


if __name__ == "__main__":
    unittest.main()
