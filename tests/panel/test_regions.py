import unittest

from primalscheme3.core.errors import CustomErrors
from primalscheme3.panel.panel_classes import Region, RegionParser


class TestRegion(unittest.TestCase):
    def test_region_initialization(self):
        region = Region("chr1", 100, 200, "region1", 0)
        self.assertEqual(region.chromname, "chr1")
        self.assertEqual(region.start, 100)
        self.assertEqual(region.stop, 200)
        self.assertEqual(region.name, "region1")
        self.assertEqual(region.score, 0)
        self.assertIsNone(region.group)

    def test_region_initialization_with_group(self):
        region = Region("chr1", 100, 200, "region1", 0, "group1")
        self.assertEqual(region.group, "group1")

    def test_region_no_group(self):
        # group not provided
        region = Region("chr1", 100, 200, "region1", 0)
        self.assertIsNone(region.group)
        # Group explicitly set to None
        region = Region("chr1", 100, 200, "region1", 0, None)
        self.assertIsNone(region.group)
        # Group empty str
        region = Region("chr1", 100, 200, "region1", 0, "")
        self.assertIsNone(region.group)

    def test_region_invalid_initialization(self):
        with self.assertRaises(ValueError):
            Region("chr1", 200, 100, "region1", 0)

    def test_region_positions(self):
        region = Region("chr1", 100, 200, "region1", 0)
        self.assertEqual(list(region.positions()), list(range(100, 200)))

    def test_region_hash(self):
        region1 = Region("chr1", 100, 200, "region1", 0)
        region2 = Region("chr1", 100, 200, "region1", 0)
        self.assertEqual(hash(region1), hash(region2))

    def test_region_equality(self):
        region1 = Region("chr1", 100, 200, "region1", 0)
        region2 = Region("chr1", 100, 200, "region1", 0)
        region3 = Region("chr1", 100, 201, "region1", 0)
        self.assertEqual(region1, region2)
        self.assertNotEqual(region1, region3)

    def test_region_to_bed(self):
        # With group
        region = Region("chr1", 100, 200, "region1", 0, "group1")
        self.assertEqual(region.to_bed(), "chr1\t100\t200\tregion1\t0\tgroup1\n")

        # Without group
        region = Region("chr1", 100, 200, "region1", 0)
        self.assertEqual(region.to_bed(), "chr1\t100\t200\tregion1\t0\t\n")


class TestRegionParser(unittest.TestCase):
    def test_from_list(self):
        bed_list = ["chr1", "100", "200", "region1", "0", "group1"]
        region = RegionParser.from_list(bed_list)
        self.assertEqual(region.chromname, "chr1")
        self.assertEqual(region.start, 100)
        self.assertEqual(region.stop, 200)
        self.assertEqual(region.name, "region1")
        self.assertEqual(region.score, 0)
        self.assertEqual(region.group, "group1")

    def test_from_list_invalid(self):
        bed_list = ["chr1", "100", "region1", "0"]
        with self.assertRaises(CustomErrors):
            RegionParser.from_list(bed_list)

    def test_from_list_invalid_integers(self):
        bed_list = ["chr1", "start", "200", "region1", "0"]
        with self.assertRaises(CustomErrors):
            RegionParser.from_list(bed_list)

    def test_from_str(self):
        # With group
        bed_str = "chr1\t100\t200\tregion1\t0\tgroup1"
        region = RegionParser.from_str(bed_str)
        self.assertEqual(region.chromname, "chr1")
        self.assertEqual(region.start, 100)
        self.assertEqual(region.stop, 200)
        self.assertEqual(region.name, "region1")
        self.assertEqual(region.score, 0)
        self.assertEqual(region.group, "group1")

        # With no group col
        bed_str = "chr1\t100\t200\tregion1\t0"
        region = RegionParser.from_str(bed_str)
        self.assertEqual(region.chromname, "chr1")
        self.assertEqual(region.start, 100)
        self.assertEqual(region.stop, 200)
        self.assertEqual(region.name, "region1")
        self.assertEqual(region.score, 0)
        self.assertIsNone(region.group)

        # With empty group col
        bed_str = "chr1\t100\t200\tregion1\t0\t"
        region = RegionParser.from_str(bed_str)
        self.assertEqual(region.chromname, "chr1")
        self.assertEqual(region.start, 100)
        self.assertEqual(region.stop, 200)
        self.assertEqual(region.name, "region1")
        self.assertEqual(region.score, 0)
        self.assertIsNone(region.group)

    def test_from_str_invalid(self):
        bed_str = "chr1\t100\tregion1\t0"
        with self.assertRaises(CustomErrors):
            RegionParser.from_str(bed_str)


if __name__ == "__main__":
    unittest.main()
