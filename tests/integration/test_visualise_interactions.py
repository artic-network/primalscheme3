import pathlib
import unittest

from primalscheme3.interaction.interaction import (
    visualise_interactions,
)


class TestVisualiseInteractions(unittest.TestCase):
    def test_visualise_interactions(self):
        bedfile = pathlib.Path("tests/test_data/test_scheme/primer.bed")
        # Read in a bedfile
        visualise_interactions(bedfile, -26)
