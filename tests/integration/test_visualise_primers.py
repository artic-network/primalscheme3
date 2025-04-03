import pathlib
import tempfile
import unittest

from Bio import SeqIO

from primalscheme3.core.msa import parse_msa
from primalscheme3.core.primer_visual import bedfile_plot_html, primer_mismatch_heatmap


class TestPrimerVisualisation(unittest.TestCase):
    def test_primer_mismatch_heatmap(self):
        # Read in a bedfile
        bedfile = pathlib.Path("tests/test_data/test_scheme/primer.bed")
        msa = pathlib.Path("tests/test_data/test_scheme/reference.fasta")

        with tempfile.TemporaryDirectory(
            dir="tests/integration", suffix="-test-mismatches"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)

            array, seqdict = parse_msa(msa)

            out_path = tempdir_path / "heatmap.html"

            with open(out_path, "w") as outfile:
                outfile.write(
                    primer_mismatch_heatmap(
                        array=array,
                        seqdict=seqdict,
                        bedfile=bedfile,
                        offline_plots=False,
                        include_seqs=True,
                    )
                )

            # Check file is written and not empty
            self.assertTrue(out_path.is_file())
            self.assertTrue(out_path.stat().st_size > 0)

    def test_bedfile_plot_html(self):
        # Read in a bedfile
        bedfile = pathlib.Path("tests/test_data/test_scheme/primer.bed")
        msa = pathlib.Path("tests/test_data/test_scheme/reference.fasta")

        ref = SeqIO.read(msa, "fasta")

        with tempfile.TemporaryDirectory(
            dir="tests/integration", suffix="-test-plot"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)

            out_path = tempdir_path / "heatmap.html"

            with open(out_path, "w") as outfile:
                outfile.write(
                    bedfile_plot_html(
                        bedfile=bedfile, ref_name=ref.id, ref_seq=str(ref.seq)
                    )
                )

            # Check file is written and not empty
            self.assertTrue(out_path.is_file())
            self.assertTrue(out_path.stat().st_size > 0)
