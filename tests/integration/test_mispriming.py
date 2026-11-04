import json
import pathlib
import random
import tempfile
import unittest
from collections import defaultdict

from primalschemers import PoolHits  # type: ignore

from primalscheme3.core.bedfiles import read_bedlines_to_bedprimerpairs
from primalscheme3.core.config import Config, MappingType
from primalscheme3.core.mismatches import MatchDB
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.scheme.scheme_main import schemecreate

MSA_PATH = pathlib.Path("./tests/core/test_mismatch.fasta").absolute()


def first_genome(path: pathlib.Path) -> str:
    seqs = path.read_text().split(">")[1:]
    return "".join(seqs[0].splitlines()[1:]).replace("-", "")


def pools_forming_products(
    primer_bed: pathlib.Path, background: pathlib.Path, config: Config
) -> list[int]:
    """Pools whose primers would amplify a product from the background"""
    matchdb = MatchDB([], [str(background)], config)
    pools = defaultdict(list)
    for pp in read_bedlines_to_bedprimerpairs(primer_bed)[0]:
        pools[pp.pool].extend(pp.find_matches(matchdb))
    return [
        pool
        for pool, hits in pools.items()
        if PoolHits().forms_product(hits, config.mismatch_product_size)
    ]


class TestBackgroundScreening(unittest.TestCase):
    def config(self, **kwargs) -> Config:
        return Config(mapping=MappingType.FIRST, use_matchdb=True, **kwargs)

    def run_scheme(self, out: pathlib.Path, config: Config) -> pathlib.Path:
        schemecreate(
            msa=[MSA_PATH],
            output_dir=out,
            pm=ProgressManager(),
            config=config,
            force=True,
            offline_plots=False,
        )
        return out / "primer.bed"

    def test_background_product_avoided(self):
        with tempfile.TemporaryDirectory() as d:
            d = pathlib.Path(d)
            plain_bed = self.run_scheme(d / "plain", self.config())

            # A background carrying the first amplicon of the plain scheme
            pp = read_bedlines_to_bedprimerpairs(plain_bed)[0][0]
            start, end = pp.fprimer.region()[0], pp.rprimer.region()[1]
            rng = random.Random(1)
            flank = "".join(rng.choice("ACGT") for _ in range(300))
            background = d / "background.fasta"
            background.write_text(
                f">decoy\n{flank}{first_genome(MSA_PATH)[start:end]}{flank}\n"
            )
            self.assertTrue(
                pools_forming_products(plain_bed, background, self.config())
            )

            config = self.config(background_genomes=[background])
            screened_bed = self.run_scheme(d / "screened", config)
            self.assertEqual(
                pools_forming_products(screened_bed, background, config), []
            )
            written = json.loads((d / "screened" / "config.json").read_text())
            self.assertEqual(written["background_genomes"], [str(background)])
            self.assertIn(str(background), written["background_data"])
            checks = written["pair_checks"]
            self.assertGreater(checks["OK"], 0)
            self.assertGreater(checks["MISPRIMING"], 0)
            self.assertIn("background.fasta", checks["mispriming_references"])


if __name__ == "__main__":
    unittest.main()
