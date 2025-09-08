import pathlib
import unittest

from primalscheme3.core.config import Config
from primalscheme3.core.msa import MSA
from primalscheme3.core.progress_tracker import ProgressManager as PM


class Test_MappingIntergration(unittest.TestCase):
    input_path = pathlib.Path("tests/test_data/test_mapping.fasta")
    config = Config()

    def test_create_msa(self):
        """
        Test that a msa can be created
        """
        msa = MSA(
            "test_mapping",
            self.input_path,
            0,
            progress_manager=PM(),
            config=self.config,
        )
        msa.digest_rs(config=self.config)

        self.assertIsNotNone(msa._mapping_array)
        mapping_list = list(msa._mapping_array)  # type: ignore

        # Assert the msa._maping_array is as expected
        self.assertEqual(len(msa._mapping_array), msa.array.shape[1])  # type: ignore

        for f in msa.fkmers:
            # Ensure all Fkmers are in the range of the mapping array
            self.assertTrue(f.end in msa._mapping_array)  # type: ignore
            # Ensure the Fkmers can be unmapped
            orginal_end = mapping_list.index(f.end)
            print(orginal_end, f.end)

        # Ensure all Rkmers are in the range of the mapping array
        for r in msa.rkmers:
            self.assertTrue(r.start in msa._mapping_array)  # type: ignore

        #


if __name__ == "__main__":
    unittest.main()
