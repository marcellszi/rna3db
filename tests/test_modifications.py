from Bio.Data import PDBData
from pathlib import Path
import unittest

from rna3db.parser import ModificationHandler


class TestModifications(unittest.TestCase):
    modification_handler = ModificationHandler(
        Path(__file__).parent / "test_data" / "modifications_cache.json"
    )

    def test_biopython_coverage(self):
        for k, v in PDBData.nucleic_letters_3to1_extended.items():
            k = k.rstrip()  # to address biopython's weird whitespace

            # we just check that we have all PDBData modifications at least
            if self.modification_handler.is_rna(k):
                self.assertTrue(
                    v == "N" or self.modification_handler.rna_letters_3to1(k) != "N"
                )


if __name__ == "__main__":
    unittest.main()
