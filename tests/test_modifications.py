import unittest
import pickle

from rna3db.parser import ModificationHandler
from Bio.Data import SCOPData

# import PDBData


class TestModifications(unittest.TestCase):
    modification_handler = ModificationHandler(
        "tests/test_data/modifications_cache.json"
    )

    # accept both A and U for selenocysteines
    # note that selenocysteines that are incorrectly encoded in SCOPData
    # shouldn't be included here
    selenocysteines = ["SEC"]

    def test_scop_equivalence(self):
        for k in SCOPData.protein_letters_3to1.keys():
            k = k.rstrip()  # to address biopython's weird whitespace

            scop_code = SCOPData.protein_letters_3to1.get(k, "X")
            scop_code = scop_code if len(scop_code) == 1 else "X"
            protein_code = self.modification_handler.protein_letters_3to1(k)

            # if it's not a protein we don't care about equivalence
            if not self.modification_handler.is_protein(k):
                continue

            # selenocysteines (U) are treated as A by AlphaFold
            # so we don't care if they are encoded as A in SCOPData
            if k in self.selenocysteines:
                self.assertTrue(protein_code == "A" or protein_code == "U")
                continue

            # these two modifications were fixed in a later version of
            # biopython so we check them manually
            if k == "4BF":
                self.assertEqual(protein_code, "F")
                continue
            if k == "PCA":
                self.assertEqual(protein_code, "Q")
                continue

            # we are happy to recover more unknowns than SCOPData
            self.assertTrue(
                protein_code == scop_code or "X" == scop_code,
            )

    def test_biopython_coverage(self):
        for k, v in SCOPData.protein_letters_3to1.items():
            k = k.rstrip()  # to address biopython's weird whitespace

            if self.modification_handler.is_protein(k):
                self.assertTrue(
                    v == "X" or self.modification_handler.protein_letters_3to1(k) != "X"
                )

            if self.modification_handler.is_rna(k):
                self.assertTrue(
                    v == "N" or self.modification_handler.rna_letters_3to1(k) != "N"
                )


if __name__ == "__main__":
    unittest.main()
