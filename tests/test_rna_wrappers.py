from rna3db.parser import Chain, Residue
import unittest


class TestRNAWrappers(unittest.TestCase):
    def test_basic_residue(self):
        res = Residue("G", "G", 1)
        self.assertEqual(res.index, 1)
        self.assertEqual(res.three_letter_code, "G")

    def test_residue_atoms(self):
        res = Residue("G", "G", 1)
        res.atoms["P"] = (0, 0, 0)

    def test_chain_multi_index(self):
        chain = Chain("A")
        chain.add_residue(Residue("A", "A", 0))
        chain.add_residue(Residue("C", "C", 1))
        chain.add_residue(Residue("G", "G", 1))
        chain.add_residue(Residue("U", "U", 2))
        self.assertEqual(len(chain), 3)
        self.assertEqual(chain.sequence, "ACU")

    def test_chain_missing_index(self):
        chain = Chain("A")
        chain.add_residue(Residue("A", "A", 0))
        chain.add_residue(Residue("C", "C", 4))
        self.assertEqual(len(chain), 5)
        self.assertEqual(chain.sequence, "ANNNC")


if __name__ == "__main__":
    unittest.main()
