import unittest
from pathlib import Path

from Bio.Data import PDBData

from rna3db.ccd.modifications import ModificationHandler


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

    def test_standard_bases(self):
        """Standard RNA bases map to themselves."""
        for code, expected in [("A", "A"), ("C", "C"), ("G", "G"), ("U", "U")]:
            with self.subTest(code=code):
                self.assertEqual(
                    self.modification_handler.rna_letters_3to1(code), expected
                )

    def test_common_modifications(self):
        """Frequently-seen RNA modifications map to the correct parent base.

        These are the modifications most commonly encountered in PDB structures.
        Regressions here would silently corrupt sequences for many structures.
        """
        known = {
            "PSU": "U",  # pseudouridine
            "5MC": "C",  # 5-methylcytidine
            "M2G": "G",  # N2-methylguanosine
            "7MG": "G",  # 7-methylguanosine
            "OMG": "G",  # 2'-O-methylguanosine
            "1MA": "A",  # 1-methyladenosine
            "OMC": "C",  # 2'-O-methylcytidine
            "OMU": "U",  # 2'-O-methyluridine
            "A2M": "A",  # 2-methyladenosine
            "4AC": "C",  # N4-acetylcytidine (found in 7pwo chain 2)
            "4SU": "U",  # 4-thiouridine
            "MIA": "A",  # 2-methylthio-N6-isopentenyladenosine
            "YYG": "G",  # wybutosine
            "GTP": "G",  # guanosine-5'-triphosphate
        }
        for code, expected in known.items():
            with self.subTest(code=code):
                self.assertTrue(
                    self.modification_handler.is_rna(code),
                    msg=f"{code} should be recognised as RNA",
                )
                self.assertEqual(
                    self.modification_handler.rna_letters_3to1(code),
                    expected,
                    msg=f"{code} should map to {expected!r}",
                )

    def test_unknown_maps_to_N(self):
        """Unknown codes return 'N'."""
        self.assertEqual(self.modification_handler.rna_letters_3to1("XYZ"), "N")
        self.assertFalse(self.modification_handler.is_rna("XYZ"))

    def test_is_rna_standard(self):
        for code in ("A", "C", "G", "U"):
            self.assertTrue(self.modification_handler.is_rna(code))


if __name__ == "__main__":
    unittest.main()
