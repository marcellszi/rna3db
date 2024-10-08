import unittest

from tempfile import NamedTemporaryFile
from pathlib import Path

from rna3db.parser import parse_file, Residue, Chain


class TestmmCIFWriter(unittest.TestCase):
    test_trna_path = Path(__file__).parent / "test_data" / "1ehz.cif"
    test_old_path = Path(__file__).parent / "test_data" / "old_1ehz_A.cif"
    test_two_chain_path = Path(__file__).parent / "test_data" / "3cgs.cif"

    def test_simple_mmcif_read(self):
        sf = parse_file(self.test_trna_path)
        self.assertEqual(list(sf.chains.keys()), ["A"])
        seq = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA"
        example_chain = Chain("A")
        for i, s in enumerate(seq):
            example_chain.add_residue(Residue(s, s, i))
        self.assertEqual(sf["A"], example_chain)

    def test_simple_mmcif_read_two_chain(self):
        sf = parse_file(self.test_two_chain_path)
        self.assertEqual(list(sf.chains.keys()), ["A", "B"])
        chain_A = Chain("A")
        for i, s in enumerate("GCGCGUAGUAGC"):
            chain_A.add_residue(Residue(s, s, i))
        chain_B = Chain("B")
        for i, s in enumerate("CGCUACUGACGCG"):
            chain_B.add_residue(Residue(s, s, i))
        self.assertEqual(sf["A"], chain_A)
        self.assertEqual(sf["B"], chain_B)

    def test_residue_eq(self):
        # exact same
        self.assertEqual(
            Residue(three_letter_code="C", one_letter_code="C", index=0, atoms={}),
            Residue(three_letter_code="C", one_letter_code="C", index=0, atoms={}),
        )

        # different index
        self.assertNotEqual(
            Residue(three_letter_code="C", one_letter_code="C", index=0, atoms={}),
            Residue(three_letter_code="C", one_letter_code="C", index=1, atoms={}),
        )

        # totally different
        self.assertNotEqual(
            Residue(three_letter_code="C", one_letter_code="C", index=0, atoms={}),
            Residue(
                three_letter_code="G",
                one_letter_code="G",
                index=1,
                atoms={"P": (0, 0, 0)},
            ),
        )

        # different atoms
        self.assertNotEqual(
            Residue(three_letter_code="C", one_letter_code="C", index=0, atoms={}),
            Residue(
                three_letter_code="C",
                one_letter_code="C",
                index=0,
                atoms={"P": (0, 0, 0)},
            ),
        )

        # exactly same with atoms
        self.assertEqual(
            Residue(
                three_letter_code="C",
                one_letter_code="C",
                index=0,
                atoms={"B": (0, 0, 1), "A": (0, 0, 0)},
            ),
            Residue(
                three_letter_code="C",
                one_letter_code="C",
                index=0,
                atoms={"A": (0, 0, 0), "B": (0, 0, 1)},
            ),
        )

        # different atom coords
        self.assertNotEqual(
            Residue(
                three_letter_code="C",
                one_letter_code="C",
                index=0,
                atoms={"A": (1, 0, 1), "B": (0, 0, 0)},
            ),
            Residue(
                three_letter_code="C",
                one_letter_code="C",
                index=0,
                atoms={"A": (0, 0, 0), "B": (0, 0, 0)},
            ),
        )

    def test_chain_eq(self):
        a = Chain("A")
        b = Chain("A")

        a.residues = [
            Residue("A", "A", 0, atoms={"A": (0, 0, 0), "B": (0, 0, 0)}),
            Residue("A", "A", 1, atoms={"A": (0, 0, 0), "B": (0, 0, 0)}),
            Residue("A", "A", 2, atoms={"A": (0, 0, 0), "B": (0, 0, 0)}),
            Residue("A", "A", 3, atoms={"A": (0, 0, 0), "B": (0, 0, 0)}),
        ]

        b.residues = [
            Residue("A", "A", 0, atoms={"A": (0, 0, 0), "B": (0, 0, 0)}),
            Residue("A", "A", 1, atoms={"A": (0, 0, 0), "B": (0, 0, 0)}),
            Residue("A", "A", 2, atoms={"A": (0, 0, 0), "B": (0, 0, 0)}),
            Residue("A", "A", 3, atoms={"A": (0, 0, 0), "B": (0, 0, 0)}),
        ]
        # exactly the same
        self.assertEqual(a, b)

        b.residues = [
            Residue("A", "A", 0, atoms={"A": (0, 0, 0), "B": (0, 0, 0)}),
            Residue("A", "A", 1, atoms={"A": (0, 0, 0), "B": (0, 0, 0)}),
            Residue("A", "A", 2, atoms={}),
            Residue("A", "A", 3, atoms={"A": (0, 0, 0), "B": (0, 0, 0)}),
        ]
        # missing atoms
        self.assertNotEqual(a, b)

        b.residues = [
            Residue("A", "A", 0, atoms={"A": (0, 0, 0), "B": (0, 0, 0)}),
            Residue("C", "C", 1, atoms={"A": (0, 0, 0), "B": (0, 0, 0)}),
            Residue("G", "G", 2, atoms={"A": (0, 0, 0), "B": (0, 0, 0)}),
            Residue("U", "U", 3, atoms={"A": (0, 0, 0), "B": (0, 0, 0)}),
        ]
        # different sequences
        self.assertNotEqual(a, b)

        b = Chain("B")
        b.residues = [
            Residue("A", "A", 0, atoms={"A": (0, 0, 0), "B": (0, 0, 0)}),
            Residue("A", "A", 1, atoms={"A": (0, 0, 0), "B": (0, 0, 0)}),
            Residue("A", "A", 2, atoms={"A": (0, 0, 0), "B": (0, 0, 0)}),
            Residue("A", "A", 3, atoms={"A": (0, 0, 0), "B": (0, 0, 0)}),
        ]
        # different author_id, otherwise same
        self.assertEqual(a, b)

    def test_full_io(self):
        sf_read = parse_file(self.test_trna_path, include_atoms=True)
        with NamedTemporaryFile() as f:
            sf_read.write_mmcif_chain(f.name + ".cif", "A")
            sf_write = parse_file(f.name + ".cif", include_atoms=True)
            for auth_id in sf_write.chains.keys():
                self.assertEqual(sf_read[auth_id], sf_write[auth_id])

    def test_full_io_two_chains(self):
        sf_read = parse_file(self.test_two_chain_path, include_atoms=True)
        with NamedTemporaryFile() as f:
            for chain in sf_read:
                sf_read.write_mmcif_chain(
                    f.name + f"_{chain.author_id}.cif", chain.author_id
                )
                sf_write = parse_file(
                    f.name + f"_{chain.author_id}.cif", include_atoms=True
                )
                for auth_id in sf_write.chains.keys():
                    self.assertEqual(sf_read[auth_id], sf_write[auth_id])

    def test_backwards_compatibility(self):
        sf_old = parse_file(self.test_old_path, include_atoms=True)
        sf_new = parse_file(self.test_trna_path, include_atoms=True)
        self.assertEqual(sf_old.chains.keys(), sf_new.chains.keys())
        for k in sf_new.chains.keys():
            self.assertEqual(sf_old[k], sf_new[k])

    def test_atom_read(self):
        sf = parse_file(self.test_trna_path, include_atoms=True)
        self.assertDictEqual(
            sf["A"][0].atoms,
            {
                "OP3": (50.193, 51.19, 50.534),
                "P": (50.626, 49.73, 50.573),
                "OP1": (49.854, 48.893, 49.562),
                "OP2": (52.137, 49.542, 50.511),
                "O5'": (50.161, 49.136, 52.023),
                "C5'": (50.216, 49.948, 53.21),
                "C4'": (50.968, 49.231, 54.309),
                "O4'": (50.45, 47.888, 54.472),
                "C3'": (52.454, 49.03, 54.074),
                "O3'": (53.203, 50.177, 54.425),
                "C2'": (52.781, 47.831, 54.957),
                "O2'": (53.018, 48.156, 56.313),
                "C1'": (51.502, 47.007, 54.836),
                "N9": (51.628, 45.992, 53.798),
                "C8": (51.064, 46.007, 52.547),
                "N7": (51.379, 44.966, 51.831),
                "C5": (52.197, 44.218, 52.658),
                "C6": (52.848, 42.992, 52.425),
                "O6": (52.826, 42.291, 51.404),
                "N1": (53.588, 42.588, 53.534),
                "C2": (53.685, 43.282, 54.716),
                "N2": (54.452, 42.733, 55.671),
                "N3": (53.077, 44.429, 54.946),
                "C4": (52.356, 44.836, 53.879),
            },
        )


if __name__ == "__main__":
    unittest.main()
