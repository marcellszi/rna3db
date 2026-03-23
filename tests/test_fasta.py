import gzip
import tempfile
import unittest

from rna3db.parsers import FASTA


class TestFASTA(unittest.TestCase):
    def test_construction(self):
        f = FASTA(["seq1", "seq2"], ["ACGU", "GGCC"])
        self.assertEqual(f.headers, ["seq1", "seq2"])
        self.assertEqual(f.sequences, ["ACGU", "GGCC"])

    def test_len(self):
        self.assertEqual(len(FASTA([], [])), 0)
        self.assertEqual(len(FASTA(["h"], ["ACGU"])), 1)

    def test_getitem(self):
        f = FASTA(["seq1", "seq2"], ["ACGU", "GGCC"])
        self.assertEqual(f[0], ("seq1", "ACGU"))
        self.assertEqual(f[1], ("seq2", "GGCC"))

    def test_iter(self):
        f = FASTA(["seq1", "seq2"], ["ACGU", "GGCC"])
        pairs = list(f)
        self.assertEqual(pairs, [("seq1", "ACGU"), ("seq2", "GGCC")])

    def test_eq(self):
        a = FASTA(["seq1"], ["ACGU"])
        b = FASTA(["seq1"], ["ACGU"])
        c = FASTA(["seq2"], ["ACGU"])
        self.assertEqual(a, b)
        self.assertNotEqual(a, c)

    def test_unpack(self):
        f = FASTA(["seq1", "seq2"], ["ACGU", "GGCC"])
        headers, sequences = f.unpack()
        self.assertEqual(headers, ["seq1", "seq2"])
        self.assertEqual(sequences, ["ACGU", "GGCC"])

    def test_eq_different_type(self):
        f = FASTA(["seq1"], ["ACGU"])
        self.assertNotEqual(f, [("seq1", "ACGU")])


class TestFASTARead(unittest.TestCase):
    def test_simple(self):
        content = ">seq1\nACGU\n>seq2\nGGCC\n"
        with tempfile.NamedTemporaryFile("w", suffix=".fa") as f:
            f.write(content)
            f.flush()
            result = FASTA.read(f.name)
        self.assertEqual(len(result), 2)
        self.assertEqual(result.headers, ["seq1", "seq2"])
        self.assertEqual(result.sequences, ["ACGU", "GGCC"])

    def test_multiline_sequence(self):
        content = ">seq1\nACGU\nGGCC\nUUAA\n"
        with tempfile.NamedTemporaryFile("w", suffix=".fa") as f:
            f.write(content)
            f.flush()
            result = FASTA.read(f.name)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.sequences[0], "ACGUGGCCUUAA")

    def test_comments_skipped(self):
        content = "# comment\n>seq1\nACGU\n# another\n>seq2\nGGCC\n"
        with tempfile.NamedTemporaryFile("w", suffix=".fa") as f:
            f.write(content)
            f.flush()
            result = FASTA.read(f.name)
        self.assertEqual(len(result), 2)
        self.assertEqual(result.sequences[0], "ACGU")

    def test_blank_lines_skipped(self):
        content = ">seq1\nACGU\n\n\n>seq2\nGGCC\n"
        with tempfile.NamedTemporaryFile("w", suffix=".fa") as f:
            f.write(content)
            f.flush()
            result = FASTA.read(f.name)
        self.assertEqual(len(result), 2)
        self.assertEqual(result.sequences[0], "ACGU")
        self.assertEqual(result.sequences[1], "GGCC")

    def test_empty_file(self):
        with tempfile.NamedTemporaryFile("w", suffix=".fa") as f:
            result = FASTA.read(f.name)
        self.assertEqual(len(result), 0)
        self.assertEqual(result, FASTA([], []))

    def test_gzip(self):
        content = b">seq1\nACGU\n>seq2\nGGCC\n"
        with tempfile.NamedTemporaryFile(suffix=".fa.gz") as f:
            with gzip.open(f.name, "wb") as gz:
                gz.write(content)
            result = FASTA.read(f.name)
        self.assertEqual(len(result), 2)
        self.assertEqual(result.headers, ["seq1", "seq2"])
        self.assertEqual(result.sequences, ["ACGU", "GGCC"])

    def test_force_gzip(self):
        """Non-.gz extension still read as gzip when force_gzip=True."""
        content = b">seq1\nACGU\n"
        with tempfile.NamedTemporaryFile(suffix=".fa") as f:
            with gzip.open(f.name, "wb") as gz:
                gz.write(content)
            result = FASTA.read(f.name, force_gzip=True)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.headers[0], "seq1")
        self.assertEqual(result.sequences[0], "ACGU")


class TestFASTAWrite(unittest.TestCase):
    def test_write(self):
        f = FASTA(["seq1", "seq2"], ["ACGU", "GGCC"])
        with tempfile.NamedTemporaryFile(suffix=".fa") as tmp:
            f.write(tmp.name)
            with open(tmp.name) as rf:
                content = rf.read()
        self.assertEqual(content, ">seq1\nACGU\n>seq2\nGGCC\n")

    def test_roundtrip(self):
        original = FASTA(
            ["1ehz_A", "3cgs_A"],
            ["GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUC", "GCGCGUAGUAGC"],
        )
        with tempfile.NamedTemporaryFile(suffix=".fa") as tmp:
            original.write(tmp.name)
            result = FASTA.read(tmp.name)
        self.assertEqual(original, result)

    def test_write_empty(self):
        f = FASTA([], [])
        with tempfile.NamedTemporaryFile(suffix=".fa") as tmp:
            f.write(tmp.name)
            with open(tmp.name) as rf:
                content = rf.read()
        self.assertEqual(content, "")


if __name__ == "__main__":
    unittest.main()
