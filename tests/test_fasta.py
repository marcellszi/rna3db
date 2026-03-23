import gzip
import tempfile
import unittest

from rna3db.parsers import fasta
from rna3db.parsers.fasta import Record


class TestRecord(unittest.TestCase):
    def test_fields(self):
        r = Record("seq1", "ACGU")
        self.assertEqual(r.header, "seq1")
        self.assertEqual(r.sequence, "ACGU")

    def test_namedtuple_unpacking(self):
        r = Record("h", "ACGU")
        header, sequence = r
        self.assertEqual(header, "h")
        self.assertEqual(sequence, "ACGU")


class TestFastaRead(unittest.TestCase):
    def test_simple(self):
        content = ">seq1\nACGU\n>seq2\nGGCC\n"
        with tempfile.NamedTemporaryFile("w", suffix=".fa") as f:
            f.write(content)
            f.flush()
            records = fasta.read(f.name)
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0], Record("seq1", "ACGU"))
        self.assertEqual(records[1], Record("seq2", "GGCC"))

    def test_multiline_sequence(self):
        content = ">seq1\nACGU\nGGCC\nUUAA\n"
        with tempfile.NamedTemporaryFile("w", suffix=".fa") as f:
            f.write(content)
            f.flush()
            records = fasta.read(f.name)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].sequence, "ACGUGGCCUUAA")

    def test_comments_skipped(self):
        content = "# comment\n>seq1\nACGU\n# another\n>seq2\nGGCC\n"
        with tempfile.NamedTemporaryFile("w", suffix=".fa") as f:
            f.write(content)
            f.flush()
            records = fasta.read(f.name)
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].sequence, "ACGU")

    def test_blank_lines_skipped(self):
        content = ">seq1\nACGU\n\n\n>seq2\nGGCC\n"
        with tempfile.NamedTemporaryFile("w", suffix=".fa") as f:
            f.write(content)
            f.flush()
            records = fasta.read(f.name)
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].sequence, "ACGU")
        self.assertEqual(records[1].sequence, "GGCC")

    def test_empty_file(self):
        with tempfile.NamedTemporaryFile("w", suffix=".fa") as f:
            records = fasta.read(f.name)
        self.assertEqual(records, [])

    def test_gzip(self):
        content = b">seq1\nACGU\n>seq2\nGGCC\n"
        with tempfile.NamedTemporaryFile(suffix=".fa.gz") as f:
            with gzip.open(f.name, "wb") as gz:
                gz.write(content)
            records = fasta.read(f.name)
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0], Record("seq1", "ACGU"))
        self.assertEqual(records[1], Record("seq2", "GGCC"))

    def test_force_gzip(self):
        """Non-.gz extension still read as gzip when force_gzip=True."""
        content = b">seq1\nACGU\n"
        with tempfile.NamedTemporaryFile(suffix=".fa") as f:
            with gzip.open(f.name, "wb") as gz:
                gz.write(content)
            records = fasta.read(f.name, force_gzip=True)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0], Record("seq1", "ACGU"))


class TestFastaWrite(unittest.TestCase):
    def test_write(self):
        records = [Record("seq1", "ACGU"), Record("seq2", "GGCC")]
        with tempfile.NamedTemporaryFile(suffix=".fa") as f:
            fasta.write(records, f.name)
            with open(f.name) as rf:
                content = rf.read()
        self.assertEqual(content, ">seq1\nACGU\n>seq2\nGGCC\n")

    def test_roundtrip(self):
        records = [
            Record("1ehz_A", "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUC"),
            Record("3cgs_A", "GCGCGUAGUAGC"),
        ]
        with tempfile.NamedTemporaryFile(suffix=".fa") as f:
            fasta.write(records, f.name)
            result = fasta.read(f.name)
        self.assertEqual(records, result)

    def test_write_empty(self):
        with tempfile.NamedTemporaryFile(suffix=".fa") as f:
            fasta.write([], f.name)
            with open(f.name) as rf:
                content = rf.read()
        self.assertEqual(content, "")


if __name__ == "__main__":
    unittest.main()
