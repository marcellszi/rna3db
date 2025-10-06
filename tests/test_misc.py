import itertools
import unittest
import string

from rna3db.utils import to_case_insensitive, to_case_sensitive


class TestCaseSensitivity(unittest.TestCase):
    def test_to_insensitive(self):
        self.assertEqual(to_case_insensitive("a"), "-a-")
        self.assertEqual(to_case_insensitive("Aa"), "A-a-")
        self.assertEqual(to_case_insensitive("aBc"), "-a-B-c-")
        self.assertEqual(to_case_insensitive("1aA"), "1-a-A")
        self.assertEqual(to_case_insensitive("AaaaA"), "A-aaa-A")

    def test_to_sensitive(self):
        self.assertEqual(to_case_sensitive("-a-"), "a")
        self.assertEqual(to_case_sensitive("-A-"), "a")

    def test_exhaustive_three(self):
        for s in map("".join, itertools.product(string.ascii_letters, repeat=3)):
            insensitive = to_case_insensitive(s)
            self.assertEqual(to_case_sensitive(insensitive), s)
            self.assertEqual(to_case_sensitive(insensitive.upper()), s)
            self.assertEqual(to_case_sensitive(insensitive.lower()), s)


if __name__ == "__main__":
    unittest.main()
