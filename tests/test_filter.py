import json
import tempfile
import unittest

from rna3db.filter import apply_filters


def _chain(sequence: str, resolution: float = 2.0) -> dict:
    return {
        "release_date": "2020-01-01",
        "structure_method": "x-ray diffraction",
        "resolution": resolution,
        "length": len(sequence),
        "sequence": sequence,
    }


# Mixed sequences that pass all default filters (no homopolymer bias, no Ns)
_MIXED_40 = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUC"  # 40 nt
_MIXED_32 = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGAC"  # 32 nt
_MIXED_31 = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGA"  # 31 nt


class TestApplyFilters(unittest.TestCase):
    def test_all_pass(self):
        data = {"A": _chain(_MIXED_40)}
        result = apply_filters(data)
        self.assertIn("A", result)

    def test_short_filtered(self):
        data = {
            "short": _chain("ACGU"),  # 4 nt < default min_length=32
            "long": _chain(_MIXED_32),  # 32 nt, passes
        }
        result = apply_filters(data, min_length=32)
        self.assertNotIn("short", result)
        self.assertIn("long", result)

    def test_short_boundary(self):
        data = {
            "exact": _chain(_MIXED_32),  # exactly at boundary, passes
            "below": _chain(_MIXED_31),  # one below, filtered
        }
        result = apply_filters(data, min_length=32)
        self.assertIn("exact", result)
        self.assertNotIn("below", result)

    def test_low_resolution_filtered(self):
        data = {
            "bad": _chain(_MIXED_40, resolution=10.0),
            "good": _chain(_MIXED_40, resolution=2.0),
        }
        result = apply_filters(data, max_resolution=9.0)
        self.assertNotIn("bad", result)
        self.assertIn("good", result)

    def test_single_ratio_filtered(self):
        # 90% G → above 0.8 cutoff
        data = {"A": _chain("G" * 36 + "ACGU")}
        result = apply_filters(data, single_ratio_cutoff=0.8)
        self.assertNotIn("A", result)

    def test_single_ratio_passes(self):
        data = {"A": _chain("GCGCGUAGUAGCGCGCGUAGUAGCGCGCGUAGUAGCGCGC")}
        result = apply_filters(data, single_ratio_cutoff=0.8)
        self.assertIn("A", result)

    def test_many_unknowns_filtered(self):
        # 50% N → above 0.3 cutoff
        data = {"A": _chain("ACGU" * 5 + "N" * 20)}
        result = apply_filters(data, max_unknown_ratio=0.3)
        self.assertNotIn("A", result)

    def test_many_unknowns_passes(self):
        # ~5% N → well below 0.3 cutoff
        data = {"A": _chain("ACGU" * 10 + "NN")}
        result = apply_filters(data, max_unknown_ratio=0.3)
        self.assertIn("A", result)

    def test_disabled_zero(self):
        """Passing 0 for any threshold disables that filter."""
        data = {
            "short": _chain("ACGU"),
            "high_res": _chain(_MIXED_40, resolution=100.0),
        }
        result = apply_filters(
            data,
            min_length=0,
            max_resolution=0,
            single_ratio_cutoff=0,
            max_unknown_ratio=0,
        )
        self.assertIn("short", result)
        self.assertIn("high_res", result)

    def test_disabled_none(self):
        """Passing None for any threshold disables that filter."""
        data = {"short": _chain("ACGU")}
        result = apply_filters(
            data,
            min_length=None,
            max_resolution=None,
            single_ratio_cutoff=None,
            max_unknown_ratio=None,
        )
        self.assertIn("short", result)

    def test_multiple_filters(self):
        data = {
            "bad": _chain("ACGU", resolution=10.0),  # too short AND high res
            "good": _chain(_MIXED_40, resolution=2.0),
        }
        result = apply_filters(data, min_length=32, max_resolution=9.0)
        self.assertNotIn("bad", result)
        self.assertIn("good", result)

    def test_filter_log(self):
        data = {
            "short": _chain("ACGU"),
            "good": _chain(_MIXED_40),
        }
        with tempfile.NamedTemporaryFile(suffix=".json") as f:
            apply_filters(data, min_length=32, filter_log_path=f.name)
            with open(f.name) as rf:
                log = json.load(rf)
        self.assertIn("is_short_sequence", log["short"])
        self.assertEqual(log["good"], [])

    def test_filter_log_records_all_chains(self):
        """Filter log includes every input chain, even those that pass."""
        data = {
            "good": _chain(_MIXED_40),
            "bad": _chain("ACGU"),
        }
        with tempfile.NamedTemporaryFile(suffix=".json") as f:
            apply_filters(data, min_length=32, filter_log_path=f.name)
            with open(f.name) as rf:
                log = json.load(rf)
        self.assertIn("good", log)
        self.assertIn("bad", log)

    def test_empty_data(self):
        result = apply_filters({})
        self.assertEqual(result, {})


if __name__ == "__main__":
    unittest.main()
