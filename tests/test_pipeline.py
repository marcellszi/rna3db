"""
Pipeline / CLI integration tests.

Fast tests (no CIF parsing) use pre-built JSON constants and run in milliseconds.
Slow tests parse real mmCIF files or run the full pipeline; they are skipped
unless the environment variable RNA3DB_SLOW_TESTS=1 is set.

    RNA3DB_SLOW_TESTS=1 python -m unittest tests.test_pipeline
"""

import io
import os
import runpy
import sys
import tempfile
import unittest

from contextlib import redirect_stderr
from pathlib import Path

from rna3db.utils import read_json, write_json


SLOW_TESTS = os.environ.get("RNA3DB_SLOW_TESTS")
slow = unittest.skipUnless(SLOW_TESTS, "Set RNA3DB_SLOW_TESTS=1 to run slow tests")

# Pre-built parse output matching rna3db parse tests/test_data/mmcifs/.
# 3cgs_A (12 nt) and 3cgs_B (13 nt) are filtered by the default min_length=32.
# 3tup_T (3.05 Å) is filtered when max_resolution < 3.05.
# All others pass the default filters.
_SAMPLE_PARSED = {
    "1ehz_A": {
        "release_date": "2000-10-02",
        "structure_method": "x-ray diffraction",
        "resolution": 1.93,
        "length": 76,
        "sequence": "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA",
    },
    "1y27_X": {
        "release_date": "2004-12-28",
        "structure_method": "x-ray diffraction",
        "resolution": 2.4,
        "length": 68,
        "sequence": "GGAUCAUAUAAUCGCGUGGAUAUGGCACGCAAGUUUCUACCGGGCACCGUAAAUGUCCGACUAUGGUC",
    },
    "3cgs_A": {
        "release_date": "2008-07-01",
        "structure_method": "x-ray diffraction",
        "resolution": 1.65,
        "length": 12,
        "sequence": "GCGCGUAGUAGC",
    },
    "3cgs_B": {
        "release_date": "2008-07-01",
        "structure_method": "x-ray diffraction",
        "resolution": 1.65,
        "length": 13,
        "sequence": "CGCUACUGACGCG",
    },
    "3tup_T": {
        "release_date": "2011-11-23",
        "structure_method": "x-ray diffraction",
        "resolution": 3.05,
        "length": 76,
        "sequence": "GCCGAGGUAGCUCAGUUGGUAGAGCAUGCGACUGAAAAUCGCAGUGUCGGCGGUUCGAUUCUGCUCCUCGGCACCA",
    },
    "5m0h_A": {
        "release_date": "2017-01-18",
        "structure_method": "x-ray diffraction",
        "resolution": 2.65,
        "length": 42,
        "sequence": "GAUAACUGAAUCGCUAAGGAUGAAAGUCUAUGCGACAUUAUC",
    },
}

# Sequence-clustered JSON (output of cluster_sequences), one cluster per chain.
# Chains that pass the default filter; all appear as query_name in test_data/tbls.
_SAMPLE_SEQ_CLUSTERED = {
    "1ehz_A": {"1ehz_A": _SAMPLE_PARSED["1ehz_A"]},
    "1y27_X": {"1y27_X": _SAMPLE_PARSED["1y27_X"]},
    "3tup_T": {"3tup_T": _SAMPLE_PARSED["3tup_T"]},
    "5m0h_A": {"5m0h_A": _SAMPLE_PARSED["5m0h_A"]},
}

# Structure-clustered JSON (output of cluster_structures on _SAMPLE_SEQ_CLUSTERED
# with tests/test_data/tbls and the default e-value cutoff of 1.0).
# 5m0h_A has no Infernal hits below 1.0, so it goes into component_0.
# 1ehz_A and 3tup_T both hit tRNA (RF00005) and share component_1.
# 1y27_X hits a different family and is placed in component_2.
_SAMPLE_CLUSTERED = {
    "component_0": {
        "5m0h_A": {"5m0h_A": _SAMPLE_PARSED["5m0h_A"]},
    },
    "component_1": {
        "3tup_T": {"3tup_T": _SAMPLE_PARSED["3tup_T"]},
        "1ehz_A": {"1ehz_A": _SAMPLE_PARSED["1ehz_A"]},
    },
    "component_2": {
        "1y27_X": {"1y27_X": _SAMPLE_PARSED["1y27_X"]},
    },
}


class _PipelineBase(unittest.TestCase):
    def _run_rna3db(self):
        with redirect_stderr(io.StringIO()):
            runpy.run_module("rna3db", run_name="__main__")


class TestCommandLineFilter(_PipelineBase):
    def _write_parsed(self, tmpdir):
        path = tmpdir + "/parsed.json"
        write_json(_SAMPLE_PARSED, path)
        return path

    def test_filter_default(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            parsed = self._write_parsed(tmpdir)
            filtered = tmpdir + "/filtered.json"
            sys.argv = ["rna3db", "filter", parsed, filtered]
            self._run_rna3db()

            result = read_json(filtered)
            self.assertIn("1ehz_A", result)
            self.assertIn("1y27_X", result)
            self.assertIn("3tup_T", result)
            self.assertIn("5m0h_A", result)
            self.assertNotIn("3cgs_A", result)  # 12 nt, too short
            self.assertNotIn("3cgs_B", result)  # 13 nt, too short

    def test_filter_permissive_min_length(self):
        """--min_length 0 disables the length filter."""
        with tempfile.TemporaryDirectory() as tmpdir:
            parsed = self._write_parsed(tmpdir)
            filtered = tmpdir + "/filtered.json"
            sys.argv = ["rna3db", "filter", parsed, filtered, "--min_length", "0"]
            self._run_rna3db()

            result = read_json(filtered)
            self.assertIn("3cgs_A", result)  # 12 nt, now passes
            self.assertIn("3cgs_B", result)  # 13 nt, now passes

    def test_filter_strict_min_length(self):
        """--min_length 75 keeps only chains at or above 75 nt."""
        with tempfile.TemporaryDirectory() as tmpdir:
            parsed = self._write_parsed(tmpdir)
            filtered = tmpdir + "/filtered.json"
            sys.argv = ["rna3db", "filter", parsed, filtered, "--min_length", "75"]
            self._run_rna3db()

            result = read_json(filtered)
            self.assertIn("1ehz_A", result)  # 76 nt, passes
            self.assertIn("3tup_T", result)  # 76 nt, passes
            self.assertNotIn("1y27_X", result)  # 68 nt, filtered
            self.assertNotIn("5m0h_A", result)  # 42 nt, filtered

    def test_filter_strict_resolution(self):
        """--max_resolution 3.0 filters 3tup_T (3.05 Å)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            parsed = self._write_parsed(tmpdir)
            filtered = tmpdir + "/filtered.json"
            sys.argv = ["rna3db", "filter", parsed, filtered, "--max_resolution", "3.0"]
            self._run_rna3db()

            result = read_json(filtered)
            self.assertNotIn("3tup_T", result)  # 3.05 Å > 3.0 cutoff
            self.assertIn("1ehz_A", result)  # 1.93 Å, passes

    def test_filter_log_path(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            parsed = self._write_parsed(tmpdir)
            filtered = tmpdir + "/filtered.json"
            log = tmpdir + "/filter_log.json"
            sys.argv = ["rna3db", "filter", parsed, filtered, "--filter_log_path", log]
            self._run_rna3db()

            self.assertTrue(Path(log).exists())
            log_data = read_json(log)
            self.assertIn("3cgs_A", log_data)
            self.assertIn("is_short_sequence", log_data["3cgs_A"])
            self.assertEqual(log_data["1ehz_A"], [])  # no filters hit

    def test_filter_all_disabled(self):
        """Setting all thresholds to 0 keeps every chain."""
        with tempfile.TemporaryDirectory() as tmpdir:
            parsed = self._write_parsed(tmpdir)
            filtered = tmpdir + "/filtered.json"
            sys.argv = [
                "rna3db",
                "filter",
                parsed,
                filtered,
                "--min_length",
                "0",
                "--max_resolution",
                "0",
                "--single_ratio_cutoff",
                "0",
                "--max_unknown_ratio",
                "0",
            ]
            self._run_rna3db()

            result = read_json(filtered)
            self.assertEqual(set(result.keys()), set(_SAMPLE_PARSED.keys()))


class TestCommandLineSplit(_PipelineBase):
    def _write_clustered(self, tmpdir):
        path = tmpdir + "/clustered.json"
        write_json(_SAMPLE_CLUSTERED, path)
        return path

    def test_split_default(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            clustered = self._write_clustered(tmpdir)
            split_out = tmpdir + "/split.json"
            sys.argv = ["rna3db", "split", clustered, split_out]
            self._run_rna3db()

            result = read_json(split_out)
            self.assertIn("train_set", result)
            self.assertIn("valid_set", result)
            self.assertIn("test_set", result)
            total = sum(len(v) for v in result.values())
            self.assertEqual(total, 3)  # all 3 components accounted for

    def test_split_custom_ratios(self):
        """Equal three-way split puts one component in each set."""
        with tempfile.TemporaryDirectory() as tmpdir:
            clustered = self._write_clustered(tmpdir)
            split_out = tmpdir + "/split.json"
            sys.argv = [
                "rna3db",
                "split",
                clustered,
                split_out,
                "--train_ratio",
                "0.3",
                "--valid_ratio",
                "0.3",
            ]
            self._run_rna3db()

            result = read_json(split_out)
            self.assertEqual(len(result["train_set"]), 1)
            self.assertEqual(len(result["valid_set"]), 1)
            self.assertEqual(len(result["test_set"]), 1)

    def test_split_force_zero_test(self):
        """--force_zero_test puts component_0 in the test set."""
        with tempfile.TemporaryDirectory() as tmpdir:
            clustered = self._write_clustered(tmpdir)
            split_out = tmpdir + "/split.json"
            sys.argv = [
                "rna3db",
                "split",
                clustered,
                split_out,
                "--train_ratio",
                "0.5",
                "--valid_ratio",
                "0.0",
                "--force_zero_test",
            ]
            self._run_rna3db()

            result = read_json(split_out)
            self.assertIn("component_0", result["test_set"])
            self.assertNotIn("component_0", result["train_set"])
            self.assertNotIn("component_0", result["valid_set"])

    def test_split_large_train_ratio(self):
        """train_ratio=1.0 puts all components into train."""
        with tempfile.TemporaryDirectory() as tmpdir:
            clustered = self._write_clustered(tmpdir)
            split_out = tmpdir + "/split.json"
            sys.argv = [
                "rna3db",
                "split",
                clustered,
                split_out,
                "--train_ratio",
                "1.0",
                "--valid_ratio",
                "0.0",
            ]
            self._run_rna3db()

            result = read_json(split_out)
            self.assertEqual(len(result["valid_set"]), 0)
            self.assertEqual(len(result["test_set"]), 0)
            self.assertEqual(len(result["train_set"]), 3)


class TestCommandLineCluster(_PipelineBase):
    tbls_path = Path(__file__).parent / "test_data" / "tbls"

    def _write_seq_clustered(self, tmpdir):
        path = tmpdir + "/seq_clustered.json"
        write_json(_SAMPLE_SEQ_CLUSTERED, path)
        return path

    def test_cluster_only_structure(self):
        """--only_structure skips mmseqs2 and clusters by Infernal hits."""
        with tempfile.TemporaryDirectory() as tmpdir:
            seq_clustered = self._write_seq_clustered(tmpdir)
            struct_clustered = tmpdir + "/struct_clustered.json"
            sys.argv = [
                "rna3db",
                "cluster",
                seq_clustered,
                struct_clustered,
                "--tbl_dir",
                str(self.tbls_path),
                "--only_structure",
            ]
            self._run_rna3db()

            result = read_json(struct_clustered)
            # 5m0h_A has no hits below the default e-value cutoff of 1.0
            self.assertIn("5m0h_A", result["component_0"])
            # 1ehz_A and 3tup_T both hit tRNA (RF00005) at very low e-values
            # and are grouped into the same component
            all_chains = {c for comp in result.values() for c in comp.keys()}
            self.assertIn("1ehz_A", all_chains)
            self.assertIn("1y27_X", all_chains)
            self.assertIn("3tup_T", all_chains)

    def test_cluster_tight_evalue_cutoff(self):
        """A very tight e-value cutoff puts all chains into component_0."""
        with tempfile.TemporaryDirectory() as tmpdir:
            seq_clustered = self._write_seq_clustered(tmpdir)
            struct_clustered = tmpdir + "/struct_clustered.json"
            sys.argv = [
                "rna3db",
                "cluster",
                seq_clustered,
                struct_clustered,
                "--tbl_dir",
                str(self.tbls_path),
                "--only_structure",
                "--structural_e_value_cutoff",
                "1e-20",
            ]
            self._run_rna3db()

            result = read_json(struct_clustered)
            # No hits pass 1e-20; all chains go to component_0
            self.assertIn("component_0", result)
            self.assertNotIn("component_1", result)


@slow
class TestCommandLineParseSlow(_PipelineBase):
    mmcif_path = Path(__file__).parent / "test_data" / "mmcifs"

    def test_parse_default(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sys.argv = [
                "rna3db",
                "parse",
                str(self.mmcif_path),
                tmpdir + "/parsed.json",
            ]
            self._run_rna3db()
            result = read_json(tmpdir + "/parsed.json")
            self.assertIn("1ehz_A", result)
            self.assertIn("1y27_X", result)
            self.assertNotIn("atoms", result["1ehz_A"])

    def test_parse_include_atoms(self):
        """--include_atoms populates per-residue atom coordinate lists."""
        with tempfile.TemporaryDirectory() as tmpdir:
            sys.argv = [
                "rna3db",
                "parse",
                str(self.mmcif_path),
                tmpdir + "/parsed.json",
                "--include_atoms",
            ]
            self._run_rna3db()
            result = read_json(tmpdir + "/parsed.json")
            self.assertIn("atoms", result["1ehz_A"])
            self.assertGreater(len(result["1ehz_A"]["atoms"]), 0)

    def test_parse_nmr_resolution(self):
        """--nmr_resolution is accepted and does not break non-NMR parsing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            sys.argv = [
                "rna3db",
                "parse",
                str(self.mmcif_path),
                tmpdir + "/parsed.json",
                "--nmr_resolution",
                "3.0",
            ]
            self._run_rna3db()
            result = read_json(tmpdir + "/parsed.json")
            self.assertIn("1ehz_A", result)

    def test_parse_cpu(self):
        """--cpu is accepted and results are unchanged."""
        with tempfile.TemporaryDirectory() as tmpdir:
            sys.argv = [
                "rna3db",
                "--cpu",
                "1",
                "parse",
                str(self.mmcif_path),
                tmpdir + "/parsed.json",
            ]
            self._run_rna3db()
            result = read_json(tmpdir + "/parsed.json")
            self.assertIn("1ehz_A", result)

    def test_parse_output_keys_are_sorted(self):
        """Parse output keys are sorted alphabetically (chain_id order)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out = tmpdir + "/parsed.json"
            sys.argv = ["rna3db", "parse", str(self.mmcif_path), out]
            self._run_rna3db()
            keys = list(read_json(out).keys())
            self.assertEqual(keys, sorted(keys))


@slow
class TestCommandPipeline(_PipelineBase):
    mmcif_path = Path(__file__).parent / "test_data" / "mmcifs"
    tbls_path = Path(__file__).parent / "test_data" / "tbls"

    def test_default_integration(self):
        """Parse → filter → cluster → split with all defaults."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # parse
            sys.argv = [
                "rna3db",
                "parse",
                str(self.mmcif_path),
                tmpdir + "/parse.json",
            ]
            self._run_rna3db()
            parse_json = read_json(tmpdir + "/parse.json")
            self.assertIn("1ehz_A", parse_json)
            self.assertIn("1y27_X", parse_json)

            # filter
            sys.argv = [
                "rna3db",
                "filter",
                tmpdir + "/parse.json",
                tmpdir + "/filter.json",
            ]
            self._run_rna3db()
            filter_json = read_json(tmpdir + "/filter.json")
            # 3cgs chains are 12/13 nt, below the default 32 nt minimum
            self.assertNotIn("3cgs_A", filter_json)
            self.assertNotIn("3cgs_B", filter_json)
            self.assertIn("1ehz_A", filter_json)

            # cluster (sequence + structure)
            sys.argv = [
                "rna3db",
                "cluster",
                tmpdir + "/filter.json",
                tmpdir + "/cluster.json",
                "--tbl_dir",
                str(self.tbls_path),
            ]
            self._run_rna3db()
            cluster_json = read_json(tmpdir + "/cluster.json")
            self.assertGreater(len(cluster_json), 0)

            # split
            sys.argv = [
                "rna3db",
                "split",
                tmpdir + "/cluster.json",
                tmpdir + "/split.json",
            ]
            self._run_rna3db()
            split_json = read_json(tmpdir + "/split.json")
            self.assertIn("train_set", split_json)
            self.assertIn("valid_set", split_json)
            self.assertIn("test_set", split_json)

    def test_run_command(self):
        """The `run` subcommand produces all four intermediate files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            sys.argv = [
                "rna3db",
                "run",
                str(self.mmcif_path),
                tmpdir,
                "--tbl_dir",
                str(self.tbls_path),
            ]
            self._run_rna3db()

            for fname in (
                "parsed.json",
                "filtered.json",
                "clustered.json",
                "split.json",
            ):
                self.assertTrue(Path(tmpdir, fname).exists(), f"{fname} not found")

            split_json = read_json(tmpdir + "/split.json")
            self.assertIn("train_set", split_json)
            self.assertIn("test_set", split_json)

    def test_integration_include_atoms(self):
        """Atoms are preserved all the way through parse → filter."""
        with tempfile.TemporaryDirectory() as tmpdir:
            sys.argv = [
                "rna3db",
                "parse",
                str(self.mmcif_path),
                tmpdir + "/parse.json",
                "--include_atoms",
            ]
            self._run_rna3db()
            parse_json = read_json(tmpdir + "/parse.json")
            self.assertIn("atoms", parse_json["1ehz_A"])

            sys.argv = [
                "rna3db",
                "filter",
                tmpdir + "/parse.json",
                tmpdir + "/filter.json",
            ]
            self._run_rna3db()
            filter_json = read_json(tmpdir + "/filter.json")
            # Atom data should survive the filter step unchanged
            self.assertIn("atoms", filter_json["1ehz_A"])

    def test_integration_two_phase_cluster(self):
        """Sequence-then-structure clustering via two separate cluster calls."""
        with tempfile.TemporaryDirectory() as tmpdir:
            sys.argv = [
                "rna3db",
                "parse",
                str(self.mmcif_path),
                tmpdir + "/parse.json",
            ]
            self._run_rna3db()

            sys.argv = [
                "rna3db",
                "filter",
                tmpdir + "/parse.json",
                tmpdir + "/filter.json",
            ]
            self._run_rna3db()

            # Phase 1: sequence clustering only
            sys.argv = [
                "rna3db",
                "cluster",
                tmpdir + "/filter.json",
                tmpdir + "/seq_clustered.json",
                "--only_sequence",
            ]
            self._run_rna3db()

            # Phase 2: structure clustering only (input is seq-clustered JSON)
            sys.argv = [
                "rna3db",
                "cluster",
                tmpdir + "/seq_clustered.json",
                tmpdir + "/cluster.json",
                "--tbl_dir",
                str(self.tbls_path),
                "--only_structure",
            ]
            self._run_rna3db()
            cluster_json = read_json(tmpdir + "/cluster.json")
            self.assertGreater(len(cluster_json), 0)


if __name__ == "__main__":
    unittest.main()
