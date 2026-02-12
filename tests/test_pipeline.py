import unittest
import tempfile
import runpy
import sys
import io

from contextlib import redirect_stderr
from pathlib import Path

from rna3db.utils import read_json


class TestCommandPipeline(unittest.TestCase):
    mmcif_path = Path(__file__).parent / "test_data" / "mmcifs"
    tbls_path = Path(__file__).parent / "test_data" / "tbls"

    def _run_rna3db(self):
        # suppress tqdm, etc.
        with redirect_stderr(io.StringIO()):
            runpy.run_module("rna3db", run_name="__main__")

    def test_default_integration(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sys.argv = [
                "rna3db",
                "parse",
                str(self.mmcif_path),
                tmpdir + "/parse.json",
            ]
            self._run_rna3db()
            parse_json = read_json(tmpdir + "/parse.json")

            sys.argv = [
                "rna3db",
                "filter",
                tmpdir + "/parse.json",
                tmpdir + "/filter.json",
            ]
            self._run_rna3db()
            filter_json = read_json(tmpdir + "/filter.json")

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

            sys.argv = [
                "rna3db",
                "split",
                tmpdir + "/cluster.json",
                tmpdir + "/split.json",
            ]
            self._run_rna3db()
            split_json = read_json(tmpdir + "/split.json")


if __name__ == "__main__":
    unittest.main()
