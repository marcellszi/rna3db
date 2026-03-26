import contextlib
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path

from rna3db.parsers import Table
from rna3db.parsers.tabular import Hit

TBL_STR = (
    "#target name         accession query name           accession mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc description of target\n"
    "#------------------- --------- -------------------- --------- --- -------- -------- -------- -------- ------ ----- ---- ---- ----- ------ --------- --- ---------------------\n"
    "mir-4850             RF03522   7lhd_A               -          cm        1       88      914      859      -    no    1 0.52   0.0   24.6     0.099 ?   mir-4850 microRNA precursor family\n"
    "Lysine               RF00168   7lhd_A               -          cm       83      168     4093     4176      +    no    1 0.48   0.0   18.1      0.18 ?   Lysine riboswitch\n"
    "tRNA                 RF00005   7osa_PSIT            -          cm        1       71        1       74      +    no    1 0.68   0.0   60.9   7.2e-14 !   tRNA"
)
TEST_DICTS = [
    {
        "target_name": "mir-4850",
        "target_accession": "RF03522",
        "query_name": "7lhd_A",
        "query_accession": None,
        "mdl": "cm",
        "mdl_from": 1,
        "mdl_to": 88,
        "seq_from": 914,
        "seq_to": 859,
        "strand": "-",
        "trunc": False,
        "pass_n": 1,
        "gc": 0.52,
        "bias": 0.0,
        "score": 24.6,
        "e_value": 0.099,
        "inc": "?",
        "description_of_target": "mir-4850 microRNA precursor family",
    },
    {
        "target_name": "Lysine",
        "target_accession": "RF00168",
        "query_name": "7lhd_A",
        "query_accession": None,
        "mdl": "cm",
        "mdl_from": 83,
        "mdl_to": 168,
        "seq_from": 4093,
        "seq_to": 4176,
        "strand": "+",
        "trunc": False,
        "pass_n": 1,
        "gc": 0.48,
        "bias": 0.0,
        "score": 18.1,
        "e_value": 0.18,
        "inc": "?",
        "description_of_target": "Lysine riboswitch",
    },
    {
        "target_name": "tRNA",
        "target_accession": "RF00005",
        "query_name": "7osa_PSIT",
        "query_accession": None,
        "mdl": "cm",
        "mdl_from": 1,
        "mdl_to": 71,
        "seq_from": 1,
        "seq_to": 74,
        "strand": "+",
        "trunc": False,
        "pass_n": 1,
        "gc": 0.68,
        "bias": 0.0,
        "score": 60.9,
        "e_value": 7.2e-14,
        "inc": "!",
        "description_of_target": "tRNA",
    },
]


class TestTabularParser(unittest.TestCase):
    def setUp(self):
        with self.tmp_txt(TBL_STR) as f:
            self.tbl = Table(f.name)

    @contextlib.contextmanager
    def tmp_txt(self, s: str):
        tmp = tempfile.NamedTemporaryFile("w")
        tmp.write(TBL_STR)
        tmp.seek(0)
        try:
            yield tmp
        finally:
            tmp.close()

    def _assert_hit(self, hit: Hit, test_dict: dict):
        for k, v in test_dict.items():
            self.assertEqual(hit.__getattribute__(k), v)

    def test_row_parse(self):
        for row_str, test_dict in zip(TBL_STR.split("\n")[2:], TEST_DICTS):
            hit = Table._parse_tbl_row(row_str)
            self._assert_hit(hit, test_dict)

    def test_parse_tbl(self):
        tbl = Table(hits=[])
        with self.tmp_txt(TBL_STR) as f:
            entries = tbl._parse_tbl(f.name)
            for hit, test_dict in zip(entries, TEST_DICTS):
                self._assert_hit(hit, test_dict)

    def test_len(self):
        empty_hit = Hit(*[None] * 18)
        self.assertEqual(len(Table(hits=[])), 0)
        self.assertEqual(len(Table(hits=[empty_hit] * 1)), 1)
        self.assertEqual(len(Table(hits=[empty_hit] * 1337)), 1337)

    def test_tophits(self):
        tbl = self.tbl.top_hits
        self.assertEqual(len(tbl), 2)
        for actual_hit, expected_hit_dict in zip(tbl, [TEST_DICTS[0], TEST_DICTS[2]]):
            self._assert_hit(actual_hit, expected_hit_dict)

    def test_get_col(self):
        actual_list = defaultdict(list)
        for d in TEST_DICTS:
            for k, v in d.items():
                actual_list[k].append(v)

        for k, v in actual_list.items():
            self.assertEqual(self.tbl.__getattribute__(k), v)


class TestTableFiltering(unittest.TestCase):
    # TBL_STR e-values: 0.099 (7lhd_A/mir-4850), 0.18 (7lhd_A/Lysine), 7.2e-14 (7osa_PSIT/tRNA)

    def setUp(self):
        self._tmp = tempfile.NamedTemporaryFile("w")
        self._tmp.write(TBL_STR)
        self._tmp.flush()
        self.tbl = Table(self._tmp.name)

    def tearDown(self):
        self._tmp.close()

    def test_filter_e_value(self):
        # cutoff 0.1 keeps 7.2e-14 and 0.099; sorts ascending
        filtered = self.tbl.filter_e_value(0.1)
        self.assertEqual(len(filtered), 2)
        self.assertEqual(filtered.e_value, [7.2e-14, 0.099])

    def test_filter_e_value_removes_all(self):
        filtered = self.tbl.filter_e_value(1e-30)
        self.assertEqual(len(filtered), 0)

    def test_filter_e_value_keeps_all(self):
        filtered = self.tbl.filter_e_value(1.0)
        self.assertEqual(len(filtered), 3)

    def test_filter_attr_by_set(self):
        filtered = self.tbl.filter_attr_by_set("query_name", {"7lhd_A"})
        self.assertEqual(len(filtered), 2)
        self.assertTrue(all(q == "7lhd_A" for q in filtered.query_name))

    def test_filter_attr_by_set_multiple(self):
        filtered = self.tbl.filter_attr_by_set("target_name", {"mir-4850", "tRNA"})
        self.assertEqual(len(filtered), 2)

    def test_filter_attr_by_value(self):
        filtered = self.tbl.filter_attr_by_value("query_name", "7osa_PSIT")
        self.assertEqual(len(filtered), 1)
        self.assertEqual(filtered.hits[0].target_name, "tRNA")

    def test_getitem(self):
        filtered = self.tbl["7lhd_A"]
        self.assertEqual(len(filtered), 2)
        self.assertTrue(all(q == "7lhd_A" for q in filtered.query_name))

    def test_reverse(self):
        rev = self.tbl.reverse
        self.assertEqual(len(rev), len(self.tbl))
        self.assertEqual(rev.hits, self.tbl.hits[::-1])

    def test_reverse_roundtrip(self):
        self.assertEqual(self.tbl.reverse.reverse.hits, self.tbl.hits)


class TestTableInit(unittest.TestCase):
    def test_init_raises_both_path_and_hits(self):
        with tempfile.NamedTemporaryFile("w") as f:
            f.write(TBL_STR)
            f.flush()
            with self.assertRaises(ValueError):
                Table(path=f.name, hits=[])

    def test_init_raises_neither(self):
        with self.assertRaises(ValueError):
            Table()

    def test_init_with_empty_hits(self):
        tbl = Table(hits=[])
        self.assertEqual(len(tbl), 0)


class TestTabularRead(unittest.TestCase):
    tbls_path = Path(__file__).parent / "test_data" / "tbls"

    def test_read_file(self):
        with tempfile.NamedTemporaryFile("w") as f:
            f.write(TBL_STR)
            f.flush()
            tbl = Table.read(f.name)
        self.assertEqual(len(tbl), 3)

    def test_read_directory(self):
        # cmscan.tbl has 11 hits, cmscan-nohits.tbl has 8 hits
        tbl = Table.read(self.tbls_path)
        self.assertEqual(len(tbl), 19)

    def test_read_directory_sorted_by_evalue(self):
        tbl = Table.read(self.tbls_path)
        e_values = tbl.e_value
        self.assertEqual(e_values, sorted(e_values))
