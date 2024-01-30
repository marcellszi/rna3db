import contextlib
import unittest
import tempfile

from collections import defaultdict

from rna3db.tabular import TabularOutput

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
            self.tbl = TabularOutput(f.name)

    @contextlib.contextmanager
    def tmp_txt(self, s: str):
        tmp = tempfile.NamedTemporaryFile("w")
        tmp.write(TBL_STR)
        tmp.seek(0)
        try:
            yield tmp
        finally:
            tmp.close()

    def _assert_hit(self, hit: TabularOutput.Hit, test_dict: dict):
        for k, v in test_dict.items():
            self.assertEqual(hit.__getattribute__(k), v)

    def test_row_parse(self):
        for row_str, test_dict in zip(TBL_STR.split("\n")[2:], TEST_DICTS):
            hit = TabularOutput._parse_tbl_row(row_str)
            self._assert_hit(hit, test_dict)

    def test_parse_tbl(self):
        tbl = TabularOutput(hits=[])
        with self.tmp_txt(TBL_STR) as f:
            entries = tbl._parse_tbl(f.name)
            for hit, test_dict in zip(entries, TEST_DICTS):
                self._assert_hit(hit, test_dict)

    def test_len(self):
        empty_hit = TabularOutput.Hit(*[None] * 18)
        self.assertEqual(len(TabularOutput(hits=[])), 0)
        self.assertEqual(len(TabularOutput(hits=[empty_hit] * 1)), 1)
        self.assertEqual(len(TabularOutput(hits=[empty_hit] * 1337)), 1337)

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
