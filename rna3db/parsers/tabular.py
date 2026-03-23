from __future__ import annotations
from collections import namedtuple
from typing import Any, Sequence
from pathlib import Path


def read(path: Path) -> "TabularOutput":
    """Read a `.tbl` file or a directory of `.tbl` files.

    Args:
        path (Path): Path to a `.tbl` file or a directory containing `.tbl` files.

    Returns:
        TabularOutput: Object containing all hits, sorted by E-value.
    """
    path = Path(path)
    if path.is_dir():
        hits = []
        for p in path.glob("*.tbl"):
            hits.extend(TabularOutput(p).hits)
        return TabularOutput(hits=sorted(hits, key=lambda x: x.e_value))
    return TabularOutput(path)


class TabularOutput:
    TBL_ROW_TYPES = {
        "target_name": str,
        "target_accession": str,
        "query_name": str,
        "query_accession": lambda x: None if x == "-" else str(x),
        "mdl": str,
        "mdl_from": int,
        "mdl_to": int,
        "seq_from": int,
        "seq_to": int,
        "strand": str,
        "trunc": lambda x: True if x == "yes" else False,
        "pass_n": int,
        "gc": float,
        "bias": float,
        "score": float,
        "e_value": float,
        "inc": str,
        "description_of_target": str,
    }

    Hit = namedtuple("Hit", list(TBL_ROW_TYPES.keys()))

    def __init__(self, path: Path = None, hits: Sequence[Hit] = None):
        """
        Args:
            path (Path, optional): Path to a single ``.tbl`` file to parse.
            hits (Sequence[Hit], optional): Pre-parsed list of Hit namedtuples.

        Raises:
            ValueError: If both or neither of ``path`` and ``hits`` are provided.
        """
        if (path is None) == (hits is None):
            raise ValueError("Invalid values for path and/or hits.")
        if path is not None:
            self.hits = self._parse_tbl(path)
        if hits is not None:
            self.hits = hits

    def __getitem__(self, query: str) -> TabularOutput:
        return self.filter_attr_by_value("query_name", query)

    def __getattribute__(self, name: str):
        if name in TabularOutput.TBL_ROW_TYPES:
            setattr(self, name, [getattr(hit, name) for hit in self.hits])
        return super().__getattribute__(name)

    def __len__(self) -> int:
        return len(self.hits)

    def __iter__(self):
        return iter(self.hits)

    def __repr__(self) -> str:
        col_width = 20
        print_cols = [
            "target_name",
            "target_accession",
            "query_name",
            "score",
            "e_value",
        ]
        max_rows = 15
        s = ""
        for col in print_cols:
            if len(col) > col_width:
                col = col[: col_width - 3] + "..."
            s += f"{col[:col_width]:{col_width}}"
        s += "\n" + "-" * len(s) + "\n"
        for hit in self.hits[:max_rows]:
            for col in print_cols:
                c = str(getattr(hit, col))
                if len(c) > col_width:
                    c = c[: col_width - 3] + "..."
                s += f"{c:{col_width}}"
            s += "\n"
        if len(self.hits) > max_rows:
            s += f"... ({len(self.hits)-max_rows} rows hidden)\n"
        return s

    @property
    def reverse(self) -> "TabularOutput":
        """Return a TabularOutput with hits in reverse order."""
        return TabularOutput(hits=self.hits[::-1])

    @property
    def top_hits(self) -> "TabularOutput":
        """Get the top hit (lowest E-value) for each query in the table.

        Note:
            If multiple hits share the same E-value, only the first encountered
            is kept. This commonly occurs when E-value == 0.0.

        Returns:
            TabularOutput: One hit per query, unsorted.
        """
        # this is ugly, but O(n)
        th = {}
        for hit in self.hits:
            if hit.query_name not in th:
                th[hit.query_name] = hit
                continue
            if hit.e_value < th[hit.query_name].e_value:
                th[hit.query_name] = hit
        return TabularOutput(hits=list(th.values()))

    def filter_e_value(self, cutoff: float) -> TabularOutput:
        """Filter hits to those with E-value <= cutoff.

        Args:
            cutoff (float): Maximum E-value to retain.

        Returns:
            TabularOutput: Filtered hits, sorted by E-value.
        """
        hits = []
        for hit in self.hits:
            if hit.e_value <= cutoff:
                hits.append(hit)
        return TabularOutput(hits=sorted(hits, key=lambda x: x.e_value))

    def filter_attr_by_set(self, attr: str, filter_set: Sequence[str]) -> TabularOutput:
        """Filter hits to those whose attribute value is in a given set.

        Args:
            attr (str): The hit attribute to filter by. Must be a key of
                ``TabularOutput.TBL_ROW_TYPES``.
            filter_set (Sequence[str]): Any object supporting ``__contains__``
                to filter by.

        Returns:
            TabularOutput: Filtered hits, sorted by E-value.

        Examples:
            >>> print(tbl.target_name)
            ['5S_rRNA', 'tRNA5', 'tRNA5', 'Cobalamin']
            >>> result = tbl.filter_attr_by_set('target_name', ['5S_rRNA', 'tRNA5'])
            >>> print(result.target_name)
            ['5S_rRNA', 'tRNA5', 'tRNA5']
        """
        hits = [hit for hit in self.hits if getattr(hit, attr) in filter_set]
        return TabularOutput(hits=sorted(hits, key=lambda x: x.e_value))

    def filter_attr_by_value(self, attr: str, val: Any) -> TabularOutput:
        """Filter hits to those whose attribute matches a single value.

        Alias for ``filter_attr_by_set(attr, [val])``.

        Args:
            attr (str): The hit attribute to filter by. Must be a key of
                ``TabularOutput.TBL_ROW_TYPES``.
            val: Value to filter by.

        Returns:
            TabularOutput: Filtered hits, sorted by E-value.

        Examples:
            >>> print(tbl.target_name)
            ['5S_rRNA', 'tRNA5', 'tRNA5', 'Cobalamin']
            >>> print(tbl.filter_attr_by_value('target_name', '5S_rRNA').target_name)
            ['5S_rRNA']
        """
        return self.filter_attr_by_set(attr, [val])

    @staticmethod
    def _parse_tbl_row(s):
        row = s.split()

        for i, field in enumerate(TabularOutput.Hit._fields):
            row[i] = TabularOutput.TBL_ROW_TYPES[field](row[i])

        # handle spaces in the last column
        row[i] = " ".join(row[i:])
        del row[i + 1 :]

        return TabularOutput.Hit(*row)

    def _parse_tbl(self, path):
        entries = []
        with open(path) as f:
            for line in f:
                if line[0] == "#":
                    continue
                if len(line.split()) == 0:
                    continue  # prob unnecessary, but w/e
                entries.append(self._parse_tbl_row(line))
        return entries
