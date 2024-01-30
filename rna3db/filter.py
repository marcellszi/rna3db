import logging
import json

from rna3db.utils import PathLike, write_json


class Filterer:
    def __init__(
        self,
        min_length: int = 32,
        max_resolution: float = 9.0,
        single_ratio_cutoff: float = 0.8,
        max_unknown_ratio: float = 0.3,
    ):
        self.min_length = min_length
        self.max_resolution = max_resolution
        self.single_ratio_cutoff = single_ratio_cutoff
        self.max_unknown_ratio = max_unknown_ratio

        self.filters = []
        if self.min_length:
            self.filters.append(self.is_short_sequence)
        if self.max_resolution:
            self.filters.append(self.is_low_resolution)
        if self.single_ratio_cutoff:
            self.filters.append(self.is_singleratio_sequence)
        if self.max_unknown_ratio:
            self.filters.append(self.sequence_has_many_unknowns)

    def is_low_resolution(self, d: dict):
        return d["resolution"] > self.max_resolution

    def is_short_sequence(self, d: dict):
        return len(d["sequence"]) < self.min_length

    def is_singleratio_sequence(self, d: dict):
        l = len(d["sequence"])
        for nt in set(d["sequence"]):
            if d["sequence"].count(nt) / l > self.single_ratio_cutoff:
                return True
        return False

    def sequence_has_many_unknowns(self, d: dict):
        ratio = d["sequence"].count("N") / len(d["sequence"])
        return ratio > self.max_unknown_ratio

    def apply_filters(self, data: dict, json_filter_log_path: PathLike = None):
        logging.info(f"Applying filters {[f.__name__ for f in self.filters]}")
        filtered_data = {}
        applied_filters = {}

        for iid, d in data.items():
            conditions = [f(d) for f in self.filters]
            if not any(conditions):
                filtered_data[iid] = d.copy()
            applied_filters[iid] = [
                self.filters[i].__name__ for i, b in enumerate(conditions) if b
            ]

        if json_filter_log_path is not None:
            with open(json_filter_log_path, "w") as f:
                json.dump(applied_filters, f, indent=4)

        return filtered_data
