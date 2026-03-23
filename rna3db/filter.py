import logging
import json

from pathlib import Path


def _is_low_resolution(d: dict, cutoff: float) -> bool:
    return d["resolution"] > cutoff


def _is_short_sequence(d: dict, min_length: int) -> bool:
    return len(d["sequence"]) < min_length


def _is_singleratio_sequence(d: dict, cutoff: float) -> bool:
    l = len(d["sequence"])
    for nt in set(d["sequence"]):
        if d["sequence"].count(nt) / l > cutoff:
            return True
    return False


def _has_many_unknowns(d: dict, cutoff: float) -> bool:
    return d["sequence"].count("N") / len(d["sequence"]) > cutoff


def apply_filters(
    data: dict,
    min_length: int = 32,
    max_resolution: float = 9.0,
    single_ratio_cutoff: float = 0.8,
    max_unknown_ratio: float = 0.3,
    filter_log_path: Path = None,
) -> dict:
    """Filter a parsed JSON of RNA chains by quality criteria.

    Each filter is only applied if its threshold is truthy (non-zero, non-None).

    Args:
        data (dict): Parsed chain data, keyed by chain ID.
        min_length (int): Remove chains shorter than this.
        max_resolution (float): Remove chains with resolution above this.
        single_ratio_cutoff (float): Remove chains where any single
            nucleotide makes up more than this fraction of residues.
        max_unknown_ratio (float): Remove chains with more than this
            fraction of unknown nucleotides (N).
        filter_log_path (Path, optional): If provided, write a JSON log
            mapping each chain ID to the filters that hit it.

    Returns:
        dict: Filtered chain data, keyed by chain ID.
    """
    # Build (name, predicate, threshold) tuples for each active filter.
    # A filter is skipped if its threshold is falsy (0 or None).
    active = []
    if min_length:
        active.append(("is_short_sequence", _is_short_sequence, min_length))
    if max_resolution:
        active.append(("is_low_resolution", _is_low_resolution, max_resolution))
    if single_ratio_cutoff:
        active.append(
            ("is_singleratio_sequence", _is_singleratio_sequence, single_ratio_cutoff)
        )
    if max_unknown_ratio:
        active.append(("has_many_unknowns", _has_many_unknowns, max_unknown_ratio))

    logging.info(f"Applying filters {[name for name, _, _ in active]}")

    filtered_data = {}
    applied_filters = {}

    for iid, d in data.items():
        hits = [name for name, pred, threshold in active if pred(d, threshold)]
        if not hits:
            filtered_data[iid] = d.copy()
        applied_filters[iid] = hits

    if filter_log_path is not None:
        with open(filter_log_path, "w") as f:
            json.dump(applied_filters, f, indent=4)

    return filtered_data
