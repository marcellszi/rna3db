"""Automatic management of the CCD modifications cache for nucleic acid residues.

Downloads the Chemical Component Dictionary (CCD) from wwPDB, processes it to
extract nucleic acid (RNA/DNA) residue mappings, and caches the result locally.
"""

import gzip
import json
import logging
import os
import urllib.request
from pathlib import Path

logger = logging.getLogger(__name__)

CCD_URL = "https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz"
# follows the same cache convention as PyTorch
# see https://github.com/pytorch/pytorch/blob/main/torch/hub.py#L183-L190
_xdg_cache = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
CACHE_DIR = Path(os.environ.get("RNA3DB_CACHE_DIR", _xdg_cache / "rna3db"))
CACHE_FILE = "modifications_cache.json"
VALID_RNA_CODES = set("ACGUT")


def get_cache_path() -> Path:
    """Return the path to the cached modifications JSON."""
    return CACHE_DIR / CACHE_FILE


def _download_ccd(dest: Path) -> Path:
    """Download components.cif.gz from wwPDB."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"Downloading CCD from {CCD_URL} ...")
    urllib.request.urlretrieve(CCD_URL, dest)
    logger.info("Download complete.")
    return dest


def _parse_components_gz(gz_path: Path) -> dict:
    """Parse all components from a gzipped CCD file into {id: cif_text}."""
    cif_strings = {}
    chem_comp_id, cif_lines = None, []
    with gzip.open(gz_path, "rt") as fp:
        for line in fp:
            if line.startswith("data_"):
                if cif_lines:
                    cif_strings[chem_comp_id] = "".join(cif_lines)
                chem_comp_id = line.split("_")[-1].rstrip()
                cif_lines = []
            cif_lines.append(line)
        if cif_lines:
            cif_strings[chem_comp_id] = "".join(cif_lines)
    return cif_strings


def _parse_cif_fields(cif_string: str):
    """Extract key fields from a single CCD component's CIF text."""
    comp_id = comp_type = one_letter = parent_id = release_status = "?"
    for line in cif_string.split("\n"):
        if line.startswith("_chem_comp.id"):
            comp_id = line.split()[-1]
        elif line.startswith("_chem_comp.type"):
            comp_type = " ".join(line.split()[1:])
        elif line.startswith("_chem_comp.one_letter_code"):
            one_letter = line.split()[-1]
        elif line.startswith("_chem_comp.mon_nstd_parent_comp_id"):
            cleaned = "".join(c for c in line.split()[-1] if c.isalnum() or c == "?")
            parent_id = (cleaned or "?")[:3].upper()
        elif line.startswith("_chem_comp.pdbx_release_status"):
            release_status = line.split()[-1]
    return comp_id, comp_type, one_letter, parent_id, release_status


def _resolve_one_letter(comp_id: str, cif_strings: dict, _seen=None) -> tuple:
    """Resolve a component's one-letter code, following parent references."""
    if _seen is None:
        _seen = set()
    if comp_id in _seen:
        return "?", "?", "?"
    _seen.add(comp_id)

    if comp_id not in cif_strings:
        return "?", "?", "?"

    cid, comp_type, one_letter, parent_id, release_status = _parse_cif_fields(
        cif_strings[comp_id]
    )

    # follow parent reference
    if parent_id != "?" and cid != parent_id:
        return _resolve_one_letter(parent_id, cif_strings, _seen)

    # if one_letter_code is itself a multi-letter code, try resolving that
    if len(one_letter) > 1 and one_letter in cif_strings:
        return _resolve_one_letter(one_letter, cif_strings, _seen)

    return one_letter, comp_type, release_status


def _generate_cache(cif_strings: dict) -> dict:
    """Process parsed CCD components into an RNA-only modifications dict."""
    data = {}

    for comp_id in cif_strings:
        # check obsolete status on the component itself, before resolving parents
        _, _, _, _, release_status = _parse_cif_fields(cif_strings[comp_id])
        if release_status == "OBS":
            continue

        one_letter, comp_type, _ = _resolve_one_letter(comp_id, cif_strings)

        # only keep nucleic acid types with valid single-letter codes
        if ("RNA" in comp_type or "DNA" in comp_type) and one_letter in VALID_RNA_CODES:
            data[comp_id] = one_letter

    return data


def generate_from_ccd() -> dict:
    """Download the CCD and generate a fresh modifications cache.

    Returns:
        dict: Mapping of 3-letter CCD codes to 1-letter RNA/DNA codes.
    """
    gz_path = CACHE_DIR / "components.cif.gz"
    try:
        _download_ccd(gz_path)
        logger.info("Parsing CCD components ...")
        cif_strings = _parse_components_gz(gz_path)
        logger.info("Generating nucleic acid modifications cache ...")
        data = _generate_cache(cif_strings)
    finally:
        # always clean up the large download
        gz_path.unlink(missing_ok=True)

    # save to cache
    cache_path = get_cache_path()
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump(data, f, indent=4)
    logger.info(f"Modifications cache saved to {cache_path} ({len(data)} entries)")

    return data


def load(json_path=None) -> dict:
    """Load the modifications cache, generating from CCD if necessary.

    Args:
        json_path: Optional explicit path to a modifications JSON file.
            Supports both the new flat format and the old {"rna": {...}} format.

    Returns:
        dict: Mapping of 3-letter CCD codes to 1-letter RNA codes.
    """
    # use explicit path if provided
    if json_path is not None:
        with open(json_path) as f:
            data = json.load(f)
        # support old format with "rna" key
        if "rna" in data and isinstance(data["rna"], dict):
            return data["rna"]
        return data

    # check for cached version
    cache_path = get_cache_path()
    if cache_path.is_file():
        with open(cache_path) as f:
            return json.load(f)

    # generate from CCD
    return generate_from_ccd()
