"""Automatic management of the CCD modifications cache for nucleic acid residues.

Downloads the Chemical Component Dictionary (CCD) from wwPDB, processes it to
extract nucleic acid (RNA/DNA) residue mappings, and caches the result locally.
"""

import json
import logging
from pathlib import Path

from rna3db.ccd import CACHE_DIR, download_ccd, parse_cif_fields, parse_components_gz

logger = logging.getLogger(__name__)

CACHE_FILE = "modifications_cache.json"
VALID_RNA_CODES = set("ACGUT")


def get_cache_path() -> Path:
    """Return the path to the cached modifications JSON."""
    return CACHE_DIR / CACHE_FILE


def _resolve_one_letter(comp_id: str, cif_strings: dict, _seen=None) -> tuple:
    """Resolve a component's one-letter code, following parent references."""
    if _seen is None:
        _seen = set()
    if comp_id in _seen:
        return "?", "?"
    _seen.add(comp_id)

    if comp_id not in cif_strings:
        return "?", "?"

    cid, comp_type, one_letter, parent_id, *_ = parse_cif_fields(cif_strings[comp_id])

    # follow parent reference
    if parent_id != "?" and cid != parent_id:
        return _resolve_one_letter(parent_id, cif_strings, _seen)

    # if one_letter_code is itself a multi-letter code, try resolving that
    if len(one_letter) > 1 and one_letter in cif_strings:
        return _resolve_one_letter(one_letter, cif_strings, _seen)

    return one_letter, comp_type


def _generate_cache(cif_strings: dict) -> dict:
    """Process parsed CCD components into an RNA-only modifications dict."""
    data = {}

    for comp_id in cif_strings:
        # check obsolete status on the component itself, before resolving parents
        _, _, _, _, release_status, _, _, _ = parse_cif_fields(cif_strings[comp_id])
        if release_status == "OBS":
            continue

        one_letter, comp_type = _resolve_one_letter(comp_id, cif_strings)

        # only keep nucleic acid types with valid single-letter codes
        if ("RNA" in comp_type or "DNA" in comp_type) and one_letter in VALID_RNA_CODES:
            data[comp_id] = one_letter

    return data


def generate_from_ccd(cif_strings: dict = None) -> dict:
    """Generate the modifications cache, optionally reusing already-parsed CCD data.

    Args:
        cif_strings: Pre-parsed CCD components (from ccd.parse_components_gz).
            If None, downloads and parses the CCD from scratch.

    Returns:
        dict: Mapping of 3-letter CCD codes to 1-letter RNA/DNA codes.
    """
    if cif_strings is None:
        gz_path = download_ccd()
        try:
            logger.info("Parsing CCD components ...")
            cif_strings = parse_components_gz(gz_path)
        finally:
            gz_path.unlink(missing_ok=True)

    logger.info("Generating nucleic acid modifications cache ...")
    data = _generate_cache(cif_strings)

    # save to cache
    cache_path = get_cache_path()
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump(data, f, indent=4)
    logger.info(f"Modifications cache saved to {cache_path} ({len(data)} entries)")

    return data


class ModificationHandler:
    def __init__(self, json_path=None):
        """Used for converting `three_letter_code`s to `one_letter_code`s, including modifications.

        On first use, automatically downloads and processes the Chemical Component
        Dictionary (CCD) from wwPDB, caching the result in ``~/.cache/rna3db/``.
        Set the ``RNA3DB_CACHE_DIR`` environment variable to override the cache location.

        Args:
            json_path (Path, optional): Explicit path to a modifications cache JSON file.
                If not provided, the cache is loaded (or generated) automatically.
        """
        self.modifications = load(Path(json_path) if json_path else None)

    def is_rna(self, three_letter_code: str) -> bool:
        """Check if `three_letter_code` is a known RNA/DNA nucleic acid residue.

        Args:
            three_letter_code (str): Three letter code to check.

        Returns:
            bool: True if `three_letter_code` is a known nucleic acid residue.
        """
        return three_letter_code in self.modifications

    def rna_letters_3to1(self, three_letter_code: str) -> str:
        """Convert RNA nucleic acid `three_letter_code` to `one_letter_code`.

        Args:
            three_letter_code (str): Three letter code to check.

        Returns:
           str: one_letter_code of RNA nucleic acid, "N" if cannot be found.
        """
        return self.modifications.get(three_letter_code, "N")


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
