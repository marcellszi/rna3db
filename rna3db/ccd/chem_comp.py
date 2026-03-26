"""Cache of standard RNA nucleotide _chem_comp fields for mmCIF output.

Extracts the CCD entries for A, C, G, U from the Chemical Component Dictionary
and caches them locally. Used by the mmCIF writer to produce valid output files.
"""

import json
import logging
from pathlib import Path

from rna3db.ccd import CACHE_DIR, download_ccd, parse_cif_fields, parse_components_gz

logger = logging.getLogger(__name__)

CACHE_FILE = "chem_comp_cache.json"
STANDARD_IDS = {"A", "C", "G", "U"}


def get_cache_path() -> Path:
    """Return the path to the cached chem_comp JSON."""
    return CACHE_DIR / CACHE_FILE


def _generate_cache(cif_strings: dict) -> list:
    """Extract _chem_comp fields for standard RNA nucleotides from parsed CCD data.

    Returns:
        list: List of dicts with id, name, formula, weight for each standard nucleotide.
    """
    data = []
    for comp_id in sorted(STANDARD_IDS):
        if comp_id not in cif_strings:
            continue
        cid, _, _, _, _, name, formula, weight = parse_cif_fields(cif_strings[comp_id])
        data.append(
            {
                "id": cid,
                "name": name,
                "formula": formula,
                "weight": weight,
            }
        )

    # N as a catch-all for unknown nucleotides (not a real CCD component)
    data.append({"id": "N", "name": "N", "formula": "?", "weight": "?"})

    return data


def generate_from_ccd(cif_strings: dict = None) -> list:
    """Generate the chem_comp cache, optionally reusing already-parsed CCD data.

    Args:
        cif_strings: Pre-parsed CCD components (from ccd.parse_components_gz).
            If None, downloads and parses the CCD from scratch.

    Returns:
        list: Standard RNA nucleotide chem_comp entries.
    """
    if cif_strings is None:
        gz_path = download_ccd()
        try:
            cif_strings = parse_components_gz(gz_path)
        finally:
            gz_path.unlink(missing_ok=True)

    data = _generate_cache(cif_strings)

    cache_path = get_cache_path()
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump(data, f, indent=4)
    logger.info(f"Chem comp cache saved to {cache_path}")

    return data


def load(cif_strings: dict = None) -> list:
    """Load the chem_comp cache, generating from CCD if necessary.

    Args:
        cif_strings: Pre-parsed CCD components to avoid re-downloading.

    Returns:
        list: List of dicts with id, name, formula, weight.
    """
    cache_path = get_cache_path()
    if cache_path.is_file():
        with open(cache_path) as f:
            return json.load(f)

    return generate_from_ccd(cif_strings)
