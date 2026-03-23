"""Regenerate the modifications cache from the Chemical Component Dictionary.

Usage:
    python scripts/generate_modifications_cache.py

Downloads the CCD from wwPDB and regenerates ~/.cache/rna3db/modifications_cache.json.
Set RNA3DB_CACHE_DIR to override the cache location.
"""

import logging
import sys

sys.path.append(".")

from rna3db.ccd.modifications import generate_from_ccd, get_cache_path

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    data = generate_from_ccd()
    print(f"Generated {len(data)} nucleic acid entries at {get_cache_path()}")
