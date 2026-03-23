"""Shared utilities for downloading and parsing the Chemical Component Dictionary (CCD).

The CCD is maintained by wwPDB at https://www.wwpdb.org/data/ccd.
"""

import gzip
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


def download_ccd() -> Path:
    """Download components.cif.gz from wwPDB into the cache directory.

    Returns:
        Path: Path to the downloaded gzipped file.
    """
    gz_path = CACHE_DIR / "components.cif.gz"
    gz_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"Downloading CCD from {CCD_URL} ...")
    urllib.request.urlretrieve(CCD_URL, gz_path)
    logger.info("Download complete.")
    return gz_path


def parse_components_gz(gz_path: Path) -> dict:
    """Parse all components from a gzipped CCD file into {id: cif_text}.

    Args:
        gz_path: Path to components.cif.gz.

    Returns:
        dict: Mapping of component ID to raw CIF text block.
    """
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


def parse_cif_fields(cif_string: str):
    """Extract key fields from a single CCD component's CIF text.

    Returns:
        tuple: (comp_id, comp_type, one_letter, parent_id, release_status,
                name, formula, formula_weight)
    """
    comp_id = comp_type = one_letter = parent_id = release_status = "?"
    name = formula = "?"
    formula_weight = "?"
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
        elif line.startswith("_chem_comp.name"):
            name = " ".join(line.split()[1:]).strip("'\"")
        elif line.startswith("_chem_comp.formula "):
            formula = " ".join(line.split()[1:]).strip("'\"")
        elif line.startswith("_chem_comp.formula_weight"):
            try:
                formula_weight = float(line.split()[-1])
            except ValueError:
                formula_weight = "?"
    return (
        comp_id,
        comp_type,
        one_letter,
        parent_id,
        release_status,
        name,
        formula,
        formula_weight,
    )
