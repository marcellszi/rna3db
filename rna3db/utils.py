import json
import re
from pathlib import Path


def to_case_insensitive(s: str) -> str:
    """Encode a string so its case can be recovered from a case-insensitive
    representation.

    Wraps each run of lowercase letters with hyphens, leaving uppercase letters
    and non-alpha characters unchanged. The result can be safely uppercased or
    lowercased (e.g. by a case-insensitive filesystem) and still round-trip
    through ``to_case_sensitive``.

    This is used to encode case-sensitive PDB chain IDs into filenames that
    survive on case-insensitive filesystems (e.g. macOS HFS+, Windows NTFS).

    Args:
        s (str): Input string, typically a PDB chain ID.

    Returns:
        str: Encoded string with lowercase runs wrapped in hyphens.

    Examples:
        >>> to_case_insensitive("A")
        'A'
        >>> to_case_insensitive("a")
        '-a-'
        >>> to_case_insensitive("Aa")
        'A-a-'
    """
    return re.sub(r"([a-z]+)", r"-\1-", s)


def to_case_sensitive(s: str) -> str:
    """Decode a string produced by ``to_case_insensitive``, recovering the
    original case.

    Splits on hyphens; even-indexed segments are uppercased and odd-indexed
    segments (originally lowercase runs) are lowercased. Works correctly
    regardless of whether the encoded string has been uppercased or lowercased
    in the interim.

    Args:
        s (str): Encoded string, as produced by ``to_case_insensitive``.

    Returns:
        str: Decoded string with original case restored.

    Examples:
        >>> to_case_sensitive("A")
        'A'
        >>> to_case_sensitive("-a-")
        'a'
        >>> to_case_sensitive("A-a-")
        'Aa'
    """
    out = []
    for i, part in enumerate(s.split("-")):
        if i % 2 == 0:
            out.append(part.upper())
        else:
            out.append(part.lower())
    return "".join(out)


def read_json(input_path: Path) -> dict:
    """Read a JSON file to a Python dictionary.

    Args:
        input_path (Path): path from which the JSON is read

    Returns:
        Python dict of read JSON
    """
    with open(input_path) as f:
        return json.load(f)


def write_json(data: dict, output_path: Path):
    """Write a Python dictionary to a JSON file.

    Args:
        data (dict): dictionary to write
        output_path (Path): path to write the JSON to
    """
    with open(output_path, "w") as f:
        json.dump(data, f, indent=4)
