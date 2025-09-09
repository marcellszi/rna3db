import json
import os
import re

from typing import Union

PathLike = Union[str, os.PathLike]


def to_case_insensitive(s: str) -> str:
    return re.sub(r"([a-z]+)", r"-\1-", s)


def to_case_sensitive(s: str) -> str:
    out = []
    for i, part in enumerate(s.split("-")):
        if i % 2 == 0:
            out.append(part.upper())
        else:
            out.append(part.lower())
    return "".join(out)


def read_json(input_path: PathLike) -> dict:
    """Read a JSON file to a Python dictionary.

    Args:
        input_path (PathLike): path from which the JSON is read

    Returns:
        Python dict of read JSON
    """
    with open(input_path) as f:
        return json.load(f)


def write_json(data: dict, output_path: PathLike):
    """Write a Python dictionary to a JSON file.

    Args:
        data (dict): dictionary to write
        output_path (PathLike): path to write the JSON to
    """
    with open(output_path, "w") as f:
        json.dump(data, f, indent=4)
