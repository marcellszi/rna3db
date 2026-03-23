from pathlib import Path
from typing import Sequence


def read(path, force_gzip=False):
    """Parse a FASTA file.

    Supports multi-line sequences.

    Args:
        path (Path): Path to input FASTA file.
        force_gzip (bool, optional): If True, will attempt to read the file as a
        gzip file.

    Returns:

    """
    if Path(path).suffix == ".gz" or force_gzip:
        import gzip

        reader = gzip.open(path, "rt")
    else:
        reader = open(path, "r")
    with reader as f:
        descriptions = []
        sequences = []
        i = -1
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                i += 1
                descriptions.append(line[1:])
                sequences.append("")
            elif line.startswith("#") or not line:
                continue
            else:
                sequences[i] += line
    return descriptions, sequences


def write(descriptions: Sequence[str], sequences: Sequence[str], output_path: Path):
    """Write to a FASTA file.

    Args:
        descriptions (Sequence): List of descriptions for each sequence.
        sequences (Sequence): List of sequences.
        output_path (Path): Path to write FASTA file to.
    """
    if len(descriptions) != len(sequences):
        raise ValueError("The length of descriptions and sequences must match.")
    with open(output_path, "w") as f:
        for k, v in zip(descriptions, sequences):
            f.write(f">{k}\n{v}\n")
