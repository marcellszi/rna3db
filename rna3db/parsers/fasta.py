from pathlib import Path
from typing import NamedTuple


class Record(NamedTuple):
    header: str
    sequence: str


def read(path: Path, force_gzip: bool = False) -> list[Record]:
    """Parse a Record file.

    Supports multi-line sequences.

    Args:
        path (Path): Path to input Record file.
        force_gzip (bool, optional): If True, will attempt to read the file as a
            gzip file.

    Returns:
        list[Record]: List of Record records.
    """
    if Path(path).suffix == ".gz" or force_gzip:
        import gzip

        reader = gzip.open(path, "rt")
    else:
        reader = open(path, "r")
    with reader as f:
        records = []
        current_header = None
        current_sequence = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header is not None:
                    records.append(Record(current_header, current_sequence))
                current_header = line[1:]
                current_sequence = ""
            elif line.startswith("#") or not line:
                continue
            else:
                current_sequence += line
        if current_header is not None:
            records.append(Record(current_header, current_sequence))
    return records


def write(records: Record[Record], output_path: Path):
    """Write Record records to a file.

    Args:
        records (Record[Record]): Record records to write.
        output_path (Path): Path to write Record file to.
    """
    with open(output_path, "w") as f:
        for record in records:
            f.write(f">{record.header}\n{record.sequence}\n")
