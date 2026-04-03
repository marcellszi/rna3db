import argparse
from collections import defaultdict
from pathlib import Path

from rna3db.parsers import fasta
from rna3db.utils import write_json

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Converts a FASTA to an RNA3DB-style JSON."
    )
    parser.add_argument("input_path", type=Path)
    parser.add_argument("output_path", type=Path)
    args = parser.parse_args()

    records = fasta.FASTA.read(args.input_path)
    data = defaultdict(dict)
    for header, sequence in records:
        data[header]["release_date"] = "1970-01-01"
        data[header]["structure_method"] = ""
        data[header]["resolution"] = 0.0
        data[header]["length"] = len(sequence)
        data[header]["sequence"] = sequence

    write_json(data, args.output_path)
