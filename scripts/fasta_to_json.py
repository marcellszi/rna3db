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

    records = fasta.read(args.input_path)
    data = defaultdict(dict)
    for r in records:
        data[r.header]["release_date"] = "1970-01-01"
        data[r.header]["structure_method"] = ""
        data[r.header]["resolution"] = 0.0
        data[r.header]["length"] = len(r.sequence)
        data[r.header]["sequence"] = r.sequence

    write_json(data, args.output_path)
