import argparse
from pathlib import Path

from rna3db.parsers import fasta
from rna3db.utils import read_json

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract only new sequences from two parse outputs, and write them to a FASTA. "
    )
    parser.add_argument("old_path", type=Path)
    parser.add_argument("new_path", type=Path)
    parser.add_argument("output_path", type=Path)
    args = parser.parse_args()

    old_parse = read_json(args.old_path)
    new_parse = read_json(args.new_path)

    records = [
        fasta.Record(header=k, sequence=new_parse[k]["sequence"])
        for k in set(new_parse.keys()) - set(old_parse.keys())
    ]

    fasta.write(records, args.output_path)
