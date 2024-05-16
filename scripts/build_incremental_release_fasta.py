from rna3db.utils import read_json
from rna3db.parser import write_fasta

from collections import defaultdict
from pathlib import Path

import argparse

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

    descriptions, sequences = [], []
    for k in set(new_parse.keys()) - set(old_parse.keys()):
        descriptions.append(k)
        sequences.append(new_parse[k]["sequence"])

    write_fasta(descriptions, sequences, args.output_path)
