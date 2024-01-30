from rna3db.utils import read_json
from pathlib import Path

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Converts a JSON produced by RNA3DB's parse command to a FASTA file."
    )
    parser.add_argument("input_path", type=Path)
    parser.add_argument("output_path", type=Path)
    args = parser.parse_args()

    with open(args.output_path, "w") as f:
        for k, v in read_json(args.input_path).items():
            f.write(f'>{k}\n{v["sequence"]}\n')
