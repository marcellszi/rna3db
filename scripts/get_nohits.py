from rna3db.tabular import read_tbls_from_dir
from rna3db.utils import read_json
from pathlib import Path

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extracts chains with no hits from a .tbl file."
    )
    parser.add_argument("input_path", type=Path)
    parser.add_argument("output_path", type=Path)
    parser.add_argument("tbls_path", type=Path)
    args = parser.parse_args()

    parse_json = read_json(args.input_path)
    all_chains = set(parse_json.keys())
    tbl = read_tbls_from_dir(args.tbls_path)
    hit_chains = set(tbl.query_name)

    with open(args.output_path, "w") as f:
        for nohit in all_chains - hit_chains:
            f.write(f'>{nohit}\n{parse_json[nohit]["sequence"]}\n')
