from rna3db.tabular import read_tbls_from_dir
from rna3db.parser import parse_fasta, write_fasta
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

    all_chains, all_sequences = parse_fasta(args.input_path)
    all_dict = {k: v for k, v in zip(all_chains, all_sequences)}
    all_chains = set(all_chains)

    tbl = read_tbls_from_dir(args.tbls_path)
    hit_chains = set(tbl.query_name)

    nohits = all_chains - hit_chains

    output_descriptions = list(nohits)
    output_sequences = [all_dict[i] for i in nohits]

    write_fasta(output_descriptions, output_sequences, args.output_path)
