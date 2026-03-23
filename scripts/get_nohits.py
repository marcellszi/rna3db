from rna3db.parsers import tabular, fasta
from pathlib import Path

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extracts chains with no hits from a .tbl file."
    )
    parser.add_argument("input_path", type=Path)
    parser.add_argument("output_path", type=Path)
    parser.add_argument("tbls_path", type=Path)
    parser.add_argument("--e_value_threshold", type=float, default=1.0)
    parser.add_argument("--length_threshold", type=int, default=64)
    args = parser.parse_args()

    all_chains, all_sequences = fasta.read(args.input_path)
    tbl = tabular.read(args.tbls_path)

    all_hits = set(tbl.query_name)
    edge_hits = set(tbl.filter_e_value(args.e_value_threshold).query_name)
    short_chains = set(
        {k for k, v in zip(all_chains, all_sequences) if len(v) < args.length_threshold}
    )

    # chains that were not hit at all
    zero_hits = set(all_chains) - all_hits

    # chains that were not hit below the e-value threshold and are shorter than
    # the length threshold
    short_bad_hits = (set(all_chains) - edge_hits) & short_chains

    # not hit at all or only bad hits
    nohits = zero_hits | short_bad_hits

    all_dict = dict(zip(all_chains, all_sequences))
    output_descriptions = list(nohits)
    output_sequences = [all_dict[i] for i in nohits]

    fasta.write(output_descriptions, output_sequences, args.output_path)
