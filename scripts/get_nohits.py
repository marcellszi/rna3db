import argparse
from pathlib import Path

from rna3db.parsers import fasta, tabular

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

    all_fasta = fasta.FASTA.read(args.input_path)
    tbl = tabular.Table.read(args.tbls_path)

    all_hits = set(tbl.query_name)
    edge_hits = set(tbl.filter_e_value(args.e_value_threshold).query_name)
    short_chains = {
        h for h, s in all_fasta if len(s) < args.length_threshold
    }

    # chains that were not hit at all
    all_headers = {h for h, _ in all_fasta}
    zero_hits = all_headers - all_hits

    # chains that were not hit below the e-value threshold and are shorter than
    # the length threshold
    short_bad_hits = (all_headers - edge_hits) & short_chains

    # not hit at all or only bad hits
    nohits = zero_hits | short_bad_hits

    all_dict = {h: s for h, s in all_fasta}
    output_headers = [h for h in nohits]
    output_sequences = [all_dict[h] for h in output_headers]

    fasta.FASTA(output_headers, output_sequences).write(args.output_path)
