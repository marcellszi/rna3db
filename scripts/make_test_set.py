import argparse
from pathlib import Path

from rna3db.tabular import read_tbls_from_dir
from rna3db.utils import read_json

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("output_path", type=Path)
    parser.add_argument("date_cutoff", type=str)
    parser.add_argument("cluster_json", type=Path)

    args = parser.parse_args()

    cluster_json = read_json(args.cluster_json)

    valid_components = set()
    num_chains = 0
    for component_k, component_v in cluster_json.items():
        curr_dates = []

        for repr_k, repr_v in component_v.items():
            for chain_k, chain_v in repr_v.items():
                curr_dates.append(chain_v["release_date"])

        if min(curr_dates) > args.date_cutoff:
            valid_components.add(component_k)
            num_chains += len(curr_dates)

    print("\n".join(sorted(valid_components)))
    print("num_chains:", num_chains)
