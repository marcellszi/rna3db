#!/usr/bin/env python
# -*- coding: utf-8 -*-
from rna3db.utils import read_json
from rna3db.parser import parse_file
from pathlib import Path
from tqdm import tqdm

import argparse


def main(args):
    data = read_json(args.input_json_path)
    num_chains = sum(
        len(i.keys()) for s in data.values() for c in s.values() for i in c.values()
    )

    pbar = tqdm(total=num_chains)

    for set_name, set_content in data.items():
        for component_name, component_content in set_content.items():
            for cluster_name, cluster_content in component_content.items():
                cluster_path = (
                    args.output_dir / set_name / component_name / cluster_name
                )
                cluster_path.mkdir(exist_ok=True, parents=True)
                for chain_name in cluster_content.keys():
                    pdb_id, author_id = chain_name.split("_")
                    output_path = cluster_path / f"{chain_name}.cif"
                    pdb_mmcif_path = args.input_mmcif_dir / f"{pdb_id}.cif"
                    if not pdb_mmcif_path.is_file():
                        print(f"WARNING: could not find {pdb_mmcif_path}")
                        continue
                    if not output_path.is_file() or not args.skip_existing:
                        sf = parse_file(pdb_mmcif_path, include_atoms=True)
                        sf.write_mmcif_chain(output_path, author_id)
                    pbar.update(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Converts a JSON produced by RNA3DB's parse command to a set of PDBx/mmCIF files."
    )
    parser.add_argument("input_json_path", type=Path)
    parser.add_argument("input_mmcif_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="don't do write again if the output file already exists",
    )

    args = parser.parse_args()

    main(args)
