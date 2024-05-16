from multiprocessing import Pool
from functools import partial
from pathlib import Path
from tqdm import tqdm

from rna3db.parser import parse_as_dict
from rna3db.filter import Filterer
from rna3db.cluster import SequenceClusterer, StructureClusterer
from rna3db.split import split
from rna3db.utils import PathLike, read_json, write_json

import argparse


def parse(
    input_dir: PathLike,
    output_path: PathLike,
    extension: str = "cif",
    cpu: int = None,
    nmr_resolution: float = None,
    include_atoms: bool = False,
):
    files = list(Path(input_dir).glob(f"*.{extension}"))
    data = {}

    f = partial(
        parse_as_dict, nmr_resolution=nmr_resolution, include_atoms=include_atoms
    )
    with Pool(processes=cpu) as p:
        with tqdm(total=len(files)) as pbar:
            for d in p.imap_unordered(f, files):
                data.update(d)
                pbar.update()

    write_json(data, output_path)


def filter(
    input_path: PathLike,
    output_path: PathLike,
    min_length: int = 32,
    max_resolution: float = 9.0,
    single_ratio_cutoff: float = 0.8,
    max_unknown_ratio: float = 0.3,
):
    data = read_json(input_path)

    filterer = Filterer(
        min_length=min_length,
        max_resolution=max_resolution,
        single_ratio_cutoff=single_ratio_cutoff,
        max_unknown_ratio=max_unknown_ratio,
    )

    filtered_data = filterer.apply_filters(
        data, json_filter_log_path=args.filter_log_path
    )

    write_json(filtered_data, output_path)


def cluster_sequence(
    input_path: PathLike,
    output_path: PathLike,
    mmseqs_binary_path: PathLike = None,
    min_seq_id: float = 0.99,
    min_coverage: float = 0.99,
    coverage_mode: int = 1,
    sensitivity: float = 7.5,
    alignment_mode: int = 3,
    max_seqs: int = 10000,
):
    seq_cluster = SequenceClusterer(
        mmseqs_binary_path,
        min_seq_id=min_seq_id,
        min_coverage=min_coverage,
        coverage_mode=coverage_mode,
        sensitivity=sensitivity,
        alignment_mode=alignment_mode,
        max_seqs=max_seqs,
    )

    cluster = seq_cluster.cluster(input_path, output_path)
    write_json(cluster, output_path)


def cluster_structure(
    input_path: PathLike,
    output_path: PathLike,
    tbl_dir: PathLike,
    e_value_cutoff: float = 1,
):
    str_cluster = StructureClusterer(e_value_cutoff)
    cluster = str_cluster.cluster(input_path, tbl_dir)
    write_json(cluster, output_path)


def main(args):
    if args.command == "parse":
        parse(
            args.input,
            args.output,
            cpu=args.cpu,
            nmr_resolution=args.nmr_resolution,
            include_atoms=args.include_atoms,
        )
    elif args.command == "filter":
        filter(
            args.input,
            args.output,
            args.min_length,
            args.max_resolution,
            args.single_ratio_cutoff,
            args.max_unknown_ratio,
        )
    elif args.command == "cluster":
        if not args.only_structure:
            cluster_sequence(
                args.input,
                args.output,
                args.mmseqs_binary_path,
                args.min_seq_id,
                args.min_seq_coverage,
                args.mmseqs_coverage_mode,
                args.mmseqs_sensitivity,
                args.mmseqs_alignment_mode,
                args.mmseqs_max_seqs,
            )
            args.input = args.output
        if not args.only_sequence:
            cluster_structure(
                args.input, args.output, args.tbl_dir, args.structural_e_value_cutoff
            )
    elif args.command == "split":
        split(
            args.input,
            args.output,
            splits=[
                args.train_ratio,
                args.valid_ratio,
                1 - args.train_ratio - args.valid_ratio,
            ],
            force_zero_last=args.force_zero_test,
        )
    else:
        raise ValueError


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RNA3DB")
    parser.add_argument(
        "--cpu", type=int, default=None, help="Number of CPUs to use when able"
    )

    # subparsers for different commands
    subparsers = parser.add_subparsers(
        dest="command", title="Available commands", required=True
    )

    # subparser for the "parse" command
    parse_parser = subparsers.add_parser("parse", help="Parse PDB and extract RNAs")
    parse_parser.add_argument(
        "input", type=Path, help="Directory containing mmCIF files to parse"
    )
    parse_parser.add_argument("output", type=Path, help="Output JSON file")
    parse_parser.add_argument(
        "--nmr_resolution",
        type=float,
        help="Resolution to use for NMR structures. By default we use float('inf').",
    )
    parse_parser.add_argument(
        "--include_atoms",
        action="store_true",
        help="Whether to include the XYZ atom coordinates in the parsed output.",
    )

    # subparser for the "filter" command
    filter_parser = subparsers.add_parser("filter", help="Filter a JSON")
    filter_parser.add_argument("input", type=Path, help="Input JSON file")
    filter_parser.add_argument("output", type=Path, help="Output JSON file")
    filter_parser.add_argument(
        "--single_ratio_cutoff",
        type=float,
        default=0.8,
        help="Filter chains where a single nucleotide makes up more than this fraction of residues",
    )
    filter_parser.add_argument(
        "--max_unknown_ratio",
        type=float,
        default=0.3,
        help="Filter chains with more than this fraction of unknown nucleotides",
    )
    filter_parser.add_argument(
        "--max_resolution",
        type=float,
        default=9.0,
        help="Filter chains over this resolution",
    )
    filter_parser.add_argument(
        "--min_length", type=int, default=32, help="Filter chains shorter than this"
    )
    filter_parser.add_argument(
        "--filter_log_path",
        type=Path,
        default=None,
        help="Path to filter log. The filter log shows which filters hit each sequence.",
    )

    # subparser for the "cluster" command
    cluster_parser = subparsers.add_parser(
        "cluster", help="Cluster RNAs by sequence and structure similarity"
    )
    cluster_parser.add_argument("input", type=Path, help="Input JSON file")
    cluster_parser.add_argument("output", type=Path, help="Output JSON file")
    cluster_parser.add_argument(
        "--tbl_dir", type=Path, help="Directory containing .tbl files"
    )
    cluster_parser.add_argument(
        "--min_seq_id", type=float, default=0.99, help="Minimum Sequence Identity"
    )
    cluster_parser.add_argument(
        "--min_seq_coverage", type=float, default=0.99, help="Minimum Sequence Coverage"
    )
    cluster_parser.add_argument(
        "--mmseqs_binary_path",
        type=Path,
        help="Path to MMseqs2 binary. May be required if RNA3DB cannot find MMseqs2's installation.",
    )
    cluster_parser.add_argument(
        "--mmseqs_coverage_mode", type=int, default=1, help="MMseqs Coverage Mode"
    )
    cluster_parser.add_argument(
        "--mmseqs_sensitivity", type=float, default=7.5, help="MMseqs Sensitivity"
    )
    cluster_parser.add_argument(
        "--mmseqs_alignment_mode", type=int, default=3, help="MMseqs Alignment Mode"
    )
    cluster_parser.add_argument(
        "--mmseqs_max_seqs", type=int, default=3, help="MMseqs max seqs"
    )
    cluster_parser.add_argument(
        "--structural_e_value_cutoff",
        type=float,
        default=1.0,
        help="Structural E-Value Cutoff used to build graph edges",
    )
    seq_struct_parser = cluster_parser.add_mutually_exclusive_group()
    seq_struct_parser.add_argument("--only_sequence", action="store_true")
    seq_struct_parser.add_argument("--only_structure", action="store_true")

    # subparser for the "split" command
    split_parser = subparsers.add_parser(
        "split", help="Split RNA data into training and test sets"
    )
    split_parser.add_argument("input", type=Path, help="Input JSON file")
    split_parser.add_argument("output", type=Path, help="Output JSON file")
    split_parser.add_argument(
        "--train_ratio",
        type=float,
        default=0.7,
        help="Ratio of data to use for the training set",
    )
    split_parser.add_argument(
        "--valid_ratio",
        type=float,
        default=0.0,
        help="Ratio of the data to use for the validation set",
    )
    split_parser.add_argument(
        "--force_zero_test",
        action="store_true",
        help="Force component zero into the test set",
    )

    args = parser.parse_args()
    main(args)
