from multiprocessing import Pool
from functools import partial
from pathlib import Path
from tqdm import tqdm
import argparse

from rna3db.parsers.structure import StructureFile
from rna3db.filter import apply_filters
from rna3db.cluster import cluster_sequences, cluster_structures
from rna3db.split import split
from rna3db.utils import read_json, write_json


def _read_as_dict(
    path: Path, nmr_resolution: float = None, include_atoms: bool = False
):
    d = {}
    try:
        sf = StructureFile(path, nmr_resolution, include_atoms)
        for chain in sf:
            chain_id = f"{sf.pdb_id}_{chain.author_id}"
            d[chain_id] = {
                "release_date": sf.release_date,
                "structure_method": sf.structure_method,
                "resolution": sf.resolution,
                "length": len(chain),
                "sequence": chain.sequence,
            }
            if include_atoms:
                d[chain_id]["atoms"] = [res.atoms for res in chain]
    except Exception as e:
        print(
            f"Unable to parse {path}. "
            f"Please report this at https://github.com/marcellszi/rna3db/issues."
        )
        print(f"Exception: {e}")
    return d


def _do_parse(args, input_path: Path, output_path: Path):
    files = list(input_path.glob("*.cif"))
    data = {}
    f = partial(
        _read_as_dict,
        nmr_resolution=args.nmr_resolution,
        include_atoms=args.include_atoms,
    )
    with Pool(processes=args.cpu) as p, tqdm(total=len(files)) as pbar:
        for d in p.imap_unordered(f, files):
            data |= d
            pbar.update()
    write_json(data, output_path)


def _do_filter(args, input_path: Path, output_path: Path):
    data = read_json(input_path)
    filtered = apply_filters(
        data,
        min_length=args.min_length,
        max_resolution=args.max_resolution,
        single_ratio_cutoff=args.single_ratio_cutoff,
        max_unknown_ratio=args.max_unknown_ratio,
        filter_log_path=args.filter_log_path,
    )
    write_json(filtered, output_path)


def _do_cluster(args, input_path: Path, output_path: Path):
    only_sequence = getattr(args, "only_sequence", False)
    only_structure = getattr(args, "only_structure", False)

    if not only_structure:
        cluster = cluster_sequences(
            input_path,
            output_path,
            mmseqs2_binary_path=args.mmseqs_binary_path,
            min_seq_id=args.min_seq_id,
            min_coverage=args.min_seq_coverage,
            coverage_mode=args.mmseqs_coverage_mode,
            sensitivity=args.mmseqs_sensitivity,
            alignment_mode=args.mmseqs_alignment_mode,
            max_seqs=args.mmseqs_max_seqs,
        )
        write_json(cluster, output_path)
        input_path = output_path

    if not only_sequence:
        cluster = cluster_structures(
            input_path, args.tbl_dir, args.structural_e_value_cutoff
        )
        write_json(cluster, output_path)


def _do_split(args, input_path: Path, output_path: Path):
    split(
        input_path,
        output_path,
        splits=[
            args.train_ratio,
            args.valid_ratio,
            1 - args.train_ratio - args.valid_ratio,
        ],
        force_zero_last=args.force_zero_test,
    )


def main(args):
    if args.command == "parse":
        _do_parse(args, args.input, args.output)
    elif args.command == "filter":
        _do_filter(args, args.input, args.output)
    elif args.command == "cluster":
        _do_cluster(args, args.input, args.output)
    elif args.command == "split":
        _do_split(args, args.input, args.output)
    elif args.command == "run":
        args.output.mkdir(parents=True, exist_ok=True)
        parsed = args.output / "parsed.json"
        filtered = args.output / "filtered.json"
        clustered = args.output / "clustered.json"
        split_out = args.output / "split.json"
        _do_parse(args, args.input, parsed)
        _do_filter(args, parsed, filtered)
        _do_cluster(args, filtered, clustered)
        _do_split(args, clustered, split_out)
    else:
        raise ValueError(f"Unknown command: {args.command}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RNA3DB")
    parser.add_argument(
        "--cpu", type=int, default=None, help="Number of CPUs to use when able"
    )

    subparsers = parser.add_subparsers(
        dest="command", title="Available commands", required=True
    )

    _parse_args = argparse.ArgumentParser(add_help=False)
    _parse_args.add_argument(
        "--nmr_resolution",
        type=float,
        default=None,
        help="Resolution to use for NMR structures. By default we use float('inf').",
    )
    _parse_args.add_argument(
        "--include_atoms",
        action="store_true",
        help="Include XYZ atom coordinates in the parsed output.",
    )

    _filter_args = argparse.ArgumentParser(add_help=False)
    _filter_args.add_argument(
        "--min_length", type=int, default=32, help="Filter chains shorter than this"
    )
    _filter_args.add_argument(
        "--max_resolution",
        type=float,
        default=9.0,
        help="Filter chains over this resolution",
    )
    _filter_args.add_argument(
        "--single_ratio_cutoff",
        type=float,
        default=0.8,
        help=(
            "Filter chains where a single nucleotide makes up more than "
            "this fraction of residues"
        ),
    )
    _filter_args.add_argument(
        "--max_unknown_ratio",
        type=float,
        default=0.3,
        help="Filter chains with more than this fraction of unknown nucleotides",
    )
    _filter_args.add_argument(
        "--filter_log_path",
        type=Path,
        default=None,
        help="Path to filter log showing which filters hit each sequence.",
    )

    _cluster_args = argparse.ArgumentParser(add_help=False)
    _cluster_args.add_argument(
        "--tbl_dir", type=Path, help="Directory containing .tbl files"
    )
    _cluster_args.add_argument(
        "--min_seq_id",
        type=float,
        default=0.99,
        help="Minimum sequence identity for a match to be retained (--min-seq-id, range 0.0-1.0).",
    )
    _cluster_args.add_argument(
        "--min_seq_coverage",
        type=float,
        default=0.99,
        help=(
            "Minimum fraction of aligned residues required for a match "
            "(-c, range 0.0-1.0). Interpreted according to --mmseqs_coverage_mode."
        ),
    )
    _cluster_args.add_argument(
        "--mmseqs_binary_path",
        type=Path,
        default=None,
        help=(
            "Path to MMseqs2 binary. May be required if RNA3DB cannot "
            "find MMseqs2's installation."
        ),
    )
    _cluster_args.add_argument(
        "--mmseqs_coverage_mode",
        type=int,
        default=1,
        help=(
            "Defines how --min_seq_coverage is applied (--cov-mode). "
            "0 = bidirectional, 1 = target coverage only, 2 = query coverage only."
        ),
    )
    _cluster_args.add_argument(
        "--mmseqs_sensitivity",
        type=float,
        default=7.5,
        help=(
            "Prefilter sensitivity (-s). Higher values find more distant homologs "
            "at the cost of speed: 1.0 fastest, 4.0 fast, 7.5 sensitive."
        ),
    )
    _cluster_args.add_argument(
        "--mmseqs_alignment_mode",
        type=int,
        default=3,
        help=(
            "Alignment information to compute (--alignment-mode). "
            "0 = automatic, 1 = score and end position, 2 = score/end/start, "
            "3 = full alignment with sequence identity, 4 = ungapped only."
        ),
    )
    _cluster_args.add_argument(
        "--mmseqs_max_seqs",
        type=int,
        default=10000,
        help=(
            "Maximum results per query passed by the prefilter (--max-seqs). "
            "Higher values increase sensitivity but slow down the search."
        ),
    )
    _cluster_args.add_argument(
        "--structural_e_value_cutoff",
        type=float,
        default=1.0,
        help="Structural E-value cutoff used to build graph edges",
    )

    _split_args = argparse.ArgumentParser(add_help=False)
    _split_args.add_argument(
        "--train_ratio",
        type=float,
        default=0.7,
        help="Ratio of data to use for the training set",
    )
    _split_args.add_argument(
        "--valid_ratio",
        type=float,
        default=0.0,
        help="Ratio of data to use for the validation set",
    )
    _split_args.add_argument(
        "--force_zero_test",
        action="store_true",
        help="Force component zero into the test set",
    )

    parse_parser = subparsers.add_parser(
        "parse", parents=[_parse_args], help="Parse mmCIF files and extract RNAs"
    )
    parse_parser.add_argument(
        "input", type=Path, help="Directory containing mmCIF files to parse"
    )
    parse_parser.add_argument("output", type=Path, help="Output JSON file")

    filter_parser = subparsers.add_parser(
        "filter", parents=[_filter_args], help="Filter a parsed JSON"
    )
    filter_parser.add_argument("input", type=Path, help="Input JSON file")
    filter_parser.add_argument("output", type=Path, help="Output JSON file")

    cluster_parser = subparsers.add_parser(
        "cluster",
        parents=[_cluster_args],
        help="Cluster RNAs by sequence and structure similarity",
    )
    cluster_parser.add_argument("input", type=Path, help="Input JSON file")
    cluster_parser.add_argument("output", type=Path, help="Output JSON file")
    seq_struct_group = cluster_parser.add_mutually_exclusive_group()
    seq_struct_group.add_argument("--only_sequence", action="store_true")
    seq_struct_group.add_argument("--only_structure", action="store_true")

    split_parser = subparsers.add_parser(
        "split", parents=[_split_args], help="Split clustered data into train/test sets"
    )
    split_parser.add_argument("input", type=Path, help="Input JSON file")
    split_parser.add_argument("output", type=Path, help="Output JSON file")

    run_parser = subparsers.add_parser(
        "run",
        parents=[_parse_args, _filter_args, _cluster_args, _split_args],
        help="Run all steps: parse, filter, cluster, split",
    )
    run_parser.add_argument(
        "input", type=Path, help="Directory containing mmCIF files to parse"
    )
    run_parser.add_argument(
        "output",
        type=Path,
        help=(
            "Output directory. Writes parsed.json, filtered.json, "
            "clustered.json, and split.json."
        ),
    )

    args = parser.parse_args()
    main(args)
