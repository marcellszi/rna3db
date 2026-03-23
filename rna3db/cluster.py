from collections import defaultdict
from typing import Sequence
from pathlib import Path

from rna3db.parsers import tabular, fasta
from rna3db.utils import read_json

import subprocess
import tempfile
import logging
import os

_MMSEQS2_CLUSTER_PREFIX = "mmseqs2"


class InfernalGraph:
    """A class for InfernalGraphs supporting family and PDB chain nodes."""

    def __init__(self):
        self.graph = {}

    def add_chain(self, node: str):
        """Adds a chain node to the graph.

        Args:
            node (str): The identifier of the chain.
        """
        if node not in self.graph:
            self.graph[node] = {"is_family": False, "neighbours": set()}

    def add_family(self, node: str):
        """Adds a family node to the graph.

        Args:
            node (str): The Rfam accession of the family.
        """
        if node not in self.graph:
            self.graph[node] = {"is_family": True, "neighbours": set()}

    def add_edge(self, node1: str, node2: str):
        """Adds an edge between two nodes in the graph.

        Args:
            node1 (str): The name of the first node.
            node2 (str): The name of the second node.
        """
        if node1 in self.graph and node2 in self.graph:
            self.graph[node1]["neighbours"].add(node2)
            self.graph[node2]["neighbours"].add(node1)
        else:
            raise ValueError

    def components(self) -> Sequence[Sequence[str]]:
        """Performs DFS to find disjoint components (chains only) of the graph.

        The returned components only contain the chains, and do not output
        families.

        Returns:
            list: A list of sets, where each set represents a connected
                component.
        """
        visited = set()
        components = []

        def dfs(node, component):
            visited.add(node)
            if not self.graph[node]["is_family"]:
                component.add(node)
            for neighbour in self.graph[node]["neighbours"]:
                if neighbour not in visited:
                    dfs(neighbour, component)

        for node in self.graph:
            if node not in visited:
                component = set()
                dfs(node, component)
                components.append(component)

        return components


def _run_mmseqs2(
    binary_path: str,
    fasta_path: Path,
    output_path: Path,
    min_seq_id: float,
    min_coverage: float,
    coverage_mode: int,
    sensitivity: float,
    alignment_mode: int,
    max_seqs: int,
):
    fasta_path = Path(fasta_path).resolve()
    curr_dir = os.getcwd()
    os.chdir(output_path)

    cmd = [
        binary_path,
        "easy-cluster",
        "--min-seq-id",
        min_seq_id,
        fasta_path,
        _MMSEQS2_CLUSTER_PREFIX,
        "mmseqs2_tmp",
        "-c",
        min_coverage,
        "--cov-mode",
        coverage_mode,
        "--max-seqs",
        max_seqs,
        "-s",
        sensitivity,
        "--alignment-mode",
        alignment_mode,
    ]
    cmd = list(map(str, cmd))

    logging.info(f'Launching subprocess {" ".join(cmd)}')
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    retcode = process.wait()
    os.chdir(curr_dir)

    if retcode:
        logging.error(f"{stderr}:\n{stdout}")
    else:
        logging.info(stdout)


def cluster_sequences(
    input_json_path: Path,
    output_json_path: Path,
    mmseqs2_binary_path: Path = None,
    min_seq_id: float = 0.99,
    min_coverage: float = 0.99,
    coverage_mode: int = 1,
    sensitivity: float = 7.5,
    alignment_mode: int = 3,
    max_seqs: int = 10000,
) -> dict:
    """Cluster sequences by similarity using MMseqs2 ``easy-cluster``.

    For a practical overview of clustering parameters and their interaction,
    see the MMseqs2 wiki tutorials:
    https://github.com/soedinglab/MMseqs2/wiki/Tutorials

    For full parameter reference, see the MMseqs2 user guide:
    https://mmseqs.com/latest/userguide.pdf

    Args:
        input_json_path (Path): Path to input JSON file.
        output_json_path (Path): Path to output JSON file. MMseqs2 scratch
            files are written to the same directory.
        mmseqs2_binary_path (Path, optional): Path to MMseqs2 binary. Inferred
            from PATH if not provided.
        min_seq_id (float): Minimum sequence identity for a match to be
            retained (``--min-seq-id``, range 0.0-1.0).
        min_coverage (float): Minimum fraction of aligned residues required
            for a match (``-c``, range 0.0-1.0). Interpreted according to
            ``coverage_mode``.
        coverage_mode (int): Defines how ``min_coverage`` is applied
            (``--cov-mode``). Common values: 0 = bidirectional (both query
            and target must meet coverage), 1 = target coverage only,
            2 = query coverage only. See user guide for modes 3-5.
        sensitivity (float): Controls the prefilter sensitivity (``-s``).
            Higher values find more distant homologs at the cost of speed;
            1.0 is fastest, 4.0 is fast, 7.5 is sensitive.
        alignment_mode (int): Determines which alignment information is
            computed (``--alignment-mode``). 0 = automatic, 1 = score and
            end position only, 2 = score/end/start position, 3 = full
            alignment with sequence identity, 4 = only ungapped alignment.
        max_seqs (int): Maximum number of results per query sequence passed
            by the prefilter (``--max-seqs``). Higher values increase
            sensitivity but slow down the search.

    Returns:
        dict: Mapping from representative sequence to cluster members.
    """
    if mmseqs2_binary_path is None:
        mmseqs2_binary_path = (
            subprocess.check_output(["which", "mmseqs"]).decode("utf-8").strip()
        )

    data = read_json(input_json_path)
    descriptions = list(data.keys())
    sequences = [v["sequence"] for v in data.values()]

    with tempfile.NamedTemporaryFile() as fasta_f:
        fasta.write(descriptions, sequences, fasta_f.name)
        _run_mmseqs2(
            binary_path=mmseqs2_binary_path,
            fasta_path=fasta_f.name,
            output_path=output_json_path.parent,
            min_seq_id=min_seq_id,
            min_coverage=min_coverage,
            coverage_mode=coverage_mode,
            sensitivity=sensitivity,
            alignment_mode=alignment_mode,
            max_seqs=max_seqs,
        )

    tsv_path = output_json_path.parent / f"{_MMSEQS2_CLUSTER_PREFIX}_cluster.tsv"
    sequence_cluster = defaultdict(dict)
    with open(tsv_path) as f:
        for line in f:
            repr_sequence, sequence = line.split()
            sequence_cluster[repr_sequence][sequence] = data[sequence]

    return sequence_cluster


def cluster_structures(
    input_json_path: Path,
    tbl_dir: Path,
    e_value_cutoff: float = 1.0,
) -> dict:
    """Cluster structures by RNA family using Infernal tabular output.

    Builds a bipartite graph of chains and Rfam families, then finds connected
    components. Chains with no Infernal hits are grouped into component_0.

    Args:
        input_json_path (Path): Path to input JSON file (sequence-clustered).
        tbl_dir (Path): Directory containing Infernal ``.tbl`` output files.
        e_value_cutoff (float): Maximum E-value for a hit to be included as a
            graph edge.

    Returns:
        dict: Mapping from component name to chains in that component.
    """
    graph = InfernalGraph()

    data = read_json(input_json_path)
    tbl = tabular.read(tbl_dir)
    tbl = tbl.filter_attr_by_set("query_name", set(data.keys()))
    tbl = tbl.filter_e_value(e_value_cutoff)

    repr_mapping = {}
    for repr_k, repr_v in data.items():
        for chain_k in repr_v.keys():
            repr_mapping[chain_k] = repr_k

    for chain_node in set(data.keys()):
        graph.add_chain(chain_node)

    for family_node in set(tbl.target_accession):
        graph.add_family(family_node)

    for hit in tbl:
        query_k = repr_mapping[hit.query_name]
        graph.add_edge(query_k, hit.target_accession)

    components = graph.components()

    components_dict = defaultdict(dict)
    hit_chains = set(tbl.query_name)
    curr = 1
    for component in sorted(components, key=len, reverse=True):
        if not any(i in hit_chains for i in component):
            name = "component_0"
        else:
            name = f"component_{curr}"
            curr += 1
        for chain in component:
            components_dict[name][chain] = data[chain]

    return components_dict
