from collections import defaultdict
from typing import Sequence
from pathlib import Path

from rna3db.tabular import read_tbls_from_dir
from rna3db.utils import PathLike, read_json
from rna3db.parser import write_fasta

import subprocess
import tempfile
import logging
import os


class InfernalGraph:
    """A class for InfernalGraphs supporting family and PDB chain nodes."""

    def __init__(self):
        self.graph = {}

    def add_chain(self, node: str):
        """
        Adds a chain node to the graph.

        Args:
            node (str): The identifier of the chain.
        """
        if node not in self.graph:
            self.graph[node] = {"is_family": False, "neighbours": set()}

    def add_family(self, node: str):
        """
        Adds a family node to the graph.

        Args:
            node (str): The Rfam accession of the family.
        """
        if node not in self.graph:
            self.graph[node] = {"is_family": True, "neighbours": set()}

    def add_edge(self, node1: str, node2: str):
        """
        Adds an edge between two nodes in the graph.

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
        """
        Performs DFS to find disjoint components (chains only) of the graph.

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


class StructureClusterer:
    def __init__(self, e_value_cutoff: float = 1):
        self.e_value_cutoff = e_value_cutoff

    def cluster(self, input_json_path: PathLike, input_tbls_dir: PathLike):
        # make new graph
        graph = InfernalGraph()

        # read required data
        data = read_json(input_json_path)
        tbl = read_tbls_from_dir(input_tbls_dir)
        tbl = tbl.filter_attr_by_set("query_name", set(data.keys()))
        tbl = tbl.filter_e_value(self.e_value_cutoff)

        # add chain nodes
        for chain_node in set(data.keys()):
            graph.add_chain(chain_node)

        # add family nodes
        for family_node in set(tbl.target_accession):
            graph.add_family(family_node)

        # add edges
        for hit in tbl:
            graph.add_edge(hit.query_name, hit.target_accession)

        # do DFS
        components = graph.components()

        # build dictionary ready to be written as JSON
        components_dict = defaultdict(dict)
        hit_chains = set(tbl.query_name)
        curr = 1
        for component in sorted(components, key=len, reverse=True):
            # force chains with no hits to go into component_0
            if not any([i in hit_chains for i in component]):
                name = "component_0"
            else:
                name = f"component_{curr}"
                curr += 1

            # add chains to components_dict
            for chain in component:
                components_dict[name][chain] = data[chain]

        return components_dict


class SequenceClusterer:
    def __init__(
        self,
        mmseqs2_binary_path: PathLike,
        min_seq_id: float = 0.99,
        min_coverage: float = 0.99,
        coverage_mode: int = 1,
        sensitivity: float = 7.5,
        alignment_mode: int = 3,
        max_seqs: int = 10000,
    ):
        """Cluster by sequence similarity using MMseqs2."""

        # attempt to infer binary path
        if mmseqs2_binary_path is None:
            self.mmseqs2_binary_path = (
                subprocess.check_output(["which", "mmseqs"]).decode("utf-8").strip()
            )
        else:
            self.mmseqs2_binary_path = mmseqs2_binary_path
        self.cluster_prefix = "mmseqs2"

        # set parameters
        self.min_seq_id = min_seq_id
        self.min_coverage = min_coverage
        self.coverage_mode = coverage_mode
        self.sensitivity = sensitivity
        self.alignment_mode = alignment_mode
        self.max_seqs = max_seqs

    def cluster(self, input_json_path: PathLike, output_json_path: PathLike):
        # read json
        data = read_json(input_json_path)

        # extract descriptions and sequence from json for intermediate FASTA
        descriptions, sequences = [], []
        for k, v in data.items():
            descriptions.append(k)
            sequences.append(v["sequence"])
        # create temp FASTA and run mmseqs2
        with tempfile.NamedTemporaryFile() as fasta_f:
            write_fasta(descriptions, sequences, fasta_f.name)
            self._mmseqs2(
                fasta_path=fasta_f.name,
                output_path=output_json_path.parent,
                min_seq_id=self.min_seq_id,
                min_coverage=self.min_coverage,
                coverage_mode=self.coverage_mode,
                sensitivity=self.sensitivity,
                alignment_mode=self.alignment_mode,
                max_seqs=self.max_seqs,
            )

        sequence_cluster = defaultdict(dict)
        with open(output_json_path.parent / f"{self.cluster_prefix}_cluster.tsv") as f:
            for line in f:
                repr_sequence, sequence = line.split()
                sequence_cluster[repr_sequence][sequence] = data[sequence]

        return sequence_cluster

    def _mmseqs2(
        self,
        fasta_path,
        output_path,
        min_seq_id: float = 0.99,
        min_coverage: float = 0.99,
        coverage_mode: int = 1,
        sensitivity: float = 7.5,
        alignment_mode: int = 3,
        max_seqs: int = 10000,
    ):
        fasta_path = Path(fasta_path)
        fasta_path = fasta_path.resolve()
        curr_dir = os.getcwd()
        os.chdir(output_path)
        tmpdir = "mmseqs2_tmp"

        cmd = [
            self.mmseqs2_binary_path,
            "easy-cluster",
            "--min-seq-id",
            min_seq_id,
            fasta_path,
            self.cluster_prefix,
            tmpdir,
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
