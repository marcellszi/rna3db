from typing import Sequence, Mapping, Tuple
from collections import defaultdict
from pathlib import Path
from Bio import PDB

from rna3db.utils import PathLike

import dataclasses
import json


def parse_as_dict(
    path: PathLike,
    modifications_cache_path: PathLike = None,
    molecule_type: str = "RNA",
    nmr_resolution=None,
):
    """Top-level API that parses an mmCIF/PDBx file as a Python dict.

    Args:
        path (PathLike): path to mmCIF file to parse
        nmr_resolution (float): resolution to use for NMR structures
            Default behaviour is to treat NMR resolution as float('inf').

    Returns:
        dict: a dictionary containing release_date, structure_method,
            resolution, length, and sequence
    """
    d = {}

    modification_handler = ModificationHandler(modifications_cache_path)

    # we want to make sure this never fails and interrupts the parse
    try:
        structure_file = StructureFile(
            path, modification_handler, molecule_type, nmr_resolution
        )
        for chain in structure_file:
            chain_id = f"{structure_file.pdb_id}_{chain.author_id}"
            d[chain_id] = {
                "release_date": structure_file.release_date,
                "structure_method": structure_file.structure_method,
                "resolution": structure_file.resolution,
                "length": len(chain),
                "sequence": chain.sequence,
            }
    except:
        # we just print an error if something went wrong, don't exit
        print(
            f"Unable to parse {path}. "
            f"Please report this at https://github.com/marcellszi/rna3db/issues."
        )

    return d


def parse_file(
    path: PathLike,
    modifications_cache_path: PathLike = None,
    molecule_type: str = "RNA",
    nmr_resolution=None,
):
    """Top-level API that parses an mmCIF/PDBx file as a StructureFile.

    Args:
        path (PathLike): path to mmCIF file to parse
        modifications_cache_path (PathLike, optinal): path to
            modifications_cache, default is None
        molecule_type (str, default="RNA"): One of "RNA" or "peptide".
            Specifies whether RNA or protein chains should be extracted from
            the file. NOTE: support for proteins is currently experimental.


    Returns:
        StructureFile: object containing data of parsed file
    """

    modification_handler = ModificationHandler(modifications_cache_path)
    return StructureFile(path, modification_handler, molecule_type, nmr_resolution)


class Residue:
    """Data class wrapping individual residues."""

    def __init__(
        self,
        three_letter_code: str,
        one_letter_code: str,
        index: int,
        atoms: Mapping[str, Tuple[float, float, float]],
    ):
        self.three_letter_code = three_letter_code
        self.one_letter_code = one_letter_code
        self.index = index
        self.atoms = atoms

    @property
    def code(self) -> str:
        return self.one_letter_code

    def __repr__(self):
        return (
            f"Residue(code={self.code}, three_letter_code={self.three_letter_code}, "
            f"index={self.index}"
        )


class Chain:
    """Data class wrapping chains. Contains a list of Residues."""

    def __init__(self, author_id: str, residues: Sequence[Residue]):
        self.author_id = author_id
        self.residues = residues

    def __iter__(self):
        if len(self) == 0:
            return None
        return iter(self.residues)

    def __getitem__(self, idx):
        if len(self) == 0:
            return None
        return self.residues[idx]

    def __len__(self):
        return len(self.residues)

    @property
    def sequence(self):
        return "".join(i.code for i in self.residues)

    def __repr__(self):
        return f"Chain(author_id={self.author_id}, len={len(self)})"

    def __str__(self):
        max_line_length = 120
        idx_steps = 50

        numbers_str = [" " for i in self.sequence]
        lbl = list(str(len(self.sequence)))
        numbers_str[-len(lbl) :] = lbl
        numbers_str[0] = "1"

        for idx in range(idx_steps, len(self.sequence), idx_steps):
            lbl = list(str(idx + 1))
            numbers_str[idx - len(lbl) : idx] = lbl

        numbers_str = "".join(numbers_str)

        S = ""
        for idx in range(0, len(self.sequence), max_line_length):
            S += numbers_str[idx : idx + max_line_length] + "\n"
            S += self.sequence[idx : idx + max_line_length] + "\n"
            S += "\n"

        return S


class ModificationHandler:
    JSON_PATH_DEFAULTS = [
        "modifications_cache.json",
        "data/modifications_cache.json",
        "../data/modifications_cache.json",
    ]

    def __init__(self, json_path: PathLike = None):
        """Used for converting `three_letter_code`s to `one_letter_code`s, including modifications.

        Args:
            json_path (PathLike, optional): Path to `modifications_cache.json`, generated from a Chemical Component
                Dictionary .cif file. If not provided, `JSON_PATH_DEFAULTS` are checked for a valid path.

        Attributes:
            JSON_PATH_DEFAULTS (Sequence[PathLike]): List of paths checked for a `modifications_cache.json` file when
                `json_path` is not initialised.
        """
        # try to see we can find the modifications_cache.json in any of the default paths
        if json_path is None:
            for path in self.JSON_PATH_DEFAULTS:
                if Path(path).is_file():
                    json_path = Path(path)
                    break

        # if we still haven't found it, we need to throw an exception
        if json_path is None:
            raise FileNotFoundError(
                "Could not find `modifications_cache.json` in JSON_PATH_DEFAULTS."
            )

        with open(json_path) as f:
            self.modifications = json.load(f)

    def is_rna(self, three_letter_code: str) -> bool:
        """Check if `three_letter_code` is RNA nucleic acid.

        Args:
            three_letter_code (str): Three letter code to check.

        Returns:
            bool: True if `three_letter_code` is typically an RNA nucleic acid, False otherwise.

        """
        return three_letter_code in self.modifications["rna"]

    def is_peptide(self, three_letter_code: str) -> bool:
        """Check if `three_letter_code` is an amino acid.

        Args:
            three_letter_code (str): Three letter code to check.

        Returns:
            bool: True if `three_letter_code` is an amino acid.

        """
        return three_letter_code in self.modifications["peptide"]

    def peptide_letters_3to1(self, three_letter_code: str) -> str:
        """Convert amino acid `three_letter_code` to `one_letter_code`.

        Args:
            three_letter_code (str): Three letter code to check.

        Returns:
           str: one_letter_code of an amino acid, "X" if cannot be found.

        """
        return self.modifications["peptide"].get(three_letter_code, "X")

    def rna_letters_3to1(self, three_letter_code: str) -> str:
        """Convert RNA nucleic acid `three_letter_code` to `one_letter_code`.

        Args:
            three_letter_code (str): Three letter code to check.

        Returns:
           str: one_letter_code of RNA nucleic acid, "N" if cannot be found.
        """
        return self.modifications["rna"].get(three_letter_code, "N")


class StructureFile:
    def __init__(
        self,
        path: PathLike,
        modification_handler: ModificationHandler,
        molecule_type: str = "RNA",
        nmr_resolution: float = 0.0,
    ):
        """A StructureFile encapsulates a parsed mmCIF or PDB file and contains convenient high-level APIs for common
        tasks.

        Note: File type is inferred from the extension. File extensions are not case-sensitive.
            Valid types: mmCIF (`.cif`, `.mmcif`)

        Args:
            path (:PathLike:): The path to the file.
            modification_handler (:ModificationHandler:):
            molecule_type (str, default="RNA"): One of "RNA" or "peptide". Specifies whether RNA or protein chains
                should be extracted from the file.

        Attributes:
            path (Path): The path of the input file.
            pdb_id (str): The PDB ID as read from the file.
            release_date (str): Date given as a `str` in ISO 8601 format.
            resolution (float): Resolution in ångströms.
            structure_method (str): Method used to resolve the structure.
            chains (Mapping[str, :Chain:])
        """
        # determine which parser to use
        path = Path(path)
        if path.suffix.lower() in [".cif", ".mmcif"]:
            file_parser = mmCIFParser
        elif path.suffix.lower() == ".pdb":
            raise NotImplementedError
            file_parser = PDBParser
        else:
            raise ValueError(
                f"The extension {self.path.suffix.lower()} is not supported."
            )

        # make the parser
        parser = file_parser(path, modification_handler, molecule_type)

        # use parser to get attributes
        self.pdb_id = parser.pdb_id
        self.release_date = parser.release_date
        self.resolution = parser.resolution
        self.structure_method = parser.structure_method
        self.chains = parser.chains

    def __repr__(self):
        return f"Structure({self.path})"

    def __getitem__(self, idx):
        return self.chains[idx]

    def __iter__(self):
        return iter(self.chains.values())

    def __repr__(self):
        return (
            f"StructureFile(pdb_id={self.pdb_id}, chains={self.chains.keys()}, "
            f"resolution={self.resolution}, release_date={self.release_date}, structure_method={self.structure_method})"
        )


class FileParser:
    def __init__(self, path, modification_handler, molecule_type, nmr_resolution=None):
        self.path = path
        self.molecule_type = molecule_type
        self.nmr_resolution = nmr_resolution
        if molecule_type == "RNA":
            self.letters_3to1 = lambda x: modification_handler.rna_letters_3to1(x)
        elif molecule_type == "peptide":
            self.letters_3to1 = lambda x: modification_handler.peptide_letters_3to1(x)
        else:
            raise ValueError('molecule_type must be one of "RNA" or "peptide".')


class mmCIFParser(FileParser):
    def __init__(self, path, modification_handler, molecule_type, nmr_resolution=None):
        super().__init__(path, modification_handler, molecule_type, nmr_resolution)
        self.parsed_info = PDB.MMCIF2Dict.MMCIF2Dict(self.path)

    @property
    def pdb_id(self):
        return self.parsed_info["_entry.id"][0].lower()

    @property
    def release_date(self):
        # prefer to get the date from the earliest revision date
        if "_pdbx_audit_revision_history.revision_date" in self.parsed_info:
            return min(self.parsed_info["_pdbx_audit_revision_history.revision_date"])
        # use deposition date if there are no revisions
        return self.parsed_info["_pdbx_database_status.recvd_initial_deposition_date"]

    @property
    def resolution(self):
        resolutions = []
        for res_key in [
            "_refine.ls_d_res_high",
            "_em_3d_reconstruction.resolution",
            "_reflns.d_resolution_high",
        ]:
            if res_key in self.parsed_info:
                if self.parsed_info[res_key][0] not in ".?":
                    resolutions.append(float(self.parsed_info[res_key][0]))

        if self.structure_method == "solution nmr" and self.nmr_resolution is not None:
            return self.nmr_resolution

        if len(resolutions) == 0:
            return float("inf")

        return max(resolutions)


    @property
    def structure_method(self):
        return ",".join(self.parsed_info["_exptl.method"]).lower()


    @dataclasses.dataclass
    class _AtomSite:
        atom_id: str
        three_letter_code: str
        author_chain_id: str
        author_seq_num: int
        mmcif_seq_num: int
        insertion_code: str
        hetatm_atom: str
        x: float
        y: float
        z: float


    @staticmethod
    def _get_atom_sites(parsed_info: PDB.MMCIF2Dict):
        return [
            mmCIFParser._AtomSite(*site)
            for site in zip(
                parsed_info["_atom_site.label_atom_id"],  # atom name
                parsed_info["_atom_site.label_comp_id"],  # residue name
                parsed_info["_atom_site.auth_asym_id"],  # author chain
                parsed_info["_atom_site.auth_seq_id"],  # author_seq_num
                parsed_info["_atom_site.label_seq_id"],  # mmcif_seq_num
                parsed_info["_atom_site.pdbx_PDB_ins_code"],  # insertion code?
                parsed_info["_atom_site.group_PDB"],  # hetatm_atom
                parsed_info["_atom_site.Cartn_x"],  # x
                parsed_info["_atom_site.Cartn_y"],  # y
                parsed_info["_atom_site.Cartn_z"],  # z
            )
        ]


    @property
    def chains(self):
        # no SEQRES chains in this file
        if "_entity_poly_seq.entity_id" not in self.parsed_info:
            return {}

        # get the mapping from entity_id to internal mmcif_chain_id
        mmcif_chain_to_entity_id = {
            mmcif_chain_id: entity_id
            for mmcif_chain_id, entity_id in zip(
                self.parsed_info["_struct_asym.id"],
                self.parsed_info["_struct_asym.entity_id"],
            )
        }

        # create mapping that maps entity_id to author_chain_id
        # (used to get "seqres" chain names)
        # NOTE: these are not unique, so there is a set of author_ids for each
        #       entity_id
        id_map = defaultdict(set)
        for author_chain_id, mmcif_chain_id in zip(
            self.parsed_info["_atom_site.auth_asym_id"],
            self.parsed_info["_atom_site.label_asym_id"],
        ):
            k = mmcif_chain_to_entity_id[mmcif_chain_id]
            id_map[k].add(author_chain_id)

        # get the chem_comp type for each mon_id
        chem_comp_type = {
            mon_id: comp_type
            for mon_id, comp_type in zip(
                self.parsed_info["_chem_comp.id"], self.parsed_info["_chem_comp.type"]
            )
        }

        # parse full chains from "seqres"
        chains_full = defaultdict(list)
        for entity_id, mon_id, idx in zip(
            self.parsed_info["_entity_poly_seq.entity_id"],
            self.parsed_info["_entity_poly_seq.mon_id"],
            self.parsed_info["_entity_poly_seq.num"],
        ):
            for author_id in id_map[entity_id]:
                chains_full[author_id].append(
                    Residue(
                        three_letter_code=mon_id,
                        one_letter_code=self.letters_3to1(mon_id),
                        index=int(idx),
                        atoms={} # placeholder until we get coordinates
                    )
                )

        # keep only chains that contain at least one 'self.molecule_type'
        chains = {}
        for author_chain_id, chain_sequence in chains_full.items():
            if any(
                [
                    self.molecule_type in chem_comp_type[i.three_letter_code]
                    for i in chain_sequence
                ]
            ):
                chains[author_chain_id] = Chain(
                    author_id=author_chain_id,
                    residues=chain_sequence,
                )


        # find starting index of relevant chains
        seq_start_num = {
            author_chain_id: min([res.index for res in chain])
            for author_chain_id, chain in chains.items()
        }

        # iterate through atom sites to get coordinates
        for site in mmCIFParser._get_atom_sites(self.parsed_info):
            # if not a relevant chain, don't care about atoms
            if site.author_chain_id not in chains:
                continue

            # TODO: verify this is fine
            # I think these are just for H_* and W?
            if site.mmcif_seq_num == ".":
                continue

            # the idx is just the mmcif_seq_num - seq_start_num
            seq_idx = int(site.mmcif_seq_num) - seq_start_num[site.author_chain_id]

            # make sure that the sites actually match, should never be a mismatch
            assert (
                site.three_letter_code
                == chains[site.author_chain_id][seq_idx].three_letter_code
            )
            assert seq_idx == chains[site.author_chain_id][seq_idx].index - 1

            chains[site.author_chain_id][seq_idx].author_index = int(site.author_seq_num)

            # add atom coordinates
            chains[site.author_chain_id][seq_idx].atoms[site.atom_id] = tuple(
                map(float, (site.x, site.y, site.z))
            )

        return chains


def parse_fasta(path):
    """Parse a FASTA file.

    Supports multi-line sequences.

    Args:
        path (PathLike): Path to input FASTA file.

    Returns:

    """
    with open(path) as f:
        descriptions = []
        sequences = []
        i = -1
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                i += 1
                descriptions.append(line[1:])
                sequences.append("")
            elif line.startswith("#") or not line:
                continue
            else:
                sequences[i] += line
    return descriptions, sequences


def write_fasta(
    descriptions: Sequence[str], sequences: Sequence[str], output_path: PathLike
):
    """Write to a FASTA file.

    Args:
        descriptions (Sequence): List of descriptions for each sequence.
        sequences (Sequence): List of sequences.
        output_path (PathLike): Path to write FASTA file to.
    """
    if len(descriptions) != len(sequences):
        raise ValueError("The length of descriptions and sequences must match.")
    with open(output_path, "w") as f:
        for k, v in zip(descriptions, sequences):
            f.write(f">{k}\n{v}\n")
