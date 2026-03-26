import dataclasses
from collections import defaultdict
from pathlib import Path
from typing import Iterator, Mapping, Sequence, Tuple

from Bio import PDB

from rna3db.ccd.chem_comp import load as load_chem_comps
from rna3db.ccd.modifications import ModificationHandler


class Residue:
    """Data class wrapping individual residues."""

    def __init__(
        self,
        three_letter_code: str,
        one_letter_code: str,
        index: int,
        atoms: dict = None,
    ):
        """
        Args:
            three_letter_code (str): Three-letter CCD residue code (e.g. ``"GTP"``).
            one_letter_code (str): Single-letter code (e.g. ``"G"``).
            index (int): Zero-based sequence index.
            atoms (dict, optional): Mapping of atom name to (x, y, z) coordinates.
        """
        self.three_letter_code = three_letter_code
        self.one_letter_code = one_letter_code
        self.index = index
        # NOTE: we need to handle dict like this, cannot use `atoms: dict = {}`
        #       in the method definition. See important warning:
        #       https://docs.python.org/3/tutorial/controlflow.html#default-argument-values
        #       (the default value is evaluated only once, causing issues with
        #       mutable dicts)
        self.atoms = atoms if atoms else {}

    @property
    def code(self) -> str:
        return self.one_letter_code

    @property
    def is_missing(self) -> bool:
        return not len(self.atoms) > 0

    def __eq__(self, other: object) -> bool:
        # NOTE: we don't care about three letter codes, only one letter
        #       this means modifications are still equal
        return (
            self.one_letter_code == other.one_letter_code
            and self.index == other.index
            and self.atoms == other.atoms
        )

    def __repr__(self):
        return (
            f"Residue(code={self.code}, three_letter_code={self.three_letter_code}, "
            f"index={self.index}, is_missing={self.is_missing})"
        )


class Chain:
    """Data class wrapping chains. Contains a list of Residues."""

    def __init__(self, author_id: str = None):
        """
        Args:
            author_id (str, optional): Author chain identifier as recorded in
                the mmCIF file.
        """
        self.author_id = author_id
        self.residues = []

    def __iter__(self):
        if len(self) == 0:
            return None
        return iter(self.residues)

    def __getitem__(self, idx):
        if len(self) == 0:
            return None
        return self.residues[idx]

    def __len__(self) -> int:
        return len(self.residues)

    def __eq__(self, other: object) -> bool:
        # NOTE: we ignore the author_id for equality checks
        if len(self) != len(other):
            return False

        for res_self, res_other in zip(self, other):
            if res_self != res_other:
                return False

        return True

    @property
    def has_atoms(self) -> bool:
        return any([not res.is_missing for res in self])

    def add_residue(self, res: Residue):
        """Add a residue to the chain.

        Note:
            Residues must be added in increasing index order. Duplicate indices
            are silently ignored (see 4x4t_G for motivation). Gaps in the
            sequence are filled with ``N`` residues.

        Args:
            res (Residue): Residue to add to the chain.

        Raises:
            ValueError: If ``res`` has an index less than the last added residue.
        """
        if len(self) == 0 or self.residues[-1].index == res.index - 1:
            self.residues.append(res)
        elif self.residues[-1].index < res.index - 1:
            for idx in range(self.residues[-1].index + 1, res.index):
                self.residues.append(Residue("N", "N", idx))
            self.residues.append(res)
        elif self.residues[-1].index == res.index:
            pass
        else:
            raise ValueError(f"Cannot add residues out of order.")

    @property
    def sequence(self) -> str:
        return "".join(i.code for i in self.residues)

    def __repr__(self) -> str:
        return f"Chain(author_id={self.author_id}, len={len(self)})"

    def __str__(self) -> str:
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


class Structure:
    def __init__(
        self,
        path: Path,
        nmr_resolution: float = None,
        include_atoms: bool = False,
        format: str = None,
    ):
        """Encapsulates a parsed structure file with high-level access to chains.

        Note:
            File format is inferred from the extension (case-insensitive) when
            ``format`` is not provided. Recognised extensions: ``.cif``,
            ``.mmcif``, ``.pdb``. Unknown extensions default to mmCIF. PDB
            format is not currently supported.

        Args:
            path (Path): Path to the structure file.
            nmr_resolution (float, optional): Resolution to assign to NMR structures.
                Default behaviour is to treat NMR resolution as float('inf').
            include_atoms (bool, optional): If True, atom coordinates are parsed and
                stored for each Residue. Default is False.
            format (str, optional): File format, one of ``"mmcif"`` or ``"pdb"``. If
                None (default), inferred from file extension; falls back to mmCIF.

        Attributes:
            pdb_id (str): The PDB ID as read from the file.
            release_date (str): Date in ISO 8601 format.
            resolution (float): Resolution in ångströms.
            structure_method (str): Method used to resolve the structure.
            chains (Mapping[str, Chain]): Mapping of author chain ID to Chain.
        """
        # determine which parser to use
        path = Path(path)
        if format is not None:
            fmt = format.lower()
        elif path.suffix.lower() in [".cif", ".mmcif"]:
            fmt = "mmcif"
        elif path.suffix.lower() == ".pdb":
            fmt = "pdb"
        else:
            fmt = "mmcif"  # default to mmCIF for unrecognised extensions

        if fmt == "pdb":
            raise NotImplementedError(
                "PDB format is not currently supported. Please use PDBx/mmCIF."
            )
        file_parser = mmCIFParser

        modification_handler = ModificationHandler()

        # make the parser
        parser = file_parser(path, modification_handler, nmr_resolution, include_atoms)

        # use parser to get attributes
        self.pdb_id = parser.pdb_id
        self.release_date = parser.release_date
        self.resolution = parser.resolution
        self.structure_method = parser.structure_method
        self.chains = parser.chains
        self.auth_asym_to_label_asym = parser.auth_asym_to_label_asym

    @classmethod
    def read(
        cls,
        path: Path,
        nmr_resolution: float = None,
        include_atoms: bool = False,
        format: str = None,
    ) -> "Structure":
        """Parse a structure file.

        Args:
            path (Path): Path to structure file to parse.
            nmr_resolution (float, optional): Resolution to assign to NMR structures.
                Default behaviour is to treat NMR resolution as float('inf').
            include_atoms (bool, optional): If True, atom coordinates are parsed and
                stored for each Residue. Default is False.
            format (str, optional): File format, one of ``"mmcif"`` or ``"pdb"``. If
                None (default), inferred from file extension; falls back to mmCIF.

        Returns:
            Structure: Parsed structure.
        """
        return cls(path, nmr_resolution, include_atoms, format=format)

    def __getitem__(self, idx: str) -> Chain:
        return self.chains[idx]

    def __iter__(self) -> Iterator[Chain]:
        return iter(self.chains.values())

    def __repr__(self) -> str:
        return (
            f"Structure(pdb_id={self.pdb_id}, chains={self.chains.keys()}, "
            f"resolution={self.resolution}, release_date={self.release_date}, "
            f"structure_method={self.structure_method})"
        )

    @staticmethod
    def _gen_mmcif_loop_str(
        name: str, headers: Sequence[str], values: Sequence[tuple]
    ) -> str:
        s = "#\nloop_\n"
        for header in headers:
            s += f"_{name}.{header}\n"

        max_widths = {k: 0 for k in headers}
        for V in values:
            for k, v in zip(headers, V):
                max_widths[k] = max(max_widths[k], len(str(v)))

        for V in values:
            row = ""
            for k, v in zip(headers, V):
                row += f"{str(v):<{max_widths[k]}} "
            s += row + "\n"

        return s

    def write(self, output_path: Path, author_id: str):
        """Write a single chain to a minimal mmCIF file.

        Args:
            output_path (Path): Path to write the mmCIF file to.
            author_id (str): Author chain identifier to write.

        Raises:
            ValueError: If no atom coordinates are available for the chain.
                Ensure the file was read with ``include_atoms=True``.
        """
        if not self[author_id].has_atoms:
            raise ValueError(
                f"Did not find any atoms for chain {author_id}. "
                f"Did you set `include_atoms=True`?"
            )
        # extract needed info
        entity_poly_seq_data = []
        atom_site_data = []
        for i, res in enumerate(self[author_id]):
            entity_poly_seq_data.append((1, res.index + 1, res.code, "n"))
            for idx, (atom_name, atom_coords) in enumerate(res.atoms.items()):
                x, y, z = atom_coords
                # fmt: off
                atom_site_data.append(
                    (
                        "ATOM",       # group_PDB
                        idx + 1,      # id
                        atom_name[0], # type_symbol (element)
                        atom_name,    # label_atom_id
                        ".",          # label_alt_id
                        res.code,     # label_comp_id
                        author_id,    # label_asym_id
                        "?",          # label_entity_id
                        i + 1,        # label_seq_id
                        "?",          # pdbx_PDB_ins_code
                        x,            # Cartn_x
                        y,            # Cartn_y
                        z,            # Cartn_z
                        1.0,          # occupancy
                        0.0,          # B_iso_or_equiv
                        "?",          # pdbx_formal_charge
                        i + 1,        # auth_seq_id
                        res.code,     # auth_comp_id
                        author_id,    # auth_asym_id
                        atom_name,    # auth_atom_id
                        1,            # pdbx_PDB_model_num
                    )
                )
                # fmt: on

        # build required strings
        header_str = (
            f"# generated by rna3db\n"
            f"#\n"
            f"data_{self.pdb_id}_{author_id}\n"
            f"_entry.id {self.pdb_id}_{author_id}\n"
            f"_pdbx_database_status.recvd_initial_deposition_date {self.release_date}\n"
            # some readers prefer one field over the other, so we write both
            f"_pdbx_audit_revision_history.revision_date {self.release_date}\n"
            f"_exptl.method '{self.structure_method.upper()}'\n"
            f"_reflns.d_resolution_high {self.resolution}\n"
            f"_entity_poly.pdbx_seq_one_letter_code_can {self[author_id].sequence}\n"
        )
        struct_asym_str = Structure._gen_mmcif_loop_str(
            "struct_asym",
            [
                "id",
                "pdbx_blank_PDB_chainid_flag",
                "pdbx_modified",
                "entity_id",
                "details",
            ],
            [(author_id, "N", "N", 1, "?")],
        )
        chem_comp_str = Structure._gen_mmcif_loop_str(
            "chem_comp",
            [
                "id",
                "type",
                "mon_nstd_flag",
                "pdbx_synonyms",
                "formula",
                "formula_weight",
            ],
            [
                (
                    comp["id"],
                    "'RNA linking'",
                    "y",
                    f'"{comp["name"]}"',
                    "?",
                    f"'{comp['formula']}'" if comp["formula"] != "?" else "?",
                    comp["weight"],
                )
                for comp in load_chem_comps()
            ],
        )
        entity_poly = Structure._gen_mmcif_loop_str(
            "entity_poly",
            [
                "entity_id",
                "type",
            ],
            [(1, "polyribonucleotide")],
        )

        entity_poly_seq_str = Structure._gen_mmcif_loop_str(
            "entity_poly_seq",
            [
                "entity_id",
                "num",
                "mon_id",
                "heter",
            ],
            entity_poly_seq_data,
        )
        atom_site_str = Structure._gen_mmcif_loop_str(
            "atom_site",
            [
                "group_PDB",
                "id",
                "type_symbol",
                "label_atom_id",
                "label_alt_id",
                "label_comp_id",
                "label_asym_id",
                "label_entity_id",
                "label_seq_id",
                "pdbx_PDB_ins_code",
                "Cartn_x",
                "Cartn_y",
                "Cartn_z",
                "occupancy",
                "B_iso_or_equiv",
                "pdbx_formal_charge",
                "auth_seq_id",
                "auth_comp_id",
                "auth_asym_id",
                "auth_atom_id",
                "pdbx_PDB_model_num",
            ],
            atom_site_data,
        )

        # write to file
        with open(output_path, "w") as f:
            f.write(header_str)
            f.write(struct_asym_str)
            f.write(chem_comp_str)
            f.write(entity_poly)
            f.write(entity_poly_seq_str)
            f.write(atom_site_str)


class mmCIFParser:
    """Low-level parser for mmCIF/PDBx files. Wraps BioPython's MMCIF2Dict."""

    def __init__(
        self,
        path: Path,
        modification_handler: ModificationHandler,
        nmr_resolution: float = None,
        include_atoms: bool = False,
    ):
        """
        Args:
            path (Path): Path to the mmCIF file.
            modification_handler (ModificationHandler): Handler for converting
                three-letter residue codes to one-letter codes.
            nmr_resolution (float, optional): Resolution to assign to NMR structures.
            include_atoms (bool, optional): Whether to parse atom coordinates.
        """
        self.path = path
        self.nmr_resolution = nmr_resolution
        self.include_atoms = include_atoms

        self.letters_3to1 = lambda x: modification_handler.rna_letters_3to1(x)

        self.parsed_info = PDB.MMCIF2Dict.MMCIF2Dict(self.path)

    @property
    def pdb_id(self) -> str:
        return self.parsed_info["_entry.id"][0].lower()

    @property
    def release_date(self) -> str:
        # prefer to get the date from the earliest revision date
        if "_pdbx_audit_revision_history.revision_date" in self.parsed_info:
            return min(self.parsed_info["_pdbx_audit_revision_history.revision_date"])
        # use deposition date if there are no revisions
        return min(
            self.parsed_info["_pdbx_database_status.recvd_initial_deposition_date"]
        )

    @property
    def resolution(self) -> float:
        resolutions = []
        for res_key in [
            "_refine.ls_d_res_high",
            "_em_3d_reconstruction.resolution",
            "_reflns.d_resolution_high",
        ]:
            if res_key in self.parsed_info:
                if self.parsed_info[res_key][0] not in ".?":
                    resolutions.append(float(self.parsed_info[res_key][0]))

        # if we have an NMR structure and we overwrite default NMR resolution
        if "solution nmr" in self.structure_method and self.nmr_resolution is not None:
            return self.nmr_resolution

        if len(resolutions) == 0:
            return float("inf")

        return max(resolutions)

    @property
    def structure_method(self) -> str:
        return ",".join(self.parsed_info["_exptl.method"]).lower()

    @dataclasses.dataclass
    class _AtomSite:
        atom_id: str
        three_letter_code: str
        author_chain_id: str
        entity_id: str
        author_seq_num: str
        mmcif_seq_num: str
        insertion_code: str
        hetatm_atom: str
        alt_id: str
        x: str
        y: str
        z: str

    @staticmethod
    def _get_atom_sites(parsed_info: PDB.MMCIF2Dict) -> list:
        # fmt: off
        return [
            mmCIFParser._AtomSite(*site)
            for site in zip(
                parsed_info["_atom_site.label_atom_id"],     # atom name
                parsed_info["_atom_site.label_comp_id"],     # residue name
                parsed_info["_atom_site.auth_asym_id"],      # author chain
                parsed_info["_atom_site.label_entity_id"],   # entity id
                parsed_info["_atom_site.auth_seq_id"],       # author_seq_num
                parsed_info["_atom_site.label_seq_id"],      # mmcif_seq_num
                parsed_info["_atom_site.pdbx_PDB_ins_code"], # insertion code
                parsed_info["_atom_site.group_PDB"],         # hetatm_atom
                parsed_info["_atom_site.label_alt_id"],      # alt conformation id
                parsed_info["_atom_site.Cartn_x"],           # x
                parsed_info["_atom_site.Cartn_y"],           # y
                parsed_info["_atom_site.Cartn_z"],           # z
            )
        ]
        # fmt: on

    @property
    def chains(self) -> Mapping[str, Chain]:
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
        self.auth_asym_to_label_asym = defaultdict(set)
        for author_chain_id, mmcif_chain_id in zip(
            self.parsed_info["_atom_site.auth_asym_id"],
            self.parsed_info["_atom_site.label_asym_id"],
        ):
            k = mmcif_chain_to_entity_id[mmcif_chain_id]
            id_map[k].add(author_chain_id)
            self.auth_asym_to_label_asym[author_chain_id].add(mmcif_chain_id)

        # parse full chains from "seqres"
        chains_full = defaultdict(Chain)
        for entity_id, mon_id, idx in zip(
            self.parsed_info["_entity_poly_seq.entity_id"],
            self.parsed_info["_entity_poly_seq.mon_id"],
            self.parsed_info["_entity_poly_seq.num"],
        ):
            for author_id in id_map[entity_id]:
                chains_full[author_id].author_id = author_id
                chains_full[author_id].add_residue(
                    Residue(
                        three_letter_code=mon_id,
                        one_letter_code=self.letters_3to1(mon_id),
                        index=int(idx) - 1,
                    )
                )

        chain_type = {}
        # we check if we have _entity_poly
        if (
            "_entity_poly.entity_id" in self.parsed_info
            and "_entity_poly.type" in self.parsed_info
        ):
            # get chain/polymer types
            for entity_id, poly_type in zip(
                self.parsed_info["_entity_poly.entity_id"],
                self.parsed_info["_entity_poly.type"],
            ):
                for author_id in id_map[entity_id]:
                    chain_type[author_id] = poly_type
        else:
            # if we don't have _entity_poly, fall back to chem_comp type for
            # each mon_id (backwards compatibility with older RNA3DB release mmCIFs)
            chem_comp_type = {
                mon_id: comp_type
                for mon_id, comp_type in zip(
                    self.parsed_info["_chem_comp.id"],
                    self.parsed_info["_chem_comp.type"],
                )
            }
            for author_id, chain_data in chains_full.items():
                # "keep" only chains that contain at least one RNA residue
                if any(
                    [
                        "RNA" in chem_comp_type[i.three_letter_code]
                        for i in chain_data.residues
                    ]
                ):
                    chain_type[author_id] = "polyribonucleotide"
                else:
                    # we just set to "other" if not an RNA
                    chain_type[author_id] = "other"

        # keep only chains of the appropriate polymer type
        chains = {}
        for author_chain_id, chain_data in chains_full.items():
            if "polyribonucleotide" in chain_type[author_chain_id]:
                chains[author_chain_id] = chain_data

        # find starting index of relevant chains
        seq_start_num = {
            author_chain_id: min([res.index for res in chain])
            for author_chain_id, chain in chains.items()
        }

        # keep track of alt conformations so we only take one
        chain_alt_id = {}

        # iterate through atom sites to get coordinates if required
        if self.include_atoms:
            for site in mmCIFParser._get_atom_sites(self.parsed_info):
                # if not a relevant chain, don't care about atoms
                if site.author_chain_id not in chains:
                    continue

                # TODO: verify this is fine
                # I think these are just for H_* and W?
                if site.mmcif_seq_num == ".":
                    continue

                # make sure we always take the same alt conformation
                if site.author_chain_id not in chain_alt_id:
                    chain_alt_id[site.author_chain_id] = site.alt_id
                if site.alt_id != chain_alt_id[site.author_chain_id]:
                    continue

                # the idx is just the mmcif_seq_num - seq_start_num
                # NOTE: we zero index everything
                seq_idx = (
                    int(site.mmcif_seq_num) - seq_start_num[site.author_chain_id] - 1
                )

                # make sure that the sites actually match, should never be a mismatch
                if (
                    site.three_letter_code
                    != chains[site.author_chain_id][seq_idx].three_letter_code
                ):
                    expected = chains[site.author_chain_id][seq_idx].three_letter_code
                    print(
                        f"WARNING: found a mismatch in "
                        f"{self.pdb_id}_{site.author_chain_id} "
                        f"({site.entity_id}) at position {seq_idx} "
                        f"(entity_poly_seq: {expected}, "
                        f"atom_site: {site.three_letter_code})"
                    )
                    continue

                # add atom coordinates
                chains[site.author_chain_id][seq_idx].atoms[site.atom_id] = tuple(
                    map(float, (site.x, site.y, site.z))
                )

        return chains
