import argparse
import logging
import json

import sys

sys.path.append(".")  # an innocent hack to get this to run from the top level

from Bio.Data import SCOPData

VALID_RNA_CODES = set("ACGUT")

# a few values are adjusted manually to make sure we are as correct as possible
# and to avoid inconsistency with AlphaFold
manual_fixes = {
    # SCOPData for a few protein cases to avoid inconsistency with AlphaFold
    **{
        k: SCOPData.protein_letters_3to1.get(k, "X")
        for k in ["0AZ", "DNM", "OTY", "TOX", "DNG", "PCA", "AEI", "PDU"]
    },
    # some selenocysteines are encoded incorrectly in SCOPData
    # we keep SCOPData's version for consistency with AlphaFold
    **{k: SCOPData.protein_letters_3to1.get(k, "X") for k in ["PSW", "SOC", "SYS"]},
    # this is from Bio.Data.PDBData from an updated version of biopython, hardcoded to avoid
    # having to upgrade the whole package
    **{"F2Y": "Y", "QM8": "L", "PCA": "Q", "4BF": "F"},
}


def parse_components(path):
    chem_comp_id, cif_lines = None, []
    with open(path) as fp:
        for line in fp:
            if line.startswith("data_"):
                if len(cif_lines) > 0:
                    yield (chem_comp_id, "".join(cif_lines))
                chem_comp_id = line.split("_")[-1].rstrip()
                cif_lines = []
            cif_lines.append(line)
        yield (chem_comp_id, "".join(cif_lines))


def parse_cif_string(cif_string):
    for line in cif_string.split("\n"):
        if line.startswith("_chem_comp.id"):
            comp_id = line.split()[-1]

        if line.startswith("_chem_comp.type"):
            comp_type = " ".join(line.split()[1:])

        if line.startswith("_chem_comp.one_letter_code"):
            comp_one_letter_code = line.split()[-1]

        if line.startswith("_chem_comp.mon_nstd_parent_comp_id"):
            line = "".join(c for c in line.split()[-1] if c.isalnum() or c == "?")
            if line == "":
                line = "?"
            parent_comp_id = line[:3].upper()

    return comp_id, comp_type, comp_one_letter_code, parent_comp_id


def parse_cif(cif_string):
    comp_id, comp_type, comp_one_letter_code, parent_comp_id = parse_cif_string(
        cif_string
    )

    # if we have parent, parse that instead
    # if parent is itself, don't re-parse
    if parent_comp_id != "?" and comp_id != parent_comp_id:
        return parse_cif(cif_strings[parent_comp_id])

    # if our one_letter_code is too long, try to parse that
    if len(comp_one_letter_code) > 1:
        return parse_cif(cif_strings[comp_one_letter_code])

    return comp_one_letter_code, comp_type


cif_strings = {}


def main(args):
    data = {"protein": {}, "rna": {}}

    for k, v in parse_components(args.components_path):
        cif_strings[k] = v
    comp_types = {}

    modifications = {}
    for k, cif_string in cif_strings.items():
        one_letter, comp_type = parse_cif(cif_string)

        if "peptide" in comp_type.lower():
            comp_type = "protein"
        elif "RNA" in comp_type or "DNA" in comp_type:
            comp_type = "rna"
        else:
            comp_type = comp_type.lower()
        comp_types[k] = comp_type

        if one_letter != "?":
            modifications[k] = one_letter

    modifications = {**modifications, **manual_fixes}

    for k, v in modifications.items():
        if comp_types[k] == "protein":
            data["protein"][k] = v
        if comp_types[k] == "rna":
            if v in VALID_RNA_CODES:
                data["rna"][k] = v

    with open(args.output_path, "w") as fp:
        fp.write(json.dumps(data, indent=4))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "components_path",
        type=str,
        help="Path for Chemical Component Dictionary .cif file",
    )
    parser.add_argument("output_path", type=str, help="Path for .json output")

    args = parser.parse_args()

    main(args)
