#!/usr/bin/env python
# -*- coding: utf-8 -*-
from rna3db.utils import read_json
from pathlib import Path

import argparse
import os

ATOMS = {'G' : "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9 C8 N7 C5 C6 O6 N1 C2 N2 N3 C4".split(), # 23
         'A' : "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9 C8 N7 C5 C6 N6 N1 C2 N3 C4".split(),    # 22
         'U' : "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1 C2 O2 N3 C4 O4 C5 C6".split(),          # 20
         'C' : "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1 C2 O2 N3 C4 N4 C5 C6".split(),          # 20
         'N' : "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1'".split(),  
         'T' : "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2'     C1' N1 C2 O2 N3 C4 N4 C5 C6 C7".split(),       # 20
    }

def format_line(atom_index,
                atom_name,
                residue_name,
                chain_id,
                residue_index,
                x,
                y,
                z,
                verbose = False ):
                                group_PDB = "ATOM"
                                atom_id = atom_index
                                type_symbol = atom_name[0]
                                label_atom_id = atom_name
                                label_alt_id = "."
                                label_comp_id = residue_name
                                label_asym_id = chain_id
                                label_entity_id = "?"
                                label_seq_id = residue_index 
                                pdbx_PDB_ins_code = "?"
                                Cartn_x = x
                                Cartn_y = y
                                Cartn_z = z
                                occupancy = 1.0
                                B_iso_or_equiv = 0.0
                                auth_seq_id = residue_index 
                                auth_asym_id = chain_id
                                pdbx_PDB_model_num = "1"

                                mmcif_format = (
                                    f"{group_PDB:<6}"
                                    f"{atom_id:>7} "
                                    f"{type_symbol:>2} "
                                    f"{label_atom_id:<4} "
                                    f"{label_alt_id:<1} "
                                    f"{label_comp_id:>3} "
                                    f"{label_asym_id:>2} "
                                    f"{label_entity_id:>1} "
                                    f"{label_seq_id:>3} "
                                    f"{pdbx_PDB_ins_code:>2} "
                                    f"{Cartn_x:>8.3f} "
                                    f"{Cartn_y:>8.3f} "
                                    f"{Cartn_z:>8.3f} "
                                    f"{occupancy} "
                                    f"{B_iso_or_equiv:>6.2f} "
                                    f"? "
                                    f"{auth_seq_id:>4} "
                                    f"{label_comp_id:>3} "
                                    f"{label_asym_id:>2} "
                                    f'{label_atom_id:>4}'
                                    f"{pdbx_PDB_model_num:>2}"
                                )

                                # Now we replace the variables in the mmCIF format
                                formatted_mmcif_line = mmcif_format.format(
                                    group_PDB=group_PDB,
                                    atom_id=atom_id,
                                    type_symbol=type_symbol,
                                    label_atom_id=label_atom_id,
                                    label_alt_id=label_alt_id,
                                    label_comp_id=label_comp_id,
                                    label_asym_id=label_asym_id,
                                    label_entity_id=label_entity_id,
                                    label_seq_id=label_seq_id,
                                    pdbx_PDB_ins_code=pdbx_PDB_ins_code,
                                    Cartn_x=Cartn_x,
                                    Cartn_y=Cartn_y,
                                    Cartn_z=Cartn_z,
                                    occupancy=occupancy,
                                    B_iso_or_equiv=B_iso_or_equiv,
                                    auth_seq_id=auth_seq_id,
                                    auth_asym_id=auth_asym_id,
                                    pdbx_PDB_model_num=pdbx_PDB_model_num
                                )

                                if verbose: print(formatted_mmcif_line)
                                return formatted_mmcif_line + '\n'

def run(data, output_path):
    log = ''
    for dataset, value in data.items():
            dataset_pdb_index = 1
            try :
                os.makedirs(output_path + os.sep + dataset)
            except FileExistsError:
                pass
            for component, value in data[dataset].items():
                for cluster, value in data[dataset][component].items():
                    for pdb, value in data[dataset][component][cluster].items():
                       
                        release_date = value['release_date']
                        resolution = value['resolution']
                        
                        structure_method = ''

                        if 'electron microscopy' in value['structure_method']:
                            structure_method = "'ELECTRON MICROSCOPY'"
                        if 'x-ray diffraction' == value['structure_method']:    
                            structure_method = "'X-RAY DIFFRACTION'"  
                        if not structure_method:
                            print('error: no method')
                            
                        pdb_id, chain_id = pdb.split('_')

                        cif_txt = f"""# generated by rna3db
#
data_{pdb_id}_{chain_id}
_entry.id {pdb_id}_{chain_id}
_pdbx_database_status.recvd_initial_deposition_date {release_date}
_exptl.method {structure_method}
_reflns.d_resolution_high {resolution}
_entity_poly.pdbx_seq_one_letter_code_can   {value['sequence']}
#
loop_
_struct_asym.id 
_struct_asym.pdbx_blank_PDB_chainid_flag 
_struct_asym.pdbx_modified 
_struct_asym.entity_id 
_struct_asym.details 
{chain_id} N N 1 ?
## 
loop_
_chem_comp.id 
_chem_comp.type 
_chem_comp.mon_nstd_flag 
_chem_comp.name 
_chem_comp.pdbx_synonyms 
_chem_comp.formula 
_chem_comp.formula_weight 
A   'RNA linking' y "ADENOSINE-5'-MONOPHOSPHATE" ? 'C10 H14 N5 O7 P'   347.221 
C   'RNA linking' y "CYTIDINE-5'-MONOPHOSPHATE"  ? 'C9 H14 N3 O8 P'    323.197 
G   'RNA linking' y "GUANOSINE-5'-MONOPHOSPHATE" ? 'C10 H14 N5 O8 P'   363.221 
U   'RNA linking' y "URIDINE-5'-MONOPHOSPHATE"   ? 'C9 H13 N2 O9 P'    324.181
T   'RNA linking' y "T"   ? ''    0
N   'RNA linking' y "N"   ? ''    0
"""

                        cif_txt_poly_seq = f"""#
# 
loop_
_entity_poly_seq.entity_id 
_entity_poly_seq.num 
_entity_poly_seq.mon_id 
_entity_poly_seq.hetero
"""
                        cif_txt_atom_site = f"""
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy 
_atom_site.B_iso_or_equiv 
_atom_site.pdbx_formal_charge 
_atom_site.auth_seq_id 
_atom_site.auth_comp_id 
_atom_site.auth_asym_id 
_atom_site.auth_atom_id 
_atom_site.pdbx_PDB_model_num 
"""

                        print(f'# {dataset_pdb_index} >{pdb_id}_{chain_id}')
                        log += f'# {dataset_pdb_index} >{pdb_id}_{chain_id} \n{value["sequence"]} {len(value["sequence"])}\n'
                        dataset_pdb_index += 1

                        seq = value['sequence']

                        atom_index = 1
                        residue_index = 0

                        for residue in value['atoms']:
                            residue_index += 1  # at the beginning of the cycle, so this will not mess up continues/break

                            # collect poly seq per earch residue 
                            cif_txt_poly_seq += f'1 {residue_index} {seq[residue_index - 1]} n\n' # 1 1 G n

                            residue_name = ''
                            # get the residue_names
                            for atom_name, xyz in residue.items():
                                if 'O6' in residue:
                                    residue_name = 'G'
                                    break
                                elif 'N6' in residue:
                                    residue_name = 'A'
                                    break
                                elif 'N4' in residue:
                                    residue_name = 'C'
                                    break
                                elif ('O4' in residue) or ('S4' in residue):
                                    residue_name = 'U'
                                    break
                                elif "C2'": # at least there is a sugar
                                    residue_name = 'N' # for gap!
                            
                            # check the index!
                            # if this is logged then you dont get missing atoms because residue_name is N and I fetch for backbone atoms 
                            # according the the dictionary at the top
                            if residue_name != seq[residue_index - 1]:
                                log += (f"Seq inconsistence {pdb_id}_{chain_id} resi {residue_index} inferResName: {residue_name} != seqres: {seq[residue_index - 1]} ats:{residue}\n")
                            residue_name = seq[residue_index - 1]
                            if False:
                                # fetch atoms for given residue!
                                for atom_name in ATOMS[residue_name]:
                                    if atom_name not in residue:
                                        log += (f'Missing atoms   {pdb_id}_{chain_id} resi {residue_index} {residue_name} {atom_name}\n')
                                        continue 
                                    xyz = ' '.join([str(f) for f in residue[atom_name]])
                                    x, y, z = residue[atom_name]

                                    cif_txt_atom_site += format_line(
                                        atom_index,
                                        atom_name,
                                        seq[residue_index - 1], #residue_name, # which name to use it? seq[residue_index - 1]
                                        chain_id,
                                        residue_index, verbose=False)
                                    atom_index += 1
                            else: # get atoms in json files, all of them
                                for atom_name in ATOMS[residue_name]:
                                    if atom_name not in residue:
                                        log += (f'Missing atoms   {pdb_id}_{chain_id} resi {residue_index} {residue_name} {atom_name}\n')
                                        continue 
                                        
                                for atom_name, xyz in residue.items():
                                    #xyz = ' '.join([str(f) for f in residue[atom_name]])
                                    x, y, z = xyz

                                    cif_txt_atom_site += format_line(
                                        atom_index,
                                        atom_name,
                                        seq[residue_index - 1], #residue_name, # which name to use it? seq[residue_index - 1]
                                        chain_id,
                                        residue_index, 
                                        x, y, z ,verbose=False)
                                    atom_index += 1                               
                                
                        directory = f'{output_path}/{dataset}/{component}/{cluster}'
                        if not os.path.exists(directory):
                          os.makedirs(directory)
                        fn = output_path + f'/{dataset}/{component}/{cluster}/{pdb_id}_{chain_id}.cif'
                        with open(fn, "w") as f:
                            print(f'save {fn}')
                            f.write(cif_txt)
                            f.write(cif_txt_poly_seq.strip())
                            f.write(cif_txt_atom_site)

    with open(f'{output_path}/rna3db.log', "w") as f:
        f.write(log)        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Converts a JSON produced by RNA3DB's parse command to a set of PDBx/mmCIF files."
    )
    parser.add_argument("input_json_file", type=Path)
    parser.add_argument("output_folder", type=Path)

    args = parser.parse_args()

    print('loading json file ...', args.input_json_file)
    data = read_json(args.input_json_file)
    output_path = str(args.output_folder)

    run(data, output_path)
