import os
import re
import subprocess
from Bio.PDB import PDBParser, Selection, PDBIO
from Bio.PDB.PDBIO import Select
import numpy as np
import sys
sys.path.append('./alphafold')
from alphafold.data import parsers, templates, mmcif_parsing, pipeline
import pickle
import string


three_to_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P',
                'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R',
                'TRP': 'W', 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

atom_types = ['N', 'CA', 'C', 'CB', 'O', 'CG', 'CG1', 'CG2', 'OG', 'OG1', 'SG', 'CD', 'CD1', 'CD2',
              'ND1', 'ND2', 'OD1', 'OD2', 'SD', 'CE', 'CE1', 'CE2', 'CE3', 'NE', 'NE1', 'NE2', 'OE1',
              'OE2', 'CH2', 'NH1', 'NH2', 'OH', 'CZ', 'CZ2', 'CZ3', 'NZ', 'OXT']

ID_TO_HHBLITS_AA = {0: 'A', 1: 'C', 2: 'D', 3: 'E', 4: 'F', 5: 'G', 6: 'H',
                    7: 'I', 8: 'K', 9: 'L', 10: 'M', 11: 'N', 12: 'P', 13: 'Q',
                    14: 'R', 15: 'S', 16: 'T', 17: 'V', 18: 'W', 19: 'Y',
                    20: 'X', 21: '-'}

AA_TO_ID_TO_HHBLITS = {v: k for k, v in ID_TO_HHBLITS_AA.items()}


def read_remark_350(pdb_path):
    pdb_id = pdb_path

    with open(pdb_id, 'r') as f:
        pdb_text = f.read()
        match_biomt1 = [m.start() for m in re.finditer(r'REMARK 350   BIOMT1', pdb_text)]
        match_biomt3 = [m.end() for m in re.finditer(r'REMARK 350   BIOMT3', pdb_text)]
        end_remark_350_block = [m.start() for m in re.finditer('\n', pdb_text[match_biomt3[-1]:])]
        transformation_blocks_indices = match_biomt1 + [match_biomt3[-1] + end_remark_350_block[0] + 1]

        transformations_list = []
        for index in range(len(transformation_blocks_indices) - 1):
            block = pdb_text[transformation_blocks_indices[index]:transformation_blocks_indices[index + 1]]
            matrix = [item.split()[4:8] for item in block.split('\n')[:-1]]
            r11, r12, r13, t1 = matrix[0]
            r21, r22, r23, t2 = matrix[1]
            r31, r32, r33, t3 = matrix[2]
            transformation = [[r11, r12, r13], [r21, r22, r23], [r31, r32, r33], [t1, t2, t3]]
            transformations_list.append(transformation)

    return transformations_list


def rot_and_trans(pdb_path, rot_tra_matrix):

    r11, r12, r13 = rot_tra_matrix[0]
    r21, r22, r23 = rot_tra_matrix[1]
    r31, r32, r33 = rot_tra_matrix[2]
    t1, t2, t3 = rot_tra_matrix[3]

    with open('pdbset.sh', 'w') as f:
        f.write(f'pdbset xyzin {pdb_path} xyzout {pdb_path.split("/")[-1].lower()} << eof\n')
        f.write(
            f'rotate {float(r11)} {float(r12)} {float(r13)} {float(r21)} {float(r22)} {float(r23)} {float(r31)} {float(r32)} {float(r33)}\n')
        f.write(f'shift {float(t1)} {float(t2)} {float(t3)}\n')
        f.write('end\n')
        f.write('eof')

    os.system(f'bash pdbset.sh')
    os.system(f'rm pdbset.sh')


def renumber_residues(structure):

    for num, res in enumerate([res for res in Selection.unfold_entities(structure, 'R')]):
        res.id = (' ', num + 1, ' ')

    return structure


def merge_chains(pdb):

    output_file = open('out.pdb', 'w')

    with open(pdb, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                output_file.write(line[:21] + 'A' + line[22:])
    output_file.close()


def remove_hydrogens(pdb):

    output_file = open('output.pdb', 'w')

    counter = 0
    with open(f'{pdb}', 'r') as f:
        for line in f.readlines():
            if line.split()[-1] in ['N', 'C', 'O', 'S']:
                counter = counter + 1
                output_file.write(line[:6] + str(counter).rjust(5) + line[11:])


def assamble(af_models_dir, pdb_id, assembly_path, pdb_db):

    isExist = os.path.exists(assembly_path)
    if not isExist:
        os.makedirs(assembly_path)
    os.chdir(assembly_path)
    os.system(f'cp {pdb_db}/{pdb_id.split("-")[-1]}.pdb {assembly_path}')
    transformations_list = read_remark_350(f'{pdb_id.split("-")[-1]}.pdb')
    chains, list_of_structures = [], []
    list_in_path = os.listdir(af_models_dir)
    for folder in list_in_path:
        if os.path.isdir(af_models_dir + '/' + folder) and folder[:-1] == pdb_id:
            chain = folder[-1]
            chains.append(chain)
            os.system(f'cp {af_models_dir}/{pdb_id}{chain}/summary/templates/BEST.pdb {assembly_path}/{pdb_id}{chain}.pdb')
            superpose_output = subprocess.Popen(
                ['superpose', f'{pdb_id}{chain}.pdb', '-s', '-all', f'{pdb_id.split("-")[-1]}.pdb', '-s', f'{chain}',
                 '-o', f'{chain}.pdb'], stdout=subprocess.PIPE).communicate()
            list_of_structures.append(PDBParser(QUIET=1).get_structure(f'{chain}', f'{chain}.pdb'))

    if len(transformations_list) == 1:
        for num, model in enumerate(list_of_structures):
            if num == 0:
                chain = list(model[0])[0]
                chain.id = chains[num]
                main_structure = model
            else:
                chain = list(model[0])[0]
                chain.id = chains[num]
                chain.detach_parent()
                main_structure[0].add(chain)
        main_structure = renumber_residues(main_structure)
        io = PDBIO()
        io.set_structure(main_structure)
        io.save(f'{pdb_id}_{"".join(chains)}.pdb')
        merge_chains(f'{pdb_id}_{"".join(chains)}.pdb')
        os.system(f'mv out.pdb template.pdb')
        remove_hydrogens('template.pdb')
        os.system('mv output.pdb template.pdb')
        for chain in chains:
            os.system(f'rm {chain}.pdb')
            os.system(f'rm {pdb_id}{chain}.pdb')
    else:
        for transformation in transformations_list[1:]:
            for chain in chains:
                rot_and_trans(f'{chain}.pdb', transformation)
            chains = chains + [ch.lower() for ch in chains]
            list_of_structures = [PDBParser(QUIET=1).get_structure(f'{ch}', f'{ch}.pdb') for ch in chains]
            for num, model in enumerate(list_of_structures):
                if num == 0:
                    chain = list(model[0])[0]
                    chain.id = chains[num]
                    main_structure = model
                else:
                    chain = list(model[0])[0]
                    chain.id = chains[num]
                    chain.detach_parent()
                    main_structure[0].add(chain)

        main_structure = renumber_residues(main_structure)
        io = PDBIO()
        io.set_structure(main_structure)
        io.save(f'{pdb_id}_{"".join(chains)}.pdb')
        merge_chains(f'{pdb_id}_{"".join(chains)}.pdb')
        # os.system(f'mv out.pdb template.pdb')
        # remove_hydrogens('template.pdb')
        # os.system('mv output.pdb template.pdb')
        for chain in chains:
            os.system(f'rm {chain}.pdb')
        for chain in chains[:len(chains) // 2]:
            os.system(f'rm {pdb_id}{chain}.pdb')

    template_with_chains = f'{assembly_path}/{pdb_id}_{"".join(chains)}.pdb'
    structure = PDBParser(QUIET=1).get_structure(pdb_id, template_with_chains)
    template_ch_list = [res for res in Selection.unfold_entities(structure, "C")]
    res_to_link = []
    for ch in template_ch_list:
        first_res, last_res = list(ch)[0], list(ch)[-1]
        res_to_link.append([first_res, last_res])
    num_gly_in_chain_list = []
    for num, item in enumerate(res_to_link):
        try:
            dist_between_res = np.linalg.norm(item[1]['CA'].coord - res_to_link[num + 1][0]['CA'].coord)
            number_of_glycines_to_add = 5 * (dist_between_res // 3.80)
            num_gly_in_chain_list.append(int(number_of_glycines_to_add))
        except:
            pass
    num_gly_in_chain_list.insert(0, 0)
    counter = 0
    offset_list = []
    for item in num_gly_in_chain_list:
        counter = counter + item
        offset_list.append(counter)

    for num, ch in enumerate(template_ch_list):
        for res in [res for res in list(ch) if res.get_parent().id == ch.id][::-1]:
            res.id = (' ', res.id[1] + offset_list[num], ' ')

    io = PDBIO()
    io.set_structure(structure)
    io.save(f'temp.pdb')

    merge_chains('temp.pdb')
    remove_hydrogens('out.pdb')
    os.system('mv output.pdb template.pdb')
    os.system('rm out.pdb')
    os.system('rm temp.pdb')


def split_chains_and_return_to_initial_res_numbering(best_pdb_path, template_pdb_path, reference_pdb):

    ref_structure = PDBParser(QUIET=1).get_structure('reference', reference_pdb)
    number_of_chains = len(list(ref_structure[0]))
    monomer_seq_length = len(list(list(ref_structure[0])[0]))

    chain_ids = string.ascii_uppercase[:number_of_chains]


    template = PDBParser(QUIET=1).get_structure('template', template_pdb_path)
    template_res_ids_list = [res.id[1] for res in Selection.unfold_entities(template, "R")]

    best = PDBParser(QUIET=1).get_structure('BEST', best_pdb_path)
    best_res_ids_to_remove = [res.id for res in Selection.unfold_entities(best, "R") if res.id[1] not in template_res_ids_list]

    for chain in best[0]:
        for id in best_res_ids_to_remove:
            chain.detach_child(id)
    io = PDBIO()
    io.set_structure(best)
    io.save(f'temp.pdb')

    best_res_ids = [res.id[1] for res in Selection.unfold_entities(best, "R")]
    best_res = [res for res in Selection.unfold_entities(best, "R")]



    to_match, to_remplace = [], []
    for item in range(0, number_of_chains):
        counter = 0
        for res in best_res[item * monomer_seq_length: item * monomer_seq_length + monomer_seq_length]:
            counter += 1
            to_match.append(res.resname + ' A' + str(res.id[1]).rjust(4))
            to_remplace.append(res.resname + chain_ids[item].rjust(2) + str(counter).rjust(4))

    with open('output.pdb', 'w') as out_f:
        with open('temp.pdb', 'r') as f:
            for line in f.readlines():
                if line[17:26] in to_match:
                    out_f.write(line[:17] + to_remplace[to_match.index(line[17:26])] + line[26:])
    out_f.close()
    os.system('rm temp.pdb')


def run_pisa(pdb_path, out_dir):


    pisa_analyse_output = subprocess.Popen(['pisa', 'temp', '-analyse', f'{pdb_path}', '--lig-exclude=all'], stdout=subprocess.PIPE).communicate()[0]
    pisa_interfaces_list_output = subprocess.Popen(['pisa', 'temp', '-list', 'interfaces'], stdout=subprocess.PIPE).communicate()[0]
    print(pisa_analyse_output.decode('utf-8'))

    # print(pisa_analyse_output.decode('utf-8'))
    # print(pisa_interfaces_list_output.decode('utf-8'))

    for num, line in enumerate(pisa_interfaces_list_output.decode('utf-8').split('\n')):
        if '## Id' in line:
            num1 = num
        if line.startswith(' ##:  serial number'):
            num2 = num

    list_of_interfaces = []
    for line in pisa_interfaces_list_output.decode('utf-8').split('\n')[num1+2:num2-1]:
        interface_num = line.split()[0]
        list_of_interfaces.append(interface_num)
    # print(list_of_interfaces)

    for serial_no in list_of_interfaces:
        pisa_interfaces_details_output = subprocess.Popen(['pisa', 'temp', '-detail', 'interfaces', f'{serial_no}'], stdout=subprocess.PIPE).communicate()[0].decode('utf-8')


        with open(f'{out_dir}/interface_{serial_no}_output.log','w') as f:
            f.write(pisa_interfaces_details_output)
        f.close()


        identification_list, energy_list = [], []
        num1_list, num2_list = [], []
        for num, line in enumerate(pisa_interfaces_details_output.split('\n')):
            if 'Interfacing Residues: Structure' in line:
                num1 = num + 4
                num1_list.append(num1)
                # print(num, line)
            if line == " -----'-'------------'--'----------------------":
                num2 = num
                num2_list.append(num2)
        lines_to_show = zip(num1_list, num2_list)
        for item in lines_to_show:
            for line in pisa_interfaces_details_output.split('\n')[item[0]:item[1]]:
                # print(line)
                chain = line[10:11]
                res_id = line[12:15]
                res_num = line[15:20].replace(" ", "")
                energy = line[39:].replace(" ", "")

                identification_list.append(res_id + ' ' + chain + res_num.rjust(4))
                energy_list.append(energy.rjust(6))

        with open(f'{out_dir}/{serial_no}_interface.pdb', 'w') as f:
            for line in open(f'{pdb_path}', 'r').readlines():
                #print(line)
                if line[17:26] in identification_list:
                    match_index = identification_list.index(line[17:26])
                    if float(energy_list[match_index]) != 0:
                        # print(line)
                        f.write(line)
