import os
import pandas as pd
import matplotlib.pyplot as plt
import pickle5 as pickle
import numpy as np
import subprocess
import re
from Bio.PDB import PDBParser, Selection


ID_TO_HHBLITS_AA_3LETTER_CODE = {0: 'ALA', 1: 'CYS', 2: 'ASP', 3: 'GLU', 4: 'PHE', 5: 'GLY', 6: 'HIS',
                                 7: 'ILE', 8: 'LYS', 9: 'LEU', 10: 'MET', 11: 'ASN', 12: 'PRO', 13: 'GLN',
                                 14: 'ARG', 15: 'SER', 16: 'THR', 17: 'VAL', 18: 'TRP', 19: 'TYR', 20: 'X',
                                 21: '-'}
atom_types = ['N', 'CA', 'C', 'CB', 'O', 'CG', 'CG1', 'CG2', 'OG', 'OG1', 'SG', 'CD', 'CD1', 'CD2',
              'ND1', 'ND2', 'OD1', 'OD2', 'SD', 'CE', 'CE1', 'CE2', 'CE3', 'NE', 'NE1', 'NE2', 'OE1',
              'OE2', 'CH2', 'NH1', 'NH2', 'OH', 'CZ', 'CZ2', 'CZ3', 'NZ', 'OXT']

atom_order = {atom_type: i for i, atom_type in enumerate(atom_types)}
order_atom = {v: k for k, v in atom_order.items()}


def superpose_pdbs(query_pdb, target_pdb, output_superposition=True):

    # WARNING: this function is only for PDBs containing only one chain and has to be executed in the same
    # query_pdb and target_pdb path

    if output_superposition:
        superpose_output = subprocess.Popen(['superpose', f'{query_pdb}', '-s', '-all', f'{target_pdb}', '-s', '-all',
                                             '-o', f'{query_pdb[:-4]}_superposed.pdb'],
                                            stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
    else:
        superpose_output = subprocess.Popen(['superpose', f'{query_pdb}', '-s', '-all', f'{target_pdb}', '-s', '-all'],
                                            stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
    for line in superpose_output.split('\n'):
        if 'r.m.s.d:' in line:
            rmsd = float(line.split()[1])
        if 'quality Q:' in line:
            quality_q = line.split()[2]
        if 'Nalign:' in line:
            nalign = line.split()[1]

    try:
        match1 = [m.start() for m in re.finditer("TEXT:Residue alignment:", superpose_output)][0]
        match2 = [m.start() for m in re.finditer("`-------------'----------'-------------'", superpose_output)][0]
        alignment_output = superpose_output[match1:match2].split('\n')[5:]
        aligned_res_list = []
        for line in alignment_output:
            if line[23:25] == '**':
                aligned_res_list.append(int(line[36:39].replace(' ', '')))
    except:

        rmsd, nalign, quality_q, aligned_res_list = None, None, None, None


    return rmsd, nalign, quality_q, aligned_res_list


def write_all_templates_in_features_pkl(features, output_path):

    for i, pdb_name in enumerate(features['template_domain_names']):
        pdb = pdb_name.decode('utf-8')
        output_pdb = open(output_path + '/' + pdb + '.pdb', 'w')
        template_domain_index = np.where(features['template_domain_names'] == pdb_name)[0][0]
        true_seq = ''
        atom_num_int = 0  # AQUI
        for index, atoms_mask in enumerate(features['template_all_atom_masks'][template_domain_index][:]):
            template_residue_masks = features['template_aatype'][template_domain_index][index]
            template_residue_masks_index = np.where(template_residue_masks == 1)[0][0]
            res_type = ID_TO_HHBLITS_AA_3LETTER_CODE[template_residue_masks_index]
            list_of_atoms_in_residue = [order_atom[i] for i, atom in enumerate(atoms_mask) if atom == 1]
            for atom in list_of_atoms_in_residue:
                atom_num_int = atom_num_int + 1
                atom_remark = 'ATOM'.ljust(6)
                atom_num = str(atom_num_int).rjust(5)
                atom_name = atom.ljust(4)
                res_name = res_type.ljust(3)
                chain_id = 'A'.rjust(1)
                res_num = str(index + 1).rjust(4)
                x_coord = str('%8.3f' % (float(str(
                    features['template_all_atom_positions'][template_domain_index][index][atom_types.index(atom)][
                        0])))).rjust(8)
                y_coord = str('%8.3f' % (float(str(
                    features['template_all_atom_positions'][template_domain_index][index][atom_types.index(atom)][
                        1])))).rjust(8)
                z_coord = str('%8.3f' % (float(str(
                    features['template_all_atom_positions'][template_domain_index][index][atom_types.index(atom)][
                        2])))).rjust(8)
                occ = str('%6.2f' % (float('1.0'))).rjust(6)
                bfact = str('%6.2f' % (float('25.0'))).ljust(6)
                atom_type = atom[0].rjust(12)
                output_pdb.write(
                    f'{atom_remark}{atom_num}  {atom_name}{res_name} {chain_id}{res_num}    {x_coord}{y_coord}{z_coord}{occ}{bfact}{atom_type}\n')

        output_pdb.close()


def plot_plddt(ranked_models_path):
    plt.clf()
    ranked_models_list = sorted([ranked for ranked in os.listdir(ranked_models_path) if ranked[:6] == 'ranked'])
    for ranked in ranked_models_list:
        plddt_list = []
        with open(ranked_models_path + '/' + ranked, 'r') as f:
            for line in f.readlines():
                if line[:4] == 'ATOM' and line[13:16] == 'CA ':
                    plddt_list.append(float(line[60:66].replace(" ", "")))
        res_list = [int(item) for item in range(1, len(plddt_list) + 1)]
        plt.plot(res_list, plddt_list, label=f'{ranked.split("/")[-1][:-4]}')
    plt.legend()
    plt.xlabel('residue number')
    plt.ylabel('pLLDT')
    plt.savefig(f'{ranked_models_path}/pLDDT.png')



def af_summary(af_job_path):

    isExist = os.path.exists(f'{af_job_path}/summary')
    if not isExist:
        os.makedirs(f'{af_job_path}/summary')

    isExist = os.path.exists(f'{af_job_path}/summary/templates')
    if not isExist:
        os.makedirs(f'{af_job_path}/summary/templates')

    plot_plddt(ranked_models_path=af_job_path)
    os.system(f'mv {af_job_path}/pLDDT.png {af_job_path}/summary')

    features_file = open(f'{af_job_path}/features.pkl', 'rb')
    features = pickle.load(features_file)

    write_all_templates_in_features_pkl(features, output_path=f'{af_job_path}/summary/templates')

    os.system(f'cp {af_job_path}/ranked_*.pdb {af_job_path}/summary/templates')


    results_dict = {}
    templates_list = [f'{template.decode("utf-8")}.pdb' for template in features['template_domain_names']]
    num = 0
    for template in templates_list:

        res_list_length = len([res for res in Selection.unfold_entities(PDBParser().get_structure('test', f'{af_job_path}/summary/templates/{template}'), 'R')])
        results_list = []
        rankeds_list = sorted([file for file in os.listdir(f'{af_job_path}/summary/templates') if 'ranked_' in file], key=lambda x: int("".join([i for i in x if i.isdigit()])))
        for ranked in rankeds_list:
            rmsd, nalign, quality_q, aligned_res_list = superpose_pdbs(query_pdb=f'{af_job_path}/summary/templates/{ranked}',
                                                                       target_pdb=f'{af_job_path}/summary/templates/{template}',
                                                                       output_superposition=False)
            results_list.append(f'{rmsd}, {nalign} ({res_list_length})')
        if template[:-4] in results_dict:
            num = num + 1
            results_dict[f'{template[:-4]}_({num})'] = results_list
        else:
            results_dict[f'{template[:-4]}'] = results_list

    #########

    with open(f'{af_job_path}/summary/templates/RMSD_RANKEDS_TEMPLATES.log', 'w') as f:
        rows = []
        for key in results_dict.keys():
            rows.append([key] + results_dict[key])

        df = pd.DataFrame(rows, columns=['template', 'ranked_0', 'ranked_1', 'ranked_2', 'ranked_3', 'ranked_4'])
        f.write(df.to_markdown())
    f.close()

    query_list = templates_list + [file for file in os.listdir(f'{af_job_path}/summary/templates') if 'ranked_' in file and not 'ranked_0.pdb' in file]
    for query in query_list:
        superpose_pdbs(query_pdb=f'{af_job_path}/summary/templates/{query}',
                       target_pdb=f'{af_job_path}/summary/templates/ranked_0.pdb',
                       output_superposition=True)
        os.system(f'mv {af_job_path}/summary/templates/{query[:-4]}_superposed.pdb {af_job_path}/summary/templates/{query}')





