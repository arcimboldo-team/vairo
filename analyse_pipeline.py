import sys
sys.path.append('./alphafold')
import pandas as pd
from alphafold.data import parsers, templates, mmcif_parsing, pipeline
import re
import os
import pickle5 as pickle
import subprocess
import numpy as np
import json
from scipy.spatial import distance
import matplotlib.pyplot as plt
from Bio.PDB import MMCIFParser, PDBParser, Selection, PDBList, Structure, PDBIO
import logging

ID_TO_HHBLITS_AA_3LETTER_CODE = {0: 'ALA', 1: 'CYS', 2: 'ASP', 3: 'GLU', 4: 'PHE', 5: 'GLY', 6: 'HIS',
                                 7: 'ILE', 8: 'LYS', 9: 'LEU', 10: 'MET', 11: 'ASN', 12: 'PRO', 13: 'GLN',
                                 14: 'ARG', 15: 'SER', 16: 'THR', 17: 'VAL', 18: 'TRP', 19: 'TYR', 20: 'X',
                                 21: '-'}
atom_types = ['N', 'CA', 'C', 'CB', 'O', 'CG', 'CG1', 'CG2', 'OG', 'OG1', 'SG', 'CD', 'CD1', 'CD2',
              'ND1', 'ND2', 'OD1', 'OD2', 'SD', 'CE', 'CE1', 'CE2', 'CE3', 'NE', 'NE1', 'NE2', 'OE1',
              'OE2', 'CH2', 'NH1', 'NH2', 'OH', 'CZ', 'CZ2', 'CZ3', 'NZ', 'OXT']

atom_order = {atom_type: i for i, atom_type in enumerate(atom_types)}
order_atom = {v: k for k, v in atom_order.items()}


def merge_pdbs(list_of_paths_of_pdbs_to_merge, merged_pdb_path):

    with open(merged_pdb_path, 'w') as f:
        counter = 0
        for pdb_path in list_of_paths_of_pdbs_to_merge:
            for line in open(pdb_path,'r').readlines():
                if line[:4] == 'ATOM':
                    counter += 1
                    f.write(line[:4] + str(counter).rjust(7) + line[11:])


def change_chain_and_apply_offset_in_single_chain(pdb_in_path, pdb_out_path, offset=None, chain=None):

    with open('pdbset.sh', 'w') as f:
        f.write(f'pdbset xyzin {pdb_in_path} xyzout {pdb_out_path} << eof\n')
        if offset:
            f.write(f'renumber increment {offset}\n')
        if chain:
            f.write(f'chain {chain}\n')
        f.write('end\n')
        f.write('eof')
    subprocess.Popen(['bash', 'pdbset.sh'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
    os.system(f'rm pdbset.sh')


def download_pdb(pdb_id, out_dir):

    pdbl = PDBList()
    logging.info(f'Downloading {pdb_id} coordinates from PDB in {out_dir}')
    pdbl.retrieve_pdb_file(pdb_id, pdir=f'{out_dir}', file_format='pdb')


def read_remark_350(pdb_path):

    # TO-DO if REMARK 350 is not found then it will run PISA to generate it

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


def rot_and_trans(pdb_path, out_pdb_path, rot_tra_matrix):

    r11, r12, r13 = rot_tra_matrix[0]
    r21, r22, r23 = rot_tra_matrix[1]
    r31, r32, r33 = rot_tra_matrix[2]
    t1, t2, t3 = rot_tra_matrix[3]

    with open('pdbset.sh', 'w') as f:
        f.write(f'pdbset xyzin {pdb_path} xyzout {out_pdb_path} << eof\n')
        f.write(
            f'rotate {float(r11)} {float(r12)} {float(r13)} {float(r21)} {float(r22)} {float(r23)} {float(r31)} {float(r32)} {float(r33)}\n')
        f.write(f'shift {float(t1)} {float(t2)} {float(t3)}\n')
        f.write('end\n')
        f.write('eof')
    subprocess.Popen(['bash', 'pdbset.sh'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
    os.system(f'rm pdbset.sh')



def run_hhsearch(query_sequence):
    with open('./query.fasta', 'w') as f:
        f.write('> query\n')
        f.write(f'{query_sequence}\n')

    out = subprocess.Popen(['hhsearch', '-i', 'query.fasta', '-o', './output.hhr', '-maxseq',
                            '1000000', '-d', '/alpha/alphauser/af2_dbs/pdb70/pdb70'],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = out.communicate()
    out_text = stdout.decode('utf-8')
    os.system('rm query.fasta')

    return out_text

def extract_single_template_features(query_sequence, hhr_path, pdb_id, chain_id):


    hhr_text = open(hhr_path, 'r').read()

    matches = re.finditer(r'No\s+[0-9]+', hhr_text)
    matches_positions = [match.start() for match in matches] + [len(hhr_text)]

    detailed_lines_list = []

    for i in range(len(matches_positions) - 1):
        detailed_lines_list.append(hhr_text[matches_positions[i]:matches_positions[i + 1]].split('\n')[:-3])

    hits_list = [detailed_lines for detailed_lines in detailed_lines_list
                 if detailed_lines[1].split('_')[0][1:] == f'{pdb_id.upper()}']
    detailed_lines = hits_list[0]

    file_id = f'{pdb_id.lower()}'
    hit = parsers._parse_hhr_hit(detailed_lines)

    template_sequence = hit.hit_sequence.replace('-', '')

    mapping = templates._build_query_to_hit_index_mapping(
        hit.query, hit.hit_sequence, hit.indices_hit, hit.indices_query,
        query_sequence)

    mmcif_string = open(f'{mmcif_db}/{file_id}.cif').read()
    parsing_result = mmcif_parsing.parse(file_id=file_id, mmcif_string=mmcif_string)

    template_features, realign_warning = templates._extract_template_features(mmcif_object=parsing_result.mmcif_object, pdb_id=file_id,
                                                             mapping=mapping, template_sequence=template_sequence,
                                                             query_sequence=query_sequence, template_chain_id=chain_id,
                                                             kalign_binary_path='kalign')
    if hit.sum_probs is None:
      template_features['template_sum_probs'] = np.array([[0]])
    else:
      template_features['template_sum_probs'] = np.array([[hit.sum_probs]])

    template_features['template_aatype'] = np.array([template_features['template_aatype']])
    template_features['template_all_atom_masks'] = np.array([template_features['template_all_atom_masks']])
    template_features['template_all_atom_positions'] = np.array([template_features['template_all_atom_positions']])
    template_features['template_domain_names'] = np.array([template_features['template_domain_names']])
    template_features['template_sequence'] = np.array([template_features['template_sequence']])

    return template_features



def features_for_single_template(query_sequence, pdb_id, chain_id, hhr_path=False):

    if not hhr_path:
        run_hhsearch(query_sequence)
        hhr_path = './output.hhr'

    template_features = extract_single_template_features(query_sequence, hhr_path, pdb_id, chain_id)
    template_sequence = template_features['template_sequence'][0].decode('utf-8')

    input_msa = f'>query\n{query_sequence}\n>template\n{template_sequence}'
    custom_msa_result = {'a3m': input_msa}
    custom_msa = parsers.parse_a3m(custom_msa_result['a3m'])
    msa_features = pipeline.make_msa_features_for_custom(custom_msa)

    input_seqs, input_descs = parsers.parse_fasta(f'>query\n{query_sequence}')
    input_sequence = input_seqs[0]
    input_description = input_descs[0]
    sequence_features = pipeline.make_sequence_features(sequence=input_sequence, description=input_description,
                                                        num_res=len(query_sequence))

    features = {**sequence_features, **msa_features, **template_features}

    return features


def write_pkl_from_features(features):

    with open('new_features.pkl', 'wb') as handle:
        pickle.dump(features, handle, protocol=pickle.HIGHEST_PROTOCOL)


def pdist(structure):

    res_list = [res for res in Selection.unfold_entities(structure, 'R')]
    CAs_coords = [res['CA'].coord for res in res_list]
    pdist = distance.pdist(CAs_coords, "euclidean")
    pdist_matrix = distance.squareform(pdist)

    return pdist_matrix


def diff_pdist(pdist1, pdist2):

    diff_pdist_matrices = np.abs(pdist1 - pdist2)

    return diff_pdist_matrices


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
            rmsd = line.split()[1]
        if 'quality Q:' in line:
            quality_q = line.split()[2]
        if 'Nalign:' in line:
            nalign = line.split()[1]

    match1 = [m.start() for m in re.finditer("TEXT:Residue alignment:", superpose_output)][0]
    match2 = [m.start() for m in re.finditer("`-------------'----------'-------------'", superpose_output)][0]

    alignment_output = superpose_output[match1:match2].split('\n')[5:]
    aligned_res_list = []
    for line in alignment_output:
        if line[23:25] == '**':
            aligned_res_list.append(line[36:39])

    return rmsd, nalign, quality_q, aligned_res_list


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
    plt.savefig('pLLDT.png')


def write_all_templates_in_features_pkl(features, output_path):

    for i, pdb_name in enumerate(features['template_domain_names']):
        pdb, chain = pdb_name.decode('utf-8').split('_')
        output_pdb = open(output_path + '/' + pdb + chain + '_template.pdb', 'w')
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
                chain_id = pdb_name.decode('utf-8').split('_')[1].rjust(1)
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

def plot_bfact_and_gaps_in_template(pdb, chain):

    cif_db = '/alpha/alphauser/af2_dbs/pdb_mmcif/mmcif_files/'
    template_num, template_bfact = zip(*[(res.id[1], res['CA'].bfactor) for res in
                                       Selection.unfold_entities(MMCIFParser(QUIET=1).get_structure('{pdb}', f'{cif_db}/{pdb}.cif'), 'R')
                                       if res.id[0] == ' ' and res.get_parent().id == chain])
    continuous_num_list = list(range(template_num[0], template_num[-1] + 1))
    gaps = []
    for res in continuous_num_list:
        if res not in template_num:
            plt.axvline(x = res, color = 'r', alpha=0.3)
            gaps.append(res)
        else:
            pass
    plt.plot(template_num, template_bfact, label=f'chain {chain}')
    plt.ylabel('B-factor ($\AA^2$)')
    plt.xlabel('res. number')
    plt.show()


def rmsd_calc(af_model_path, template_model_path):

    # WARNING: this function is only for models that share the same residue numeration!!

    os.system(f'cp {af_model_path} ./tmp1.pdb')
    os.system(f'cp {template_model_path} ./tmp2.pdb')
    superpose_output = subprocess.Popen(['superpose', 'tmp1.pdb', '-s', '-all', 'tmp2.pdb', '-s', '-all',
                                         '-o', 'tmp3.pdb'], stdout=subprocess.PIPE).communicate()

    af_structure = PDBParser(QUIET=1).get_structure(f' ', './tmp3.pdb')

    template_structure = PDBParser(QUIET=1).get_structure(f' ', template_model_path)
    template_num_list, template_res_list = zip(
        *[(res.id[1], res) for res in Selection.unfold_entities(template_structure, 'R')])

    rmsd_list = []
    for num in template_num_list:
        try:
            af_CA_coord = [res for res in Selection.unfold_entities(af_structure, 'R') if res.id[1] == num][0][
                'CA'].coord
            template_CA_coord = \
            [res for res in Selection.unfold_entities(template_structure, 'R') if res.id[1] == num][0]['CA'].coord
            rmsd_list.append(round(np.linalg.norm(af_CA_coord - template_CA_coord), 2))
        except:
            pass
    rmsd = round(sum(rmsd_list) / len(rmsd_list), 2)
    os.system('rm tmp[0-9].pdb')

    return rmsd


def af_summary(af_job_path):

    cwd = os.getcwd()
    output_dict = {}

    isExist = os.path.exists(f'{af_job_path}/summary')
    if not isExist:
        os.makedirs(f'{af_job_path}/summary')

    # PLDDT PLOT:
    plot_plddt(af_job_path)
    os.system(f'mv pLLDT.png {af_job_path}/summary')

    ranked_models_list = sorted([ranked for ranked in os.listdir(af_job_path) if ranked[:6] == 'ranked'])
    output_dict['ranked_models'] = np.array(ranked_models_list)

    # INTERDISTANCES IN TEMPLATES:
    output_dict['dist_in_templates'] = {}
    output_dict['dist_diff_between_template_and_rankeds'] = {}
    output_dict['rmsd_between_template_and_rankeds'] = {}

    isExist = os.path.exists(f'{af_job_path}/summary/templates')
    if not isExist:
        os.makedirs(f'{af_job_path}/summary/templates')
    
    #features = pd.read_pickle(f'{af_job_path}/features.pkl')
    
    file = open(f'{af_job_path}/features.pkl', 'rb')
    features = pickle.load(file)
    file.close()
    
    write_all_templates_in_features_pkl(features, f'{af_job_path}/summary/templates')

    os.chdir(f'{af_job_path}/summary/templates')
    os.system(f'cp {af_job_path}/ranked_*.pdb {af_job_path}/summary/templates')
    templates_list = [template for template in os.listdir(f'{af_job_path}/summary/templates') if '_template.pdb' in template]
    for template in templates_list:
        template_structure = PDBParser(QUIET=1).get_structure(f'{template[:-4]}', f'{af_job_path}/summary/templates/{template}')
        templat_CAs_list = [res['CA'].coord for res in Selection.unfold_entities(template_structure, "R")]
        template_res_id_list = [res.id[1] for res in Selection.unfold_entities(template_structure, "R")]
        template_pdist = distance.pdist(templat_CAs_list, "euclidean")
        template_pdist_matrix = distance.squareform(template_pdist)
        output_dict['dist_in_templates'][template.split('_')[0]] = np.array(template_pdist_matrix)

        # INTERDISTANCE DIFERENCES BETWEEN TEMPLATE AND RANKED MODELS:
        diff_pdist_matrices_list = []
        rmsd_between_template_and_ranked_list, nalign_between_template_and_ranked_list = [], []
        for num, model in enumerate(ranked_models_list):
            rmsd, nalign, quality_q, aligned_res_list = superpose_pdbs(model, template, output_superposition=True)
            rmsd_between_template_and_ranked_list.append(round(float(rmsd), 2))
            nalign_between_template_and_ranked_list.append(int(nalign))
            ranked_structure = PDBParser(QUIET=1).get_structure(f'{model}', model)
            ranked_CAs_list = [res['CA'].coord for res in Selection.unfold_entities(ranked_structure, "R") if
                               res.id[1] in template_res_id_list]
            ranked_pdist = distance.pdist(ranked_CAs_list, "euclidean")
            ranked_pdist_matrix = distance.squareform(ranked_pdist)
            diff_pdist_matrices = np.abs(ranked_pdist_matrix - template_pdist_matrix)
            diff_pdist_matrices_list.append(diff_pdist_matrices)
        output_dict['dist_diff_between_template_and_rankeds'][template.split('_')[0]] = np.array(diff_pdist_matrices_list)

        # RMSD BETWEEN TEMPLATE AND RANKED MODELS:
        rmsd_all_CAs_list = []
        for model in ranked_models_list:
            rmsd_all_CAs_list.append(rmsd_calc(model, f'{template}'))
        output_dict['rmsd_between_template_and_rankeds'][template.split('_')[0]] = np.array(rmsd_all_CAs_list)

        with open(f'{af_job_path}/summary/templates/rmsd_between_templates_and_rankeds.log', 'w') as f:
            rows = []
            for item in output_dict['rmsd_between_template_and_rankeds']:
                rows.append([item] + output_dict['rmsd_between_template_and_rankeds'][item].tolist())
            df = pd.DataFrame(rows, columns=['template', 'ranked_0', 'ranked_1', 'ranked_2', 'ranked_3', 'ranked_4'])
            f.write(df.to_markdown())
        f.close()

    if len(templates_list) == 1: # TODO: think about output if list of templates is biger than 1
        template_id = templates_list[0].split('_')[0]
        min_rmsd = output_dict['rmsd_between_template_and_rankeds'][template_id].min()
        ranked_index = output_dict['rmsd_between_template_and_rankeds'][template_id].tolist().index(min_rmsd)
        os.system(f'cp {af_job_path}/summary/templates/{output_dict["ranked_models"][ranked_index][:-4]}_superposed.pdb {af_job_path}/summary/templates/BEST.pdb')
        for i in range(0,5):
            os.system(f'mv ranked_{i}_superposed.pdb ranked_{i}.pdb')

    # INTERDISTANCES IN RANKED MODELS:
    interdistances_in_ranked_models = []
    for num, model in enumerate(ranked_models_list):
        structure = PDBParser(QUIET=1).get_structure(f'{model}', f'{af_job_path}/{model}')
        pdist_matrix = pdist(structure)
        interdistances_in_ranked_models.append(pdist_matrix)
    output_dict['dist_in_rankeds'] = np.array(interdistances_in_ranked_models)

    # RANKED SUPERPOSITION:
    os.system(f'cp {af_job_path}/ranked_0.pdb {af_job_path}/summary')
    os.chdir(f'{af_job_path}')
    with open(f'{af_job_path}/summary/ranked_superpositions.log', 'w') as f:
        superpose_results_list = []
        for ranked in ranked_models_list[1:]:
            rmsd, nalign, quality_q, aligned_res_list = superpose_pdbs(f'{ranked}', 'ranked_0.pdb')
            superpose_results_list.append(['ranked_0', f'{ranked[:-4]}', rmsd, nalign, quality_q])
            os.system(f'mv {ranked[:-4]}_superposed.pdb {af_job_path}/summary/{ranked}')
        df = pd.DataFrame(superpose_results_list,
                          columns=['target model', 'query model', 'core rmsd', 'Naligned', 'quality Q'])
        f.write(df.to_markdown())
    f.close()

    with open(f'{af_job_path}/summary/output_summary.pkl', 'wb') as file:
        pickle.dump(output_dict, file, protocol=pickle.HIGHEST_PROTOCOL)
    os.chdir(cwd)

af_summary(af_job_path='/localdata/albert/ATZR_PAPER/OTHER_JOBS/1ixc-3fxq')