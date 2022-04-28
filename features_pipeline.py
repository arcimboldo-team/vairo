import pandas as pd
import numpy as np
import sys
import analyse_pipeline
sys.path.append('./alphafold')
from alphafold.data import parsers, templates, mmcif_parsing, pipeline
import subprocess
import os
import re
import pickle
from Bio.PDB import PDBParser, PDBIO, Selection
import logging
import string

three_to_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P',
                'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R',
                'TRP': 'W', 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

ID_TO_HHBLITS_AA_3LETTER_CODE = {0: 'ALA', 1: 'CYS', 2: 'ASP', 3: 'GLU', 4: 'PHE', 5: 'GLY', 6: 'HIS',
                                 7: 'ILE', 8: 'LYS', 9: 'LEU', 10: 'MET', 11: 'ASN', 12: 'PRO', 13: 'GLN',
                                 14: 'ARG', 15: 'SER', 16: 'THR', 17: 'VAL', 18: 'TRP', 19: 'TYR', 20: 'X',
                                 21: '-'}

ID_TO_HHBLITS_AA = {0: 'A', 1: 'C', 2: 'D', 3: 'E', 4: 'F', 5: 'G', 6: 'H',
                    7: 'I', 8: 'K', 9: 'L', 10: 'M', 11: 'N', 12: 'P', 13: 'Q',
                    14: 'R', 15: 'S', 16: 'T', 17: 'V', 18: 'W', 19: 'Y',
                    20: 'X', 21: '-'}

atom_types = ['N', 'CA', 'C', 'CB', 'O', 'CG', 'CG1', 'CG2', 'OG', 'OG1', 'SG', 'CD', 'CD1', 'CD2',
              'ND1', 'ND2', 'OD1', 'OD2', 'SD', 'CE', 'CE1', 'CE2', 'CE3', 'NE', 'NE1', 'NE2', 'OE1',
              'OE2', 'CH2', 'NH1', 'NH2', 'OH', 'CZ', 'CZ2', 'CZ3', 'NZ', 'OXT']

atom_order = {atom_type: i for i, atom_type in enumerate(atom_types)}
order_atom = {v: k for k, v in atom_order.items()}

AA_TO_ID_TO_HHBLITS = {v: k for k, v in ID_TO_HHBLITS_AA.items()}


def extract_query_sequence(query_fasta_path):

    with open(query_fasta_path, 'r') as f:
        fasta_lines = f.readlines()
        fasta_name, query_sequence = fasta_lines[0][1:-1], fasta_lines[1].split('\n')[0]

    logging.info(f'Query sequence read:\n{query_sequence}')

    return fasta_name, query_sequence


def run_hhsearch(query_sequence, pdb70_db, output_dir):
    with open('./query.fasta', 'w') as f:
        f.write('> query\n')
        f.write(f'{query_sequence}\n')

    out = subprocess.Popen(['hhsearch', '-i', 'query.fasta', '-o', f'{output_dir}/output.hhr', '-maxseq',
                            '1000000', '-d', pdb70_db],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = out.communicate()
    out_text = stdout.decode('utf-8')
    os.system('rm query.fasta')

    return out_text


def extract_template_features_from_pdb(query_sequence, hhr_path, pdb_id, chain_id, mmcif_db):

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

    template_features, realign_warning = templates._extract_template_features(mmcif_object=parsing_result.mmcif_object,
                                                                              pdb_id=file_id,
                                                                              mapping=mapping,
                                                                              template_sequence=template_sequence,
                                                                              query_sequence=query_sequence,
                                                                              template_chain_id=chain_id,
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

def write_pkl_from_features(features, out_path):

    with open(out_path, 'wb') as handle:
        pickle.dump(features, handle, protocol=pickle.HIGHEST_PROTOCOL)


def template_features_for_pdb_aligned_with_query_sequence(query_sequence, pdb_path, chain_ID, template_name):

    # WARNING: input PDB must be aligned to the MSA part in features #

    seq_length = len(query_sequence)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('test', pdb_path)

    template_sequence = '-' * (seq_length)
    template_res_list = [res for res in Selection.unfold_entities(structure, "R")
                         if res.get_parent().id == chain_ID and res.id[0] != 'W']
    for res in template_res_list:
        template_sequence = template_sequence[:res.id[1]] + three_to_one[res.resname] + template_sequence[res.id[1]:]
    template_sequence = np.array([template_sequence[:seq_length + 1]])[0]

    atom_masks = []
    for i, res in enumerate(template_sequence):
        if res == '-':
            a37_in_res = [0] * 37
            atom_masks.append(a37_in_res)
        else:
            a37_in_res = [0] * 37
            list_of_atoms_in_res = [atom.id for atom in [resi for resi in template_res_list if resi.id[1] == i][0]]
            for atom in list_of_atoms_in_res:
                index = atom_types.index(atom)
                a37_in_res[index] = 1.
            atom_masks.append(a37_in_res)
    template_all_atom_masks = np.array([atom_masks[1:]])

    template_container = []
    for i, res in enumerate(template_all_atom_masks[0]):
        res_container = []
        for j, atom in enumerate(res):
            if atom == 1.:
                resi = [res for res in Selection.unfold_entities(structure, "R")
                        if res.get_parent().id == chain_ID and res.id[0] != 'W' and res.id[1] == (i + 1)][0]
                res_container.append(resi[atom_types[j]].coord)
            else:
                res_container.append(np.array([0.] * 3))
        template_container.append(res_container)
    template_all_atom_positions = np.array([template_container])

    template_domain_names = np.array([(f'{template_name}_A').encode('ascii')])

    template_aatype_container = []
    for res in template_sequence[1:]:
        aa_container = [0] * 22
        aa_container[AA_TO_ID_TO_HHBLITS[res]] = 1
        template_aatype_container.append(aa_container)
    template_aatype = np.array([template_aatype_container])

    template_sum_probs = np.array([100.])

    template_sequence_to_add = np.array([template_sequence[1:].encode('ascii')])
    template_all_atom_masks_to_add = template_all_atom_masks
    template_all_atom_positions_to_add = template_all_atom_positions
    template_domain_names_to_add = template_domain_names
    template_aatype_to_add = template_aatype
    template_sum_probs_to_add = template_sum_probs

    template_features = {}

    template_features['template_sequence'] = template_sequence_to_add
    template_features['template_all_atom_masks'] = template_all_atom_masks_to_add
    template_features['template_all_atom_positions'] = template_all_atom_positions_to_add
    template_features['template_domain_names'] = template_domain_names_to_add
    template_features['template_aatype'] = template_aatype_to_add
    template_features['template_sum_probs'] = np.array([template_sum_probs_to_add])

    logging.info(f'Template features are ready from {pdb_path}.')

    return template_features


def convert_template_to_polyala(pdb_in_path, pdb_out_path, list_of_res_ranges):

    ala_res_list = []

    logging.info(f'Template will be converted to polyalanine model in the following ranges: {list_of_res_ranges}.')

    if len(list_of_res_ranges) > 1:
        for item in list_of_res_ranges:
            ala_res_list.extend(list(range(int(item[0]),int(item[1]))))
    else:
        ala_res_list = list(range(int(list_of_res_ranges[0][0]), int(list_of_res_ranges[0][1])))

    ala_atoms_list = ['N', 'CA', 'C', 'CB', 'O']

    with open(f'{pdb_out_path}', 'w') as f_out:
        with open(f'{pdb_in_path}') as f_in:
            lines = f_in.readlines()
            num = 0
            for line in lines:
                if line[:4] == 'ATOM':
                    if int(line[22:26].replace(' ', '')) in ala_res_list:
                        if line[13:16].replace(' ', '') in ala_atoms_list:
                            num = num + 1
                            f_out.write(line[:7] + str(num).rjust(4) + line[11:17] + 'ALA' + line[20:])
                    else:
                        num = num + 1
                        f_out.write(line[:7] + str(num).rjust(4) + line[11:])


class Features:

    def __init__(self):

        self.features = None

    def features_from_features_pkl(self, features_pkl_path):

        self.features = pd.read_pickle(features_pkl_path)

        return self.features

    def features_from_pdb_id(self, query_sequence, pdb_id, chain_id, out_dir, mmcif_db, pdb70_db, hhr_path=None):

        if hhr_path == None:
            logging.info('Running HHSEARCH...')
            run_hhsearch(query_sequence, pdb70_db, out_dir)
            hhr_path = f'{out_dir}/output.hhr'
            logging.info(f'...HHSEARCH finished. output.hhr saved at {out_dir}/output.hhr')
        else:
            logging.info(f'Reading HHR provided by the user ({hhr_path}).')

        template_features = extract_template_features_from_pdb(query_sequence, hhr_path, pdb_id, chain_id, mmcif_db)
        logging.info(f'Building template features for PDB {pdb_id}_{chain_id}.')
        template_sequence = template_features['template_sequence'][0].decode('utf-8')
        logging.info(f'Template sequence ({pdb_id}_{chain_id}):\n{template_sequence}')

        logging.info(f'Building MSA features.')
        template_sequence_copy, query_sequence_copy = list(template_sequence), list(query_sequence)
        counter = 0
        for num, item in enumerate(template_features['template_sequence'][0].decode('utf-8')):
            if item == '-':
                pass
            else:
                if template_sequence_copy[num] == query_sequence_copy[num]:
                    pass
                else:
                    counter = counter + 1
        if counter == 0:
            logging.info('NOTE: Template and query sequences are the same: MSA will only contain the query sequence.')
            input_msa = f'>query\n{query_sequence}'
        else:
            input_msa = f'>query\n{query_sequence}\n>template\n{template_sequence}'

        custom_msa_result = {'a3m': input_msa}
        custom_msa = parsers.parse_a3m(custom_msa_result['a3m'])
        msa_features = pipeline.make_msa_features_for_custom(custom_msa)

        logging.info(f'Building sequence features.')
        input_seqs, input_descs = parsers.parse_fasta(f'>query\n{query_sequence}')
        input_sequence = input_seqs[0]
        input_description = input_descs[0]
        sequence_features = pipeline.make_sequence_features(sequence=input_sequence, description=input_description,
                                                            num_res=len(query_sequence))
        logging.info('Input features are ready for the inference!')
        self.features = {**sequence_features, **msa_features, **template_features}

        return self.features

    def features_for_query_sequence_and_experimental_assembly_in_pdb(self, query_subunit_sequence, pdb_id, mmcif_db,
                                                                     pdb70_db, out_dir, hhr_path=False):

        query_seq_length = len(query_subunit_sequence)

        analyse_pipeline.download_pdb(pdb_id, out_dir)
        os.system(f'cp {out_dir}/pdb{pdb_id}.ent {out_dir}/{pdb_id}.pdb')

        structure = PDBParser(QUIET=1).get_structure(f'{pdb_id}', f'{out_dir}/{pdb_id}.pdb')
        chain_list = [chain.id for chain in Selection.unfold_entities(structure, 'C')]

        if not hhr_path:
            logging.info('Running HHSEARCH...')
            run_hhsearch(query_subunit_sequence, pdb70_db, out_dir)
            hhr_path = f'{out_dir}/output.hhr'
            logging.info(f'...HHSEARCH finished. output.hhr saved at {out_dir}/output.hhr')
        else:
            logging.info(f'Reading HHR provided by the user ({hhr_path}).')

        transformations_list = analyse_pipeline.read_remark_350(f'{out_dir}/pdb{pdb_id.lower()}.ent')
        logging.info(f'REMARK 350 for PDB {pdb_id} contains {len(transformations_list)} --> Assembly will contain'
                     f' {len(transformations_list) * len(chain_list)} subunits.')

        new_chain_list = list(string.ascii_uppercase)[:len(transformations_list) * len(chain_list)]

        counter = 0 # TODO Deal with those cases where there is more than one assembly in ASU!
        for chain in chain_list:
            self.features = extract_template_features_from_pdb(query_subunit_sequence, hhr_path, pdb_id, chain,
                                                               mmcif_db)
            self.write_all_templates_in_features(output_path=out_dir)
            for num, transformation in enumerate(transformations_list):
                counter += 1
                logging.info(f'Rot and trans ({transformation}) is applied on {pdb_id}{chain}_template.pdb')
                analyse_pipeline.rot_and_trans(f'{out_dir}/{pdb_id}{chain}_template.pdb',
                                               f'{out_dir}/{counter}.pdb',
                                               transformation)

                if counter == 1:
                    analyse_pipeline.change_chain_and_apply_offset_in_single_chain(f'{out_dir}/{counter}.pdb',
                                                                                   f'{out_dir}/{new_chain_list[counter-1]}.pdb',
                                                                                   offset=None, chain='A')
                else:
                    offset = query_seq_length*(counter-1)+50*(counter-1) # 50 glycines as offset !!!
                    analyse_pipeline.change_chain_and_apply_offset_in_single_chain(f'{out_dir}/{counter}.pdb',
                                                                                   f'{out_dir}/{new_chain_list[counter - 1]}.pdb',
                                                                                   offset=offset, chain='A')
                os.system(f'rm {out_dir}/{counter}.pdb')

        list_of_paths_of_pdbs_to_merge = [f'{out_dir}/{ch}.pdb' for ch in new_chain_list]
        analyse_pipeline.merge_pdbs(list_of_paths_of_pdbs_to_merge, merged_pdb_path=f'{out_dir}/template.pdb')
        os.system(f'rm {out_dir}/[A-Z].pdb')
        logging.info(f'Assembly template is ready for the inference ({out_dir}/template.pdb)')

        assembly_sequence_with_linkers = ''
        for num in range(len(new_chain_list)):
            assembly_sequence_with_linkers = assembly_sequence_with_linkers + query_subunit_sequence + 50*'G'
        assembly_sequence_with_linkers = assembly_sequence_with_linkers[:-50]
        logging.info('The assembly sequence with 50 Gly residues between subunits is ready for the inference:')
        logging.info(assembly_sequence_with_linkers)

        template_features = template_features_for_pdb_aligned_with_query_sequence(
            query_sequence=assembly_sequence_with_linkers,
            pdb_path=f'{out_dir}/template.pdb',
            chain_ID='A',
            template_name=pdb_id)

        template_sequence = template_features['template_sequence'][0].decode('utf-8')

        logging.info(f'Building MSA features.')
        template_sequence_copy, query_sequence_copy = list(template_sequence), list(assembly_sequence_with_linkers)
        counter = 0
        for num, item in enumerate(template_features['template_sequence'][0].decode('utf-8')):
            if item == '-':
                pass
            else:
                if template_sequence_copy[num] == query_sequence_copy[num]:
                    pass
                else:
                    counter = counter + 1
        if counter == 0:
            logging.info('NOTE: Template and query sequences are the same: MSA will only contain the query sequence.')
            input_msa = f'>query\n{assembly_sequence_with_linkers}'
        else:
            input_msa = f'>query\n{assembly_sequence_with_linkers}\n>template\n{template_sequence}'

        custom_msa_result = {'a3m': input_msa}
        custom_msa = parsers.parse_a3m(custom_msa_result['a3m'])
        msa_features = pipeline.make_msa_features_for_custom(custom_msa)

        logging.info(f'Building sequence features.')
        input_seqs, input_descs = parsers.parse_fasta(f'>query\n{assembly_sequence_with_linkers}')
        input_sequence = input_seqs[0]
        input_description = input_descs[0]
        sequence_features = pipeline.make_sequence_features(sequence=input_sequence, description=input_description,
                                                            num_res=len(assembly_sequence_with_linkers))
        logging.info('Input features are ready for the inference!')
        self.features = {**sequence_features, **msa_features, **template_features}

        return self.features

    def features_for_custom_pdb(self, query_sequence, pdb_path):

        '''poly_ala can be set to False, all or list of ranges'''

        template_features = template_features_for_pdb_aligned_with_query_sequence(
            query_sequence=query_sequence,
            pdb_path=pdb_path,
            chain_ID='A',
            template_name=pdb_path.split('/')[-1][:4])

        template_sequence = template_features['template_sequence'][0].decode('utf-8')

        logging.info(f'Building MSA features.')
        template_sequence_copy, query_sequence_copy = list(template_sequence), list(query_sequence)
        counter = 0
        for num, item in enumerate(template_features['template_sequence'][0].decode('utf-8')):
            if item == '-':
                pass
            else:
                if template_sequence_copy[num] == query_sequence_copy[num]:
                    pass
                else:
                    counter = counter + 1
        if counter == 0 :
            logging.info('NOTE: Template and query sequences are the same: MSA will only contain the query sequence.')
            input_msa = f'>query\n{query_sequence}'
        else:
            input_msa = f'>query\n{query_sequence}\n>template\n{template_sequence}'

        custom_msa_result = {'a3m': input_msa}
        custom_msa = parsers.parse_a3m(custom_msa_result['a3m'])
        msa_features = pipeline.make_msa_features_for_custom(custom_msa)

        logging.info(f'Building sequence features.')
        input_seqs, input_descs = parsers.parse_fasta(f'>query\n{query_sequence}')
        input_sequence = input_seqs[0]
        input_description = input_descs[0]
        sequence_features = pipeline.make_sequence_features(sequence=input_sequence, description=input_description,
                                                            num_res=len(query_sequence))
        logging.info('Input features are ready for the inference!')
        self.features = {**sequence_features, **msa_features, **template_features}

        return self.features

    def write_all_templates_in_features(self, output_path):

        for i, pdb_name in enumerate(self.features['template_domain_names']):
            pdb, chain = pdb_name.decode('utf-8').split('_')
            output_pdb = open(output_path + '/' + pdb + chain + '_template.pdb', 'w')
            template_domain_index = np.where(self.features['template_domain_names'] == pdb_name)[0][0]
            true_seq = ''
            atom_num_int = 0  # AQUI
            for index, atoms_mask in enumerate(self.features['template_all_atom_masks'][template_domain_index][:]):
                template_residue_masks = self.features['template_aatype'][template_domain_index][index]
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
                        self.features['template_all_atom_positions'][template_domain_index][index][
                            atom_types.index(atom)][
                            0])))).rjust(8)
                    y_coord = str('%8.3f' % (float(str(
                        self.features['template_all_atom_positions'][template_domain_index][index][
                            atom_types.index(atom)][
                            1])))).rjust(8)
                    z_coord = str('%8.3f' % (float(str(
                        self.features['template_all_atom_positions'][template_domain_index][index][
                            atom_types.index(atom)][
                            2])))).rjust(8)
                    occ = str('%6.2f' % (float('1.0'))).rjust(6)
                    bfact = str('%6.2f' % (float('25.0'))).ljust(6)
                    atom_type = atom[0].rjust(12)
                    output_pdb.write(
                        f'{atom_remark}{atom_num}  {atom_name}{res_name} {chain_id}{res_num}    {x_coord}{y_coord}{z_coord}{occ}{bfact}{atom_type}\n')

            output_pdb.close()

    def append_row_in_msa(self, sequence):

        sequence_array = np.array([AA_TO_ID_TO_HHBLITS[res] for res in sequence])
        self.features['msa'] = np.vstack([self.features['msa'], sequence_array])
        self.features['msa_uniprot_accession_identifiers'] = np.hstack(
            [self.features['msa_uniprot_accession_identifiers'], ''])
        self.features['deletion_matrix_int'] = np.vstack(
            [self.features['deletion_matrix_int'], np.zeros(self.features['msa'].shape[1])])
        self.features['msa_species_identifiers'] = np.hstack([self.features['msa_species_identifiers'], ''])
        self.features['num_alignments'] = np.full(self.features['num_alignments'].shape, len(self.features['msa']))

        return self.features

    def delete_rows_in_msa(self, list_to_remove):

        self.features['msa'] = np.delete(self.features['msa'], list_to_remove, axis=0)
        self.features['deletion_matrix_int'] = np.delete(self.features['deletion_matrix_int'], list_to_remove, axis=0)
        self.features['msa_uniprot_accession_identifiers'] = np.delete(self.features['msa_uniprot_accession_identifiers'],
                                                                  list_to_remove, axis=0)
        self.features['msa_species_identifiers'] = np.delete(self.features['msa_species_identifiers'], list_to_remove, axis=0)
        self.features['num_alignments'] = np.full(self.features['num_alignments'].shape, len(self.features['msa']))

        return self.features

    def replace_columns_by_gaps_in_msa(self, list_to_remove, specific_row=False):

        for col in list_to_remove:
            if specific_row == False:
                self.features['msa'][:, col] = 21
                self.features['deletion_matrix_int'][:, col] = 0
            else:
                self.features['msa'][specific_row, col] = 21
                self.features['deletion_matrix_int'][specific_row, col] = 0

        return self.features


def append_template_features_from_pdb(features, features_for_append): # TESTING!!!

    features['template_sequence'] = np.hstack([features['template_sequence'], features_for_append['template_sequence']])
    features['template_all_atom_masks'] = np.vstack([features['template_all_atom_masks'], features_for_append['template_all_atom_masks']])
    features['template_all_atom_positions'] = np.vstack([features['template_all_atom_positions'], features_for_append['template_all_atom_positions']])
    features['template_domain_names'] = np.hstack([features['template_domain_names'], features_for_append['template_domain_names']])
    features['template_aatype'] = np.vstack([features['template_aatype'], features_for_append['template_aatype']])
    features['template_sum_probs'] = np.vstack([features['template_sum_probs'], features_for_append['template_sum_probs']])

    return features


# append_template_features_from_pdb()
