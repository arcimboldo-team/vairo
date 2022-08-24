from typing import Dict, List, Set
from Bio.PDB import PDBParser, Selection
import os
import re
from ALPHAFOLD.alphafold.data import parsers, pipeline, templates, mmcif_parsing, pipeline, msa_identifiers
from ALPHAFOLD.alphafold.common import residue_constants
import numpy as np
import pickle
import logging


three_to_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P',
                'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R',
                'TRP': 'W', 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
AA_TO_ID_TO_HHBLITS = {v: k for k, v in residue_constants.ID_TO_HHBLITS_AA.items()}
ID_TO_HHBLITS_AA_3LETTER_CODE = {0: 'ALA', 1: 'CYS', 2: 'ASP', 3: 'GLU', 4: 'PHE', 5: 'GLY', 6: 'HIS',
                                 7: 'ILE', 8: 'LYS', 9: 'LEU', 10: 'MET', 11: 'ASN', 12: 'PRO', 13: 'GLN',
                                 14: 'ARG', 15: 'SER', 16: 'THR', 17: 'VAL', 18: 'TRP', 19: 'TYR', 20: 'X',
                                 21: '-'}
atom_types = ['N', 'CA', 'C', 'CB', 'O', 'CG', 'CG1', 'CG2', 'OG', 'OG1', 'SG', 'CD', 'CD1', 'CD2',
              'ND1', 'ND2', 'OD1', 'OD2', 'SD', 'CE', 'CE1', 'CE2', 'CE3', 'NE', 'NE1', 'NE2', 'OE1',
              'OE2', 'CH2', 'NH1', 'NH2', 'OH', 'CZ', 'CZ2', 'CZ3', 'NZ', 'OXT']
atom_order = {atom_type: i for i, atom_type in enumerate(atom_types)}
order_atom = {v: k for k, v in atom_order.items()}


def empty_msa_features(query_sequence):

    msa = {'a3m': f'>query\n{query_sequence}'}
    custom_msa = parsers.parse_a3m(msa['a3m'])

    msas = [custom_msa]  # ACT: it is needed in order to introduce MSA inside a list in the code

    int_msa = []
    deletion_matrix = []
    uniprot_accession_ids = []
    species_ids = []
    seen_sequences = set()
    for msa_index, msa in enumerate(msas):
        if not msa:
            raise ValueError(f'MSA {msa_index} must contain at least one sequence.')
        for sequence_index, sequence in enumerate(msa.sequences):
            if sequence in seen_sequences:
                continue
            seen_sequences.add(sequence)
            int_msa.append(
                [residue_constants.HHBLITS_AA_TO_ID[res] for res in sequence])
            deletion_matrix.append(msa.deletion_matrix[sequence_index])
            identifiers = msa_identifiers.get_identifiers(
                msa.descriptions[sequence_index])
            uniprot_accession_ids.append(
                identifiers.uniprot_accession_id.encode('utf-8'))
            species_ids.append(identifiers.species_id.encode('utf-8'))

    num_res = len(msas[0].sequences[0])
    num_alignments = len(int_msa)
    features = {}
    features['deletion_matrix_int'] = np.array(deletion_matrix, dtype=np.int32)
    features['msa'] = np.array(int_msa, dtype=np.int32)
    features['num_alignments'] = np.array(
        [num_alignments] * num_res, dtype=np.int32)
    features['msa_uniprot_accession_identifiers'] = np.array(
        uniprot_accession_ids, dtype=np.object_)
    features['msa_species_identifiers'] = np.array(species_ids, dtype=np.object_)
    return features


def empty_template_features(query_sequence):

    ln = (len(query_sequence) if isinstance(query_sequence, str) else sum(len(s) for s in query_sequence))
    output_templates_sequence = "A" * ln

    templates_all_atom_positions = np.zeros((ln, residue_constants.atom_type_num, 3))
    templates_all_atom_masks = np.zeros((ln, residue_constants.atom_type_num))
    templates_aatype = residue_constants.sequence_to_onehot(output_templates_sequence, residue_constants.HHBLITS_AA_TO_ID)
    template_sum_probs = f'None'
    template_features = {
        "template_all_atom_positions": np.tile(templates_all_atom_positions[None], [0, 1, 1, 1]),
        "template_all_atom_masks": np.tile(templates_all_atom_masks[None], [0, 1, 1]),
        "template_sequence": [f"None".encode()] * 0,
        "template_aatype": np.tile(np.array(templates_aatype)[None], [0, 1, 1]),
        "template_domain_names": [f"None".encode()] * 0,
        "template_sum_probs": np.tile(template_sum_probs, [0, 1])
    }
    return template_features

def extract_template_features_from_pdb(query_sequence, hhr_path, pdb_id, chain_id, mmcif_db): # TODO: template names must be in lowercase

    hhr_text = open(hhr_path, 'r').read()

    matches = re.finditer(r'No\s+\d+', hhr_text)
    matches_positions = [match.start() for match in matches] + [len(hhr_text)]

    detailed_lines_list = []
    for i in range(len(matches_positions) - 1):
        detailed_lines_list.append(hhr_text[matches_positions[i]:matches_positions[i + 1]].split('\n')[:-3])

    hits_list = [detailed_lines for detailed_lines in detailed_lines_list if pdb_id in detailed_lines[1]]

    detailed_lines = hits_list[0]

    file_id = f'{pdb_id.lower()}'
    hit = parsers._parse_hhr_hit(detailed_lines)

    template_sequence = hit.hit_sequence.replace('-', '')
    mapping = templates._build_query_to_hit_index_mapping(
            hit.query, hit.hit_sequence, hit.indices_hit, hit.indices_query,
            query_sequence)

    mmcif_string = open(f'{mmcif_db}/{pdb_id}.cif').read()
    parsing_result = mmcif_parsing.parse(file_id=file_id, mmcif_string=mmcif_string)
    
    template_features, _ = templates._extract_template_features(
            mmcif_object=parsing_result.mmcif_object,
            pdb_id=file_id,
            mapping=mapping,
            template_sequence=template_sequence,
            query_sequence=query_sequence,
            template_chain_id=chain_id,
            kalign_binary_path='kalign')
    template_features['template_sum_probs'] = np.array([[hit.sum_probs]])
    template_features['template_aatype'] = np.array([template_features['template_aatype']])
    template_features['template_all_atom_masks'] = np.array([template_features['template_all_atom_masks']])
    template_features['template_all_atom_positions'] = np.array([template_features['template_all_atom_positions']])
    template_features['template_domain_names'] = np.array([template_features['template_domain_names']])
    template_features['template_sequence'] = np.array([template_features['template_sequence']])

    return template_features

def extract_template_features_from_aligned_pdb_and_sequence(query_sequence: str, pdb_path: str, pdb_id: str, chain_id: str):

    # WARNING: input PDB must be aligned to the MSA part in features #

    seq_length = len(query_sequence)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_path)

    template_sequence = '-' * (seq_length)
    template_res_list = [res for res in Selection.unfold_entities(structure, "R")
                            if res.get_parent().id == chain_id and res.id[0] != 'W']
    for res in template_res_list:
        template_sequence = template_sequence[:res.id[1]] + three_to_one[res.resname] + template_sequence[
                                                                                        res.id[1]:]
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
                        if res.get_parent().id == chain_id and res.id[0] != 'W' and res.id[1] == (i + 1)][0]
                res_container.append(resi[atom_types[j]].coord)
            else:
                res_container.append(np.array([0.] * 3))
        template_container.append(res_container)
    template_all_atom_positions = np.array([template_container])

    template_domain_names = np.array([(f'{pdb_id}_{chain_id}').encode('ascii')])

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

    return template_features

def print_features(pkl_in_path: str):
    with open(f"{pkl_in_path}", "rb") as input_file:
        features = pickle.load(input_file)

    for key in features.keys():
        logging.info(f'{key} {features[key].shape}')

    logging.info('\n')
    logging.info('MSA:')
    for num, name in enumerate(features['msa_uniprot_accession_identifiers']):
        logging.info(name.decode('utf-8'))
        logging.info('\n')
        logging.info(''.join([residue_constants.ID_TO_HHBLITS_AA[res] for res in features['msa'][num].tolist()]))
        logging.info('\n')

    logging.info('TEMPLATES:')
    for num, seq in enumerate(features['template_sequence']):
        logging.info(f'{features["template_domain_names"][num].decode("utf-8")}:\n')
        for i in range(4):
            logging.info('\t'+''.join(np.array_split(list(seq.decode('utf-8')),4)[i].tolist()))
        logging.info('\n')

class Features:

    def __init__(self, query_sequence: str):
        self.query_sequence: str
        self.sequence_features: pipeline.FeatureDict
        self.msa_features: Dict
        self.template_features: Dict

        self.query_sequence = query_sequence
        self.sequence_features = pipeline.make_sequence_features(sequence=self.query_sequence,
                                                                 description='Query',
                                                                 num_res=len(self.query_sequence))
        self.msa_features = empty_msa_features(query_sequence=self.query_sequence)
        self.template_features = empty_template_features(query_sequence=self.query_sequence)        

    def append_new_template_features(self, new_template_features: Dict, custom_sum_prob: int = None) -> Dict:

        self.template_features['template_all_atom_positions'] = np.vstack([self.template_features['template_all_atom_positions'], new_template_features['template_all_atom_positions']])
        self.template_features['template_all_atom_masks'] = np.vstack([self.template_features['template_all_atom_masks'], new_template_features['template_all_atom_masks']])
        self.template_features['template_aatype'] = np.vstack([self.template_features['template_aatype'], new_template_features['template_aatype']])
        self.template_features['template_sequence'] = np.hstack([self.template_features['template_sequence'], new_template_features['template_sequence']])
        self.template_features['template_domain_names'] = np.hstack([self.template_features['template_domain_names'], new_template_features['template_domain_names']])
        if not custom_sum_prob:
            self.template_features['template_sum_probs'] = np.vstack([self.template_features['template_sum_probs'], new_template_features['template_sum_probs']])
        else:
            self.template_features['template_sum_probs'] = np.vstack([self.template_features['template_sum_probs'], custom_sum_prob])

        return self.template_features

    def append_row_in_msa(self, sequence: str, msa_uniprot_accession_identifiers):

        sequence_array = np.array([AA_TO_ID_TO_HHBLITS[res] for res in sequence])
        self.msa_features['msa'] = np.vstack([self.msa_features['msa'], sequence_array])
        self.msa_features['msa_uniprot_accession_identifiers'] = np.hstack([self.msa_features['msa_uniprot_accession_identifiers'], msa_uniprot_accession_identifiers.encode()])
        self.msa_features['deletion_matrix_int'] = np.vstack([self.msa_features['deletion_matrix_int'], np.zeros(self.msa_features['msa'].shape[1])])
        self.msa_features['msa_species_identifiers'] = np.hstack([self.msa_features['msa_species_identifiers'], ''])
        self.msa_features['num_alignments'] = np.full(self.msa_features['num_alignments'].shape, len(self.msa_features['msa']))

    def write_all_templates_in_features(self, output_dir: str) -> Dict:

        templates_dict = {}

        for pdb_name in self.template_features['template_domain_names']:
            pdb = pdb_name.decode('utf-8')
            pdb_path = os.path.join(output_dir,f'{pdb}_template.pdb')
            templates_dict[pdb] = pdb_path
            with open(pdb_path, 'w') as output_pdb:
                template_domain_index = np.where(self.template_features['template_domain_names'] == pdb_name)[0][0]
                atom_num_int = 0
                for index, atoms_mask in enumerate(self.template_features['template_all_atom_masks'][template_domain_index][:]):
                    template_residue_masks = self.template_features['template_aatype'][template_domain_index][index]
                    template_residue_masks_index = np.where(template_residue_masks == 1)[0][0]
                    res_type = ID_TO_HHBLITS_AA_3LETTER_CODE[template_residue_masks_index]
                    list_of_atoms_in_residue = [order_atom[i] for i, atom in enumerate(atoms_mask) if atom == 1]
                    for atom in list_of_atoms_in_residue:
                        atom_num_int = atom_num_int + 1
                        atom_remark = 'ATOM'.ljust(6)
                        atom_num = str(atom_num_int).rjust(5)
                        atom_name = atom.ljust(4)
                        res_name = res_type.ljust(3)
                        res_num = str(index + 1).rjust(4)
                        x_coord = str('%8.3f' % (float(str(
                            self.template_features['template_all_atom_positions'][template_domain_index][index][
                                atom_types.index(atom)][
                                0])))).rjust(8)
                        y_coord = str('%8.3f' % (float(str(
                            self.template_features['template_all_atom_positions'][template_domain_index][index][
                                atom_types.index(atom)][
                                1])))).rjust(8)
                        z_coord = str('%8.3f' % (float(str(
                            self.template_features['template_all_atom_positions'][template_domain_index][index][
                                atom_types.index(atom)][
                                2])))).rjust(8)
                        occ = str('%6.2f' % (float('1.0'))).rjust(6)
                        bfact = str('%6.2f' % (float('25.0'))).ljust(6)
                        atom_type = atom[0].rjust(12)
                        output_pdb.write(
                            f'{atom_remark}{atom_num}  {atom_name}{res_name} A{res_num}    {x_coord}{y_coord}{z_coord}{occ}{bfact}{atom_type}\n')

        return templates_dict

    def complete_msa_from_template_features(self, template_features):

        msa_from_templates_list = [(''.join(residue_constants.ID_TO_HHBLITS_AA[res] for res in self.msa_features['msa'][0]), 'Query')]
        for num, seq in enumerate(template_features['template_sequence']):
            msa = ''.join([f'>{seq[1]}\n{seq[0]}\n' for seq in msa_from_templates_list])
            if seq.decode('utf-8') not in msa:
                msa_from_templates_list.append((seq.decode('utf-8'), template_features['template_domain_names'][num].decode('utf-8')))

        for seq in msa_from_templates_list[1:]:
            self.append_row_in_msa(sequence=seq[0], msa_uniprot_accession_identifiers=seq[1])

    def write_pkl(self, output_dir: str):

        logging.info(f'Writting all input features in {output_dir}')

        merged_features = self.merge_features()
        with open(output_dir, 'wb') as handle:
            pickle.dump(merged_features, handle, protocol=pickle.HIGHEST_PROTOCOL)

    def merge_features(self) -> Set:

        logging.info(f'Merging sequence, msa and template features!')

        return {**self.sequence_features, **self.msa_features, **self.template_features}