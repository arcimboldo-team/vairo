import os
import re
import numpy as np
import pickle
import logging
from typing import Dict, List, Union, Any
from Bio.PDB import PDBParser, Selection
from alphafold.data import parsers, pipeline, templates, mmcif_parsing, pipeline, msa_identifiers
from alphafold.common import residue_constants
from libs import bioutils, utils

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

        self.template_features['template_all_atom_positions'] = np.vstack(
            [self.template_features['template_all_atom_positions'],
             new_template_features['template_all_atom_positions']])
        self.template_features['template_all_atom_masks'] = np.vstack(
            [self.template_features['template_all_atom_masks'], new_template_features['template_all_atom_masks']])
        self.template_features['template_aatype'] = np.vstack(
            [self.template_features['template_aatype'], new_template_features['template_aatype']])
        self.template_features['template_sequence'] = np.hstack(
            [self.template_features['template_sequence'], new_template_features['template_sequence']])
        self.template_features['template_domain_names'] = np.hstack(
            [self.template_features['template_domain_names'], new_template_features['template_domain_names']])
        if not custom_sum_prob:
            self.template_features['template_sum_probs'] = np.vstack(
                [self.template_features['template_sum_probs'], new_template_features['template_sum_probs']])
        else:
            self.template_features['template_sum_probs'] = np.vstack(
                [self.template_features['template_sum_probs'], custom_sum_prob])

        return self.template_features

    def append_row_in_msa(self, sequence: str, sequence_id: str):

        sequence_array = np.array([AA_TO_ID_TO_HHBLITS[res] for res in sequence])
        self.msa_features['msa'] = np.vstack([self.msa_features['msa'], sequence_array])
        self.msa_features['accession_ids'] = np.hstack(
            [self.msa_features['accession_ids'], sequence_id.encode()])
        self.msa_features['deletion_matrix_int'] = np.vstack(
            [self.msa_features['deletion_matrix_int'], np.zeros(self.msa_features['msa'].shape[1])])
        self.msa_features['msa_species_identifiers'] = np.hstack([self.msa_features['msa_species_identifiers'], ''])
        self.msa_features['num_alignments'] = np.full(self.msa_features['num_alignments'].shape,
                                                      len(self.msa_features['msa']))

    def complete_msa_from_template_features(self, template_features):

        msa_from_templates_list = [
            (''.join(residue_constants.ID_TO_HHBLITS_AA[res] for res in self.msa_features['msa'][0]), 'Query')]
        for num, seq in enumerate(template_features['template_sequence']):
            msa = ''.join([f'>{seq[1]}\n{seq[0]}\n' for seq in msa_from_templates_list])
            if seq.decode('utf-8') not in msa:
                msa_from_templates_list.append(
                    (seq.decode('utf-8'), template_features['template_domain_names'][num].decode('utf-8')))

        for seq in msa_from_templates_list[1:]:
            self.append_row_in_msa(sequence=seq[0], sequence_id=seq[1])

    def write_all_templates_in_features(self, output_dir: str, chain='A', print_number=True) -> Dict:

        return write_templates_in_features(self.template_features, output_dir, chain, print_number)

    def write_pkl(self, pkl_path: str):

        logging.info(f'Writing all input features in {pkl_path}')

        merged_features = self.merge_features()
        with open(pkl_path, 'wb') as f_out:
            pickle.dump(merged_features, f_out, protocol=pickle.HIGHEST_PROTOCOL)

    def get_names_templates(self) -> List[str]:

        return [x.decode() for x in self.template_features['template_domain_names']]

    def get_names_msa(self) -> List[str]:

        return [x.decode() for x in self.msa_features['accession_ids']]

    def get_msa_by_name(self, name: str) -> Union[str, None]:
        index = np.flatnonzero(
            np.core.defchararray.find(name.encode(), self.msa_features['accession_ids']) != -1)
        
        if len(index) > 1:
            return (''.join(
                [residue_constants.ID_TO_HHBLITS_AA[res] for res in self.msa_features['msa'][index[-1]].tolist()]))
        
        return None

    def get_sequence_by_name(self, name: str) -> str:

        index = np.where(self.template_features['template_domain_names'] == name.encode())
        return self.template_features['template_sequence'][index][0].decode()

    def merge_features(self) -> Dict[Any, Any]:

        logging.info(f'Merging sequence, msa and template features!')

        return {**self.sequence_features, **self.msa_features, **self.template_features}

    def slicing_features(self, mosaic: int, overlap: int) -> List:
        # This function will generate as many features
        # as required per size. It will return a list with
        # the path of all the generated features

        sequence = (''.join([residue_constants.ID_TO_HHBLITS_AA[res] for res in self.msa_features['msa'][0].tolist()]))
        chunk_list = utils.chunk_string(length=len(sequence), number_partitions=mosaic, overlap=overlap)
        features_list = []
        for start_min, start_max in chunk_list:
            new_features = Features(query_sequence=sequence[start_min:start_max])
            for i in range(1, len(self.msa_features['accession_ids'])):
                sequence = (
                    ''.join([residue_constants.ID_TO_HHBLITS_AA[res] for res in self.msa_features['msa'][i].tolist()]))
                new_features.append_row_in_msa(sequence=sequence[start_min:start_max],
                                               sequence_id=self.msa_features['accession_ids'][i].decode("utf-8"))
            for i in range(0, len(self.template_features['template_sequence'])):
                template_dict = {
                    'template_all_atom_positions': np.array(
                        [self.template_features['template_all_atom_positions'][i][start_min:start_max]]),
                    'template_all_atom_masks': np.array(
                        [self.template_features['template_all_atom_masks'][i][start_min:start_max]]),
                    'template_aatype': np.array([self.template_features['template_aatype'][i][start_min:start_max]]),
                    'template_sequence': np.array(
                        [self.template_features['template_sequence'][i][start_min:start_max]]),
                    'template_domain_names': np.array([self.template_features['template_domain_names'][i]]),
                    'template_sum_probs': np.array([self.template_features['template_sum_probs'][i]])
                }
                new_features.append_new_template_features(template_dict)
            features_list.append(new_features)
        logging.info(
            f'Features has been sliced in {len(features_list)} partitions with the following sizes: {chunk_list}')
        return features_list


def empty_msa_features(query_sequence):
    msa = {'a3m': f'>query\n{query_sequence}'}
    custom_msa = parsers.parse_a3m(msa['a3m'])

    msas = [custom_msa]  # ACT: it is needed in order to introduce MSA inside a list in the code

    int_msa = []
    deletion_matrix = []
    accession_ids = []
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
            identifiers = msa_identifiers.get_identifiers(msa.descriptions[sequence_index])
            accession_ids.append(str('').encode('utf-8'))
            species_ids.append(identifiers.species_id.encode('utf-8'))

    num_res = len(msas[0].sequences[0])
    num_alignments = len(int_msa)
    features = {'deletion_matrix_int': np.array(deletion_matrix, dtype=np.int32),
                'msa': np.array(int_msa, dtype=np.int32), 'num_alignments': np.array(
            [num_alignments] * num_res, dtype=np.int32), 'accession_ids': np.array(
            accession_ids, dtype=np.object_), 'msa_species_identifiers': np.array(species_ids, dtype=np.object_)}
    return features


def empty_template_features(query_sequence):
    ln = (len(query_sequence) if isinstance(query_sequence, str) else sum(len(s) for s in query_sequence))
    output_templates_sequence = "A" * ln

    templates_all_atom_positions = np.zeros((ln, residue_constants.atom_type_num, 3))
    templates_all_atom_masks = np.zeros((ln, residue_constants.atom_type_num))
    templates_aatype = residue_constants.sequence_to_onehot(output_templates_sequence,
                                                            residue_constants.HHBLITS_AA_TO_ID)
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


def extract_template_features_from_pdb(query_sequence, hhr_path, cif_path, chain_id) -> List[str]:
    pdb_id = utils.get_file_name(cif_path)
    hhr_text = open(hhr_path, 'r').read()

    matches = re.finditer(r'No\s+\d+', hhr_text)
    matches_positions = [match.start() for match in matches] + [len(hhr_text)]

    detailed_lines_list = []
    for i in range(len(matches_positions) - 1):
        detailed_lines_list.append(hhr_text[matches_positions[i]:matches_positions[i + 1]].split('\n')[:-3])

    hits_list = [detailed_lines for detailed_lines in detailed_lines_list if
                 pdb_id[:10] + ':' + chain_id in detailed_lines[1]]

    if not hits_list:
        logging.info(f'No hits in the alignment of the chain {chain_id}. Skipping chain.')
        return None, None
    detailed_lines = hits_list[0]

    file_id = f'{pdb_id.lower()}'
    hit = parsers._parse_hhr_hit(detailed_lines)

    match = re.findall(r'No 1.*[\r\n]+.*\n+(.*\n)', hhr_text)
    identities = re.findall(r'Identities=+([0-9]+)', match[0])[0]
    aligned_columns = re.findall(r'Aligned_cols=+([0-9]+)', match[0])[0]
    total_columns = len(hit.hit_sequence)
    evalue = re.findall(r'E-value=+(.*?) ', match[0])[0]
    logging.info(f'Alignment results:')
    logging.info(
        f'Aligned columns: {aligned_columns} ({len(hit.hit_sequence)}), Evalue: {evalue}, Identities: {identities}')

    template_sequence = hit.hit_sequence.replace('-', '')

    mapping = templates._build_query_to_hit_index_mapping(
        hit.query, hit.hit_sequence, hit.indices_hit, hit.indices_query,
        query_sequence)

    mmcif_string = open(cif_path).read()
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

    return template_features, mapping, identities, aligned_columns, total_columns, evalue


def extract_template_features_from_aligned_pdb_and_sequence(query_sequence: str, pdb_path: str, pdb_id: str,
                                                            chain_id: str):
    # WARNING: input PDB must be aligned to the MSA part in features #

    seq_length = len(query_sequence)

    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure(pdb_id, pdb_path)
    except Exception as e:
        raise Exception(f'The template {pdb_id} could not be aligned and inserted to the features file')

    template_sequence = '-' * seq_length
    template_res_list = [res for res in Selection.unfold_entities(structure, "R")
                         if res.get_parent().id == chain_id and res.id[0] != 'W']

    for res in template_res_list:
        if res.resname != 'X' and res.resname != '-':
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

    template_domain_names = np.array([pdb_id.encode('ascii')])

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

    template_features = {'template_sequence': template_sequence_to_add,
                         'template_all_atom_masks': template_all_atom_masks_to_add,
                         'template_all_atom_positions': template_all_atom_positions_to_add,
                         'template_domain_names': template_domain_names_to_add,
                         'template_aatype': template_aatype_to_add,
                         'template_sum_probs': np.array([template_sum_probs_to_add])}

    return template_features


def write_templates_in_features(template_features: Dict, output_dir: str, chain='A', print_number=True) -> Dict:
    templates_dict = {}

    for pdb_name in template_features['template_domain_names']:
        pdb = pdb_name.decode('utf-8')
        number = '1' if print_number else ''
        pdb_path = os.path.join(output_dir, f'{pdb}{number}.pdb')
        templates_dict[utils.get_file_name(pdb_path)] = pdb_path
        with open(pdb_path, 'w') as output_pdb:
            template_domain_index = np.where(template_features['template_domain_names'] == pdb_name)[0][0]
            atom_num_int = 0
            for index, atoms_mask in enumerate(template_features['template_all_atom_masks'][template_domain_index][:]):
                template_residue_masks = template_features['template_aatype'][template_domain_index][index]
                template_residue_masks_index = np.where(template_residue_masks == 1)[0][0]
                res_type = ID_TO_HHBLITS_AA_3LETTER_CODE[template_residue_masks_index]
                list_of_atoms_in_residue = [order_atom[i] for i, atom in enumerate(atoms_mask) if atom == 1]
                for atom in list_of_atoms_in_residue:
                    atom_num_int = atom_num_int + 1
                    atom_remark = 'ATOM'
                    atom_num = str(atom_num_int)
                    atom_name = atom
                    res_name = res_type
                    res_num = str(index + 1)
                    x_coord = str('%8.3f' % (float(str(
                        template_features['template_all_atom_positions'][template_domain_index][index][
                            atom_types.index(atom)][
                            0]))))
                    y_coord = str('%8.3f' % (float(str(
                        template_features['template_all_atom_positions'][template_domain_index][index][
                            atom_types.index(atom)][
                            1]))))
                    z_coord = str('%8.3f' % (float(str(
                        template_features['template_all_atom_positions'][template_domain_index][index][
                            atom_types.index(atom)][
                            2]))))
                    occ = '1.0'
                    bfact = '25.0'
                    atom_type = atom[0]
                    atom_line = bioutils.get_atom_line(remark=atom_remark, num=int(atom_num), name=atom_name,
                                                       res=res_name, chain=chain, resseq=res_num, x=float(x_coord),
                                                       y=float(y_coord), z=float(z_coord), occ=occ, bfact=bfact,
                                                       atype=atom_type)
                    output_pdb.write(atom_line)

    return templates_dict


def print_features_from_file(pkl_in_path: str):
    with open(f"{pkl_in_path}", "rb") as input_file:
        features_dict = pickle.load(input_file)

    for key in features_dict.keys():
        try:
            logging.info(f'{key} {features_dict[key].shape}')
        except Exception as e:
            pass

    logging.info('\n')
    logging.info('MSA:')
    for num, name in enumerate(features_dict['msa']):
        logging.info(num)
        logging.info('\n')
        logging.info(''.join([residue_constants.ID_TO_HHBLITS_AA[res] for res in features_dict['msa'][num].tolist()]))
        logging.info('\n')

    logging.info('TEMPLATES:')
    for num, seq in enumerate(features_dict['template_sequence']):
        logging.info(f'{features_dict["template_domain_names"][num].decode("utf-8")}:\n')
        for i in range(4):
            logging.info('\t' + ''.join(np.array_split(list(seq.decode('utf-8')), 4)[i].tolist()))
        logging.info('\n')


def create_features_from_file(pkl_in_path: str) -> Features:
    # Read features.pkl and generate a features class

    with open(f"{pkl_in_path}", "rb") as input_file:
        features_dict = pickle.load(input_file)
    new_features = Features(query_sequence=features_dict['sequence'][0].decode('utf-8'))
    for i in range(1, len(features_dict['msa'])):
        sequence = (''.join([residue_constants.ID_TO_HHBLITS_AA[res] for res in features_dict['msa'][i].tolist()]))
        new_features.append_row_in_msa(sequence=sequence,
                                       sequence_id=str(i))


    for i in range(0, len(features_dict['template_sequence'])):
        template_dict = {
            'template_all_atom_positions': np.array([features_dict['template_all_atom_positions'][i]]),
            'template_all_atom_masks': np.array([features_dict['template_all_atom_masks'][i]]),
            'template_aatype': np.array([features_dict['template_aatype'][i]]),
            'template_sequence': np.array([features_dict['template_sequence'][i]]),
            'template_domain_names': np.array([features_dict['template_domain_names'][i]]),
            'template_sum_probs': np.array([features_dict['template_sum_probs'][i]])
        }
        new_features.append_new_template_features(template_dict)

    return new_features
