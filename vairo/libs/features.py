import logging
import os
import pickle
import re
import tempfile
from typing import Dict, List, Union, Any
import numpy as np
from Bio.PDB import PDBParser, Selection
from alphafold.common import residue_constants
from alphafold.data import parsers, templates, mmcif_parsing, pipeline, msa_identifiers
from libs import bioutils, utils
from libs.global_variables import ATOM_TYPES, ID_TO_HHBLITS_AA_3LETTER_CODE, ORDER_ATOM


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
        self.extra_info = {'msa_coverage': [],
                           'templates_coverage': [],
                           'num_msa': 0,
                           'num_templates:': 0
                           }

    def modify_query_sequence(self, new_query_sequence: str):
        self.query_sequence = new_query_sequence
        self.sequence_features = pipeline.make_sequence_features(sequence=self.query_sequence,
                                                                 description='Query',
                                                                 num_res=len(self.query_sequence))

    def append_new_template_features(self, new_template_features: Dict, custom_sum_prob: int = None) -> Dict:
        self.template_features['template_all_atom_positions'] = np.vstack(
            [self.template_features['template_all_atom_positions'],
             new_template_features['template_all_atom_positions']])
        self.template_features['template_all_atom_masks'] = np.vstack(
            [self.template_features['template_all_atom_masks'],
             new_template_features['template_all_atom_masks']])
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

    def append_row_in_msa_from_features(self, new_msa_features: Dict):
        self.msa_features['msa'] = np.vstack([self.msa_features['msa'], new_msa_features['msa']])
        self.msa_features['accession_ids'] = np.hstack(
            [self.msa_features['accession_ids'], new_msa_features['accession_ids']])
        self.msa_features['deletion_matrix_int'] = np.vstack(
            [self.msa_features['deletion_matrix_int'], new_msa_features['deletion_matrix_int']])
        self.msa_features['msa_species_identifiers'] = np.hstack(
            [self.msa_features['msa_species_identifiers'], new_msa_features['msa_species_identifiers']])
        self.msa_features['num_alignments'] = np.full(self.msa_features['num_alignments'].shape,
                                                      len(self.msa_features['msa']))

    def delete_templates(self, rows_to_delete: List):
        self.template_features['template_all_atom_positions'] = np.delete(
            self.template_features['template_all_atom_positions'], rows_to_delete, axis=0)
        self.template_features['template_all_atom_masks'] = np.delete(self.template_features['template_all_atom_masks'],
                                                                      rows_to_delete, axis=0)
        self.template_features['template_aatype'] = np.delete(self.template_features['template_aatype'], rows_to_delete,
                                                              axis=0)
        self.template_features['template_sequence'] = np.delete(self.template_features['template_sequence'],
                                                                rows_to_delete, axis=0)
        self.template_features['template_domain_names'] = np.delete(self.template_features['template_domain_names'],
                                                                    np.array(rows_to_delete), axis=0)
        self.template_features['template_sum_probs'] = np.delete(self.template_features['template_sum_probs'],
                                                                 rows_to_delete, axis=0)

    def delete_msas(self, rows_to_delete: List):
        self.msa_features['msa'] = np.delete(self.msa_features['msa'], rows_to_delete, axis=0)
        self.msa_features['accession_ids'] = np.delete(self.msa_features['accession_ids'], rows_to_delete, axis=0)
        self.msa_features['deletion_matrix_int'] = np.delete(self.msa_features['deletion_matrix_int'], rows_to_delete,
                                                             axis=0)
        self.msa_features['msa_species_identifiers'] = np.delete(self.msa_features['msa_species_identifiers'],
                                                                 rows_to_delete, axis=0)
        self.msa_features['num_alignments'] = np.full(self.msa_features['num_alignments'].shape,
                                                      len(self.msa_features['msa']))

    def append_row_in_msa(self, sequence_in: str, sequence_id: str, position: int = None):
        seq_length = len(self.msa_features['msa'][0])
        if position:
            sequence_in = '-' * (position - 1) + sequence_in + '-' * (seq_length - (position - 1) - len(sequence_in))
        sequence_in = sequence_in[:seq_length]
        sequence_array = np.array([residue_constants.HHBLITS_AA_TO_ID[res] for res in sequence_in])
        self.msa_features['msa'] = np.vstack([self.msa_features['msa'], sequence_array])
        self.msa_features['accession_ids'] = np.hstack([self.msa_features['accession_ids'], sequence_id.encode()])
        self.msa_features['deletion_matrix_int'] = np.vstack(
            [self.msa_features['deletion_matrix_int'], np.zeros(self.msa_features['msa'].shape[1])])
        self.msa_features['msa_species_identifiers'] = np.hstack([self.msa_features['msa_species_identifiers'], ''])
        self.msa_features['num_alignments'] = np.full(self.msa_features['num_alignments'].shape,
                                                      len(self.msa_features['msa']))

    def write_all_templates_in_features(self, output_dir: str, chain='A', print_number=True) -> Dict:
        return write_templates_in_features(self.template_features, output_dir, chain, print_number)

    def write_pkl(self, pkl_path: str):
        logging.error(f'Writing all input information to {pkl_path}')
        merged_features = self.merge_features()
        with open(pkl_path, 'wb') as f_out:
            pickle.dump(merged_features, f_out, protocol=pickle.HIGHEST_PROTOCOL)

    def get_names_templates(self) -> List[str]:
        return [x.decode() for x in self.template_features['template_domain_names']]

    def get_names_msa(self) -> List[str]:
        return [x.decode() for x in self.msa_features['accession_ids']]

    def get_msa_by_name(self, name: str) -> Union[str, None]:
        index = np.where(self.msa_features['accession_ids'] == name.encode())[0]
        if index:
            return bioutils.convert_msa_sequence(self.msa_features['msa'][index[0]].tolist())
        return None

    def get_msa_sequences(self) -> List[str]:
        results_list = []
        for msa_seq in self.msa_features['msa']:
            results_list.append(bioutils.convert_msa_sequence(msa_seq))
        return results_list

    def get_index_by_name(self, name: str) -> Union[int, None]:
        index = np.where(self.template_features['template_domain_names'] == name.encode())
        if index[0].size != 0:
            return index[0][0]
        else:
            return None

    def get_index_msa_by_name(self, name: str) -> Union[int, None]:
        index = np.where(self.msa_features['accession_ids'] == name.encode())
        if index[0].size != 0:
            return index[0][0]
        else:
            return None

    def get_msa_length(self) -> int:
        return len(self.msa_features['msa'])

    def get_templates_length(self) -> int:
        return len(self.template_features['template_sequence'])

    def get_sequence_by_name(self, name: str) -> str:
        index = self.get_index_by_name(name)
        return self.template_features['template_sequence'][index].decode()

    def merge_features(self) -> Dict[Any, Any]:
        logging.info(f'Merging sequence, msa and template features!')
        return {**self.sequence_features, **self.msa_features, **self.template_features, **self.extra_info}

    def get_template_by_index(self, index: int) -> Dict:
        template_dict = {
            'template_all_atom_positions': np.array([self.template_features['template_all_atom_positions'][index]]),
            'template_all_atom_masks': np.array([self.template_features['template_all_atom_masks'][index]]),
            'template_aatype': np.array([self.template_features['template_aatype'][index]]),
            'template_sequence': np.array([self.template_features['template_sequence'][index]]),
            'template_domain_names': np.array([self.template_features['template_domain_names'][index]]),
            'template_sum_probs': np.array([self.template_features['template_sum_probs'][index]])
        }
        return template_dict

    def get_max_id(self) -> int:
        starting_list = [int(x) for x in self.get_names_msa() if x.isdigit()]
        starting_id = 1
        if starting_list:
            starting_id = max(starting_list) + 1
        return starting_id

    def set_msa_features(self, new_msa: Dict, start: int = 1, finish: int = -1,
                         delete_positions: List[int] = []) -> int:
        coverage_msa = []
        self.delete_residues_msa(delete_positions=delete_positions)
        for i in range(start, len(new_msa['msa'])):
            coverage_msa.append(len([residue for residue in new_msa['msa'][i] if residue != 21]))
        if finish == -1:
            coverage_msa = list(range(0, len(new_msa['msa']) - start))
        else:
            arr = np.array(coverage_msa)
            coverage_msa = arr.argsort()[-finish:][::-1]
            coverage_msa = np.sort(coverage_msa)
        msa_dict = create_empty_msa_list(len(coverage_msa))
        starting_id = self.get_max_id()
        for i, num in enumerate(coverage_msa):
            msa_dict['msa'][i] = new_msa['msa'][num + start]
            msa_dict['accession_ids'][i] = str(starting_id).encode()
            msa_dict['deletion_matrix_int'][i] = new_msa['deletion_matrix_int'][num + start]
            msa_dict['msa_species_identifiers'][i] = new_msa['msa_species_identifiers'][num + start]
            msa_dict['num_alignments'][i] = np.zeros(new_msa['num_alignments'].shape)
            starting_id += 1

        if len(msa_dict['msa']) > 0:
            self.append_row_in_msa_from_features(msa_dict)
        return len(msa_dict['msa'])

    def set_template_features(self, new_templates: Dict, finish: int = -1) -> int:
        finish = len(new_templates['template_sequence']) if finish == -1 else finish
        template_dict = create_empty_template_list(finish)
        for i in range(finish):
            template_name = new_templates['template_domain_names'][i].decode()
            index = self.get_index_by_name(name=template_name)
            if index is not None:
                j = 1
                while True:
                    new_name = f'{template_name}_{j}'
                    index2 = self.get_index_by_name(name=new_name)
                    if index2 is None:
                        template_name = new_name
                        break
                    j += 1
            template_dict['template_all_atom_positions'][i] = new_templates['template_all_atom_positions'][i]
            template_dict['template_all_atom_masks'][i] = new_templates['template_all_atom_masks'][i]
            template_dict['template_aatype'][i] = new_templates['template_aatype'][i]
            template_dict['template_sequence'][i] = new_templates['template_sequence'][i]
            template_dict['template_domain_names'][i] = template_name.encode()
            template_dict['template_sum_probs'][i] = new_templates['template_sum_probs'][i]

        if len(template_dict['template_all_atom_positions']) > 0:
            self.append_new_template_features(template_dict)
        return len(template_dict['template_all_atom_positions'])

    def cut_expand_features(self, query_sequence: str, modifications_list: List[int]):
        new_features = Features(query_sequence=query_sequence)
        msa_dict = create_empty_msa_list(self.get_msa_length() - 1)
        ext_len = len(query_sequence)
        # We skip the first one, so that's the -1, because it is the query sequence one, and we have created
        # another new_features with a fake query_sequence. This will be skipped in set_msa_features
        for i in range(self.get_msa_length() - 1):
            msa_dict['msa'][i] = np.full(ext_len, 21)
            msa_dict['deletion_matrix_int'][i] = np.full(ext_len, 0)
            msa_dict['msa_species_identifiers'][i] = self.msa_features['msa_species_identifiers'][i + 1]
            for j, mod in enumerate(modifications_list):
                if mod is not None and mod <= len(self.query_sequence):
                    msa_dict['msa'][i][j] = self.msa_features['msa'][i + 1][mod - 1]
                    msa_dict['deletion_matrix_int'][i][j] = self.msa_features['deletion_matrix_int'][i + 1][mod - 1]
            msa_dict['accession_ids'][i] = self.msa_features['accession_ids'][i + 1]
            msa_dict['num_alignments'][i] = np.full(self.msa_features['num_alignments'].shape,
                                                    len(self.msa_features['msa']))

        if len(msa_dict['msa']) > 0:
            new_features.append_row_in_msa_from_features(msa_dict)

        template_dict = create_empty_template_list((self.get_templates_length()))
        for i in range(self.get_templates_length()):
            template_dict['template_all_atom_positions'][i] = np.zeros(
                (ext_len, residue_constants.atom_type_num, 3))
            template_dict['template_all_atom_masks'][i] = np.zeros(
                (ext_len, residue_constants.atom_type_num))
            template_dict['template_aatype'][i] = residue_constants.sequence_to_onehot('A' * ext_len,
                                                                                       residue_constants.HHBLITS_AA_TO_ID)
            template_dict['template_domain_names'][i] = self.template_features['template_domain_names'][i]
            template_dict['template_sum_probs'][i] = self.template_features['template_sum_probs'][i]
            template_sequence = list(('-' * ext_len))
            for j, mod in enumerate(modifications_list):
                if mod is not None and mod <= len(self.query_sequence):
                    template_dict['template_all_atom_positions'][i][j] = \
                        self.template_features['template_all_atom_positions'][i][mod - 1]
                    template_dict['template_all_atom_masks'][i][j] = \
                        self.template_features['template_all_atom_masks'][i][mod - 1]
                    template_dict['template_aatype'][i][j] = \
                        self.template_features['template_aatype'][i][mod - 1]
                    template_sequence[j] = self.template_features['template_sequence'][i].decode()[mod - 1]
            template_dict['template_sequence'][i] = ''.join(template_sequence).encode()

        if len(template_dict['template_all_atom_positions']) > 0:
            new_features.append_new_template_features(template_dict)

        return new_features

    def slice_features(self, ini: int, end: int) -> List:
        new_features = Features(query_sequence=self.query_sequence[ini:end])
        msa_dict = create_empty_msa_list(self.get_msa_length() - 1)
        for i in range(self.get_msa_length() - 1):
            msa_dict['msa'][i] = self.msa_features['msa'][i + 1][ini:end]
            msa_dict['accession_ids'][i] = str(i).encode()
            msa_dict['deletion_matrix_int'][i] = self.msa_features['deletion_matrix_int'][i + 1][
                                                 ini:end]
            msa_dict['msa_species_identifiers'][i] = self.msa_features['msa_species_identifiers'][i + 1]
            msa_dict['num_alignments'][i] = np.zeros(self.msa_features['num_alignments'].shape)
        if len(msa_dict['msa']) > 0:
            new_features.append_row_in_msa_from_features(msa_dict)

        template_dict = create_empty_template_list((self.get_templates_length()))
        for i in range(self.get_templates_length()):
            template_dict['template_all_atom_positions'][i] = self.template_features['template_all_atom_positions'][
                                                                  i][ini:end]
            template_dict['template_all_atom_masks'][i] = self.template_features['template_all_atom_masks'][i][
                                                          ini:end]
            template_dict['template_aatype'][i] = self.template_features['template_aatype'][i][ini:end]
            template_dict['template_sequence'][i] = self.template_features['template_sequence'][i][
                                                    ini:end]
            template_dict['template_domain_names'][i] = self.template_features['template_domain_names'][i]
            template_dict['template_sum_probs'][i] = self.template_features['template_sum_probs'][i]
        if len(template_dict['template_all_atom_positions']) > 0:
            new_features.append_new_template_features(template_dict)
        return new_features

    def slicing_features(self, chunk_list: List) -> List:
        # This function will generate as many features
        # as required per size. It will return a list with
        # the path of all the generated features
        features_list = []
        for start_min, start_max in chunk_list:
            features_list.append(self.slice_features(start_min, start_max))
        if len(chunk_list) > 1:
            logging.error(
                f'Query sequence and the input information has been cut into {len(features_list)} partitions with the following sizes:')
            for start_min, start_max in chunk_list:
                logging.error(f'      - {start_min}-{start_max}')
        else:
            logging.error(f'Query sequence has the following size: {chunk_list[0][0]}-{chunk_list[0][1]}')
        return features_list

    def select_msa_templates(self, sequence_assembled, minimum_percentage: float = 0.50):
        # Trim the templates that has a 50% percentage of it in the glycines part
        # Trim the msa sequences that has a 50% percentage of it in the glycines part

        delete_msa = []
        for i in range(self.get_msa_length()):
            sequence_in = self.msa_features['msa'][i]
            res_num, perc = sequence_assembled.get_percentage_sequence(sequence_in)
            length_sequence = len(sequence_in[sequence_in != 21]) * minimum_percentage
            if sum(res_num) < length_sequence:
                delete_msa.append(i)
        if delete_msa:
            logging.error(f'{len(delete_msa)} sequences filtered from the MSA due to not enough sequence coverage')
            self.delete_msas(delete_msa)
        delete_templates = []
        for i in range(self.get_templates_length()):
            sequence_in = self.template_features['template_sequence'][i].decode()
            res_num, perc = sequence_assembled.get_percentage_sequence(sequence_in)
            length_sequence = len(sequence_in[sequence_in != 21]) * minimum_percentage
            if sum(res_num) < length_sequence:
                logging.error(
                    f'Template {self.template_features["template_domain_names"][i].decode()} has been filtered:')
                logging.error(f'    Not enough sequence coverage')
                delete_templates.append(i)

        if delete_templates:
            self.delete_templates(delete_templates)

        self.delete_linkers_regions(sequence_assembled)


    def delete_linkers_regions(self, sequence_assembled):
        delete_positions = sequence_assembled.get_list_linker_numbering()
        self.delete_residues_msa(delete_positions=delete_positions, starting=1)


    def delete_by_id(self, id_list: List[str]):
        # Given a list of ids, check if it belongs to a msa or a templates. Delete them.
        delete_templates_list = []
        delete_msa_list = []

        for id_delete in id_list:
            index = self.get_index_msa_by_name(id_delete)
            if index is None:
                index = self.get_index_by_name(id_delete)
                if index is None:
                    continue
                else:
                    delete_templates_list.append(index)
            else:
                delete_msa_list.append(index)

        if delete_msa_list:
            self.delete_msas(delete_msa_list)

        if delete_templates_list:
            self.delete_templates(delete_templates_list)

    def delete_by_range(self, min_identity: float, max_identity: float):
        # Given a minimum identity value, and a maximum identity value, delete them the outlayers.
        delete_msa = []
        for i in range(self.get_msa_length()):
            sequence_in = self.msa_features['msa'][i]
            identity = bioutils.sequence_identity(self.query_sequence, sequence_in)
            if min_identity < identity > max_identity:
                delete_msa.append(i)
        if delete_msa:
            self.delete_msas(delete_msa)
        delete_templates = []
        for i in range(self.get_templates_length()):
            sequence_in = self.template_features['template_sequence'][i].decode()
            identity = bioutils.sequence_identity(self.query_sequence, sequence_in)
            if min_identity < identity > max_identity:
                delete_templates.append(i)
        if delete_templates:
            self.delete_templates(delete_templates)

    def set_extra_info(self):
        self.extra_info['num_templates'] = self.get_templates_length()
        if self.get_templates_length() > 0:
            seq_templates = [seq.decode() for seq in self.template_features['template_sequence']]
            self.extra_info['templates_coverage'] = bioutils.calculate_coverage_scaled(query_seq=self.query_sequence,
                                                                                       sequences=seq_templates)
        else:
            self.extra_info['templates_coverage'] = [0] * len(self.query_sequence)

        self.extra_info['num_msa'] = self.get_msa_length() - 1
        if self.get_msa_length() - 1 > 0:
            seq_msa = [''.join(residue_constants.ID_TO_HHBLITS_AA[res] for res in msa.tolist()) for msa in
                       self.msa_features['msa'][1:]]
            self.extra_info['msa_coverage'] = bioutils.calculate_coverage_scaled(query_seq=self.query_sequence,
                                                                                 sequences=seq_msa)
        else:
            self.extra_info['msa_coverage'] = [0] * len(self.query_sequence)

    def delete_residues_msa(self, delete_positions: List[int], starting: int = 0):
        # Delete the specifics residues in the msa.
        if delete_positions:
            for i in range(starting, self.get_msa_length()):
                for delete in delete_positions:
                    if delete <= self.msa_features['msa'][i].size:
                        self.msa_features['msa'][i][delete - 1] = 21
                        self.msa_features['deletion_matrix_int'][i][delete - 1] = 0
                    else:
                        break

    def replace_sequence_template(self, sequence_in: str):
        # Replace the sequence for the new sequence
        if sequence_in:
            seq_bytes = sequence_in.encode()
            seq_len = len(sequence_in)
            for i in range(self.get_templates_length()):
                target = self.template_features['template_sequence'][i]
                self.template_features['template_sequence'][i] = (
                        seq_bytes + target[seq_len:]
                )
                for index in range(seq_len):
                    aa_container = [0] * 22
                    aa_container[residue_constants.HHBLITS_AA_TO_ID[sequence_in[index]]] = 1
                    self.template_features['template_aatype'][i][index] = aa_container


    def mutate_residues(self, mutation_dict: dict):
        if mutation_dict:
            for i in range(self.get_msa_length()):
                for res, numbering in mutation_dict.items():
                    for value in numbering:
                        self.msa_features['msa'][i, value - 1] = residue_constants.HHBLITS_AA_TO_ID[res]

            for i in range(self.get_templates_length()):
                sequence_in = list(self.template_features['template_sequence'][i].decode())
                for res, numbering in mutation_dict.items():
                    for value in numbering:
                        sequence_in[value - 1] = res
                        aa_container = [0] * 22
                        aa_container[residue_constants.HHBLITS_AA_TO_ID[res]] = 1
                        self.template_features['template_aatype'][i][value - 1] = aa_container
                self.template_features['template_sequence'][i] = "".join(sequence_in).encode()


def create_empty_msa_list(length: int) -> Dict:
    msa_dict = {
        'msa': [None] * length,
        'accession_ids': [None] * length,
        'deletion_matrix_int': [None] * length,
        'msa_species_identifiers': [None] * length,
        'num_alignments': [None] * length
    }
    return msa_dict


def create_empty_template_list(length: int) -> Dict:
    template_dict = {
        'template_all_atom_positions': [None] * length,
        'template_all_atom_masks': [None] * length,
        'template_aatype': [None] * length,
        'template_sequence': [None] * length,
        'template_domain_names': [None] * length,
        'template_sum_probs': [None] * length
    }
    return template_dict


def empty_msa_features(query_sequence):
    # Generate an empty msa, containing one element in the msa (the sequence)
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
        for sequence_index, sequence_in in enumerate(msa.sequences):
            if sequence_in in seen_sequences:
                continue
            seen_sequences.add(sequence_in)
            int_msa.append(
                [residue_constants.HHBLITS_AA_TO_ID[res] for res in sequence_in])
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


def extract_template_features_from_pdb(query_sequence: str, hhr_path: str, cif_path: str, sequence_id: str,
                                       chain_id: str) -> List[str]:
    pdb_id = utils.get_file_name(cif_path)
    hhr_text = open(hhr_path, 'r').read()
    matches = re.finditer(r'No\s+\d+', hhr_text)
    matches_positions = [match.start() for match in matches] + [len(hhr_text)]

    detailed_lines_list = []
    for i in range(len(matches_positions) - 1):
        detailed_lines_list.append(hhr_text[matches_positions[i]:matches_positions[i + 1]].split('\n')[:-3])

    hits_list = [detailed_lines for detailed_lines in detailed_lines_list if
                 sequence_id + ':' + chain_id in detailed_lines[1]]

    if not hits_list:
        logging.error(f'No hits in the alignment of the chain {chain_id}.')
        return None, None, None, 0, 0, 0
    detailed_lines = hits_list[0]

    try:
        hit = parsers._parse_hhr_hit(detailed_lines)
    except:
        return None, None, None, 0, 0, 0

    file_id = f'{pdb_id.lower()}'
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

    match = re.findall(r'No 1.*[\r\n]+.*\n+(.*\n)', hhr_text)
    identities = re.findall(r'Identities=+([0-9]+)', match[0])[0]
    aligned_columns = re.findall(r'Aligned_cols=+([0-9]+)', match[0])[0]
    total_columns = len(parsing_result.mmcif_object.chain_to_seqres[chain_id])
    evalue = re.findall(r'E-value=+(.*?) ', match[0])[0]

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
            template_sequence = template_sequence[:res.id[1]] + residue_constants.restype_3to1[
                res.resname] + template_sequence[
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
                index = ATOM_TYPES.index(atom)
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
                res_container.append(resi[ATOM_TYPES[j]].coord)
            else:
                res_container.append(np.array([0.] * 3))
        template_container.append(res_container)
    template_all_atom_positions = np.array([template_container])

    template_domain_names = np.array([pdb_id.encode('ascii')])

    template_aatype_container = []
    for res in template_sequence[1:]:
        aa_container = [0] * 22
        aa_container[residue_constants.HHBLITS_AA_TO_ID[res]] = 1
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


def write_template_in_features(template_features: Dict, template_code: str, output_path: str, chain='A') -> bool:
    with open(output_path, 'w') as output_pdb:
        template_domain_index = np.where(template_features['template_domain_names'] == template_code)[0][0]
        atom_num_int = 0
        for index, atoms_mask in enumerate(template_features['template_all_atom_masks'][template_domain_index][:]):
            template_residue_masks = template_features['template_aatype'][template_domain_index][index]
            template_residue_masks_index = np.where(template_residue_masks == 1)[0][0]
            res_type = ID_TO_HHBLITS_AA_3LETTER_CODE[template_residue_masks_index]
            list_of_atoms_in_residue = [ORDER_ATOM[i] for i, atom in enumerate(atoms_mask) if atom == 1]
            for atom in list_of_atoms_in_residue:
                atom_num_int = atom_num_int + 1
                atom_remark = 'ATOM'
                atom_num = str(atom_num_int)
                atom_name = atom
                res_name = res_type
                res_num = str(index + 1)
                x_coord = str('%8.3f' % (float(str(
                    template_features['template_all_atom_positions'][template_domain_index][index][
                        ATOM_TYPES.index(atom)][
                        0]))))
                y_coord = str('%8.3f' % (float(str(
                    template_features['template_all_atom_positions'][template_domain_index][index][
                        ATOM_TYPES.index(atom)][
                        1]))))
                z_coord = str('%8.3f' % (float(str(
                    template_features['template_all_atom_positions'][template_domain_index][index][
                        ATOM_TYPES.index(atom)][
                        2]))))
                occ = '1.0'
                bfact = '25.0'
                atom_type = atom[0]
                atom_line = bioutils.get_atom_line(remark=atom_remark, num=int(atom_num), name=atom_name,
                                                   res=res_name, chain=chain, resseq=res_num, x=float(x_coord),
                                                   y=float(y_coord), z=float(z_coord), occ=occ, bfact=bfact,
                                                   atype=atom_type)
                output_pdb.write(atom_line)

    if os.stat(output_path).st_size == 0:
        return False
    else:
        return True


def write_templates_in_features(template_features: Dict, output_dir: str, chain='A', print_number=True) -> Dict:
    templates_dict = {}
    for pdb_name in template_features['template_domain_names']:
        pdb = pdb_name.decode('utf-8')
        number = '1' if print_number else ''
        pdb_path = os.path.join(output_dir, f'{pdb}{number}.pdb')
        written = write_template_in_features(template_features=template_features, template_code=pdb_name,
                                             output_path=pdb_path,
                                             chain=chain)
        if written:
            templates_dict[utils.get_file_name(pdb_path)] = pdb_path

    return templates_dict


def print_features_from_file(pkl_in_path: str):
    with open(f"{pkl_in_path}", "rb") as input_file:
        features_dict = pickle.load(input_file)
    for key in features_dict.keys():
        try:
            logging.error(f'{key} {features_dict[key].shape}')
        except Exception as e:
            pass
    logging.error('\n')
    logging.error('MSA:')
    for num, name in enumerate(features_dict['msa']):
        logging.error(f'> {num}')
        logging.error(''.join([residue_constants.ID_TO_HHBLITS_AA[res] for res in features_dict['msa'][num].tolist()]))

    logging.error('TEMPLATES:')
    for num, seq in enumerate(features_dict['template_sequence']):
        logging.error(f'{features_dict["template_domain_names"][num].decode("utf-8")}:\n')
        print(features_dict["template_aatype"][num])
        for i in range(4):
            logging.error('\t' + ''.join(np.array_split(list(seq.decode('utf-8')), 4)[i].tolist()))
        logging.error('\n')

    logging.error('INFORMATION:')
    keys = ['num_msa', 'num_templates', 'msa_coverage', 'templates_coverage']
    for key in keys:
        if key in features_dict:
            logging.error(f"{key}: {features_dict[key]}")


def create_features_from_file(pkl_in_path: str) -> Features:
    # Read features.pkl and generate a feature class
    with open(f'{pkl_in_path}', 'rb') as input_file:
        features_dict = pickle.load(input_file)
    new_features = Features(query_sequence=features_dict['sequence'][0].decode('utf-8'))
    new_features.set_msa_features(features_dict)
    new_features.set_template_features(features_dict)
    return new_features


def delete_seq_from_msa(pkl_in_path: str, delete_list: List[str], pkl_out_path: str = None):
    feature = create_features_from_file(pkl_in_path=pkl_in_path)
    feature.delete_msas(delete_list)
    if pkl_out_path is None:
        pkl_out_path = pkl_in_path
    feature.write_pkl(pkl_out_path)


def extract_features_info(pkl_in_path: str, regions_list: List):
    feature = create_features_from_file(pkl_in_path=pkl_in_path)
    features_info_dict = {'msa': {}, 'templates': {}}
    msa_sequences = feature.get_msa_sequences()
    region_query = ''
    for k, msa_seq in enumerate(msa_sequences[1:], start=1):
        identity, region_query, region_msa = bioutils.sequence_identity_regions(feature.query_sequence, msa_seq,
                                                                                regions_list)
        global_identity = bioutils.sequence_identity(feature.query_sequence, msa_seq)
        coverage = (sum(1 for res in msa_seq if res != '-') / len(msa_seq) * 100)
        if identity != 0:
            features_info_dict['msa'][feature.msa_features['accession_ids'][k].decode()] = {
                'global_identity': round(global_identity, 2),
                'identity': round(identity, 2),
                'coverage': round(coverage, 2),
                'seq': msa_seq,
                'seq_query': region_query,
                'seq_msa': region_msa
            }
    for i, template_seq in enumerate(feature.template_features['template_sequence']):
        identity, region_query, region_msa = bioutils.sequence_identity_regions(feature.query_sequence,
                                                                                template_seq.decode(),
                                                                                regions_list)
        global_identity = bioutils.sequence_identity(feature.query_sequence, template_seq.decode())
        coverage = (sum(1 for res in template_seq.decode() if res != '-') / len(template_seq.decode()) * 100)
        pdb_name = feature.template_features['template_domain_names'][i]
        with tempfile.NamedTemporaryFile() as temp_file:
            write_template_in_features(template_features=feature.template_features, template_code=pdb_name,
                                       output_path=temp_file.name)
            pdb_info = temp_file.read().decode()

        features_info_dict['templates'][feature.template_features['template_domain_names'][i].decode()] = {
            'identity': round(identity, 2),
            'coverage': round(coverage, 2),
            'global_identity': round(global_identity, 2),
            'seq': template_seq.decode(),
            'seq_query': region_query,
            'seq_msa': region_msa,
            'pdb': pdb_info
        }

    return features_info_dict, region_query, feature.query_sequence
