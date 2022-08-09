from audioop import rms
import copy
import os
import logging
import shutil
import string
from ALPHAFOLD.alphafold.data import pipeline
from libs import arcimboldo_air, bioutils, features, hhsearch, utils
from typing import Dict, List

class Template:

    def __init__ (self, parameters_list: List, output_dir: str):
        self.pdb_path: str
        self.pdb_id: str
        self.template_path: str
        self.chain: str = ''
        self.polyala_res_list: List = []
        self.add_to_msa: bool = False
        self.add_to_templates: bool = False
        self.sum_prob: bool = False
        self.aligned: bool = False
        self.template_features: Dict = None
        
        self.pdb_path = self.__check_pdb(parameters_list['pdb'], output_dir)
        self.pdb_id = utils.get_path_name(self.pdb_path)
        self.chain = parameters_list.get('chain', self.chain)
        self.polyala_res_list = parameters_list.get('polyala_res_list', self.polyala_res_list)
        self.add_to_msa = parameters_list.get('add_to_msa', self.add_to_msa)
        self.add_to_templates = parameters_list.get('add_to_templates', self.add_to_templates)
        self.sum_prob = parameters_list.get('sum_prob', self.sum_prob)
        self.aligned = parameters_list.get('aligned', self.aligned)
        self.template_path = f'{output_dir}/{self.pdb_id}_template.pdb'
    
    def __check_pdb(self, pdb: str, output_dir: str) -> str:
        if not os.path.exists(pdb) and output_dir != '':
            if os.path.exists(f'{output_dir}/{pdb}.pdb'):
                logging.info(f'{output_dir}/{pdb}.pdb already exists in {output_dir}.')
            else:
                bioutils.download_pdb(pdb_id=pdb, output_dir=output_dir)
            pdb = f'{output_dir}/{pdb}.pdb'           
        return pdb
    
    def generate_features(self, a_air):

        output_hhr = f'{a_air.output_dir}/output.hhr'
        pdb70_path = f'{a_air.output_dir}/pdb70'

        if a_air.num_of_copies == 1:
            print(f'Looking for alignment between {self.pdb_id} and query sequence using hhsearch:')

            if not self.aligned:
                bioutils.pdb2mmcif(output_dir=a_air.output_dir, pdb_in_path=self.pdb_path, cif_out_path=f'{a_air.output_dir}/{self.pdb_id}.cif')
                #bioutils.pdb2cif(pdb_id=self.pdb_id, pdb_in_path=self.pdb_path, cif_out_path=f'{a_air.output_dir}/{self.pdb_id}.cif')
                hhsearch.generate_hhsearch_db(pdb_id=self.pdb_id, template_cif_path=f'{a_air.output_dir}/{self.pdb_id}.cif', output_dir=a_air.output_dir)
                hhsearch.run_hhsearch(fasta_path=a_air.fasta_path, pdb70_db=pdb70_path, output_path=output_hhr)

                template_features = features.extract_template_features_from_pdb(
                                            query_sequence=a_air.query_sequence,
                                            hhr_path = output_hhr,
                                            pdb_id=self.pdb_id,
                                            chain_id=self.chain,
                                            mmcif_db=a_air.output_dir)

            else:
                shutil.copyfile(self.pdb_path, self.template_path)

                template_features = features.extract_template_features_from_aligned_pdb_and_sequence(
                    query_sequence=a_air.query_sequence,
                    pdb_path=self.template_path,
                    chain_id=self.chain,
                    pdb_id=self.pdb_id)

            self.template_features = copy.deepcopy(template_features)

        else:
            if not self.aligned:
                bioutils.pdb2mmcif(output_dir=a_air.output_dir, pdb_in_path=self.pdb_path, cif_out_path=f'{a_air.output_dir}/{self.pdb_id}.cif')
                #bioutils.pdb2cif(pdb_in_path=self.pdb_path, cif_out_path=f'{a_air.output_dir}/{self.pdb_id}.cif')
                hhsearch.generate_hhsearch_db(pdb_id=self.pdb_id, template_cif_path=f'{a_air.output_dir}/{self.pdb_id}.cif', output_dir=a_air.output_dir)

                hhsearch.run_hhsearch(fasta_path=a_air.fasta_path, pdb70_db=pdb70_path,
                            output_path=output_hhr)

                chain_list, transformations_list = bioutils.read_remark_350(pdb_path=f'{a_air.output_dir}/{self.pdb_id}.cif')
                new_chain_list = list(string.ascii_uppercase)[:len(transformations_list) * len(chain_list)]
        
                print('Assembly can be build using chain(s)', *chain_list, 'by applying the following transformations:')
                for matrix in transformations_list:
                    print(*matrix)
                
                if len(new_chain_list) != a_air.num_of_copies:
                    raise Exception(f'Assembly description from REMARK 350 contains {len(new_chain_list)} subunits. Please, try to '
                        f'generate a new REMARK 350 (manually or using e.g. PISA) for considering a new assembly with '
                        f'{a_air.num_of_copies} subunits')

                query_seq_length = len(a_air.query_sequence)
                counter = 0
                g = features.Features(query_sequence=a_air.query_sequence)

                for chain in chain_list:
                    template_features = features.extract_template_features_from_pdb(
                        query_sequence=a_air.query_sequence,
                        hhr_path=output_hhr,
                        pdb_id=self.pdb_id,
                        chain_id=chain,
                        mmcif_db=a_air.output_dir)
                                            
                    g.append_new_template_features(new_template_features=template_features, custom_sum_prob=self.sum_prob)
                    g.write_all_templates_in_features(output_path=a_air.output_dir)
                    for transformation in transformations_list:
                        counter += 1
                        bioutils.rot_and_trans(pdb_path=f'{a_air.output_dir}/{self.pdb_id}{chain}_template.pdb', 
                                    out_pdb_path=f'{a_air.output_dir}/{counter}.pdb',
                                    rot_tra_matrix=transformation)

                        if counter == 1:
                            bioutils.change_chain_and_apply_offset_in_single_chain(pdb_in_path=f'{a_air.output_dir}/{counter}.pdb',
                                                                        pdb_out_path=f'{a_air.output_dir}/{new_chain_list[counter - 1]}.pdb',
                                                                        offset=None,
                                                                        chain='A')
                        else:
                            offset = query_seq_length * (counter - 1) + 50 * (counter - 1)  # 50 glycines as offset !!!
                            bioutils.change_chain_and_apply_offset_in_single_chain(pdb_in_path=f'{a_air.output_dir}/{counter}.pdb',
                                                                        pdb_out_path=f'{a_air.output_dir}/{new_chain_list[counter - 1]}.pdb',
                                                                        offset=offset,
                                                                        chain='A')
                list_of_paths_of_pdbs_to_merge = [f'{a_air.output_dir}/{ch}.pdb' for ch in new_chain_list]
                bioutils.merge_pdbs(list_of_paths_of_pdbs_to_merge=list_of_paths_of_pdbs_to_merge,
                        merged_pdb_path=self.template_path)

            else:
                shutil.copyfile(self.pdb_path, self.template_path)

            template_features = features.extract_template_features_from_aligned_pdb_and_sequence(
                query_sequence=a_air.query_sequence_assembled,
                pdb_path=self.template_path,
                chain_id='A',
                pdb_id=self.pdb_id)
            
            self.template_features = copy.deepcopy(template_features)
        
    def __repr__(self):
        return f' \
        pdb_path: {self.pdb_path} \n \
        pdb_id: {self.pdb_id} \n \
        template_path: {self.template_path} \n \
        chain: {self.chain} \n \
        polyala_res_list: {self.polyala_res_list} \n \
        add_to_msa: {self.add_to_msa} \n \
        add_to_templates: {self.add_to_templates} \n \
        sum_prob: {self.sum_prob} \n \
        aligned: {self.aligned} \n \
        template_features: {self.template_features}'