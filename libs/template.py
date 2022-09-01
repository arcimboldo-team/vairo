import copy
import os
import logging
from pickle import NONE
import re
import shutil
from Bio.PDB import PDBParser, Selection

from libs import arcimboldo_air, bioutils, features, hhsearch, utils, change_res
from typing import Dict, List


class Template:

    def __init__ (self, parameters_dict: Dict, output_dir: str, input_dir: str, num_of_copies: int):

        self.pdb_path: str
        self.pdb_id: str
        self.template_path: str
        self.chain: str = ''
        self.generate_multimer: bool = True if num_of_copies > 1 else False
        self.change_res: List[change_res.ChangeResidues] = []
        self.add_to_msa: bool = False
        self.add_to_templates: bool = False
        self.sum_prob: bool = False
        self.aligned: bool = False
        self.template_features: Dict = None
        self.merged_templates_list: List = None
        
        self.pdb_path = self.check_pdb(utils.get_mandatory_value(parameters_dict, 'pdb'), input_dir)
        self.pdb_id = utils.get_file_name(self.pdb_path)
        self.add_to_msa = parameters_dict.get('add_to_msa', self.add_to_msa)
        self.add_to_templates = parameters_dict.get('add_to_templates', self.add_to_templates)
        self.sum_prob = parameters_dict.get('sum_prob', self.sum_prob)
        self.aligned = parameters_dict.get('aligned', self.aligned)
        self.template_path = f'{output_dir}/{self.pdb_id}_template.pdb'
        self.chain = parameters_dict.get('chain', self.chain)

        self.generate_multimer = parameters_dict.get('generate_multimer', self.generate_multimer)

        structure = bioutils.get_structure(pdb_path=self.pdb_path)
        chains = [chain.get_id() for chain in structure.get_chains()]

        if num_of_copies == 1:
            if self.chain == '':
                if len(chains) == 1:
                    self.chain = chains.pop()
                else:
                    raise Exception('There is more than one chain in the structure. Select one in the configuration file.')
            else:
                if not self.chain in chains:
                    raise Exception('The choosen chain does not exist in the structure.')
        else:
            self.chain = chains

        for paramaters_change_res in parameters_dict.get('change_res', self.change_res):
            self.change_res.append(change_res.ChangeResidues(paramaters_change_res, chains))

    def check_pdb(self, pdb: str, output_dir: str) -> str:
        #Check if pdb is a path, and if it doesn't exist, download it.
        #If the pdb is a path, copy it to our input folder

        if not os.path.exists(pdb):
            bioutils.download_pdb(pdb_id=pdb, output_dir=output_dir)
            pdb = f'{output_dir}/{pdb}.pdb'
        else:
            pdb_aux = f'{output_dir}/{os.path.basename(pdb)}'
            if pdb != pdb_aux:
                shutil.copy2(pdb, pdb_aux)
                pdb = pdb_aux                

        return pdb
    
    def generate_features(self, a_air):
        #If aligned, skip to the end, it is just necessary to create extract the features
        #If not aligned:
        #   - Convert pdb to cif, this is necessary to run hhsearch.
        #   - Create a fake hhsearch database with just our .cif in it. 
        #   - Run hhsearch, align the .cif with the given sequence.
        #   - For each chain, extract the features from the alignment result, create the
        #       the templates for each chain, each template has the chain 'A'.
        #   - Change the specified residues in the input in all the templates.
        #   - Generate offset.
        #   - Apply the generated offset to all the templates.
        #   - Build the new template merging all the templates.
        #   - Create features for the new template.

        output_hhr = f'{a_air.run_dir}/output.hhr'
        pdb70_path = f'{a_air.run_dir}/pdb70'

        if not self.aligned:
            bioutils.pdb2mmcif(output_dir=a_air.run_dir, pdb_in_path=self.pdb_path, cif_out_path=f'{a_air.run_dir}/{self.pdb_id}.cif')
            hhsearch.generate_hhsearch_db(template_cif_path=f'{a_air.run_dir}/{self.pdb_id}.cif', output_dir=a_air.run_dir)
            hhsearch.run_hhsearch(fasta_path=a_air.fasta_path, pdb70_db=pdb70_path,output_path=output_hhr)

            query_seq_length = len(a_air.query_sequence)
            template_dict = {}
            structure = bioutils.get_structure(pdb_path=self.pdb_path)
            chain_list = bioutils.get_chains(structure)

            for chain in chain_list:
                template_features = features.extract_template_features_from_pdb(
                    query_sequence=a_air.query_sequence,
                    hhr_path=output_hhr,
                    pdb_id=self.pdb_id,
                    chain_id=chain,
                    mmcif_db=a_air.run_dir)

                g = features.Features(query_sequence=a_air.query_sequence)
                g.append_new_template_features(new_template_features=template_features, custom_sum_prob=self.sum_prob)
                aux_dict = g.write_all_templates_in_features(output_dir=a_air.run_dir)
                template_path = list(aux_dict.values())[0]
                template_dict[chain] = [template_path]

            if self.generate_multimer:
                template_dict = self.generate_multimer_chains(template_dict)
            template_list = [val for sublist in template_dict.values() for val in sublist]

            for change_residues in self.change_res:
                for chain, paths_list in template_dict.items():
                    if chain in change_residues.change_dict.keys():
                        for path in paths_list:
                            change_residues.change_residues(pdb_in_path=path, pdb_out_path=path, real_chain=chain)
                                
            merge_list = []
            
            composition_list = self.sort_templates_into_positions(template_list, a_air.num_of_copies)

            for i, pdb_path in enumerate(composition_list):
                if pdb_path is not None:
                    offset = query_seq_length * (i) + a_air.glycines * (i)
                    new_pdb_path = f'{a_air.run_dir}/{utils.get_file_name(pdb_path)}_{offset}.pdb'
                    bioutils.change_chain(pdb_in_path=pdb_path,
                                    pdb_out_path=new_pdb_path,
                                    offset=offset, chain='A')
                    merge_list.append(new_pdb_path)

            bioutils.merge_pdbs(list_of_paths_of_pdbs_to_merge=utils.sort_by_digit(merge_list),
                    merged_pdb_path=self.template_path)

        else:
            shutil.copy2(self.pdb_path, self.template_path)

        template_features = features.extract_template_features_from_aligned_pdb_and_sequence(
            query_sequence=a_air.query_sequence_assembled,
            pdb_path=self.template_path,
            pdb_id=self.pdb_id,
            chain_id='A')
            
        self.template_features = copy.deepcopy(template_features)
        

    def sort_templates_into_positions(template_list: List, num_of_copies: int):




    def generate_multimer_chains(self, template_dict: Dict) -> Dict:
        #Read remark to get the transformations and the new chains
        #Apply transformations to generate the new ones
        #Rename chains with A1, A2...
        #Store a dict with the relation between old chains and new chains
        # Dict -> A: [path_to_A1, path_to_A2]

        chain_list, transformations_list = bioutils.read_remark_350(self.pdb_path)
        multimer_dict = {}
        
        logging.info('Assembly can be build using chain(s) '+ str(chain_list) + ' by applying the following transformations:')
        for matrix in transformations_list:
            logging.info(str(matrix))

        for chain in chain_list:
            pdb_path = template_dict[chain][0]
            multimer_new_chains = []
            for i, transformation in enumerate(transformations_list):
                new_pdb_path = utils.replace_last_number(text=pdb_path, value=i+1)
                print(new_pdb_path)
                bioutils.change_chain(pdb_in_path=pdb_path, 
                            pdb_out_path=new_pdb_path,
                            rot_tra_matrix=transformation)
                multimer_new_chains.append(new_pdb_path)
            multimer_dict[chain] = multimer_new_chains

        return multimer_dict

    def __repr__(self):
        #Print class

        return f' \
        pdb_path: {self.pdb_path} \n \
        pdb_id: {self.pdb_id} \n \
        template_path: {self.template_path} \n \
        chain: {self.chain} \n \
        polyala_res_dict: {self.polyala_res_dict} \n \
        add_to_msa: {self.add_to_msa} \n \
        add_to_templates: {self.add_to_templates} \n \
        sum_prob: {self.sum_prob} \n \
        aligned: {self.aligned} \n \
        template_features: {self.template_features}'