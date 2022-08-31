import copy
from operator import length_hint
import os
import logging
import shutil
import string
from Bio.PDB import PDBParser, Selection

from libs import arcimboldo_air, bioutils, features, hhsearch, utils, change_res
from typing import Dict, List
import pandas as pd


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

    def choose_best_offset(self, target, merged_list: List) -> Dict:

        self.merged_list = merged_list
        if self is target:
            return

        results_algorithm = []
        for x, query_pdb in enumerate(merged_list):
            reference_algorithm = []
            for y, target_pdb in enumerate(target.merged_list):
                reference_algorithm.append((x, y, bioutils.pdist(query_pdb=query_pdb, target_pdb=target_pdb)))   
            results_algorithm.append(reference_algorithm)

        best_offset_list = bioutils.calculate_best_offset(results_algorithm)
        return_offset_dict = {}
        for x,y,_ in best_offset_list:
            return_offset_dict[merged_list[x]] = y
        
        return return_offset_dict
    
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
                template_dict[chain] = template_path

                for change_residues in self.change_res:
                    change_residues.change_residues(template_path, template_path, chain)

            template_list = sorted(list(template_dict.values()))
            merge_list = []
            offset_dict = self.choose_best_offset(a_air.templates[0], template_list)
            
            for i, pdb_path in enumerate(template_list):
                if bool(offset_dict):
                    offset_value = offset_dict[pdb_path]
                else:
                    offset_value = i
                offset = query_seq_length * (offset_value) + a_air.glycines * (offset_value)
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
        
    def create_multimer(self, a_air):
        #Read remark to get the transformations and the new chains
        #Split chains and apply transformations to generate the new ones
        #Rename chains with the new_chains_list
        #Store a dict with the relation between old chains and new chains
        # Dict -> A: [A, B]
        #Change template path to the new multimer path

        chain_list, transformations_list = bioutils.read_remark_350(self.pdb_path)
        new_chain_list = list(string.ascii_uppercase)[:len(transformations_list) * len(chain_list)]
        len_transformations_list = len(transformations_list)
        merge_pdbs_list = []
        multimer_dict = {}
        
        logging.info('Assembly can be build using chain(s) '+ str(chain_list) + ' by applying the following transformations:')
        for matrix in transformations_list:
            logging.info(str(matrix))

        new_pdb_path = f'{a_air.run_dir}/{self.pdb_id}.pdb'

        for i, chain in enumerate(chain_list):
            pdb_path_name = f'{a_air.run_dir}/{self.pdb_id}_{chain}'
            pdb_path = f'{a_air.run_dir}/{self.pdb_id}_{chain}.pdb'
            bioutils.chain_splitter(pdb_in_path=self.pdb_path, pdb_out_path=pdb_path, chain=chain)
            multimer_new_chains = []
            for j, transformation in enumerate(transformations_list):
                bioutils.change_chain(pdb_in_path=pdb_path, 
                            pdb_out_path=f'{pdb_path_name}_{j+1}.pdb',
                            rot_tra_matrix=transformation,
                            chain=new_chain_list[j+len_transformations_list*i])
                merge_pdbs_list.append(f'{pdb_path_name}_{j+1}.pdb')
                multimer_new_chains.append(new_chain_list[j+len_transformations_list*i])
            multimer_dict[chain] = multimer_new_chains
        
        bioutils.merge_pdbs(list_of_paths_of_pdbs_to_merge=sorted(merge_pdbs_list),
                    merged_pdb_path=new_pdb_path) 
        self.pdb_path = new_pdb_path
        self.update_chains(multimer_dict)

    def update_chains(self, new_chains_dict: Dict):
        # Update the chains after changing them generating the monomer
        # new_chains_dict: A: [A, B], B: [C,D]

        for change_residues in self.change_res:
            change_residues.update_new_chains(new_chains_dict)

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