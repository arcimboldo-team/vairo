import copy
import os
import logging
from pickle import NONE
import re
import shutil
from Bio.PDB import PDBParser, Selection
import pandas as pd

from libs import bioutils, features, hhsearch, match_restrictions, utils, change_res
from typing import Dict, List, Tuple


class Template:

    def __init__ (self, parameters_dict: Dict, output_dir: str, input_dir: str, num_of_copies: int):

        self.pdb_path: str
        self.pdb_id: str
        self.template_path: str
        self.chain: str = ''
        self.generate_multimer: bool = True if num_of_copies > 1 else False
        self.change_res_list: List[change_res.ChangeResidues] = []
        self.add_to_msa: bool = False
        self.add_to_templates: bool = False
        self.sum_prob: bool = False
        self.aligned: bool = False
        self.template_features_dict: Dict = None
        self.match_restrict_list: List[match_restrictions.MatchRestrictions] = []
        self.results_path_position: List = []
        self.reference: str = None
        self.automatch: bool = True
        
        self.pdb_path = self.check_pdb(utils.get_mandatory_value(parameters_dict, 'pdb'), input_dir)
        self.pdb_id = utils.get_file_name(self.pdb_path)
        self.add_to_msa = parameters_dict.get('add_to_msa', self.add_to_msa)
        self.add_to_templates = parameters_dict.get('add_to_templates', self.add_to_templates)
        self.sum_prob = parameters_dict.get('sum_prob', self.sum_prob)
        self.aligned = parameters_dict.get('aligned', self.aligned)
        self.template_path = f'{output_dir}/{self.pdb_id}_template.pdb'
        self.chain = parameters_dict.get('chain', self.chain)
        self.reference = parameters_dict.get('reference', self.reference)
        self.generate_multimer = parameters_dict.get('generate_multimer', self.generate_multimer)
        self.automatch = parameters_dict.get('automatch', self.automatch)

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

        for paramaters_change_res in parameters_dict.get('change_res', self.change_res_list):
            change_res_dict = {}
            
            resname = utils.get_mandatory_value(paramaters_change_res, 'resname')
            del paramaters_change_res['resname']

            for values in paramaters_change_res.items():
                if len(values) != 2:
                    raise Exception('Wrong change residues format')
                chain, change = values 
                change_list = utils.expand_residues(change) 
                change_chain_list = chains if chain.lower() == 'all' else [chain]
                change_res_dict.update({key: list(set(change_list)) for key in change_chain_list})                       
            self.change_res_list.append(change_res.ChangeResidues(change_dict=change_res_dict, resname=resname))

        for parameters_match_dict in parameters_dict.get('match', self.match_restrict_list):
            self.match_restrict_list.append(match_restrictions.MatchRestrictions(parameters_match_dict))

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

    def get_reference_list(self) -> List:
        #Get all the references from another template
        #Those references can be in the match class or in the
        #template itself

        return_references_list = [match.reference for match in self.match_restrict_list]
        return_references_list.append(self.reference)
        return list(filter(None, return_references_list))
    
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
            extracted_chain_dict = {}
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
                aux_dict = g.write_all_templates_in_features(output_dir=a_air.run_dir, chain=chain)
                extracted_chain_path = list(aux_dict.values())[0]
                extracted_chain_dict[chain] = [extracted_chain_path]
                
            if self.generate_multimer:
                extracted_chain_dict = self.generate_multimer_chains(extracted_chain_dict)

            extracted_chain_list = [val for sublist in extracted_chain_dict.values() for val in sublist]

            for change_residues in self.change_res_list:
                for chain, paths_list in extracted_chain_dict.items():
                    if chain in change_residues.change_dict.keys():
                        for path in paths_list:
                            change_residues.change_residues(pdb_in_path=path, pdb_out_path=path)
                                
            merge_list = []
            self.sort_chains_into_positions(extracted_chain_dict, a_air)
            a_air.append_line_in_templates(self.results_path_position)

            for i, pdb_path in enumerate(self.results_path_position):
                if pdb_path is not None:
                    offset = query_seq_length * (i) + a_air.glycines * (i)
                    new_pdb_path = f'{a_air.run_dir}/{offset}.pdb'
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
            
        self.template_features_dict = copy.deepcopy(template_features)

    def sort_chains_into_positions(self, chain_dict: Dict, a_air):
        
        composition_path_list = [None] * a_air.num_of_copies

        new_target_path_list = []
        deleted_positions = []

        for i, match in enumerate(self.match_restrict_list):
            new_pdb = chain_dict[match.chain][0]
            path_pdb = chain_dict[match.chain][0]
            if match.chain is None or match.chain not in chain_dict:
                logging.Logger.info('Restriction could not be applied')
                continue
            if match.residues is not None:
                name = utils.get_file_name(path_pdb)
                new_pdb = os.path.join(os.path.dirname(path_pdb), f'{i}_{name}.pdb')
                match.residues.change_residues(path_pdb, new_pdb)
            if match.position is not None and match.position < len(composition_path_list):
                composition_path_list[match.position] = new_pdb
                deleted_positions.append(match.position)
                continue
            if match.reference is not None and match.reference_chain is not None:
                positions = utils.get_positions_by_chain(match.reference.results_path_position, match.reference_chain)
                for position in positions:
                    if composition_path_list[position] is None: 
                        composition_path_list[position] = new_pdb
                        deleted_positions.append(position)
                        break
                continue
            if new_pdb not in composition_path_list:
                new_target_path_list.append(new_pdb)
        
        for chain, paths in chain_dict.items():
            number_of_paths = len(paths)
            number_of_chains = len(utils.get_paths_by_chain(new_target_path_list+composition_path_list, chain))
            for i in range(number_of_chains, number_of_paths):
                new_target_path_list.append(paths[i])
        
        reference = self.reference if self.reference is not None else None
        reference = a_air.reference if reference is None else reference

        if reference != self:
            new_target_path_list = self.choose_best_offset(reference, deleted_positions, new_target_path_list)

        for i, path in enumerate(new_target_path_list):
            if composition_path_list[i] is None:
                composition_path_list[i] = path

        self.results_path_position = composition_path_list

    def set_reference_templates(self, a_air):
        #Change pdb_id str to the Template reference
    
        if self.reference is not None:
            self.reference = a_air.get_template_by_id(self.reference)
        for match in self.match_restrict_list:
            if match.reference is not None:
                new_reference = a_air.get_template_by_id(match.reference)
                match.set_reference(new_reference)

    def choose_best_offset(self, reference, deleted_positions: List, pdb_list: List) -> Dict:
        
        results_pdist = []
        results_pdist.append(['reference'] + [utils.get_file_name(file) for i, file in enumerate(reference.results_path_position) if not i in deleted_positions ])

        results_algorithm = []

        for x, query_pdb in enumerate(pdb_list):
            reference_pdist_list = []
            reference_algorithm = []
            for y, target_pdb in enumerate(reference.results_path_position):
                if not y in deleted_positions:
                    reference_algorithm.append((x, y, bioutils.pdist(query_pdb=query_pdb, target_pdb=target_pdb)))   
                    reference_pdist_list.append(bioutils.pdist(query_pdb=query_pdb, target_pdb=target_pdb))
            results_algorithm.append(reference_algorithm)
            results_pdist.append([utils.get_file_name(query_pdb)] + reference_pdist_list)

        best_offset_list = bioutils.calculate_auto_offset(results_algorithm)

        return_offset_list = [None] * (len(best_offset_list)+len(deleted_positions))

        for x,y,_ in best_offset_list:
            return_offset_list[y] = pdb_list[x]

        output_txt = f'/Users/pep/work/test/arcimboldo_air/2/output/best_offset.txt'

        with open(output_txt, 'w') as f_in:

            f_in.write('Calculated with distances:\n')

            print(results_pdist[1:])
            print(results_pdist[0])
            rows = []
            for x in results_pdist[1:]:
                rows.append(x)
            df = pd.DataFrame(rows, columns=results_pdist[0])
            f_in.write(df.to_markdown())        
        
            f_in.write('\n\n')
            f_in.write('Best combination is (using pdist):')
            f_in.write(f'{str(return_offset_list)}')


        return return_offset_list

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
                bioutils.change_chain(pdb_in_path=pdb_path, 
                            pdb_out_path=new_pdb_path,
                            rot_tra_matrix=transformation)
                multimer_new_chains.append(new_pdb_path)
            multimer_dict[chain] = multimer_new_chains

        self.chain = chain_list
        return multimer_dict

    def __repr__(self):
        #Print class

        return f' \
        pdb_path: {self.pdb_path} \n \
        pdb_id: {self.pdb_id} \n \
        template_path: {self.template_path} \n \
        chain: {self.chain} \n \
        add_to_msa: {self.add_to_msa} \n \
        add_to_templates: {self.add_to_templates} \n \
        sum_prob: {self.sum_prob} \n \
        aligned: {self.aligned} \n'