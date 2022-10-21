import copy
import os
import logging
import shutil
import collections
import pandas as pd
from libs import bioutils, features, hhsearch, match_restrictions, utils, change_res
from typing import Dict, List
from libs.sequence import SequenceAssembled

class Template:

    def __init__ (self, parameters_dict: Dict, output_dir: str, input_dir: str, num_of_copies: int):

        self.pdb_path: str
        self.pdb_id: str
        self.template_path: str
        self.chains: List = []
        self.generate_multimer: bool = True if num_of_copies > 1 else False
        self.change_res_list: List[change_res.ChangeResidues] = []
        self.add_to_msa: bool = False
        self.add_to_templates: bool = True
        self.sum_prob: bool = False
        self.aligned: bool = False
        self.legacy: bool = False
        self.template_features: features.Features = None
        self.match_restrict_list: List[match_restrictions.MatchRestrictions] = []
        self.results_path_position: List = [None] * num_of_copies
        self.hhr_path: str = None
        self.reference: str = None
        
        self.pdb_path = bioutils.check_pdb(utils.get_mandatory_value(parameters_dict, 'pdb'), input_dir)
        self.pdb_id = utils.get_file_name(self.pdb_path)
        self.add_to_msa = parameters_dict.get('add_to_msa', self.add_to_msa)
        self.add_to_templates = parameters_dict.get('add_to_templates', self.add_to_templates)
        self.sum_prob = parameters_dict.get('sum_prob', self.sum_prob)
        self.legacy = parameters_dict.get('legacy', self.legacy)
        if self.legacy:
            self.aligned = True
        self.aligned = parameters_dict.get('aligned', self.aligned)
        self.template_path = f'{output_dir}/{self.pdb_id}_template.pdb'
        self.reference = parameters_dict.get('reference', self.reference)
        self.generate_multimer = parameters_dict.get('generate_multimer', self.generate_multimer)
        structure = bioutils.get_structure(pdb_path=self.pdb_path)
        self.chains = [chain.get_id() for chain in structure.get_chains()]
        self.hhr_path = f'{output_dir}/{self.pdb_id}_output.hhr'

        for paramaters_change_res in parameters_dict.get('change_res', self.change_res_list):
            change_res_dict = {}
            resname = utils.get_mandatory_value(paramaters_change_res, 'resname')
            del paramaters_change_res['resname']
            for values in paramaters_change_res.items():
                if len(values) != 2:
                    raise Exception('Wrong change residues format')
                chain, change = values 
                change_list = utils.expand_residues(change) 
                change_chain_list = self.chains if chain.lower() == 'all' else [chain]
                change_res_dict.update({key: list(set(change_list)) for key in change_chain_list})                       
            self.change_res_list.append(change_res.ChangeResidues(chain_res_dict=change_res_dict, resname=resname))

        for parameters_match_dict in parameters_dict.get('match', self.match_restrict_list):
            self.match_restrict_list.append(match_restrictions.MatchRestrictions(parameters_match_dict))

    def get_reference_list(self) -> List:
        #Get all the references from another template
        #Those references can be in the match class or in the
        #template itself
        
        return_references_list = [match.reference for match in self.match_restrict_list]
        return_references_list.append(self.reference)
        return list(filter(None, return_references_list))
    
    def generate_features(self, output_dir: str, alignment_dict: Dict, global_reference, sequence_assembled: SequenceAssembled):
        #   - Generate offset.
        #   - Apply the generated offset to all the templates.
        #   - Build the new template merging all the templates.
        #   - Create features for the new template.

        logging.info(f'Generating features of template {self.pdb_id}')

        if not self.legacy:
            merge_list = []
            self.results_path_position = self.sort_chains_into_positions(
                                                alignment_dict=alignment_dict, 
                                                sequence_name_list=sequence_assembled.get_list_name(),
                                                global_reference=global_reference)
            for i, pdb_path in enumerate(self.results_path_position):
                if pdb_path is not None:
                    offset = sequence_assembled.get_starting_length(i)
                    new_pdb_path = os.path.join(output_dir, f'{self.pdb_id}_{offset}.pdb')
                    bioutils.change_chain(pdb_in_path=pdb_path,
                                    pdb_out_path=new_pdb_path,
                                    offset=offset, chain='A')
                    merge_list.append(new_pdb_path)

            bioutils.merge_pdbs(list_of_paths_of_pdbs_to_merge=utils.sort_by_digit(merge_list),
                    merged_pdb_path=self.template_path)
        else:
            shutil.copy2(self.pdb_path, self.template_path)
            aux_path = os.path.join(output_dir, f'{utils.get_file_name(self.pdb_path)}_splitted.pdb')
            positions = bioutils.split_chains_assembly(
                        pdb_in_path=self.template_path, 
                        pdb_out_path=aux_path, 
                        sequence_assembled=sequence_assembled)
            chain_dict = bioutils.chain_splitter(aux_path)
            for i, pos in enumerate(list(positions.keys())):
                self.results_path_position[i] = chain_dict[pos] if pos in chain_dict else None

        template_features = features.extract_template_features_from_aligned_pdb_and_sequence(
            query_sequence=sequence_assembled.sequence_assembled,
            pdb_path=self.template_path,
            pdb_id=self.pdb_id,
            chain_id='A')
            
        self.template_features = copy.deepcopy(template_features)

        logging.info(f'Positions of chains in the template {self.pdb_id}: {self.results_path_position}')
        return self.results_path_position

    def apply_changes(self, chain_dict: Dict):
        ##Apply changes in the pdb, change residues.

        for change_residues in self.change_res_list:
            for chain, paths_list in chain_dict.items():
                if chain in change_residues.chain_res_dict.keys():
                    for path in paths_list:
                        change_residues.change_residues(pdb_in_path=path, pdb_out_path=path)

    def align(self, output_dir, fasta_path) -> Dict:
        #If aligned, skip to the end, it is just necessary to extract the features
        #If not aligned:
        #   - Convert pdb to cif, this is necessary to run hhsearch.
        #   - Create a fake hhsearch database with just our .cif in it. 
        #   - Run hhsearch, align the .cif with the given sequence.
        #   - For each chain, extract the features from the alignment result, create the
        #       the templates for each chain, each template has the chain 'A'.
        #   - Change the specified residues in the input in all the templates.

        pdb70_path = os.path.join(output_dir, 'pdb70')
        query_sequence = bioutils.extract_sequence(fasta_path)
        if not self.aligned:
            extracted_chain_dict = {}
            bioutils.pdb2mmcif(output_dir=output_dir, pdb_in_path=self.pdb_path, cif_out_path=os.path.join(output_dir,f'{self.pdb_id}.cif'))
            hhsearch.generate_hhsearch_db(template_cif_path=os.path.join(output_dir,f'{self.pdb_id}.cif'), output_dir=output_dir)
            hhsearch.run_hhsearch(fasta_path=fasta_path, pdb70_db=pdb70_path,output_path=self.hhr_path)

            for chain in self.chains:
                template_features, mapping = features.extract_template_features_from_pdb(
                    query_sequence=query_sequence,
                    hhr_path=self.hhr_path,
                    pdb_id=self.pdb_id,
                    chain_id=chain,
                    mmcif_db=output_dir)

                if template_features is not None:
                    self.mapping_has_changed(chain=chain, mapping=mapping)
                    g = features.Features(query_sequence=query_sequence)
                    g.append_new_template_features(new_template_features=template_features, custom_sum_prob=self.sum_prob)
                    aux_dict = g.write_all_templates_in_features(output_dir=output_dir, chain=chain)
                    extracted_chain_path = list(aux_dict.values())[0]
                    extracted_chain_dict[chain] = [extracted_chain_path]
        else:
            aux_path = os.path.join(output_dir, os.path.basename(self.pdb_path))
            shutil.copy2(self.pdb_path, aux_path)
            chain_dict = bioutils.chain_splitter(aux_path)
            extracted_chain_dict = {k: [v] for k, v in chain_dict.items()}
        
        if self.generate_multimer:
            extracted_chain_dict = bioutils.generate_multimer_chains(self.pdb_path, extracted_chain_dict)

        self.apply_changes(chain_dict=extracted_chain_dict)

        return extracted_chain_dict

    def mapping_has_changed(self, chain: str, mapping: Dict):
        #It is necessary to update the mapping that geneartes the alignment
        #as the residues numbering has changed.

        structure = bioutils.get_structure(self.pdb_path)        
        residues_list = list(structure[0][chain].get_residues())
        idres_list = list([bioutils.get_resseq(res) for res in residues_list])
        mapping_keys = list(map(lambda x: x+1, list(mapping.keys())))
        mapping_values = list(map(lambda x: x+1, list(mapping.values())))
        mapping = dict(zip(mapping_keys, mapping_values))
        if idres_list != mapping_keys and len(idres_list) == len(mapping_keys):
            for match in self.match_restrict_list:
                if match.residues is not None:
                    match.residues.apply_mapping(chain, mapping)
            for res in self.change_res_list:
                res.apply_mapping(chain, mapping)

    def sort_chains_into_positions(self, alignment_dict: Dict, sequence_name_list: List[str], global_reference) -> List[str]:
        
        composition_path_list = [None] * len(sequence_name_list)
        new_target_code_list = []
        deleted_positions = []

        new_dict = collections.defaultdict(list)
        for _, chain_dict in alignment_dict.items():
            for chain, paths in chain_dict.items():
                codes = [utils.get_chain_and_number(path) for path in paths]
                [new_dict[chain].append(f'{code[0]}{code[1]}') for code in codes]
        chain_dict = {chain: sorted(list(set(values))) for chain, values in new_dict.items()}

        for i, match in enumerate(self.match_restrict_list):
            if match.chain is None or match.chain not in new_dict:
                logging.info('Restriction could not be applied')
                continue
            try:
                code_pdb = chain_dict[match.chain].pop(0)
            except:
                logging.info('Not enough chains in the pdb.')

            if match.residues is not None:
                paths = utils.get_paths_in_alignment(align_dict=alignment_dict, code=code_pdb)
                for path in paths:
                    match.residues.delete_residues_inverse(path, path)

            if match.position != '' and match.position != 'None':
                if (int(match.position) < len(composition_path_list)):
                    composition_path_list[match.position] = utils.select_path_from_code(align_dict=alignment_dict,
                                                                                        code=code_pdb,
                                                                                        position=match.position,
                                                                                        sequence_name_list=sequence_name_list)
                    deleted_positions.append(match.position)
                    continue
                logging.info(f'Position exceed the length of the sequence, selecting a random position for chain {match.chain}')
            elif match.position == 'None':
                continue
            elif match.reference is not None and match.reference_chain is not None:
                positions = utils.get_positions_by_chain(match.reference.results_path_position, match.reference_chain)
                for position in positions:
                    if composition_path_list[position] is None: 
                        composition_path_list[position] = utils.select_path_from_code(align_dict=alignment_dict,
                                                                                    code=code_pdb,
                                                                                    position=match.position,
                                                                                    sequence_name_list=sequence_name_list)
                        deleted_positions.append(position)
                        break
                continue

            new_target_code_list.append(code_pdb)

        for chain, paths in chain_dict.items():
            [new_target_code_list.append(paths[i]) for i in range(len(paths))]
                
        reference = self.reference if self.reference is not None else None
        reference = global_reference if reference is None else reference

        if new_target_code_list:
            if reference != self:
                new_target_path_list = self.choose_best_offset(reference=reference, 
                                                            deleted_positions=deleted_positions,
                                                            align_dict=alignment_dict, 
                                                            code_list=new_target_code_list, 
                                                            name_list=sequence_name_list)
                for i, path in enumerate(new_target_path_list):
                    if composition_path_list[i] is None:
                        composition_path_list[i] = path
            else:
                for code in new_target_code_list:
                    for i in range(len(composition_path_list)):
                        if composition_path_list[i] is None:
                            composition_path_list[i] = utils.select_path_from_code(align_dict=alignment_dict,
                                                                                    code=code,
                                                                                    position=i,
                                                                                    sequence_name_list=sequence_name_list)
                            break

        return composition_path_list

    def set_reference_templates(self, a_air):
        #Change pdb_id str to the Template reference
    
        if self.reference is not None:
            self.reference = a_air.get_template_by_id(self.reference)
        for match in self.match_restrict_list:
            if match.reference is not None:
                new_reference = a_air.get_template_by_id(match.reference)
                match.set_reference(new_reference)

    def choose_best_offset(self, reference, deleted_positions: List[int], align_dict: Dict, code_list: List[str], name_list: List[str]) -> Dict:
        
        results_algorithm = []

        for x, code_query_pdb in enumerate(code_list):
            reference_pdist_list = []
            reference_algorithm = []
            for y, target_pdb in enumerate(reference.results_path_position):
                if not y in deleted_positions:
                    query_pdb = utils.select_path_from_code(align_dict=align_dict,
                                                            code=code_query_pdb,
                                                            position=y,
                                                            sequence_name_list=name_list)
                    if query_pdb is not None:
                        reference_algorithm.append((x, y, bioutils.pdist(query_pdb=query_pdb, target_pdb=target_pdb)))   
                        reference_pdist_list.append(bioutils.pdist(query_pdb=query_pdb, target_pdb=target_pdb))
            results_algorithm.append(reference_algorithm)

        return_offset_list = [None] * (len(reference.results_path_position))
        best_offset_list = bioutils.calculate_auto_offset(results_algorithm, len(return_offset_list)-len(deleted_positions))
        for x,y,_ in best_offset_list:
            return_offset_list[y] = utils.select_path_from_code(align_dict=align_dict,
                                                                code=code_list[x],
                                                                position=y,
                                                                sequence_name_list=name_list)

        return return_offset_list

    def __repr__(self):
        #Print class

        return f' \
        pdb_path: {self.pdb_path} \n \
        pdb_id: {self.pdb_id} \n \
        template_path: {self.template_path} \n \
        chain: {self.chains} \n \
        add_to_msa: {self.add_to_msa} \n \
        add_to_templates: {self.add_to_templates} \n \
        sum_prob: {self.sum_prob} \n \
        aligned: {self.aligned} \n'