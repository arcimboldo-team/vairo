import logging
import os
import shutil
from typing import Dict, List, Optional, Union
from libs import bioutils, features, hhsearch, utils, alphafold_classes, template_chains
from libs import structures, sequence, template_modifications


class Template:

    def __init__(self, parameters_dict: Dict, output_dir: str, num_of_copies: int, new_name=None):
        self.pdb_path: str
        self.pdb_id: str
        self.template_path: str
        self.template_chains_dir: str
        self.template_originalseq_path: str
        self.generate_multimer: bool
        self.modifications_struct: template_modifications.TemplateModifications = template_modifications.TemplateModifications()
        self.add_to_msa: bool
        self.add_to_templates: bool
        self.sum_prob: bool
        self.aligned: bool
        self.legacy: bool
        self.strict: bool
        self.template_features: Optional[dict] = None
        self.results_path_position: List = [None] * num_of_copies
        self.template_chains_struct: template_chains.TemplateChainsList = template_chains.TemplateChainsList()
        self.selected_positions: bool = False


        pdb_path = utils.get_input_value(name='pdb', section='template', input_dict=parameters_dict)
        if new_name is not None:
            pdb_out_path = os.path.join(output_dir, f'{new_name}.pdb')
        else:
            pdb_out_path = os.path.join(output_dir, f'{utils.get_file_name(pdb_path)}.pdb')
        bioutils.check_pdb(pdb_path, pdb_out_path)
        if utils.check_ranked(os.path.basename(pdb_out_path)):
            raise Exception(f'Template {pdb_out_path} has a protected name (ranked). Please change the name before '
                            f'continuing, as it can cause issues with the VAIRO output.')

        self.pdb_id = utils.get_file_name(pdb_out_path)
        self.template_chains_dir = os.path.join(output_dir, f'{self.pdb_id}_template_data')
        self.pdb_path = os.path.join(self.template_chains_dir, f'{self.pdb_id}.pdb')
        utils.create_dir(self.template_chains_dir, delete_if_exists=True)
        shutil.move(pdb_out_path, self.pdb_path)
        
        self.template_originalseq_path = f'{os.path.join(self.template_chains_dir, self.pdb_id)}_template_originalseq.pdb'
        self.template_path = f'{os.path.join(self.template_chains_dir, self.pdb_id)}_template.pdb'


        self.add_to_msa = utils.get_input_value(name='add_to_msa', section='template', input_dict=parameters_dict)
        self.add_to_templates = utils.get_input_value(name='add_to_templates', section='template',
                                                      input_dict=parameters_dict)
        self.sum_prob = utils.get_input_value(name='sum_prob', section='template', input_dict=parameters_dict)
        self.legacy = utils.get_input_value(name='legacy', section='template', input_dict=parameters_dict)
        self.strict = utils.get_input_value(name='strict', section='template', input_dict=parameters_dict)
        self.aligned = utils.get_input_value(name='aligned', section='template', input_dict=parameters_dict,
                                             override_default=self.legacy)
        self.generate_multimer = utils.get_input_value(name='generate_multimer', section='template',
                                                       input_dict=parameters_dict)

        for parameter_modification in utils.get_input_value(name='modifications', section='template',
                                                            input_dict=parameters_dict):

            position = utils.get_input_value(name='position', section='modifications',
                                             input_dict=parameter_modification)

            if position != -1:
                position = position - 1
                self.selected_positions = True

            maintain_residues = utils.expand_residues(
                utils.get_input_value(name='maintain_residues', section='modifications',
                                      input_dict=parameter_modification))
            delete_residues = utils.expand_residues(
                utils.get_input_value(name='delete_residues', section='modifications',
                                      input_dict=parameter_modification))

            when = utils.get_input_value(name='when', section='modifications', input_dict=parameter_modification)

            replace_list = []
            parameter_replace = utils.get_input_value(name='mutations', section='modifications', input_dict=parameter_modification)
            if parameter_replace:
                for replace in parameter_replace:
                    residues = utils.expand_residues(utils.get_input_value(name='numbering_residues', section='mutations', input_dict=replace))
                    mutate_with = utils.get_input_value(name='mutate_with', section='mutations', input_dict=replace)
                    replace_list.append(template_modifications.ResidueMutate(mutate_residues_number=residues, mutate_with=mutate_with))

            chains = utils.get_input_value(name='chain', section='modifications', input_dict=parameter_modification)
            if chains.lower() == 'all' or self.legacy:
                chains = ['all']
            else:
                chains = chains.replace(" ", "").split(',')

            self.modifications_struct.append_modification(chains=chains, maintain_residues=maintain_residues,
                                                          delete_residues=delete_residues, position=position,
                                                          mutations=replace_list, when=when)

        cryst_card = bioutils.extract_cryst_card_pdb(pdb_in_path=self.pdb_path)
        bioutils.remove_hetatm(self.pdb_path, self.pdb_path)
        bioutils.remove_hydrogens(self.pdb_path, self.pdb_path)

        if cryst_card is not None:
            bioutils.add_cryst_card_pdb(pdb_in_path=self.pdb_path, cryst_card=cryst_card)

    def generate_features(self, global_reference, sequence_assembled: sequence.SequenceAssembled):
        #   - Generate offset.
        #   - Apply the generated offset to all the templates.
        #   - Build the new template merging all the templates.
        #   - Create features for the new template.

        logging.info(f'Generating features of template {self.pdb_id}')
        self.template_chains_struct.apply_changes()

        merge_list = []
        self.results_path_position = self.sort_chains_into_positions(
            sequence_name_list=sequence_assembled.get_list_name(),
            global_reference=global_reference)
        for i, pdb_path in enumerate(self.results_path_position):
            if pdb_path is not None:
                offset = sequence_assembled.get_starting_length(i)
                new_pdb_path = os.path.join(self.template_chains_dir, f'{self.pdb_id}_{offset}.pdb')
                bioutils.change_chain(pdb_in_path=pdb_path,
                                      pdb_out_path=new_pdb_path,
                                      offset=offset, chain='A')
                merge_list.append(new_pdb_path)
        bioutils.merge_pdbs(list_of_paths_of_pdbs_to_merge=utils.sort_by_digit(merge_list),
                            merged_pdb_path=self.template_path)


        aux_path_list = []
        chain_name = 'A'
        for path in self.results_path_position:
            if path is not None:
                tchain = self.template_chains_struct.get_template_chain(path)
                tchain = tchain.path_before_changes
                bioutils.change_chain(pdb_in_path=tchain,
                                      pdb_out_path=tchain,
                                      chain=chain_name)
                aux_path_list.append(tchain)
            chain_name = chr(ord(chain_name) + 1)
        bioutils.merge_pdbs(list_of_paths_of_pdbs_to_merge=aux_path_list,
                            merged_pdb_path=self.template_originalseq_path)
        self.template_features = features.extract_template_features_from_aligned_pdb_and_sequence(
            query_sequence=sequence_assembled.sequence_assembled,
            pdb_path=self.template_path,
            pdb_id=self.pdb_id,
            chain_id='A')

        logging.error(
            f'Positions of chains in the template {self.pdb_id}: {" | ".join([str(element) for element in self.results_path_position])}')

    def generate_chains(self, sequence_assembled: sequence.SequenceAssembled):
        if not self.legacy:
            template_chains_aux = bioutils.split_pdb_in_chains(output_dir=self.template_chains_dir, pdb_path=self.pdb_path)
            for chain, template_chain_path in template_chains_aux.items():
                if self.generate_multimer:
                    try:
                        chain_dict = bioutils.generate_multimer_chains(self.pdb_path, {chain: template_chain_path})
                        path_list = chain_dict[chain]
                    except:
                        path_list = [template_chain_path]
                else:
                    path_list = [template_chain_path]
                if not self.selected_positions:
                    modifications_list = template_modifications.TemplateModifications(self.modifications_struct.get_modifications_by_chain(chain=chain))
                    for sequence_in in sequence_assembled.sequence_list:
                        for new_path in path_list:
                            self.template_chains_struct.new_chain_sequence(path=new_path, sequence=sequence_in,
                                                                           modifications_list=modifications_list)
                else:
                    modifications_list = self.modifications_struct.get_modifications_position_by_chain(chain=chain)
                    for i, (modification, chain_path) in enumerate(zip(modifications_list[:len(path_list)], path_list)):
                        modification_pos_list = template_modifications.TemplateModifications(
                            self.modifications_struct.get_modifications_by_chain_and_position(
                                chain=chain, position=modification.position))
                        self.template_chains_struct.new_chain_sequence(path=chain_path,
                                                                       sequence=sequence_assembled.sequence_list_expanded[
                                                                           modification.position],
                                                                       modifications_list=modification_pos_list)

        else:
            aux_path = os.path.join(self.template_chains_dir, f'{utils.get_file_name(self.pdb_path)}_split.pdb')
            positions = bioutils.split_chains_assembly(
                pdb_in_path=self.pdb_path,
                pdb_out_path=aux_path,
                sequence_assembled=sequence_assembled)
            chain_dict = bioutils.split_pdb_in_chains(pdb_path=aux_path)
            for i, pos in enumerate(list(positions.keys())):
                if pos in chain_dict:
                    modification_pos_list = template_modifications.TemplateModifications(
                        self.modifications_struct.get_modifications_by_chain_and_position(chain=pos, position=i))
                    modification_pos_list.append_modification(chains=[pos], position=i)
                    self.template_chains_struct.new_chain_sequence(path=chain_dict[pos],
                                                                   sequence=sequence_assembled.get_sequence_name(i),
                                                                   modifications_list=modification_pos_list)

    def align(self, databases: alphafold_classes.AlphaFoldPaths):
        template_database_dir = os.path.join(self.template_chains_dir, 'databases')
        utils.create_dir(template_database_dir, delete_if_exists=True)
        for chain_struct in self.template_chains_struct.template_chains_list:
            alignment_chain_dir = os.path.join(chain_struct.sequence.alignment_dir, f'{chain_struct.chain}{chain_struct.code}')
            utils.create_dir(alignment_chain_dir, delete_if_exists=True)
            hhr_path = self.template_chains_struct.set_same_alignment(chain_struct)
            extracted_chain, alignment_chain = hhsearch.run_hh(output_dir=alignment_chain_dir,
                                                               database_dir=template_database_dir,
                                                               chain_in_path=chain_struct.path,
                                                               query_sequence_path=chain_struct.sequence.fasta_path,
                                                               databases=databases,
                                                               temp_name=self.pdb_id,
                                                               hhr_path=hhr_path)
            if extracted_chain:
                chain_struct.set_alignment(alignment_chain)
                chain_struct.set_extracted_chain(extracted_chain)

    def sort_chains_into_positions(self, sequence_name_list: List[str], global_reference) -> List[str]:
        # Given a sequence list and if there is any global reference:
        # Sort all template chains in the corresponding positions.
        # If the user has set up any match, only the chains specified in the match will be set.
        # Otherwise, there is an algorithm that will sort the chains into the positions,
        # taking into account the pdist between the reference and the chain.
        # If the evalues are high, the program will stop.

        composition_path_list = [None] * len(sequence_name_list)
        deleted_positions = []

        for chain_match in self.template_chains_struct.get_chains_with_matches_pos():
            position = chain_match.check_position()
            if int(position) < len(composition_path_list):
                composition_path_list[position] = chain_match.path
                deleted_positions.append(position)
                chain_match.check_alignment(stop=self.strict)

        if not any(composition_path_list):
            new_targets_list = self.template_chains_struct.get_chains_not_in_list(composition_path_list)
            if new_targets_list and len(deleted_positions) < len(sequence_name_list):
                results_targets_list = self.choose_best_offset(reference=global_reference,
                                                               deleted_positions=deleted_positions,
                                                               template_chains_list=new_targets_list,
                                                               name_list=sequence_name_list)
                for i, element in enumerate(results_targets_list):
                    if composition_path_list[i] is None:
                        composition_path_list[i] = element

        if not any(composition_path_list):
            raise Exception(
                f'Not possible to meet the requisites for the template {self.pdb_id}. No chains have good alignments')

        if self.template_chains_struct.get_number_chains() != sum(x is not None for x in composition_path_list) \
                and not all(composition_path_list):
            logging.error(f'Not all chains have been selected in the template {self.pdb_id}')

        return composition_path_list

    def choose_best_offset(self, reference, deleted_positions: List[int],
                           template_chains_list: List[template_chains.TemplateChain],
                           name_list: List[str]) -> List[Optional[str]]:

        results_algorithm = []
        code_value = 0
        codes_dict = {}
        for template_chain in template_chains_list:
            chain_code = template_chain.get_chain_code()
            chain_code = f'{chain_code[0]}{chain_code[1]}'
            if chain_code not in codes_dict:
                codes_dict[chain_code] = code_value
                code_value += 1
            x = codes_dict[chain_code]
            reference_algorithm = []
            for y, target_pdb in enumerate(reference.results_path_position):
                if y not in deleted_positions and name_list[y] == template_chain.sequence.name:
                    alignment = template_chain.check_alignment(stop=False)
                    if not self.strict or (self.strict and alignment):
                        reference_algorithm.append(
                            (x, y, bioutils.pdist(query_pdb=template_chain.path, target_pdb=target_pdb), alignment,
                             template_chain.path))

            if reference_algorithm:
                results_algorithm.append(reference_algorithm)

        return_offset_list = [None] * (len(reference.results_path_position))
        best_offset_list = bioutils.calculate_auto_offset(results_algorithm,
                                                          len(return_offset_list) - len(deleted_positions))
        for x, y, _, _, path in best_offset_list:
            return_offset_list[y] = path

        return return_offset_list

    def get_old_sequence(self, sequence_list: List[sequence.Sequence], glycines: int) -> str:
        old_sequence = []
        for i, path in enumerate(self.results_path_position):
            if path is not None:
                seq = self.template_chains_struct.get_old_sequence(path)
            else:
                seq = ''
            while len(seq) < sequence_list[i].length:
                seq += '-'
            if i != len(sequence_list) - 1:
                seq += '-' * glycines
            old_sequence.append(seq)
        return ''.join(old_sequence)

    def get_changes(self) -> List:
        # Get the changes that have been done to the templates.
        # Return all those residues that have been changed.
        chains_changed = []
        fasta_changed = []
        chains_deleted = []
        old_sequence_changed = []
        for path in self.results_path_position:
            if path is not None:
                changed, fasta, deleted, old_sequence = self.template_chains_struct.get_changes(path)
            else:
                changed, fasta, deleted, old_sequence = None, None, None, None
            chains_changed.append(changed)
            fasta_changed.append(fasta)
            chains_deleted.append(deleted)
            old_sequence_changed.append(old_sequence)
        return chains_changed, fasta_changed, chains_deleted, old_sequence_changed

    def get_results_alignment(self) -> List[Union[None, structures.Alignment]]:
        # Return the alignments corresponding to the positions.
        return [self.template_chains_struct.get_alignment_by_path(path) if path is not None else None for path
                in self.results_path_position]

    def get_chain_by_position(self) -> str:
        # Return the chain used in the result path position. If there is any, if not, return None
        return [utils.get_chain_and_number(path)[0] if path is not None else None for path
                in self.results_path_position]

    def __repr__(self):
        # Print class
        return f' \
        pdb_path: {self.pdb_path} \n \
        pdb_id: {self.pdb_id} \n \
        template_path: {self.template_path} \n \
        add_to_msa: {self.add_to_msa} \n \
        add_to_templates: {self.add_to_templates} \n \
        sum_prob: {self.sum_prob} \n \
        aligned: {self.aligned} \n'
