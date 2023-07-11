import logging
import os
import shutil
from typing import Dict, List, Optional, Union
from libs import bioutils, features, hhsearch, match_restrictions, utils, change_res, alphafold_classes, template_chains
from libs import structures, sequence


class Template:

    def __init__(self, parameters_dict: Dict, output_dir: str, num_of_copies: int, new_name=None):
        self.pdb_path: str
        self.pdb_id: str
        self.cif_path: str
        self.template_path: str
        self.template_sequence: str
        self.generate_multimer: bool
        self.change_res_struct: change_res.ChangeResiduesList = change_res.ChangeResiduesList()
        self.match_restrict_struct: match_restrictions.MatchRestrictionsList = match_restrictions.MatchRestrictionsList()
        self.add_to_msa: bool
        self.add_to_templates: bool
        self.sum_prob: bool
        self.aligned: bool
        self.legacy: bool
        self.strict: bool
        self.template_features: Optional[features.Features] = None
        self.results_path_position: List = [None] * num_of_copies
        self.reference: Optional[str] = None
        self.template_chains_struct: template_chains.TemplateChainsList = template_chains.TemplateChainsList()
        self.alignment_database: List[structures.AlignmentDatabase] = []

        pdb_path = utils.get_input_value(name='pdb', section='template', input_dict=parameters_dict)
        if new_name is not None:
            pdb_out_path = os.path.join(output_dir, f'{new_name}.pdb')
        else:
            pdb_out_path = os.path.join(output_dir, f'{utils.get_file_name(pdb_path)}.pdb')
        self.pdb_path = bioutils.check_pdb(pdb_path, pdb_out_path)
        if utils.check_ranked(os.path.basename(self.pdb_path)):
            raise Exception(f'Template {self.pdb_path} has a protected name (ranked). Please change the name before '
                            f'continuing, as it can cause issues with the ARCIMBOLDO_AIR output.')

        self.pdb_id = utils.get_file_name(self.pdb_path)
        self.add_to_msa = utils.get_input_value(name='add_to_msa', section='template', input_dict=parameters_dict)
        self.add_to_templates = utils.get_input_value(name='add_to_templates', section='template',
                                                      input_dict=parameters_dict)
        self.sum_prob = utils.get_input_value(name='sum_prob', section='template', input_dict=parameters_dict)
        self.legacy = utils.get_input_value(name='legacy', section='template', input_dict=parameters_dict)
        self.strict = utils.get_input_value(name='strict', section='template', input_dict=parameters_dict)
        self.aligned = utils.get_input_value(name='aligned', section='template', input_dict=parameters_dict, override_default=self.legacy)
        self.template_path = f'{os.path.join(output_dir, self.pdb_id)}_template.pdb'
        self.reference = utils.get_input_value(name='reference', section='template', input_dict=parameters_dict)
        self.generate_multimer = utils.get_input_value(name='generate_multimer', section='template',
                                                       input_dict=parameters_dict)
        self.generate_multimer = False if self.legacy else self.generate_multimer

        for paramaters_change_res in utils.get_input_value(name='change_res', section='template',
                                                           input_dict=parameters_dict):
            change_res_dict = {}
            resname = utils.get_input_value(name='resname', section='change_res', input_dict=paramaters_change_res)
            fasta_path = utils.get_input_value(name='fasta_path', section='change_res',
                                               input_dict=paramaters_change_res)
            when = utils.get_input_value(name='when', section='change_res', input_dict=paramaters_change_res)
            for chain, change in paramaters_change_res.items():
                if chain.lower() == 'all' or len(chain) == 1:
                    change_list = utils.expand_residues(change)
                    change_chain_list = bioutils.get_chains(self.pdb_path) if chain.lower() == 'all' else [chain]
                    change_res_dict.update({key: list(set(change_list)) for key in change_chain_list})

            self.change_res_struct.append_change(chain_res_dict=change_res_dict,
                                                 resname=resname,
                                                 fasta_path=fasta_path,
                                                 when=when)

        for parameters_match_dict in utils.get_input_value(name='match', section='template',
                                                           input_dict=parameters_dict):
            chain_match = utils.get_input_value(name='chain', section='match', input_dict=parameters_match_dict)
            residues_match_list = utils.get_input_value(name='residues', section='match',
                                                        input_dict=parameters_match_dict)
            residues = None
            if residues_match_list is not None:
                change_list = utils.expand_residues(residues_match_list)
                new_dict = {chain_match: change_list}
                residues = change_res.ChangeResidues(chain_res_dict=new_dict)
            position = utils.get_input_value(name='position', section='match', input_dict=parameters_match_dict)
            if position != -1:
                position = position - 1
            reference = utils.get_input_value(name='reference', section='match', input_dict=parameters_match_dict)
            reference_chain = utils.get_input_value(name='reference_chain', section='match',
                                                    input_dict=parameters_match_dict)
            self.match_restrict_struct.append_match(chain=chain_match, position=position, residues=residues,
                                                    reference=reference, reference_chain=reference_chain)

        cryst_card = bioutils.extract_cryst_card_pdb(pdb_in_path=self.pdb_path)
        bioutils.remove_hetatm(self.pdb_path, self.pdb_path)
        bioutils.remove_hydrogens(self.pdb_path, self.pdb_path)
        tmp_dir = os.path.join(output_dir, 'tmp')
        utils.create_dir(tmp_dir)
        aux_chain_dict = bioutils.split_pdb_in_chains(output_dir=tmp_dir, pdb_in_path=self.pdb_path)
        self.apply_changes(chain_dict=aux_chain_dict, when='before_alignment')
        aux_list = utils.dict_values_to_list(aux_chain_dict)
        bioutils.merge_pdbs(list_of_paths_of_pdbs_to_merge=aux_list, merged_pdb_path=self.pdb_path)
        if cryst_card is not None:
            bioutils.add_cryst_card_pdb(pdb_in_path=self.pdb_path, cryst_card=cryst_card)
        self.cif_path = bioutils.pdb2mmcif(pdb_in_path=self.pdb_path,
                                           cif_out_path=os.path.join(output_dir, f'{self.pdb_id}.cif'))
        self.template_sequence = bioutils.extract_sequence_from_file(file_path=self.pdb_path)

    def get_templates_references(self) -> List:
        # Get all the references from another template
        # Those references can be in the match class or in the
        # template itself
        return_references_list = self.match_restrict_struct.get_reference_list()
        return_references_list.append(self.reference)
        return list(filter(None, return_references_list))

    def generate_features(self, output_dir: str, global_reference, sequence_assembled: sequence.SequenceAssembled):
        #   - Generate offset.
        #   - Apply the generated offset to all the templates.
        #   - Build the new template merging all the templates.
        #   - Create features for the new template.

        logging.debug(f'Generating features of template {self.pdb_id}')

        if not self.legacy:
            merge_list = []
            self.results_path_position = self.sort_chains_into_positions(
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
            aux_path = os.path.join(output_dir, f'{utils.get_file_name(self.pdb_path)}_split.pdb')
            positions = bioutils.split_chains_assembly(
                pdb_in_path=self.template_path,
                pdb_out_path=aux_path,
                sequence_assembled=sequence_assembled)
            chain_dict = bioutils.chain_splitter(aux_path)
            for i, pos in enumerate(list(positions.keys())):
                if pos in chain_dict:
                    self.results_path_position[i] = chain_dict[pos]
                    self.template_chains_struct.from_dict_to_struct(chain_dict={pos: chain_dict[pos]},
                                                                    alignment_dict={},
                                                                    sequence=sequence_assembled.get_sequence_name(i),
                                                                    change_res_list=self.change_res_struct,
                                                                    match_restrict_list=self.match_restrict_struct)
                else:
                    self.results_path_position[i] = None   

        self.template_features = features.extract_template_features_from_aligned_pdb_and_sequence(
            query_sequence=sequence_assembled.sequence_assembled,
            pdb_path=self.template_path,
            pdb_id=self.pdb_id,
            chain_id='A')

        logging.error(f'Positions of chains in the template {self.pdb_id}: {" | ".join([str(element) for element in self.results_path_position])}')

    def apply_changes(self, chain_dict: Dict, when: str):
        # Apply changes in the pdb, change residues.
        for change_residues in self.change_res_struct.change_residues_list:
            if change_residues.when == when:
                for chain, path in chain_dict.items():
                    if chain in change_residues.chain_res_dict.keys():
                        change_residues.change_residues(pdb_in_path=path, pdb_out_path=path)

    def generate_database(self, output_dir: str, databases: alphafold_classes.AlphaFoldPaths):
        for key, seq in self.template_sequence.items():
            aux_key = key.replace('>', '')
            sequence_name = aux_key.split(':')[0]
            sequence_chain = aux_key.split(':')[1]
            fasta_path = os.path.join(output_dir, f'{self.pdb_id}_{sequence_chain}.fasta')
            bioutils.write_sequence(sequence_name=f'{sequence_name}:{sequence_chain}', sequence_amino=seq,
                                    sequence_path=fasta_path)
            new_database = hhsearch.create_database_from_pdb(fasta_path=fasta_path,
                                                             databases=databases,
                                                             output_dir=output_dir)

            database = structures.AlignmentDatabase(fasta_path=fasta_path, chain=sequence_chain,
                                                    database_path=new_database)
            self.alignment_database.append(database)

    def align(self, output_dir: str, sequence_in: sequence.Sequence, databases: alphafold_classes.AlphaFoldPaths):
        # If aligned, skip to the end, it is just necessary to extract the features
        # If not aligned:
        #   - Convert pdb to cif, this is necessary to run hhsearch.
        #   - Run hhsearch, align the .cif with the given sequence.
        #   - For each chain, extract the features from the alignment result, create
        #       the templates for each chain, each template has the chain 'A'.
        #   - Change the specified residues in the input in all the templates.

        query_sequence = bioutils.extract_sequence(sequence_in.fasta_path)
        extracted_chain_dict = {}
        alignment_chain_dict = {}
        if not self.aligned:
            a3m_path = None
            for database in self.alignment_database:
                # Sometimes with A3M doesn't work, so first we will try with fasta.
                hhr_path = os.path.join(output_dir, f'{utils.get_file_name(database.fasta_path)}.hhr')
                hhsearch.run_hhsearch(a3m_path=sequence_in.fasta_path,
                                      database_path=database.database_path,
                                      output_path=hhr_path)
                template_features, mapping, identities, aligned_columns, total_columns, evalue = \
                    features.extract_template_features_from_pdb(
                        query_sequence=query_sequence,
                        hhr_path=hhr_path,
                        cif_path=self.cif_path,
                        chain_id=database.chain
                    )

                if int(aligned_columns) != int(total_columns):
                    hhr_path2 = os.path.join(output_dir, f'{utils.get_file_name(database.fasta_path)}2.hhr')
                    if not a3m_path:
                        a3m_path = hhsearch.create_a3m(fasta_path=sequence_in.fasta_path,
                                                       databases=databases,
                                                       output_dir=output_dir)
                    hhsearch.run_hhsearch(a3m_path=a3m_path,
                                          database_path=database.database_path,
                                          output_path=hhr_path2)

                    template_features2, mapping2, identities2, aligned_columns2, total_columns2, evalue2 = \
                        features.extract_template_features_from_pdb(
                            query_sequence=query_sequence,
                            hhr_path=hhr_path2,
                            cif_path=self.cif_path,
                            chain_id=database.chain
                        )
                    if int(aligned_columns) < int(aligned_columns2):
                        os.remove(hhr_path)
                        shutil.move(hhr_path2, hhr_path)
                        template_features, mapping, identities, aligned_columns, total_columns, evalue = \
                            template_features2, mapping2, identities2, aligned_columns2, total_columns2, evalue2


                logging.info(f'Alignment results for pdb {self.pdb_id} and chain {database.chain}, with sequence {utils.get_file_name(sequence_in.fasta_path)}:')
                logging.info(
                    f'Aligned columns: {aligned_columns} ({total_columns}), Evalue: {evalue}, Identities: {identities}')
                if template_features is not None:
                    g = features.Features(query_sequence=query_sequence)
                    g.append_new_template_features(new_template_features=template_features,
                                                   custom_sum_prob=self.sum_prob)
                    aux_dict = g.write_all_templates_in_features(output_dir=output_dir, chain=database.chain)
                    extracted_chain_dict[database.chain] = list(aux_dict.values())[0]
                alignment_chain_dict[database.chain] = structures.Alignment(hhr_path=hhr_path, identities=identities,
                                                                            aligned_columns=aligned_columns,
                                                                            total_columns=total_columns, evalue=evalue,
                                                                            database=database, mapping=mapping)
        else:
            extracted_chain_dict = bioutils.split_pdb_in_chains(output_dir=output_dir, pdb_in_path=self.pdb_path)

        if not extracted_chain_dict:
            return extracted_chain_dict

        if self.generate_multimer:
            try:
                extracted_chain_dict = bioutils.generate_multimer_chains(self.pdb_path, extracted_chain_dict)
            except Exception as e:
                logging.error(f'Not possible to generate multimer for {self.pdb_path}')


        self.template_chains_struct.from_dict_to_struct(chain_dict=extracted_chain_dict,
                                                        alignment_dict=alignment_chain_dict,
                                                        sequence=sequence_in.name,
                                                        change_res_list=self.change_res_struct,
                                                        match_restrict_list=self.match_restrict_struct)

    def sort_chains_into_positions(self, sequence_name_list: List[str], global_reference) -> List[str]:
        # Given the sequences list and if there is any global reference:
        # Sort all template chains in the corresponding positions.
        # If the user has set up any match, only the chains specified in the match will be set.
        # Otherwise, there is an algorithm that will sort the chains into the positions,
        # taking into account the pdist between the reference and the chain.
        # If the evalues are high, the program will stop.

        composition_path_list = [None] * len(sequence_name_list)
        deleted_positions = []

        for chain_match in self.template_chains_struct.get_chains_with_matches_ref():
            positions = utils.get_positions_by_chain(chain_match.match.reference.results_path_position,
                                                     chain_match.match.reference_chain)
            for position in positions:
                if composition_path_list[position] is None and \
                        chain_match.sequence == sequence_name_list[chain_match.match.position]:
                    composition_path_list[position] = chain_match.path
                    deleted_positions.append(position)
                    chain_match.check_alignment(stop=self.strict)
                    break

        for chain_match in self.template_chains_struct.get_chains_with_matches_pos():
            position = chain_match.match.position
            if int(position) < len(composition_path_list):
                if chain_match.sequence == sequence_name_list[position]:
                    composition_path_list[position] = chain_match.path
                    deleted_positions.append(position)
                    chain_match.check_alignment(stop=self.strict)

        reference = self.reference if self.reference is not None else None
        reference = global_reference if reference is None else reference

        if not any(composition_path_list):
            new_targets_list = self.template_chains_struct.get_chains_not_in_list(composition_path_list)
            if new_targets_list and len(deleted_positions) < len(sequence_name_list):
                results_targets_list = self.choose_best_offset(reference=reference,
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
                if y not in deleted_positions and name_list[y] == template_chain.sequence:
                    alignment = template_chain.check_alignment(stop=False)
                    if not self.strict or (self.strict and alignment):
                        reference_algorithm.append(
                            (x, y, bioutils.pdist(query_pdb=template_chain.path, target_pdb=target_pdb), alignment, template_chain.path))

            if reference_algorithm:
                results_algorithm.append(reference_algorithm)

        return_offset_list = [None] * (len(reference.results_path_position))
        best_offset_list = bioutils.calculate_auto_offset(results_algorithm,
                                                          len(return_offset_list) - len(deleted_positions))
        for x, y, _, _, path in best_offset_list:
            return_offset_list[y] = path

        return return_offset_list

    def set_reference_templates(self, a_air):
        # Change pdb_id str to the Template reference
        if self.reference is not None:
            self.reference = a_air.get_template_by_id(self.reference)
        for match in self.match_restrict_struct.match_restrict_list:
            if match.reference is not None:
                new_reference = a_air.get_template_by_id(match.reference)
                match.set_reference(new_reference)

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
        # Return all those residues that has been changed.
        chains_changed = []
        fasta_changed = []
        chains_deleted = []
        for path in self.results_path_position:
            if path is not None:
                changed, fasta, deleted = self.template_chains_struct.get_changes(path)
            else:
                changed, fasta, deleted = None, None, None
            chains_changed.append(changed)
            fasta_changed.append(fasta)
            chains_deleted.append(deleted)
        return chains_changed, fasta_changed, chains_deleted

    def get_results_alignment(self) -> List[Union[None, structures.Alignment]]:
        # Return the alignments corresponding to the positions.
        return [self.template_chains_struct.get_alignment_by_path(path) if path is not None else None for path
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
