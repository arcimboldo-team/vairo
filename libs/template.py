import collections
import copy
import logging
import os
import shutil
from typing import Dict, List, Optional, Tuple, Union

from libs import bioutils, features, hhsearch, match_restrictions, utils, change_res
from libs import structures, sequence


class Template:

    def __init__(self, parameters_dict: Dict, output_dir: str, input_dir: str, num_of_copies: int, new_name=None):
        self.pdb_path: str
        self.pdb_id: str
        self.cif_path: str
        self.template_path: str
        self.generate_multimer: bool = True if num_of_copies > 1 else False
        self.change_res_list: List[change_res.ChangeResidues] = []
        self.add_to_msa: bool = False
        self.add_to_templates: bool = True
        self.sum_prob: bool = False
        self.aligned: bool = False
        self.legacy: bool = False
        self.strict: bool = True
        self.template_features: Optional[features.Features] = None
        self.match_restrict_list: List[match_restrictions.MatchRestrictions] = []
        self.results_path_position: List = [None] * num_of_copies
        self.reference: Optional[str] = None
        self.alignment_dict = {}
        self.alignments: List[structures.Alignment] = []
        self.alignment_database: List[structures.AlignmentDatabase] = []
        self.changes: List[structures.ChangesChains] = []

        self.pdb_path = bioutils.check_pdb(utils.get_mandatory_value(parameters_dict, 'pdb'), input_dir)
        if new_name is not None:
            self.pdb_path = shutil.move(self.pdb_path, os.path.join(os.path.dirname(self.pdb_path), f'{new_name}.pdb'))
        self.pdb_id = utils.get_file_name(self.pdb_path)
        bioutils.remove_atoms_types(self.pdb_path, self.pdb_path)

        self.add_to_msa = parameters_dict.get('add_to_msa', self.add_to_msa)
        self.add_to_templates = parameters_dict.get('add_to_templates', self.add_to_templates)
        self.sum_prob = parameters_dict.get('sum_prob', self.sum_prob)
        self.legacy = parameters_dict.get('legacy', self.legacy)
        self.strict = parameters_dict.get('strict', self.strict)
        self.aligned = self.legacy
        self.aligned = parameters_dict.get('aligned', self.aligned)
        self.template_path = f'{output_dir}/{self.pdb_id}_template.pdb'
        self.reference = parameters_dict.get('reference', self.reference)
        self.generate_multimer = parameters_dict.get('generate_multimer', self.generate_multimer)

        for paramaters_change_res in parameters_dict.get('change_res', self.change_res_list):
            change_res_dict = {}
            resname = paramaters_change_res.get('resname', None)
            fasta_path = paramaters_change_res.get('fasta_path', None)
            when = paramaters_change_res.get('when', 'after_alignment')
            try:
                del paramaters_change_res['resname']
            except:
                pass
            try:
                del paramaters_change_res['fasta_path']
            except:
                pass
            try:
                del paramaters_change_res['when']
            except:
                pass
            for values in paramaters_change_res.items():
                if len(values) != 2:
                    raise Exception('Wrong change residues format')
                chain, change = values
                change_list = utils.expand_residues(change)
                change_chain_list = bioutils.get_chains(self.pdb_path) if chain.lower() == 'all' else [chain]
                change_res_dict.update({key: list(set(change_list)) for key in change_chain_list})
            self.change_res_list.append(
                change_res.ChangeResidues(chain_res_dict=change_res_dict, resname=resname, fasta_path=fasta_path,
                                          when=when))

        tmp_dir = os.path.join(output_dir, 'tmp')
        utils.create_dir(tmp_dir)
        cryst_card = bioutils.extract_cryst_card_pdb(pdb_in_path=self.pdb_path)
        aux_chain_dict = bioutils.split_pdb_in_chains(output_dir=tmp_dir, pdb_in_path=self.pdb_path)
        self.apply_changes(chain_dict=aux_chain_dict, when='before_alignment')
        aux_list = utils.dict_values_to_list(aux_chain_dict)
        aux_list = utils.remove_list_layer(input_list=aux_list)
        bioutils.merge_pdbs(list_of_paths_of_pdbs_to_merge=aux_list, merged_pdb_path=self.pdb_path)
        if cryst_card is not None:
            bioutils.add_cryst_card_pdb(pdb_in_path=self.pdb_path, cryst_card=cryst_card)

        self.cif_path = bioutils.pdb2mmcif(output_dir=output_dir, pdb_in_path=self.pdb_path,
                                           cif_out_path=os.path.join(output_dir, f'{self.pdb_id}.cif'))

        for parameters_match_dict in parameters_dict.get('match', self.match_restrict_list):
            self.match_restrict_list.append(match_restrictions.MatchRestrictions(parameters_match_dict))

    def get_reference_list(self) -> List:
        # Get all the references from another template
        # Those references can be in the match class or in the
        # template itself

        return_references_list = [match.reference for match in self.match_restrict_list]
        return_references_list.append(self.reference)
        return list(filter(None, return_references_list))

    def generate_features(self, output_dir: str, global_reference, sequence_assembled: sequence.SequenceAssembled):
        #   - Generate offset.
        #   - Apply the generated offset to all the templates.
        #   - Build the new template merging all the templates.
        #   - Create features for the new template.

        logging.info(f'Generating features of template {self.pdb_id}')

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
                self.results_path_position[i] = chain_dict[pos] if pos in chain_dict else None

        template_features = features.extract_template_features_from_aligned_pdb_and_sequence(
            query_sequence=sequence_assembled.sequence_assembled,
            pdb_path=self.template_path,
            pdb_id=self.pdb_id,
            chain_id='A')

        self.template_features = copy.deepcopy(template_features)

        logging.info(f'Positions of chains in the template {self.pdb_id}: {self.results_path_position}')

    def apply_changes(self, chain_dict: Dict, when: str):
        # Apply changes in the pdb, change residues.
        for change_residues in self.change_res_list:
            if change_residues.when == when:
                for chain, paths_list in chain_dict.items():
                    if chain in change_residues.chain_res_dict.keys():
                        for path in paths_list:
                            change_residues.change_residues(pdb_in_path=path, pdb_out_path=path)

    def generate_database(self, output_dir: str, database_path: str):
        template_sequence = bioutils.extract_sequence_from_file(self.cif_path)
        for extracted_sequence in template_sequence:
            sequence_chain = extracted_sequence[extracted_sequence.find(':') + 1:extracted_sequence.find('\n')]
            sequence_name = extracted_sequence.splitlines()[0].replace('>', '')
            spit_sequence = extracted_sequence.splitlines()[1]
            fasta_path = os.path.join(output_dir, f'{self.pdb_id}_{sequence_chain}.fasta')
            bioutils.write_sequence(sequence_name=sequence_name, sequence_amino=spit_sequence, sequence_path=fasta_path)
            new_database = hhsearch.create_database_from_pdb(fasta_path=fasta_path,
                                                             database_path=database_path,
                                                             output_dir=output_dir)

            database = structures.AlignmentDatabase(fasta_path=fasta_path, chain=sequence_chain,
                                                    database_path=new_database)
            self.alignment_database.append(database)

    def align(self, output_dir: str, sequence_in: sequence.Sequence):
        # If aligned, skip to the end, it is just necessary to extract the features
        # If not aligned:
        #   - Convert pdb to cif, this is necessary to run hhsearch.
        #   - Create a fake hhsearch database with just our .cif in it. 
        #   - Run hhsearch, align the .cif with the given sequence.
        #   - For each chain, extract the features from the alignment result, create
        #       the templates for each chain, each template has the chain 'A'.
        #   - Change the specified residues in the input in all the templates.

        query_sequence = bioutils.extract_sequence(sequence_in.fasta_path)

        extracted_chain_dict = {}

        if not self.aligned:
            for database in self.alignment_database:
                hhr_path = os.path.join(output_dir, f'{utils.get_file_name(database.fasta_path)}.hhr')

                hhsearch.run_hhsearch(fasta_path=sequence_in.fasta_path,
                                      database_path=database.database_path,
                                      output_path=hhr_path)
                template_features, mapping, identities, aligned_columns, total_columns, evalue = \
                    features.extract_template_features_from_pdb(
                        query_sequence=query_sequence,
                        hhr_path=hhr_path,
                        cif_path=self.cif_path,
                        chain_id=database.chain
                    )

                if template_features is not None:
                    self.mapping_has_changed(chain=database.chain, mapping=mapping)
                    g = features.Features(query_sequence=query_sequence)
                    g.append_new_template_features(new_template_features=template_features,
                                                   custom_sum_prob=self.sum_prob)
                    aux_dict = g.write_all_templates_in_features(output_dir=output_dir, chain=database.chain)
                    extracted_chain_path = list(aux_dict.values())[0]
                    extracted_chain_dict[database.chain] = [extracted_chain_path]

                self.alignments.append(
                    structures.Alignment(hhr_path=hhr_path, identities=identities, aligned_columns=aligned_columns,
                                         total_columns=total_columns, evalue=evalue, database=database,
                                         extracted_path=extracted_chain_path))
        else:
            extracted_chain_dict = bioutils.split_pdb_in_chains(output_dir=output_dir, pdb_in_path=self.pdb_path)

        if not extracted_chain_dict:
            return extracted_chain_dict

        if self.generate_multimer:
            try:
                extracted_chain_dict = bioutils.generate_multimer_chains(self.pdb_path, extracted_chain_dict)
            except Exception as e:
                logging.info(f'Not possible to generate multimer for {self.pdb_path}')


        self.apply_changes(chain_dict=extracted_chain_dict, when='after_alignment')
        self.alignment_dict[sequence_in.name] = extracted_chain_dict

    def mapping_has_changed(self, chain: str, mapping: Dict):
        # It is necessary to update the mapping that generates the alignment
        # as the residues numbering has changed.

        structure = bioutils.get_structure(self.pdb_path)
        residues_list = list(structure[0][chain].get_residues())
        idres_list = list([bioutils.get_resseq(res) for res in residues_list])
        mapping_keys = list(map(lambda x: x + 1, list(mapping.keys())))
        mapping_values = list(map(lambda x: x + 1, list(mapping.values())))
        mapping = dict(zip(mapping_keys, mapping_values))
        if idres_list != mapping_keys and len(idres_list) == len(mapping_keys):
            for match in self.match_restrict_list:
                if match.residues is not None:
                    match.residues.apply_mapping(chain, mapping)
            for res in self.change_res_list:
                res.apply_mapping(chain, mapping)

    def sort_chains_into_positions(self, sequence_name_list: List[str], global_reference) \
            -> List[Tuple[str, None]]:
        # Given the sequences list and if there is any global reference:
        # Sort all template chains in the corresponding positions.
        # If the user has set up any match, only the chains specified in the match will be set.
        # Otherwise, there is an algorithm that will sort the chains into the positions,
        # taking into account the pdist between the reference and the chain.
        # If the evalues are high, the program will stop.

        composition_path_list = [None] * len(sequence_name_list)
        new_target_code_list = []
        deleted_positions = []

        new_dict = collections.defaultdict(list)
        for _, chain_dict in self.alignment_dict.items():
            for chain, paths in chain_dict.items():
                codes = [utils.get_chain_and_number(path) for path in paths]
                [new_dict[chain].append(f'{code[0]}{code[1]}') for code in codes]
        chain_dict = {chain: sorted(list(set(values))) for chain, values in new_dict.items()}

        for i, match in enumerate(self.match_restrict_list):
            if match.chain is None or match.chain not in new_dict:
                logging.info('Restriction could not be applied')
                continue

            code_pdb = chain_dict[match.chain].pop(0)
            chain_dict[match.chain].append(code_pdb)

            if match.reference is not None and match.reference_chain is not None:
                positions = utils.get_positions_by_chain(match.reference.results_path_position,
                                                         match.reference_chain)
                for position in positions:
                    if composition_path_list[position] is None:
                        composition_path_list[position] = utils.select_path_from_code(align_dict=self.alignment_dict,
                                                                                      code=code_pdb,
                                                                                      position=match.position,
                                                                                      sequence_name_list=sequence_name_list)
                        deleted_positions.append(position)
                        self.check_alignment(pdb_path=composition_path_list[position], stop=self.strict)
                        break
                continue

            if match.position != -1:
                if int(match.position) < len(composition_path_list):
                    path = utils.select_path_from_code(align_dict=self.alignment_dict,
                                                       code=code_pdb,
                                                       position=match.position,
                                                       sequence_name_list=sequence_name_list)
                    if match.residues is not None:
                        new_path = os.path.join(os.path.dirname(path), f'{i}_{utils.get_file_name(path)}.pdb')
                        match.residues.delete_residues_inverse(path, new_path)
                        path = new_path

                    composition_path_list[match.position] = path
                    self.check_alignment(pdb_path=path, stop=self.strict)

                    deleted_positions.append(match.position)
                    continue
                logging.info(
                    f'Position exceed the length of the sequence, selecting a random position for chain {match.chain}')
                new_target_code_list.append(code_pdb)
            else:
                new_target_code_list.append(code_pdb)
                continue

        if not self.match_restrict_list:
            for chain, paths in chain_dict.items():
                [new_target_code_list.append(paths[i]) for i in range(len(paths))]

        reference = self.reference if self.reference is not None else None
        reference = global_reference if reference is None else reference

        if new_target_code_list:
            new_target_path_list = self.choose_best_offset(reference=reference,
                                                           deleted_positions=deleted_positions,
                                                           code_list=new_target_code_list,
                                                           name_list=sequence_name_list)

            for i, path in enumerate(new_target_path_list):
                if composition_path_list[i] is None:
                    composition_path_list[i] = path

            if len(new_target_code_list) != sum(x is not None for x in new_target_path_list):
                logging.info(
                    f'Not all chains have been selected in the template {self.pdb_id}. Probably there are chains with bad alignment.')

            if not any(new_target_path_list):
                raise Exception(
                    f'Not possible to meet the requisites for the template {self.pdb_id}. No chains have good alignments')

        return composition_path_list

    def choose_best_offset(self, reference, deleted_positions: List[int], code_list: List[str],
                           name_list: List[str]) -> List[Optional[str]]:
        results_algorithm = []

        for x, code_query_pdb in enumerate(code_list):
            reference_algorithm = []
            for y, target_pdb in enumerate(reference.results_path_position):
                if y not in deleted_positions:
                    query_pdb = utils.select_path_from_code(align_dict=self.alignment_dict,
                                                            code=code_query_pdb,
                                                            position=y,
                                                            sequence_name_list=name_list)
                    alignment = self.check_alignment(query_pdb) if self.strict else True
                    reference_algorithm.append(
                        (x, y, bioutils.pdist(query_pdb=query_pdb, target_pdb=target_pdb), alignment))
            results_algorithm.append(reference_algorithm)

        return_offset_list = [None] * (len(reference.results_path_position))
        best_offset_list = bioutils.calculate_auto_offset(results_algorithm,
                                                          len(return_offset_list) - len(deleted_positions))
        for x, y, _, _ in best_offset_list:
            return_offset_list[y] = utils.select_path_from_code(align_dict=self.alignment_dict,
                                                                code=code_list[x],
                                                                position=y,
                                                                sequence_name_list=name_list)

        return return_offset_list

    def set_reference_templates(self, a_air):
        # Change pdb_id str to the Template reference

        if self.reference is not None:
            self.reference = a_air.get_template_by_id(self.reference)
        for match in self.match_restrict_list:
            if match.reference is not None:
                new_reference = a_air.get_template_by_id(match.reference)
                match.set_reference(new_reference)

    def get_changes(self):
        # Get the changes that have been done to the templates.
        # Return all those residues that has been changed.
        chains_changed = [None] * len(self.results_path_position)
        fasta_changed = [None] * len(self.results_path_position)
        chains_deleted = [None] * len(self.results_path_position)
        for i, path in enumerate(self.results_path_position):
            if path is not None:
                chain, _ = utils.get_chain_and_number(path)
                changed_residues = []
                changed_fasta = []
                for change in self.change_res_list:
                    if chain in change.chain_res_dict:
                        if change.resname is not None:
                            changed_residues.extend(change.chain_res_dict[chain])
                        elif change.fasta_path is not None:
                            changed_fasta.extend(change.chain_res_dict[chain])
                if changed_residues:
                    chains_changed[i] = changed_residues
                if changed_fasta:
                    fasta_changed[i] = changed_fasta
                for change in self.match_restrict_list:
                    if change.residues is not None and chain in change.residues.chain_res_dict:
                        chains_deleted.extend(change.residues.chain_res_dict[chain])

        return chains_changed, fasta_changed, chains_deleted

    def get_alignment_by_pdb(self, pdb_path: str) -> structures.Alignment:
        # Search for the alignment that has the same name as the pdb_path
        for alignment in self.alignments:
            if os.path.dirname(pdb_path) == os.path.dirname(alignment.extracted_path):
                chain1, _ = utils.get_chain_and_number(alignment.extracted_path)
                chain2, _ = utils.get_chain_and_number(pdb_path)
                if chain1 == chain2:
                    return alignment
        return None

    def get_results_alignment(self) -> List[Union[None, structures.Alignment]]:
        # Return the alignments corresponding to the positions.

        return_alignments = []
        for path in self.results_path_position:
            alignment = None
            if path is not None:
                alignment = self.get_alignment_by_pdb(path)
            return_alignments.append(alignment)
        return return_alignments

    def check_alignment(self, pdb_path: str, stop: bool = False) -> bool:
        # Check if the alignment is OK.
        # Give the pdb alignment path and check the evalues.
        alignment = self.get_alignment_by_pdb(pdb_path)
        if alignment:
            if float(alignment.evalue) > 0.01:
                if not stop:
                    return False
                else:
                    raise Exception(
                        f'Match could not be done. Poor alignment in the template {self.pdb_id}. Stopping the run.')
            return True
        return True

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
