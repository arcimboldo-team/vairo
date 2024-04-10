import copy
import logging
import os
import shutil
from typing import Union, List, Dict

from libs import utils, bioutils, structures, template_modifications


class TemplateChain:
    def __init__(self,
                 chain: str,
                 path: str,
                 path_before_changes: str,
                 code: str,
                 sequence: str,
                 modification: template_modifications.ChainModifications = None,
                 alignment: Union[None, structures.Alignment] = None,
                 deleted_residues: List[int] = [],
                 changed_residues: List[int] = [],
                 fasta_residues: List[int] = [],
                 sequence_before_changes: str = ''
                 ):
        self.chain = chain
        self.path = path
        self.path_before_changes = path_before_changes
        self.code = code
        self.sequence = sequence
        self.modification = modification
        self.alignment = alignment
        self.deleted_residues = deleted_residues
        self.changed_residues = changed_residues
        self.fasta_residues = fasta_residues
        self.sequence_before_changes = sequence_before_changes

    def get_chain_code(self) -> List:
        return self.chain, self.code

    def check_alignment(self, stop: bool):
        # Check if it has been a good alignment. If stop, then throw an error.
        if self.alignment:
            if float(self.alignment.evalue) > 0.01:
                if not stop:
                    return False
                else:
                    raise Exception(
                        f'Match could not be done. Poor alignment in the template {utils.get_file_name(self.path)}. '
                        f'Stopping the run.')
            return True
        return True

    def __repr__(self):
        # Print class
        return f'pdb_path: {self.path}'


class TemplateChainsList:
    def __init__(self):
        self.template_chains_list: List[TemplateChain] = []

    def get_template_chain(self, pdb_path: str) -> TemplateChain:
        chain1, code1 = utils.get_chain_and_number(pdb_path)
        pdb_dirname = os.path.dirname(pdb_path)
        for template_chain in self.template_chains_list:
            # if there is an alignment, we check the directory too, because it can be from different alignments
            if pdb_dirname == os.path.dirname(
                    template_chain.path) and chain1 == template_chain.chain and code1 == template_chain.code:
                return template_chain
        return None

    def get_old_sequence(self, pdb_path: str) -> str:
        template_chain = self.get_template_chain(pdb_path)
        if template_chain is not None:
            return template_chain.sequence_before_changes
        else:
            return None

    def get_changes(self, pdb_path: str) -> List:
        template_chain = self.get_template_chain(pdb_path)
        if template_chain is not None:
            return template_chain.changed_residues, template_chain.fasta_residues, template_chain.deleted_residues, template_chain.sequence_before_changes
        else:
            return None

    def get_number_chains(self) -> int:
        return len({(chain_template.chain, chain_template.code) for chain_template in self.template_chains_list})

    def get_alignment_by_path(self, pdb_path: str):
        # Search for the alignment that has the same name as the pdb_path
        template_chain = self.get_template_chain(pdb_path)
        if template_chain is not None and template_chain.alignment:
            return template_chain.alignment
        return None

    def get_chains_not_in_list(self, input_list: List[str]) -> List[TemplateChain]:
        # It will create a List of TemplateChains where the chain and the code is not in the input list.
        # This can be used to eliminate duplicates in case there is more than one database.
        input_chains = [(utils.get_chain_and_number(path)[0], utils.get_chain_and_number(path)[1]) for path in
                        input_list if path is not None]
        return [template_chain for template_chain in self.template_chains_list if
                (template_chain.chain, template_chain.code) not in input_chains]

    def get_chains_with_matches_ref(self) -> List[TemplateChain]:
        # Get all the chains that have a match with references
        return [template_chain for template_chain in self.template_chains_list if
                template_chain.match is not None and template_chain.match.check_references()]

    def get_chains_with_matches_pos(self) -> List[TemplateChain]:
        # Get all the chains that have a match with a determinate position.
        return [template_chain for template_chain in self.template_chains_list if
                template_chain.modification is not None and template_chain.modification.check_position()]

    def new_chain_sequence(self, path: str, sequence: str):
        chain, number = utils.get_chain_and_number(path)


        self.template_chains_list.append(TemplateChain(chain=chain,
                                                       path=path,
                                                       path_before_changes=path_before_changes,
                                                       code=number,
                                                       sequence=sequence,
                                                       modification=match,
                                                       alignment=alignment,
                                                       deleted_residues=deleted_residues,
                                                       changed_residues=changed_residues,
                                                       fasta_residues=fasta_residues,
                                                       sequence_before_changes=sequence_before_changes))

        chain_aux_path = tempfile.NamedTemporaryFile(delete=False)
        temp_aux = template_modifications.TemplateModifications()
        temp_aux.append_chain_modifications(modification_pos_list)
        temp_aux.modify_template(pdb_in_path=chain_path, pdb_out_path=chain_aux_path.name,
                                 type_modify=['mutate', 'delete'], when='before_alignment')
        if not self.aligned:
            alignment_chain_dir = os.path.join(sequence_in.alignment_dir, chain)
            utils.create_dir(alignment_chain_dir)
            extracted_chain, alignment_chain = hhsearch.run_hh(output_dir=alignment_chain_dir,
                                                               database_dir=template_database_dir,
                                                               chain_in_path=chain_aux_path.name,
                                                               query_sequence_path=sequence_in.fasta_path,
                                                               databases=databases,
                                                               name=self.pdb_id)
            alignment_dict = {chain: alignment_chain}
        else:
            alignment_dict = {}
            extracted_chain = chain_aux_path
        if extracted_chain:
            self.template_chains_struct.from_dict_to_struct(chain=chain, path=extracted_chain,
                                                            alignment_dict=alignment_dict,
                                                            sequence=sequence_in.name,
                                                            modifications_list=self.modifications_struct,
                                                            modify_chain=modification,
                                                            generate_multimer=self.generate_multimer,
                                                            pdb_path=self.pdb_path)



    def from_dict_to_struct(self, chain: str, path: str, alignment_dict: Dict, sequence: str,
                            modifications_list: template_modifications.TemplateModifications,
                            generate_multimer: bool = False, pdb_path: str = '',
                            modify_chain: template_modifications.ChainModifications = None):
        # Given a dict, with all the information of a Chain (there can be more than one chain in case of multimer)
        # Read all the information, and create as many TemplateChains as paths and append them to the list.

        if generate_multimer:
            try:
                chain_dict = bioutils.generate_multimer_chains(pdb_path, {chain: path})
                path_list = chain_dict[chain]
            except Exception as e:
                logging.error(f'Not possible to generate multimer for {pdb_path}')
        if modify_chain is None:
            modify_chain_copy = copy.deepcopy(modifications_list)
        else:
            modify_chain_copy = copy.deepcopy(modify_chain)

        position = -1
        if isinstance(modify_chain_copy, template_modifications.TemplateModifications):
            match_list = modify_chain_copy.get_modifications_by_chain(chain=chain)
        else:
            position = modify_chain_copy.position
            match_list = [modify_chain_copy]
            filtered_list = list(filter(lambda x: x not in [y.path for y in self.template_chains_list], path_list))
            if filtered_list:
                path_list = [filtered_list[0]]
            else:
                path_list = [path_list[0]]

        alignment = alignment_dict.get(chain)
        for path in path_list:
            change_res_copy = copy.deepcopy(modifications_list)
            change_list = change_res_copy.get_modifications_by_chain_and_position(chain=chain, position=position)
            match = None
            # If there is an alignment, there might be a mapping. It is necessary to apply that mapping
            # to change_res (match can have a change_res too) because the residue numbering has changed during
            # the alignment
            if alignment:
                structure = bioutils.get_structure(path)
                residues_list = list(structure[0][chain].get_residues())
                idres_list = list([bioutils.get_resseq(res) for res in residues_list])
                mapping_keys = list(map(lambda x: x + 1, list(alignment.mapping.keys())))
                mapping_values = list(map(lambda x: x + 1, list(alignment.mapping.values())))
                mapping = dict(zip(mapping_keys, mapping_values))
                if idres_list != mapping_keys and len(idres_list) == len(mapping_keys):
                    for match_l in match_list:
                        match_l.apply_mapping(mapping)
                    for change in change_list:
                        change.apply_mapping(mapping)

            # Store the sequence before changing the residues, as if we want to add it in the MSA, it would not
            # make sense
            sequence_before_changes = list(bioutils.extract_sequence_msa_from_pdb(path).values())[0]
            path_before_changes = os.path.join(os.path.dirname(path), f'{utils.get_file_name(path)}_originalseq.pdb')
            if os.path.exists(path_before_changes):
                os.remove(path_before_changes)
            shutil.copy2(path, path_before_changes)

            # We just use one match per chain. As that chain just can be in one position. All the delete residues
            # have to be in a one-line sentence also.

            for match in match_list:
                change = template_modifications.TemplateModifications()
                change.append_chain_modification(match)
                change.modify_template(path, path, type_modify=['delete', 'mutate'], when='after_alignment')
                match = match_list.pop(0)
                break

            change = template_modifications.TemplateModifications()
            change.append_chain_modifications(change_list)
            change.modify_template(pdb_in_path=path, pdb_out_path=path, type_modify=['delete', 'mutate'],
                                   when='after_alignment')

            deleted_residues = []
            if match:
                deleted_residues = match.get_deleted_residues()
            changed_residues, fasta_residues = change_res_copy.get_residues_changed_by_chain(chain)
            chain, number = utils.get_chain_and_number(path)
            if not self.get_template_chain(path):
                self.template_chains_list.append(TemplateChain(chain=chain,
                                                               path=path,
                                                               path_before_changes=path_before_changes,
                                                               code=number,
                                                               sequence=sequence,
                                                               modification=match,
                                                               alignment=alignment,
                                                               deleted_residues=deleted_residues,
                                                               changed_residues=changed_residues,
                                                               fasta_residues=fasta_residues,
                                                               sequence_before_changes=sequence_before_changes))