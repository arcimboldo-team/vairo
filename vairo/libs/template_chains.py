import copy
import logging
import os
import shutil
from typing import Union, List, Dict

from libs import match_restrictions, utils, bioutils, structures


class TemplateChain:
    def __init__(self,
                 chain: str,
                 path: str,
                 path_before_changes: str,
                 code: str,
                 sequence: str,
                 match: match_restrictions.MatchRestrictions = None,
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
        self.match = match
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
            #if there is alignment, we check the directory too, because it can be from different alignments
            if pdb_dirname == os.path.dirname(template_chain.path) and chain1 == template_chain.chain and code1 == template_chain.code:
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
        # Get all the chains that has a match with references
        return [template_chain for template_chain in self.template_chains_list if
                template_chain.match is not None and template_chain.match.check_references()]

    def get_chains_with_matches_pos(self) -> List[TemplateChain]:
        # Get all the chains that has a match with a determinate position.
        return [template_chain for template_chain in self.template_chains_list if
                template_chain.match is not None and template_chain.match.check_position()]

    def from_dict_to_struct(self, chain_dict: Dict, alignment_dict: Dict, sequence: str, change_res_list,
                            match_restrict_list: match_restrictions.MatchRestrictionsList, generate_multimer: bool = False, pdb_path: str = ''):
        # Given a dict, with all the information of a Chain (there can be more than one chain in case of multimer)
        # Read all the information, and create as many TemplateChains as paths and append them to the list.
        
        if not chain_dict:
            return chain_dict

        if generate_multimer:
            try:
                chain_dict = bioutils.generate_multimer_chains(pdb_path, chain_dict)
            except Exception as e:
                logging.error(f'Not possible to generate multimer for {pdb_path}')
        for chain, paths in chain_dict.items():
            path_list = []
            if isinstance(paths, list):
                for path in paths:
                    path_list.append(path)
            else:
                path_list.append(paths)

            match_restrict_copy = copy.deepcopy(match_restrict_list)
            if isinstance(match_restrict_copy, match_restrictions.MatchRestrictionsList):
                match_list = match_restrict_copy.get_matches_by_chain(chain=chain)
            else:
                match_list = [match_restrict_copy]
                filtered_list = list(filter(lambda x: x not in [y.path for y in self.template_chains_list], path_list))
                if filtered_list:
                    path_list = [filtered_list[0]]
                else:
                    path_list = [path_list[0]]

            alignment = alignment_dict.get(chain)

            for path in path_list:
                change_res_copy = copy.deepcopy(change_res_list)
                change_list = change_res_copy.get_changes_by_chain(chain=chain, when='after_alignment')
                match = None
                # If there is an alignment, there might be a mapping. It is necessary to apply that mapping
                # to change_res (match can have a change_res too) because the residues numbering has changed during
                # the alignment
                if alignment:
                    structure = bioutils.get_structure(path)
                    residues_list = list(structure[0][chain].get_residues())
                    idres_list = list([bioutils.get_resseq(res) for res in residues_list])
                    mapping_keys = list(map(lambda x: x + 1, list(alignment.mapping.keys())))
                    mapping_values = list(map(lambda x: x + 1, list(alignment.mapping.values())))
                    mapping = dict(zip(mapping_keys, mapping_values))
                    if idres_list != mapping_keys and len(idres_list) == len(mapping_keys):
                        for match in match_list:
                            if match.residues is not None:
                                match.residues.apply_mapping(chain, mapping)
                        for res in change_list:
                            res.apply_mapping(chain, mapping)

                # We just use one match per chain. As that chain just can be in one position. All the delete residues
                # has to be in a one line sentence also.
                for match in match_list:
                    if match.residues:
                        match.residues.delete_residues_inverse(path, path)
                    match = match_list.pop(0)
                    break

                # Store the sequence before changing the residues, as if we want to add it in the MSA, it would not
                # make sense
                sequence_before_changes = list(bioutils.extract_sequence_msa_from_pdb(path).values())[0]
                path_before_changes = os.path.join(os.path.dirname(path), f'{utils.get_file_name(path)}_originalseq.pdb')
                if os.path.exists(path_before_changes):
                    os.remove(path_before_changes)
                shutil.copy2(path, path_before_changes)
                for change in change_list:
                    change.change_residues(path, path)

                deleted_residues = []
                if match:
                    deleted_residues = match.get_deleted_residues(chain=chain)
                changed_residues, fasta_residues = change_res_copy.get_residues_changed_by_chain(chain)
                chain, number = utils.get_chain_and_number(path)

                if not self.get_template_chain(path):
                    self.template_chains_list.append(TemplateChain(chain=chain,
                                                                   path=path,
                                                                   path_before_changes=path_before_changes,
                                                                   code=number,
                                                                   sequence=sequence,
                                                                   match=match,
                                                                   alignment=alignment,
                                                                   deleted_residues=deleted_residues,
                                                                   changed_residues=changed_residues,
                                                                   fasta_residues=fasta_residues,
                                                                   sequence_before_changes=sequence_before_changes))
