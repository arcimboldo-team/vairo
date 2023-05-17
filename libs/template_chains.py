import copy
import os
from typing import Union, List, Dict

from libs import match_restrictions, utils, bioutils
from libs.structures import Alignment


class TemplateChain:
    def __init__(self,
                 chain: str,
                 path: str,
                 code: str,
                 sequence_before_changes: str,
                 sequence: str,
                 match: match_restrictions.MatchRestrictions = None,
                 alignment: Union[None, Alignment] = None,
                 deleted_residues: List[int] = [],
                 changed_residues: List[int] = []
                 ):
        self.chain = chain
        self.path = path
        self.code = code
        self.sequence = sequence
        self.match = match
        self.alignment = alignment
        self.deleted_residues = deleted_residues
        self.changed_residues = changed_residues
        self.sequence_before_changes = sequence_before_changes

    def check_alignment(self, stop: bool):
        if self.alignment:
            if float(self.alignment.evalue) > 0.01:
                if not stop:
                    return False
                else:
                    raise Exception(
                        f'Match could not be done. Poor alignment in the template {self.pdb_id}. Stopping the run.')
            return True
        return True


class TemplateChainsList:
    def __init__(self):
        self.template_chains_list: List[TemplateChain] = []

    def get_alignment_by_path(self, pdb_path: str) -> Alignment:
        # Search for the alignment that has the same name as the pdb_path
        for template_chain in self.template_chains_list:
            if template_chain.alignment:
                for alignment in self.alignments:
                    if os.path.dirname(pdb_path) == os.path.dirname(alignment.extracted_path):
                        chain1, _ = utils.get_chain_and_number(alignment.extracted_path)
                        chain2, _ = utils.get_chain_and_number(pdb_path)
                        if chain1 == chain2:
                            return alignment
        return None

    def get_chains_not_in_list(self, input_list: List[str]) -> List[TemplateChain]:
        return_list = [
            template_chain
            for template_chain in self.template_chains_list
            if all(
                utils.get_chain_and_number(path)[0] == template_chain.chain
                and utils.get_chain_and_number(path)[1] == template_chain.code
                for path in input_list
                if path is not None
            )
        ]
        return return_list

    def get_chains_with_matches_ref(self) -> List[TemplateChain]:
        return [template_chain for template_chain in self.template_chains_list if template_chain.match is not None and template_chain.match.check_references()]

    def get_chains_with_matches_pos(self) -> List[TemplateChain]:
        return [template_chain for template_chain in self.template_chains_list if template_chain.match is not None and template_chain.match.check_position()]

    def from_dict_to_struct(self, chain_dict: Dict, alignment_dict: Dict, sequence: str, change_res_list,
                            match_restrict_list: match_restrictions.MatchRestrictionsList):
        for chain, paths in chain_dict.items():
            path_list = []
            if isinstance(paths, list):
                for path in paths:
                    path_list.append(path)
            else:
                path_list.append(paths)
            match_restrict_copy = copy.deepcopy(match_restrict_list)
            match_list = match_restrict_copy.get_matches_by_chain(chain=chain)

            alignment = alignment_dict.get(chain)

            for path in path_list:
                sequence_before_changes = bioutils.extract_sequence_from_file(path)
                change_res_copy = copy.deepcopy(change_res_list)
                change_list = change_res_copy.get_changes_by_chain(chain=chain, when='after_alignment')
                match = None
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

                for change in change_list:
                    change.change_residues(path, path)
                for match in match_list:
                    if match.residues:
                        match.residues.delete_residues_inverse(path, path)
                    match = match_list.pop(0)
                    break

                deleted_residues = match_restrict_copy.get_residues_deleted_by_chain(chain)
                changed_residues = change_res_copy.get_residues_changed_by_chain(chain)
                chain, number = utils.get_chain_and_number(path)
                self.template_chains_list.append(TemplateChain(chain=chain,
                                                 path=path,
                                                 code=number,
                                                 sequence=sequence,
                                                 match=match,
                                                 alignment=alignment,
                                                 deleted_residues=deleted_residues,
                                                 changed_residues=changed_residues,
                                                 sequence_before_changes=sequence_before_changes))