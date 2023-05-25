from typing import Union, List
from libs import change_res

class MatchRestrictions:
    def __init__(self, chain: str, position: int, residues,
                 reference: str, reference_chain: str):
        self.chain: str
        self.position: int = -1
        self.mask_region: bool
        self.residues: change_res.ChangeResidues = None
        self.reference
        self.reference_chain: Union[str, None]

        self.chain = chain
        self.position = position
        self.residues = residues
        self.reference = reference
        self.reference_chain = reference_chain

    def set_reference(self, new_reference):
        # Change the reference from pdb_id to the real Template
        self.reference = new_reference

    def check_references(self) -> bool:
        return self.reference is not None and self.reference_chain is not None

    def check_position(self) -> bool:
        return self.position != -1

    def get_deleted_residues(self, chain: str) -> List[int]:
        if self.residues:
            return self.residues.chain_res_dict[chain]


class MatchRestrictionsList:
    def __init__(self):
        self.match_restrict_list: List[MatchRestrictions] = []

    def append_match(self, chain: str, position: int, residues: change_res.ChangeResidues,
                     reference: str, reference_chain: str):
        self.match_restrict_list.append(MatchRestrictions(chain=chain, position=position, residues=residues,
                                                          reference=reference, reference_chain=reference_chain))

    def get_reference_list(self) -> List[str]:
        # Get all the references from another template
        # Those references can be in the match class or in the
        # template itself
        return [match.reference for match in self.match_restrict_list]

    def get_matches_by_chain(self, chain: str) -> List[MatchRestrictions]:
        # Return all the matches for a specific chain
        return [match for match in self.match_restrict_list if match.chain == chain]
