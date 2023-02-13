from typing import Dict, Union
from libs import change_res, utils


class MatchRestrictions:

    def __init__(self, parameters_dict: Dict):

        self.chain: str
        self.position: int = -1
        self.residues: Union[change_res.ChangeResidues, None] = None
        self.reference = None
        self.reference_chain: Union[str, None] = None

        self.chain = utils.get_mandatory_value(parameters_dict, 'chain')
        residues_list = parameters_dict.get('residues', None)
        if residues_list is not None:
            change_list = utils.expand_residues(residues_list)
            new_dict = {self.chain: change_list}
            self.residues = change_res.ChangeResidues(chain_res_dict=new_dict)

        self.position = parameters_dict.get('position', self.position)
        if self.position != -1:
            self.position = self.position - 1
        self.reference = parameters_dict.get('reference', self.reference)
        self.reference_chain = parameters_dict.get('reference_chain', self.reference_chain)

    def set_reference(self, new_reference):
        # Change the reference from pdb_id to the real Template

        self.reference = new_reference
