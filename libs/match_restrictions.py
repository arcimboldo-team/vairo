from typing import Dict, List
from libs import change_res, utils


class MatchRestrictions:

    def __init__ (self, parameters_dict: Dict):
        
        self.chain: str
        self.position: str = None
        self.residues: change_res.ChangeRes = None
        self.reference = None
        self.reference_chain: str = None

        self.chain = utils.get_mandatory_value(parameters_dict, 'chain')
        self.residues = parameters_dict.get('residues', self.residues)
        if self.residues is not None:
            change_list = utils.expand_residues(self.residues)
            new_dict = {self.chain: change_list}
            self.residues = change_res.ChangeResidues(change_dict=new_dict, resname='', inverted=True)
        self.position =  parameters_dict.get('position', self.position)
        self.reference =  parameters_dict.get('reference', self.reference)
        self.reference_chain =  parameters_dict.get('reference_chain', self.reference_chain)

    def set_reference(self, new_reference):
        #Change the reference from pdb_id to the real Template

        self.reference = new_reference