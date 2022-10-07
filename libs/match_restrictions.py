from curses.ascii import isdigit
import logging
from typing import Dict
from libs import change_res, utils


class MatchRestrictions:

    def __init__ (self, parameters_dict: Dict):
        
        self.chain: str
        self.position: str = ''
        self.residues: change_res.ChangeRes = None
        self.reference = None
        self.reference_chain: str = None

        self.chain = utils.get_mandatory_value(parameters_dict, 'chain')
        self.residues = parameters_dict.get('residues', self.residues)
        if self.residues is not None:
            change_list = utils.expand_residues(self.residues)
            new_dict = {self.chain: change_list}
            self.residues = change_res.ChangeResidues(chain_res_dict=new_dict)
        self.position =  parameters_dict.get('position', self.position)
        if self.position != '' and self.position != 'None' and not str(self.position).isdigit():
            raise Exception(f'Error setting position in match {self.chain}: {self.position}')
        if str(self.position).isdigit():
            self.position = self.position - 1
        self.reference =  parameters_dict.get('reference', self.reference)
        self.reference_chain =  parameters_dict.get('reference_chain', self.reference_chain)

    def set_reference(self, new_reference):
        #Change the reference from pdb_id to the real Template

        self.reference = new_reference