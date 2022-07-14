import os
import logging
from libs import utils
from typing import List

class Template:

    def __init__ (self, parameters_list: List, output_dir=""):
        self.pdb: str
        self.chain: str = ""
        self.polyala_res_list: List = []
        self.custom: bool = True
        self.add_to_msa: bool = False
        self.add_to_templates: bool = False
        self.sum_prob: bool = False
        self.aligned: bool = False
        
        self.custom = parameters_list.get('custom', self.custom)
        self.pdb = self.__check_pdb(parameters_list['pdb'], output_dir)
        self.chain = parameters_list.get('chain', self.chain)
        self.polyala_res_list = parameters_list.get('polyala_res_list', self.polyala_res_list)
        self.add_to_msa = parameters_list.get('add_to_msa', self.add_to_msa)
        self.add_to_templates = parameters_list.get('add_to_templates', self.add_to_templates)
        self.sum_prob = parameters_list.get('sum_prob', self.sum_prob)
        self.aligned = parameters_list.get('aligned', self.aligned)
    
    def __check_pdb(self, pdb, output_dir):
        if not os.path.exists(pdb) and output_dir != "":
            if os.path.exists(f'{output_dir}/{pdb}.pdb'):
                logging.info(f'{output_dir}/{pdb}.pdb already exists in {output_dir}.')
            else:
                utils.download_pdb(pdb_id=pdb, output_dir=output_dir)
            self.custom = False
            pdb = f'{output_dir}/{pdb}.pdb'            
        return pdb