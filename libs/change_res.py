import logging
from typing import Dict, List
from libs import bioutils, utils
from Bio.PDB import Select, PDBIO


class ChangeResidues:

    def __init__ (self, parameters_dict: Dict, chain_list: List):
        #Read parameters and create ChangeResidues class
        #The change is a mandatory value, can be a dict or a list
        #If change is a list, the specific residues will be changed from
        #all the chains specified in the chain_list.
        #If change is a dict, only the chain will be changed.

        self.change_dict: Dict = {}
        self.resname: str = 'ALA'
        
        change = utils.get_mandatory_value(parameters_dict, 'change')
        self.resname = parameters_dict.get('resname', self.resname)

        if isinstance(change, dict):
            for key in change.keys():
                if key not in chain_list:
                    raise Exception('Has not been possible to convert template to polyala. '
                                    f'Chain: {key} does not exist. Available chains: {chain_list}.')
        elif isinstance(change, list):
            if len(chain_list) > 1:
                for chain in chain_list:
                  self.change_dict = {chain_list[chain]: change}  
            else:
                self.change_dict = {chain_list[0]: change}
        else:
            raise Exception('Has not been possible to convert template to polyala.')

        for key, value in change.items():
            change_list = []
            for res in value:
                res_list = str(res).replace(' ', '').split('-')
                if len(res_list) == 2:
                    res_list = list(range(int(res_list[0]), int(res_list[1])+1))
                elif len(res_list) > 2:
                    raise Exception('Has not been possible to change residues.')
                change_list.extend(map(int,res_list))
            self.change_dict[key] = list(set(change_list))
            
        logging.info(f'The following residues are going to be converted to {self.resname}: {self.change_dict}')

    
    def update_new_chains(self, update_dict: Dict):
        # Update the chains after changing them generating the monomer
        # new_chains_dict: A: [A, B], B: [C,D]

        change_dict = {}
        for chain, changes_list in self.change_dict.items():
            if chain in update_dict:
                for new_chain in update_dict[chain]:
                    change_dict[new_chain] = changes_list

        self.change_dict = change_dict

    def change_residues(self, pdb_in_path: str, pdb_out_path: str, real_chain:str = None):
        # Change residues from the pdb_in and write the pdb in pdb_out
        # If the chains have been already modified, real_chain can be specified
        # to show which was the real chain.

        structure = bioutils.get_structure(pdb_in_path)
        chains_struct = bioutils.get_chains(structure)
        if real_chain is not None:
            fake_chain = chains_struct[0]
            chains_struct=[real_chain]
        chains_change = list(self.change_dict.keys())
        chains_inter = set(chains_struct).intersection(chains_change)

        res_atoms_list = ['N', 'CA', 'C', 'CB', 'O']
        atoms_del_list = []

        for chain in chains_inter:
            chain2 = chain if real_chain is None else fake_chain
            for res in structure[0][chain2]:
                if bioutils.get_resseq(res) in self.change_dict[chain]:
                    for atom in res:
                        res.resname = self.resname
                        if not atom.name in res_atoms_list:
                            atoms_del_list.append(atom.get_serial_number())

        class Atom_select(Select):
                def accept_atom(self, atom):
                    if atom.get_serial_number() in atoms_del_list:
                        return 0
                    else:
                        return 1
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_out_path, select=Atom_select(), preserve_atom_numbering = True)