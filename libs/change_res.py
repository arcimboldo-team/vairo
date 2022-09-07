import logging
from typing import Dict, List
from libs import bioutils, utils
from Bio.PDB import Select, PDBIO


class ChangeResidues:

    def __init__ (self, change_dict: Dict, resname: str, inverted: bool = False):
        #Read parameters and create ChangeResidues class
        #The change is a mandatory value, can be a dict, a list or a int
        #If change is a list, the specific residues will be changed from
        #all the chains specified in the chain_list.
        #If change is a dict, only the chain will be changed.

        self.change_dict: Dict = {}
        self.resname: str
        self.inverted: bool = False

        self.resname = resname
        self.change_dict = change_dict
        self.inverted = inverted
        
        if resname != '':
            logging.info(f'The following residues are going to be converted to {self.resname}: {self.change_dict}')
    
    def change_residues(self, pdb_in_path: str, pdb_out_path: str):
        # Change residues from the pdb_in and write the pdb in pdb_out

        structure = bioutils.get_structure(pdb_in_path)
        chains_struct = bioutils.get_chains(structure)
        chains_change = list(self.change_dict.keys())
        chains_inter = set(chains_struct).intersection(chains_change)

        res_atoms_list = ['N', 'CA', 'C', 'CB', 'O']
        atoms_del_list = []

        for chain in chains_inter:
            for res in structure[0][chain]:
                if self.inverted:
                    if bioutils.get_resseq(res) not in self.change_dict[chain]:
                        for atom in res:
                            atoms_del_list.append(atom.get_serial_number())
                else:               
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