import logging
from typing import Dict, List, Union
from libs import bioutils, features, utils
from Bio.PDB import Select, PDBIO
from alphafold.common import residue_constants


class ChangeResidues:

    def __init__(self, chain_res_dict: Dict, resname: str = None, chain_bfactors_dict: Dict = None, fasta_path: str = None, when: str = 'after_alignment'):
        # Read parameters and create ChangeResidues class
        # The change is a mandatory value, can be a dict, a list or an int
        # If change is a list, the specific residues will be changed from
        # all the chains specified in the chain_list.
        # If change is a dict, only the chain will be changed.

        self.chain_res_dict: Dict
        self.chain_bfactors_dict: Union[Dict, None] = None
        self.when: str = 'after_alignment'
        self.resname: Union[str, None] = None
        self.sequence: Union[str, None] = None

        self.resname = resname
        self.chain_res_dict = chain_res_dict
        self.chain_bfactors_dict = chain_bfactors_dict
        if fasta_path is not None:
            self.sequence = bioutils.extract_sequence(fasta_path=fasta_path)
        self.when = when

        if self.sequence is not None:
            logging.info(f'The following residues are going to be converted to {self.sequence}: {self.chain_res_dict}')
  
        if self.resname is not None:
            logging.info(f'The following residues are going to be converted to {self.resname}: {self.chain_res_dict}')

    def apply_mapping(self, chain: str, mapping: Dict):
        # Change residues numbering by the ones in mapping
        if chain in self.chain_res_dict:
            residues = self.chain_res_dict[chain]
            results = [utils.get_key_for_value(res, mapping) for res in residues]
            self.chain_res_dict[chain] = [x for x in results if x is not None]

    def delete_residues(self, pdb_in_path: str, pdb_out_path: str):
        self.__change_residues(pdb_in_path, pdb_out_path, 'delete')

    def delete_residues_inverse(self, pdb_in_path: str, pdb_out_path: str):
        self.__change_residues(pdb_in_path, pdb_out_path, 'delete_inverse')

    def change_bfactors(self, pdb_in_path: str, pdb_out_path: str):
        self.__change_residues(pdb_in_path, pdb_out_path, 'change_bfactors')

    def change_residues(self, pdb_in_path: str, pdb_out_path: str):
        self.__change_residues(pdb_in_path, pdb_out_path, 'change')

    def __change_residues(self, pdb_in_path: str, pdb_out_path: str, type: str):
        # Chainge residues of chains specified in chain_res_dict

        structure = bioutils.get_structure(pdb_in_path)
        chains_struct = bioutils.get_chains(pdb_in_path)
        chains_change = list(self.chain_res_dict.keys())
        chains_inter = set(chains_struct).intersection(chains_change)
        atoms_del_list = []

        for chain in chains_inter:
            for res in structure[0][chain]:
                if type == 'delete_inverse':
                    if bioutils.get_resseq(res) not in self.chain_res_dict[chain]:
                        for atom in res:
                            atoms_del_list.append(atom.get_serial_number())
                if type == 'delete':
                    if bioutils.get_resseq(res) in self.chain_res_dict[chain]:
                        for atom in res:
                            atoms_del_list.append(atom.get_serial_number())
                if type == 'change':
                    if bioutils.get_resseq(res) in self.chain_res_dict[chain]:
                        if self.resname is not None:
                            name = self.resname
                        elif self.sequence is not None:
                            name = utils.get_key_for_value(value=self.sequence[bioutils.get_resseq(res)-1], search_dict=features.three_to_one)
                        for atom in res:
                            res.resname = name
                            if not atom.name in residue_constants.residue_atoms[res.resname]:
                                atoms_del_list.append(atom.get_serial_number())
                if type == 'change_bfactors':
                    res_bfactor_index = self.chain_res_dict[chain].index(bioutils.get_resseq(res))
                    bfactor = self.chain_bfactors_dict[chain][res_bfactor_index]
                    for atom in res:
                        atom.set_bfactor(bfactor)

        class AtomSelect(Select):
            def accept_atom(self, atom):
                if atom.get_serial_number() in atoms_del_list:
                    return 0
                else:
                    return 1

        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_out_path, select=AtomSelect(), preserve_atom_numbering=True)
