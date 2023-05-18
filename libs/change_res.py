from typing import Dict, Union, List
from Bio.PDB import Select, PDBIO
from alphafold.common import residue_constants
from libs import bioutils, utils


class ChangeResidues:

    def __init__(self, chain_res_dict: Dict, resname: str = None, chain_bfactors_dict: Dict = None,
                 fasta_path: str = None, when: str = 'after_alignment'):
        # Read parameters and create ChangeResidues class
        # The change is a mandatory value, can be a dict, a list or an int
        # If change is a list, the specific residues will be changed from
        # all the chains specified in the chain_list.
        # If change is a dict, only the chain will be changed.

        self.chain_res_dict: Dict
        self.chain_group_res_dict: Dict = None
        self.chain_bfactors_dict: Union[Dict, None] = None
        self.when: str = 'after_alignment'
        self.resname: Union[str, None] = None
        self.sequence: Union[str, None] = None
        self.fasta_path: Union[str, None] = None

        self.resname = resname
        self.chain_res_dict = chain_res_dict
        self.chain_bfactors_dict = chain_bfactors_dict
        if fasta_path is not None:
            self.sequence = bioutils.extract_sequence(fasta_path=fasta_path)
            self.fasta_path = fasta_path
        self.when = when
        self.group_change_res()

    def apply_mapping(self, chain: str, mapping: Dict):
        # Change residues numbering by the ones in mapping
        residues = self.chain_res_dict.get(chain)
        if residues:
            results = [utils.get_key_by_value(res, mapping) for res in residues]
            self.chain_res_dict[chain] = [x[0] for x in results if x]

    def group_change_res(self):
        self.chain_group_res_dict = {}
        for key, value in self.chain_res_dict.items():
            grouped_list = utils.get_consecutive_numbers(value)
            for i, group_range in enumerate(grouped_list):
                grouped_list[i] = f'{group_range[0]}-{group_range[1]}' if group_range[0] != group_range[1] else str(
                    group_range[0])
            self.chain_group_res_dict[key] = grouped_list

    def delete_residues(self, pdb_in_path: str, pdb_out_path: str):
        self.__change_residues(pdb_in_path, pdb_out_path, 'delete')

    def delete_residues_inverse(self, pdb_in_path: str, pdb_out_path: str):
        self.__change_residues(pdb_in_path, pdb_out_path, 'delete_inverse')

    def change_bfactors(self, pdb_in_path: str, pdb_out_path: str):
        self.__change_residues(pdb_in_path, pdb_out_path, 'change_bfactors')

    def change_residues(self, pdb_in_path: str, pdb_out_path: str):
        self.__change_residues(pdb_in_path, pdb_out_path, 'change')

    def __change_residues(self, pdb_in_path: str, pdb_out_path: str, type: str):
        # Change residues of chains specified in chain_res_dict

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
                            name = utils.get_key_by_value(value=self.sequence[bioutils.get_resseq(res) - 1],
                                                          search_dict=residue_constants.restype_3to1)[0]

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
        io.save(pdb_out_path, select=AtomSelect())


class ChangeResiduesList:
    def __init__(self):
        self.change_residues_list: List[ChangeResidues] = []

    def get_residues_changed_by_chain(self, chain: str) -> List:
        # Return all the changes for a specific chain.
        # In the dict, there will be the residue name as key
        # and all the residues to change in a list
        fasta = set()
        resname = set()
        for change in self.change_residues_list:
            residues = change.chain_res_dict.get(chain)
            if residues:
                if change.fasta_path:
                    fasta.update(residues)
                elif change.resname:
                    resname.update(residues)
        return list(resname), list(fasta)

    def get_changes_by_chain(self, chain: str, when: str = '') -> List[ChangeResidues]:
        return [change for change in self.change_residues_list if
                change.chain_res_dict.get(chain) and (when == '' or change.when == when)]

    def append_change(self, chain_res_dict: Dict, resname: str, fasta_path: str, when: str):
        self.change_residues_list.append(ChangeResidues(chain_res_dict=chain_res_dict,
                                                        resname=resname,
                                                        fasta_path=fasta_path,
                                                        when=when))
