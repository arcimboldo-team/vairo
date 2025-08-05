from typing import List, Dict
import os

from libs import global_variables, bioutils, utils
from Bio.PDB import Select, PDBIO
from alphafold.common import residue_constants


class ResidueMutate:
    def __init__(self, mutate_residues_number: List[int], mutate_with: str):
        self.mutate_residues_number: List[int]
        self.mutate_residues_number_print: str
        self.mutate_with: str
        self.isSequence: bool = False

        self.mutate_residues_number = mutate_residues_number
        self.mutate_with = mutate_with
        if mutate_with not in global_variables.ID_TO_HHBLITS_AA_3LETTER_CODE.values():
            self.isSequence = True
            if os.path.exists(mutate_with):
                self.mutate_with = bioutils.extract_sequence(fasta_path=mutate_with)
        self.mutate_residues_number_print = utils.print_consecutive_numbers(mutate_residues_number)


class ChainModifications:
    def __init__(self, chain: str, bfactors: List[float] = [], position: int = -1,
                 maintain_residues: List[int] = [], delete_residues: List[int] = [],
                 mutations: List[ResidueMutate] = [], when: str = 'after_alignment'):
        # It will define the changes to apply to a template.
        # A list with all the chains that those changes will be applied.
        # A list of bfactors for each chain.
        # A position where, if it is just one chain, will be located in the query sequence
        # Accepted residues for each chain, the other residues will be deleted.
        # Deleted residues for each chain, the other will be accepted.
        self.chain: str
        self.bfactors: List[float]
        self.position: int = -1
        self.maintain_residues: List[int] = []
        self.delete_residues: List[int] = []
        self.mutations: List[ResidueMutate] = []
        self.when: str = 'after_alignment'

        self.chain = chain
        self.bfactors = bfactors
        self.position = position
        self.maintain_residues = maintain_residues
        self.delete_residues = delete_residues
        self.mutations = mutations
        self.when = when

    def apply_mapping(self, mapping: Dict):
        def convert(residues, mapp):
            results = [utils.get_key_by_value(res, mapp) for res in residues]
            return [x[0] for x in results if x]

        self.maintain_residues = convert(self.maintain_residues, mapping)
        self.delete_residues = convert(self.delete_residues, mapping)
        for mutation in self.mutations:
            mutation.mutate_residues_number = convert(mutation.mutate_residues_number, mapping)

    def get_change(self, resseq: int, when: str = '') -> str:
        name = None
        for mutation in self.mutations:
            if resseq in mutation.mutate_residues_number:
                if mutation.isSequence and len(mutation.mutate_with) > resseq - 1 and (self.when == when or when == ''):
                    name = utils.get_key_by_value(value=mutation.mutate_with[resseq - 1],
                                                  search_dict=residue_constants.restype_3to1)[0]
                else:
                    name = mutation.mutate_with
                break
        return name

    def get_deleted_residues(self) -> List[int]:
        return_list = []
        if self.delete_residues:
            return_list.extend(self.delete_residues)
        if self.maintain_residues:
            return_list.extend(list((1001 - num for num in self.maintain_residues)))
        return return_list


class TemplateModifications:
    def __init__(self, modifications_list: List[ChainModifications] = []):
        self.modifications_list: List[ChainModifications] = []
        if modifications_list:
            self.modifications_list = modifications_list

    def append_modification(self, chains: List[str], position: int = -1, maintain_residues: List[List[int]] = [],
                            delete_residues: List[List[int]] = [], mutations: List[ResidueMutate] = [],
                            bfactors: List[List[float]] = [], when: str = 'after_alignment'):

        for i, chain in enumerate(chains):
            chain_class = ChainModifications(chain=chain,
                                             position=position,
                                             when=when,
                                             bfactors=bfactors[i] if bfactors and all(
                                                 isinstance(item, list) for item in bfactors) else bfactors,
                                             maintain_residues=maintain_residues[i] if maintain_residues and all(
                                                 isinstance(item, list) for item in
                                                 maintain_residues) else maintain_residues,
                                             delete_residues=delete_residues[i] if delete_residues and all(
                                                 isinstance(item, list) for item in
                                                 delete_residues) else delete_residues,
                                             mutations=mutations)

            self.modifications_list.append(chain_class)

    def append_chain_modification(self, chain_mod: ChainModifications):
        self.modifications_list.append(chain_mod)

    def append_chain_modifications(self, chain_mod: List[ChainModifications]):
        self.modifications_list.extend(chain_mod)

    def apply_mapping(self, mapping: Dict, chain: str = 'all'):
        # Change residues numbering by the ones in mapping
        for modification in self.modifications_list:
            if chain in [modification.chain, 'all'] or modification.chain == 'all':
                modification.apply_mapping(mapping)

    def check_position(self, chain: str = 'all'):
        for modification in self.modifications_list:
            if (chain in [modification.chain, 'all'] or modification.chain == 'all') and modification.position != -1:
                return modification.position
        return -1

    def get_modifications_by_chain(self, chain: str = 'all', when: str = '') -> List[ChainModifications]:
        return [modification for modification in self.modifications_list if
                (chain in [modification.chain, 'all'] or modification.chain == 'all') and (
                            when == '' or modification.when == when)]

    def get_modifications_by_chain_and_position(self, position: int, chain: str = 'all') -> List[ChainModifications]:
        # Return all the matches for a specific chain and position
        return [modification for modification in self.modifications_list if
                (chain in [modification.chain, 'all'] or modification.chain == 'all') and (
                        modification.position == -1 or modification.position == position)]

    def get_modifications_position_by_chain(self, chain: str = 'all') -> List[ChainModifications]:
        # Return all the matches for a specific chain
        return [modification for modification in self.modifications_list if
                (chain in [modification.chain, 'all'] or modification.chain == 'all') and modification.position != -1]

    def get_residues_changed_by_chain(self, chain: str = 'all') -> List:
        # Return all the changes for a specific chain.
        # In the dict, there will be the residue name as a key
        # and all the residues to change in a list
        fasta = set()
        resname = set()
        modification_chain = self.get_modifications_by_chain(chain=chain)
        for modification in modification_chain:
            for mutation in modification.mutations:
                if mutation.isSequence:
                    fasta.update(mutation.mutate_residues_number)
                else:
                    resname.update(mutation.mutate_residues_number)
        return list(resname), list(fasta)

    def get_deleted_residues(self, chain: str = 'all') -> List[int]:
        delete_list = []
        modification_chain = self.get_modifications_by_chain(chain=chain)
        for modification in modification_chain:
            delete_list.extend(modification.delete_residues)
        return delete_list

    def modify_template(self, pdb_in_path: str, pdb_out_path: str, type_modify: List[str], when: str = ''):
        # Change residues of chains specified in chain_res_dict
        structure = bioutils.get_structure(pdb_in_path)
        chains_struct = bioutils.get_chains(pdb_in_path)
        atoms_del_list = []
        res_del_dict = {}
        for chain in chains_struct:
            modification_chain = self.get_modifications_by_chain(chain=chain, when=when)
            res_del_dict[chain] = []
            for i, res in enumerate(structure[0][chain].get_residues()):
                resseq = bioutils.get_resseq(res)
                for modify in modification_chain:
                    if 'delete' in type_modify:
                        if (
                                modify.maintain_residues and resseq not in modify.maintain_residues) or resseq in modify.delete_residues:
                            res_del_dict[chain].append(res.id)

                    if 'mutate' in type_modify:
                        change_name = modify.get_change(resseq, when)
                        if change_name is not None:
                            for atom in res:
                                res.resname = change_name
                                if not atom.name in residue_constants.residue_atoms[res.resname]:
                                    atoms_del_list.append(atom.get_serial_number())

                    if 'bfactors' in type_modify:
                        if modify.bfactors:
                            for atom in res:
                                atom.set_bfactor(modify.bfactors[i])

        for key_chain, residue_list in res_del_dict.items():
            chain = structure[0][key_chain]
            for id in residue_list:
                try:
                    chain.detach_child(id)
                except:
                    pass

        class AtomSelect(Select):
            def accept_atom(self, atom):
                return not atom.get_serial_number() in atoms_del_list

        io = PDBIO()
        io.set_structure(structure)
        if atoms_del_list:
            io.save(pdb_out_path, select=AtomSelect())
        else:
            io.save(pdb_out_path)
