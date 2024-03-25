from typing import List, Dict
import os

from vairo.libs import global_variables, bioutils, utils
from Bio.PDB import Select, PDBIO
from alphafold.common import residue_constants


class ResidueReplace:
    def __init__(self, residues: List[int], by: str, when: str = ''):
        residues: List[int]
        by: str
        when: str
        isSequence: bool = False

        residues = residues
        by = by
        when = when
        if by not in global_variables.ID_TO_HHBLITS_AA_3LETTER_CODE.values():
            isSequence = True
            if os.path.exists(by):
                by = bioutils.extract_sequence(fasta_path=by)


class ChainModifications:
    def __init__(self, chain: str, bfactors: List[float] = [], position: int = None,
                 accepted_residues: List[int] = [], deleted_residues: List[int] = [],
                 mutations: List[ResidueReplace] = []):
        # It will define the changes to apply to a template.
        # A list with all the chains that those changes will be applied.
        # A list of bfactors for each chain.
        # A position where, if it is just one chain, will be located in the query sequence
        # Accepted residues for each chain, the other residues will be deleted.
        # Deleted residues for each chain, the other will be accepted.
        chain: str
        bfactors: List[float]
        position: int = -1
        accepted_residues: List[int] = []
        deleted_residues: List[int] = []
        mutations: List[ResidueReplace] = []

        self.chain = chain
        self.bfactors = bfactors
        self.position = position
        self.accepted_residues = accepted_residues
        self.deleted_residues = deleted_residues
        self.mutations = mutations

    def apply_mapping(self, mapping: Dict):
        def convert(residues, mapp):
            results = [utils.get_key_by_value(res, mapp) for res in residues]
            return [x[0] for x in results if x]

        self.accepted_residues = convert(self.accepted_residues, mapping)
        self.deleted_residues = convert(self.deleted_residues, mapping)
        for mutation in self.mutations:
            mutation.residues = convert(mutation.residues, mapping)

    def get_change(self, resseq: int, when: str = '') -> str:
        name = None
        for mutation in self.mutations:
            if resseq in mutation.residues:
                if mutation.isSequence and len(mutation.by) > resseq - 1 and when == when:
                    name = utils.get_key_by_value(value=mutation.by[resseq - 1],
                                                  search_dict=residue_constants.restype_3to1)[0]
                else:
                    name = mutation.by
                break
        return name

    def check_position(self) -> bool:
        return self.position != -1

    def set_position(self, position):
        self.position = position

    def get_deleted_residues(self) -> List[int]:
        return_list = []
        if self.deleted_residues:
            return_list.extend(self.deleted_residues)
        if self.accepted_residues:
            return_list.extend(list((1001 - num for num in self.accepted_residues)))
        return return_list


class TemplateModifications:
    def __init__(self):
        self.modifications_list: List[ChainModifications] = []

    def append_change(self, chains: List[str], position: int = None, accepted_residues: List[List[int]] = [],
                      deleted_residues: List[List[int]] = [], mutations: List[ResidueReplace] = [],
                      bfactors: List[List[float]] = []):

        for i, chain in enumerate(chains):
            chain_class = ChainModifications(chain=chain,
                                             position=position,
                                             bfactors=bfactors[i] if all(
                                                 isinstance(item, list) for item in bfactors) else bfactors,
                                             accepted_residues=accepted_residues[i] if all(
                                                 isinstance(item, list) for item in
                                                 accepted_residues) else accepted_residues,
                                             deleted_residues=deleted_residues[i] if all(
                                                 isinstance(item, list) for item in
                                                 deleted_residues) else deleted_residues,
                                             mutations=mutations)

            self.modifications_list.append(chain_class)

    def apply_mapping(self, chain: str, mapping: Dict):
        # Change residues numbering by the ones in mapping
        for modification in self.modifications_list:
            if modification.chain == chain:
                modification.apply_mapping(mapping)

    def get_modifications_by_chain(self, chain: str, when: str = '') -> List[ChainModifications]:
        return [modification for modification in self.modifications_list if
                modification.chain == chain and (when == '' or modification.when == when)]

    def get_modifications_position_by_chain_position(self, chain: str) -> List[ChainModifications]:
        # Return all the matches for a specific chain
        return [modification for modification in self.modifications_list if
                modification.chain == chain and modification.check_position()]

    def get_residues_changed_by_chain(self, chain: str) -> List:
        # Return all the changes for a specific chain.
        # In the dict, there will be the residue name as key
        # and all the residues to change in a list
        fasta = set()
        resname = set()
        modification_chain = self.get_modifications_by_chain(chain=chain)
        for modification in modification_chain:
            for mutation in modification.mutations:
                if mutation.isSequence:
                    fasta.update(mutation.residues)
                elif mutation.resname:
                    resname.update(mutation.residues)
        return list(resname), list(fasta)

    def get_deleted_residues(self, chain: str) -> List[int]:
        delete_list = []
        modification_chain = self.get_modifications_by_chain(chain=chain)
        for modification in modification_chain:
            delete_list.expand(modification.deleted_residues())
        return delete_list

    def modify_template(self, pdb_in_path: str, pdb_out_path: str, when: str = '', type: str = ''):
        # Change residues of chains specified in chain_res_dict
        structure = bioutils.get_structure(pdb_in_path)
        chains_struct = bioutils.get_chains(pdb_in_path)
        chains_change = list(self.chain_res_dict.keys())
        chains_inter = set(chains_struct).intersection(chains_change)
        atoms_del_list = []
        res_del_dict = {}
        for chain in chains_inter:
            modification_chain = self.get_modifications_by_chain(chain=chain)
            res_del_dict[chain] = []
            for i, res in enumerate(structure[0][chain].get_residues()):
                resseq = bioutils.get_resseq(res)
                for modify in modification_chain:
                    if (modify.accepted_residues and resseq not in modify.accepted_residues) or res in modify.deleted_residues:
                        res_del_dict[chain].append(res.id)
                    change_name = modify.get_change(resseq, when)
                    if change_name is not None:
                        for atom in res:
                            res.resname = change_name
                            if not atom.name in residue_constants.residue_atoms[res.resname]:
                                atoms_del_list.append(atom.get_serial_number())

                    if modify.bfactors:
                        for atom in res:
                            atom.set_bfactor(modify.bfactors[i])

        for key_chain, residue_list in res_del_dict.items():
            chain = structure[0][key_chain]
            [chain.detach_child(id) for id in residue_list]

        class AtomSelect(Select):
            def accept_atom(self, atom):
                return not atom.get_serial_number() in atoms_del_list

        io = PDBIO()
        io.set_structure(structure)
        if atoms_del_list:
            io.save(pdb_out_path, select=AtomSelect())
        else:
            io.save(pdb_out_path)
