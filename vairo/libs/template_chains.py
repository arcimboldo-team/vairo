import os
import shutil
from typing import List
from libs import utils, bioutils, structures, template_modifications, hhsearch


class TemplateChain:
    def __init__(self,
                 chain: str,
                 path: str,
                 code: str,
                 sequence: str,
                 modifications_list: template_modifications.TemplateModifications = None,
                 ):
        self.chain = chain
        self.path = path
        self.code = code
        self.sequence = sequence
        self.modifications_list = modifications_list
        self.path_before_changes = None
        self.alignment: structures.Alignment = None
        self.deleted_residues = []
        self.changed_residues = []
        self.fasta_residues = []
        self.sequence_before_changes = ''
        self.sequence_before_alignment = ''

    def apply_changes(self, when: str):
        if when == 'before_alignment':
            self.modifications_list.modify_template(pdb_in_path=self.path, pdb_out_path=self.path,
                                                    type_modify=['mutate', 'delete'], when='before_alignment')
            self.sequence_before_alignment = list(bioutils.extract_sequence_msa_from_pdb(self.path).values())[0]
        else:
            self.sequence_before_changes = list(bioutils.extract_sequence_msa_from_pdb(self.path).values())[0]
            self.path_before_changes = os.path.join(os.path.dirname(self.path),
                                                    f'{utils.get_file_name(self.path)}_originalseq.pdb')
            if os.path.exists(self.path_before_changes):
                os.remove(self.path_before_changes)
            shutil.copy2(self.path, self.path_before_changes)
            self.modifications_list.modify_template(pdb_in_path=self.path, pdb_out_path=self.path,
                                                    type_modify=['mutate', 'delete'], when='after_alignment')
            self.deleted_residues = self.modifications_list.get_deleted_residues()
            self.changed_residues, self.fasta_residues = self.modifications_list.get_residues_changed_by_chain()

    def get_chain_code(self) -> List:
        return self.chain, self.code

    def check_alignment(self, stop: bool):
        # Check if it has been a good alignment. If stop, then throw an error.
        if self.alignment:
            if float(self.alignment.evalue) > 0.01:
                if not stop:
                    return False
                else:
                    raise Exception(
                        f'Match could not be done. Poor alignment in the template {utils.get_file_name(self.path)}. '
                        f'Stopping the run.')
            return True
        return True

    def set_alignment(self, alignment: structures.Alignment):
        structure = bioutils.get_structure(self.path)
        residues_list = list(structure[0][self.chain].get_residues())
        idres_list = list([bioutils.get_resseq(res) for res in residues_list])
        mapping_keys = list(map(lambda x: x + 1, list(alignment.mapping.keys())))
        mapping_values = list(map(lambda x: x + 1, list(alignment.mapping.values())))
        mapping = dict(zip(mapping_keys, mapping_values))
        if idres_list != mapping_keys and len(idres_list) == len(mapping_keys):
            self.modifications_list.apply_mapping(mapping=mapping)
        self.alignment = alignment

    def check_position(self):
        return self.modifications_list.check_position()

    
    def set_extracted_chain(self, extracted_path: str):
        self.path = extracted_path

    def __repr__(self):
        # Print class
        return f'pdb_path: {self.path}'


class TemplateChainsList:
    def __init__(self):
        self.template_chains_list: List[TemplateChain] = []

    def get_template_chain(self, pdb_path: str) -> TemplateChain:
        chain1, code1 = utils.get_chain_and_number(pdb_path)
        pdb_dirname = os.path.dirname(pdb_path)
        for template_chain in self.template_chains_list:
            # if there is an alignment, we check the directory too, because it can be from different alignments
            if pdb_dirname == os.path.dirname(
                    template_chain.path) and chain1 == template_chain.chain and code1 == template_chain.code:
                return template_chain
        return None

    def get_old_sequence(self, pdb_path: str) -> str:
        template_chain = self.get_template_chain(pdb_path)
        if template_chain is not None:
            return template_chain.sequence_before_changes
        else:
            return None

    def get_changes(self, pdb_path: str) -> List:
        template_chain = self.get_template_chain(pdb_path)
        if template_chain is not None:
            return template_chain.changed_residues, template_chain.fasta_residues, template_chain.deleted_residues, template_chain.sequence_before_changes
        else:
            return None

    def get_number_chains(self) -> int:
        return len({(chain_template.chain, chain_template.code) for chain_template in self.template_chains_list})

    def get_alignment_by_path(self, pdb_path: str):
        # Search for the alignment that has the same name as the pdb_path
        temp_chain = self.get_template_chain(pdb_path)
        if temp_chain is not None and temp_chain.alignment:
            return temp_chain.alignment
        return None

    def get_chains_not_in_list(self, input_list: List[str]) -> List[TemplateChain]:
        # It will create a List of TemplateChains where the chain and the code is not in the input list.
        # This can be used to eliminate duplicates in case there is more than one database.
        input_chains = [(utils.get_chain_and_number(path)[0], utils.get_chain_and_number(path)[1]) for path in
                        input_list if path is not None]
        return [temp_chain for temp_chain in self.template_chains_list if
                (temp_chain.chain, temp_chain.code) not in input_chains]

    def get_chains_with_matches_pos(self) -> List[TemplateChain]:
        # Get all the chains that have a match with a determinate position.(
        return [temp_chain for temp_chain in self.template_chains_list if temp_chain.check_position() != -1]

    def apply_changes(self, when: str = 'after_alignment'):
        for temp_chain in self.template_chains_list:
            temp_chain.apply_changes(when)

    def new_chain_sequence(self, path: str,
                           sequence: str,
                           modifications_list: template_modifications.TemplateModifications):
        chain, number = utils.get_chain_and_number(path)
        template_chain_struct = TemplateChain(chain=chain, path=path, code=number, sequence=sequence,
                                              modifications_list=modifications_list)
        template_chain_struct.apply_changes('before_alignment')
        self.template_chains_list.append(template_chain_struct)

    def set_same_alignment(self, chain_struct: TemplateChain):
        # If some chain struct share the same alignment, as they have the same sequence,
        # copy the same information, otherwise, the alignment should have to be done.
        for chain_s in self.template_chains_list:
            if chain_s.alignment and chain_s.sequence_before_alignment == chain_struct.sequence_before_alignment and chain_s.sequence == chain_struct.sequence:
                return chain_s.alignment.hhr_path
        return None
