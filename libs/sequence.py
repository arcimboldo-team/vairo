import os
import shutil
from typing import Dict, List
from libs import bioutils, utils


class Sequence:
    def __init__(self, parameters_dict: Dict, input_dir: str):
        self.fasta_path: str
        self.sequence: str
        self.name: str
        self.length: int
        self.num_of_copies: int = 1
        self.positions: List[int] = []

        fasta_path = utils.get_mandatory_value(input_load=parameters_dict, value='fasta_path')
        positions = parameters_dict.get('positions', '')
        if positions == '':
            self.num_of_copies = parameters_dict.get('num_of_copies', self.num_of_copies)
            self.positions = [-1] * self.num_of_copies
        else:
            positions_list = positions.replace(' ', '').split(',')
            for position in positions_list:
                position = int(position) - 1 if int(position) != -1 else int(position)
                self.positions.append(position)
            self.num_of_copies = len(self.positions)

        if self.num_of_copies == 0:
            raise Exception(f'Set num_of_copies or positions for sequence {fasta_path}')

        if not os.path.exists(fasta_path):
            raise Exception(f'{fasta_path} does not exist')
        else:
            self.fasta_path = os.path.join(input_dir, os.path.basename(fasta_path))
            try:
                shutil.copy2(fasta_path, self.fasta_path)
            except shutil.SameFileError:
                pass
            
            self.name = utils.get_file_name(self.fasta_path)
            self.sequence = bioutils.extract_sequence(self.fasta_path)
            self.length = len(self.sequence)


class SequenceAssembled:
    def __init__(self, sequence_list: List[Sequence], glycines: int):

        self.sequence_assembled: str = ''
        self.sequence_list: List[Sequence] = []
        self.sequence_list_expanded: List[Sequence] = []
        self.length: str
        self.glycines: int = glycines
        self.total_copies: int = 0

        self.total_copies = sum([sequence.num_of_copies for sequence in sequence_list])
        positions_to_fill = []
        self.sequence_list_expanded = [None] * self.total_copies
        self.sequence_list = sequence_list

        for sequence in sequence_list:
            for position in sequence.positions:
                if position == -1:
                    positions_to_fill.append(sequence)
                else:
                    if self.sequence_list_expanded[position] is None:
                        self.sequence_list_expanded[position] = sequence
                    else:
                        raise Exception('Wrong sequence requirements. Review sequence positions')

        for i, position in enumerate(self.sequence_list_expanded):
            if position is None:
                sequence = positions_to_fill.pop(0)
                self.sequence_list_expanded[i] = sequence

            self.sequence_assembled += self.sequence_list_expanded[i].sequence + self.glycines * 'G'

        self.sequence_assembled = self.sequence_assembled[:-self.glycines]
        self.length = len(self.sequence_assembled)

    def get_sequence_length(self, i: int) -> int:

        return len(self.sequence_list_expanded[i].sequence)

    def get_sequence_name(self, i: int) -> str:

        return self.sequence_list_expanded[i].name

    def get_list_name(self) -> List[str]:

        return [sequence.name for sequence in self.sequence_list_expanded]

    def get_starting_length(self, i: int) -> int:

        offset = 0
        for j in range(i):
            offset += len(self.sequence_list_expanded[j].sequence) + self.glycines
        return offset
