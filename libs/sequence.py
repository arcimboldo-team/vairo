import os
import shutil
import collections
from typing import Dict, List, Tuple
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
            positions_list = str(positions).replace(' ', '').split(',')
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
            
            self.name = parameters_dict.get('name', utils.get_file_name(self.fasta_path))
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
        #Get the starting position of the assembled sequence.
        offset = 0
        for j in range(i):
            offset += len(self.sequence_list_expanded[j].sequence) + self.glycines
        return offset

    def get_real_residue_number(self, i: int, residue: int) -> int:
        #Given a position (i) and a residue, get the residue number without being splitted in chains
        init = self.get_starting_length(i)
        if residue+init <= self.get_starting_length(i)+self.get_sequence_length(i):
            return residue+init
        return None

    def partition(self, number_partitions: int, overlap: int) -> List[Tuple[int, int]]:
        # Slice string in chunks of size
        length = len(self.sequence_assembled)
        if number_partitions == 1:
            return [(0, length)]
        elif len(self.sequence_list_expanded) == 1:
            reminder = length % number_partitions
            chunk_list = []
            size = int((length - reminder) / number_partitions)
            for chunk in range(0, length - reminder, size):
                chunk_list.append((chunk, size + chunk + overlap))
            last_element = chunk_list[-1]
            chunk_list[-1] = (last_element[0], last_element[1] + reminder - overlap)
            return chunk_list
        else:
            length_list = [sequence.length for sequence in self.sequence_list_expanded]
            aprox_length = length/number_partitions
            actual_partition = 0
            partitions = collections.defaultdict(list)
            for i, element in enumerate(length_list):
                if number_partitions-actual_partition == len(length_list)-i and not partitions[actual_partition]:
                    partitions[actual_partition].append(element)
                    actual_partition += 1
                elif (number_partitions-1)-actual_partition == 0:
                    partitions[actual_partition].append(element)
                else:
                    length_partition = sum([length for length in partitions[actual_partition]])
                    if (length_partition+element) > aprox_length*1.2:
                        actual_partition += 1
                    partitions[actual_partition].append(element)

            starting_position = 0
            count = 0
            chunk_list = []
            for i in range(0, number_partitions):
                count += len(partitions[i])
                starting_length = self.get_starting_length(count-1)
                sequence_length = self.get_sequence_length(count-1)
                end_position = starting_length + int(sequence_length/2)
                if i == 0:
                    chunk_list.append((int(starting_position), int(end_position + overlap/2)))
                elif i == number_partitions-1:
                    chunk_list.append((int(starting_position - overlap/2), int(length)))
                else:
                    chunk_list.append((int(starting_position - overlap/2), int(end_position + overlap/2)))
                starting_position =  end_position
            return chunk_list

            