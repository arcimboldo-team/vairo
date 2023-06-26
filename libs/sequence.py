import collections
import os
import shutil
from typing import Dict, List, Tuple
from libs import bioutils, utils
from alphafold.common import residue_constants


class Sequence:
    def __init__(self, parameters_dict: Dict, input_dir: str):
        self.fasta_path: str
        self.fasta_mutated_path: str
        self.sequence: str
        self.sequence_mutated: str
        self.name: str
        self.length: int
        self.num_of_copies: int
        self.positions: List[int] = []
        self.mutations_dict: Dict = {}

        fasta_path = utils.get_input_value(name='fasta_path', section='sequence', input_dict=parameters_dict)
        positions = utils.get_input_value(name='positions', section='sequence', input_dict=parameters_dict)
        if positions is None:
            self.num_of_copies = utils.get_input_value(name='num_of_copies', section='sequence',
                                                       input_dict=parameters_dict)
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

            self.name = utils.get_input_value(name='name', section='sequence', input_dict=parameters_dict)
            if self.name is None:
                self.name = utils.get_file_name(self.fasta_path)
            self.sequence = bioutils.extract_sequence(self.fasta_path)
            self.length = len(self.sequence)

        self.sequence_mutated = list(self.sequence)
        mutations = utils.get_input_value(name='mutations', section='sequence', input_dict=parameters_dict)
        if mutations:
            for mutation in mutations:
                key = list(mutation.keys())[0]
                values = utils.expand_residues(list(mutation.values())[0])
                if key not in list(residue_constants.restype_1to3.keys()):
                    raise Exception(
                        f'Mutation residues {"".join(values)} in {key} could not be possible. Residue {key} does not exist')
                self.mutations_dict.setdefault(key, []).extend(values)
                for value in values:
                    if value <= len(self.sequence):
                        self.sequence_mutated[value - 1] = key
        self.sequence_mutated = ''.join(self.sequence_mutated)

        mutated_name = f'{self.name}_mutated'
        self.fasta_mutated_path = os.path.join(input_dir, f'{mutated_name}.fasta')
        bioutils.write_sequence(sequence_name=mutated_name, sequence_amino=self.sequence_mutated, sequence_path=self.fasta_mutated_path)


class SequenceAssembled:
    def __init__(self, sequence_list: List[Sequence], glycines: int):

        self.sequence_assembled: str = ''
        self.sequence_mutated_assembled: str = ''
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
            self.sequence_mutated_assembled += self.sequence_list_expanded[i].sequence_mutated + self.glycines * 'G'
        self.sequence_assembled = self.sequence_assembled[:-self.glycines]
        self.sequence_mutated_assembled = self.sequence_mutated_assembled[:-self.glycines]
        self.length = len(self.sequence_assembled)

    def get_mutated_residues_list(self) -> List[int]:
        _, changes_dict = bioutils.compare_sequences(self.sequence_assembled, self.sequence_mutated_assembled)
        return list(changes_dict.keys())

    def get_mutated_residues_dict(self) -> Dict:
        _, changes_dict = bioutils.compare_sequences(self.sequence_assembled, self.sequence_mutated_assembled)
        return changes_dict

    def get_sequence_length(self, i: int) -> int:
        return len(self.sequence_list_expanded[i].sequence)

    def get_sequence_name(self, i: int) -> str:
        return self.sequence_list_expanded[i].name

    def get_list_name(self) -> List[str]:
        return [sequence.name for sequence in self.sequence_list_expanded]

    def get_starting_length(self, i: int) -> int:
        # Get the starting position of the assembled sequence.
        offset = 0
        for j in range(i):
            offset += len(self.sequence_list_expanded[j].sequence) + self.glycines
        return offset

    def get_finishing_length(self, i: int) -> int:
        # Return the starting length plus de sequence length, so the number the sequence it finishes
        return self.get_starting_length(i) + self.get_sequence_length(i) - 1

    def get_position_by_residue_number(self, res_num: int) -> int:
        # Get the number of a residue. Return the position of the sequence it belongs
        for i in range(0, self.total_copies):
            if res_num - 1 <= self.get_finishing_length(i):
                return i

    def get_real_residue_number(self, i: int, residue: int) -> int:
        # Given a position (i) and a residue, get the residue number without being split in chains
        init = self.get_starting_length(i)
        if residue + init <= self.get_starting_length(i) + self.get_sequence_length(i):
            return residue + init
        return None

    def get_range_residues(self, position_ini, position_end) -> List[int]:
        return [self.get_starting_length(position_ini), self.get_finishing_length(position_end)]

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
            aprox_length = length / number_partitions
            actual_partition = 0
            partitions = collections.defaultdict(list)
            for i, element in enumerate(length_list):
                if number_partitions - actual_partition == len(length_list) - i and not partitions[actual_partition]:
                    partitions[actual_partition].append(element)
                    actual_partition += 1
                elif (number_partitions - 1) - actual_partition == 0:
                    partitions[actual_partition].append(element)
                else:
                    length_partition = sum([length for length in partitions[actual_partition]])
                    if (length_partition + element) > aprox_length * 1.2:
                        actual_partition += 1
                    partitions[actual_partition].append(element)

            starting_position = 0
            count = 0
            chunk_list = []
            for i in range(0, number_partitions):
                count += len(partitions[i])
                starting_length = self.get_starting_length(count - 1)
                sequence_length = self.get_sequence_length(count - 1)
                end_position = starting_length + int(sequence_length / 2)
                if i == 0:
                    chunk_list.append((int(starting_position), int(end_position + overlap / 2)))
                elif i == number_partitions - 1:
                    chunk_list.append((int(starting_position - overlap / 2), int(length)))
                else:
                    chunk_list.append((int(starting_position - overlap / 2), int(end_position + overlap / 2)))
                starting_position = end_position
            return chunk_list

    def get_percentages(self, path_in: str) -> List[int]:
        structure = bioutils.get_structure(path_in)
        result_list = [0] * self.total_copies
        for residue in structure[0].get_residues():
            i = self.get_position_by_residue_number(bioutils.get_resseq(residue))
            result_list[i] += 1
        for i, num in enumerate(result_list):
            result_list[i] = result_list[i] / self.get_sequence_length(i)
        return result_list
