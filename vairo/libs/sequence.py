import collections
import os
import shutil
import logging
from typing import Dict, List, Tuple
from libs import bioutils, utils
from alphafold.common import residue_constants


class Sequence:
    def __init__(self, parameters_dict: Dict, input_dir: str, run_dir: str, predicted: bool):
        self.fasta_path: str
        self.fasta_mutated_path: str
        self.sequence: str
        self.sequence_predicted: str
        self.sequence_mutated: str
        self.name: str
        self.length: int
        self.num_of_copies: int
        self.positions: List[int] = []
        self.mutations_dict: Dict = {}
        self.alignment_dir: str
        self.predict_region: List[int] = []

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
            self.mutations_dict = utils.read_mutations_dict(input_mutations=mutations)
            for residue, numbering in self.mutations_dict.items():
                for value in numbering:
                    if value <= len(self.sequence):
                        self.sequence_mutated[value - 1] = residue

        self.sequence_mutated = ''.join(self.sequence_mutated)

        mutated_name = f'{self.name}_mutated'
        self.fasta_mutated_path = os.path.join(input_dir, f'{mutated_name}.fasta')
        bioutils.write_sequence(sequence_name=mutated_name, sequence_amino=self.sequence_mutated,
                                sequence_path=self.fasta_mutated_path)

        self.alignment_dir = os.path.join(run_dir, self.name)
        utils.create_dir(self.alignment_dir, delete_if_exists=True)

        predict_aux = utils.get_input_value(name='predict_region', section='sequence', input_dict=parameters_dict)
        if predict_aux and predicted:
            predict_values = list(map(int, str(predict_aux).replace(' ', '').split('-')))
            if len(predict_values) == 2:
                self.predict_region = [predict_values[0], predict_values[1]]
                self.sequence_mutated = self.sequence_mutated[self.predict_region[0]-1:self.predict_region[1]]
                self.sequence = self.sequence[self.predict_region[0]-1:self.predict_region[1]]
                self.length = len(self.sequence)
            else:
                raise Exception(f'predict_values input paramter should have an starting value and an ending value')


class SequenceAssembled:
    def __init__(self, sequence_list: List[Sequence], glycines: int):
        self.sequence_assembled: str = ''
        self.sequence_mutated_assembled: str = ''
        self.sequence_list: List[Sequence] = []
        self.sequence_list_expanded: List[Sequence] = []
        self.length: str
        self.glycines: int = glycines
        self.total_copies: int = 0
        self.mutated: bool = False

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
                self.sequence_list_expanded[i] = positions_to_fill.pop(0)

            sequence_part = self.sequence_list_expanded[i].sequence
            sequence_mutated_part = self.sequence_list_expanded[i].sequence_mutated

            self.sequence_assembled += sequence_part + 'G' * self.glycines
            self.sequence_mutated_assembled += sequence_mutated_part + 'G' * self.glycines

        self.sequence_assembled = self.sequence_assembled[:-self.glycines]
        self.sequence_mutated_assembled = self.sequence_mutated_assembled[:-self.glycines]
        self.length = len(self.sequence_assembled)

        if self.sequence_mutated_assembled != self.sequence_assembled:
            self.mutated = True

        if self.total_copies > 1:
            logging.error(f'Merging {self.total_copies} sequences into one, each separated by {self.glycines} glycines')
        logging.error(f'Total size of the query sequence is {self.length} amino acids')

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
            if self.get_starting_length(i) <= res_num - 1 <= self.get_finishing_length(i):
                return i
        return None

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
        if number_partitions == 1:
            return [(0, self.length)]
        elif len(self.sequence_list_expanded) == 1:
            reminder = self.length % number_partitions
            chunk_list = []
            size = int((self.length - reminder) / number_partitions)
            for chunk in range(0, self.length - reminder, size):
                chunk_list.append((chunk, size + chunk + overlap))
            chunk_list[-1] = (chunk_list[-1][0], self.length)
            return chunk_list
        else:
            length_list = [sequence.length for sequence in self.sequence_list_expanded]
            aprox_length = self.length / number_partitions
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
                    chunk_list.append((int(starting_position - overlap / 2), int(self.length)))
                else:
                    chunk_list.append((int(starting_position - overlap / 2), int(end_position + overlap / 2)))
                starting_position = end_position
            return chunk_list

    def get_percentages(self, path_in: str) -> List[float]:
        structure = bioutils.get_structure(path_in)
        result_list = [0] * self.total_copies
        for residue in structure[0].get_residues():
            pos = self.get_position_by_residue_number(bioutils.get_resseq(residue))
            if pos is not None:
                result_list[pos] += 1
        for i, num in enumerate(result_list):
            result_list[i] = result_list[i] / self.get_sequence_length(i)
        return result_list

    def get_percentage_sequence(self, sequence_in: str) -> List[float]:
        result_list = [0] * self.total_copies
        perc_list = [0] * self.total_copies
        for i, residue in enumerate(sequence_in, start=1):
            if residue != 21 and residue != '-':
                pos = self.get_position_by_residue_number(i)
                if pos is not None:
                    result_list[pos] += 1
        for i, num in enumerate(result_list):
            perc_list[i] = result_list[i] / self.get_sequence_length(i)
        return result_list, perc_list

    def get_region_starting_shifts(self) -> List[int]:
        return [(seq.predict_region[0] - 1) if seq.predict_region else 0 for seq in self.sequence_list_expanded]
    
    def get_list_linker_numbering(self) -> List[int]:
        # Return a list with the numbering of all the residues that are linkers
        linkers_list = []
        for i in range(0, self.total_copies):
            ini = self.get_finishing_length(i) + 2
            end = self.get_finishing_length(i) + 1 + self.glycines
            linkers_list.extend(utils.expand_residues(f'{ini}-{end}'))
        return linkers_list
                