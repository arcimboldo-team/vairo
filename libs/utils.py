import base64
import copy
import errno
import glob
import json
import os
import logging
import io
import re
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union
from sklearn import preprocessing
from itertools import groupby, takewhile
from operator import itemgetter


def print_msg_box(msg, indent=1, title=None):
    lines = msg.split('\n')
    space = " " * indent

    width = max(map(len, lines))
    box = f'╔{"═" * (width + indent * 2)}╗\n'  # upper_border
    if title:
        box += f'║{space}{title:<{width}}{space}║\n'  # title
        box += f'║{space}{"-" * len(title):<{width}}{space}║\n'  # underscore
    box += ''.join([f'║{space}{line:<{width}}{space}║\n' for line in lines])
    box += f'╚{"═" * (width + indent * 2)}╝'  # lower_border
    logging.info('\n')
    logging.info(box)
    logging.info('\n')


def print_matrix(matrix: List):
    print('\n'.join('\t'.join(map(str, row)) for row in matrix))


def normalize_list(input_list: List):
    # Normalize list of values

    return preprocessing.normalize(input_list)[0]


def get_mandatory_value(input_load: dict, value: str) -> str:
    # Read value, Raise exception if value could not be found

    read_value = input_load.get(value)
    if read_value is None:
        raise Exception(f'{value} is mandatory')
    return read_value


def get_file_extension(path: str) -> str:
    # Get the extension of the file
    _, extension = os.path.splitext(path)
    return extension


def get_file_name(path: str) -> str:
    # Get the name of the file, without path or extension

    return os.path.splitext(os.path.basename(path))[0]


def get_readme() -> str:
    # Get README.md file

    return os.path.join(os.path.dirname(get_parent_folder(str(Path(__file__)))), 'README.md')


def get_main_path() -> Path:
    # Get the path of the main.py

    return Path(__file__).parent.parent.absolute()


def get_parent_folder(dir_path: str) -> Path:
    return Path(dir_path).parent.absolute()


def get_working_dir() -> str:
    # Get working directory

    return os.getcwd()


def chunk_string(length: int, number_partitions: int, overlap: int) -> List[Tuple[int, int]]:
    # Slice string in chunks of size
    if number_partitions == 1:
        return [(0, length)]
    else:
        reminder = length % number_partitions
        chunk_list = []
        size = int((length - reminder) / number_partitions)
        for chunk in range(0, length - reminder, size):
            chunk_list.append((chunk, size + chunk + overlap))
        last_element = chunk_list[-1]
        chunk_list[-1] = (last_element[0], last_element[1] + reminder - overlap)
        return chunk_list


def dict_values_to_list(input_dict: Dict):
    # Given a Dict, return all the values from the dict in a list

    return [value for value in input_dict.values()]


def get_key_for_value(value: str, search_dict: Dict) -> Union[List[str], None]:
    # Given a value, get the list of all keys that contains that value
    try:
        return list(search_dict.keys())[list(search_dict.values()).index(value)]
    except Exception as e:
        return None


def get_positions_by_chain(path_list: List[str], chain: str) -> List[int]:
    # Give a list of paths, return all the positions in the list
    # that contain the chain

    return [path_list.index(path) for path in get_paths_by_chain(path_list, chain)]


def get_paths_by_chain(path_list: List[str], search_chain: str) -> List[str]:
    # Return all the paths that contain the chain

    return_list = []
    for path in path_list:
        if path is not None and get_chain_and_number(path)[0] == search_chain:
            return_list.append(path)
    return return_list


def get_consecutive_numbers(number_list: List[int]) -> List[Tuple[int, int]]:
    # Given an integer list, return ranges of consecutive numbers

    result_list = []
    for k, g in groupby(enumerate(number_list), lambda x: x[0] - x[1]):
        group = (map(itemgetter(1), g))
        group = list(map(int, group))
        result_list.append((group[0], group[-1]))
    return result_list


def get_chain_and_number(path_pdb: str) -> Tuple[str, int]:
    # Given a path: ../../template_A1.pdb return A and 1
    # Return CHAIN and NUMBER
    name = get_file_name(path_pdb)
    code = name.split('_')[-1]
    return code[0], int(code[1:])


def select_paths_in_dict(chain_dict: Dict, code: str) -> str:
    # Search for the files in all the dict that
    # finish with code

    for _, paths in chain_dict.items():
        for path in paths:
            split_code = get_chain_and_number(path)
            if f'{split_code[0]}{split_code[1]}' == code:
                return path


def get_paths_in_alignment(align_dict: Dict, code: str) -> List[str]:
    # Search for the files in all to align dict that
    # finish with code
    return_list = []
    for _, chain_dict in align_dict.items():
        return_list.append(select_paths_in_dict(chain_dict=chain_dict, code=code))
    return return_list


def select_path_from_code(align_dict: Dict, code: str, position: int, sequence_name_list: List[str]) -> str:
    sequence_name = sequence_name_list[position]
    return select_paths_in_dict(chain_dict=align_dict[sequence_name], code=code)


def replace_last_number(text: str, value: int) -> str:
    # Replace the last number of text by the value

    return re.sub(r'\d+.pdb', str(value), str(text)) + '.pdb'


def expand_residues(res: str) -> List:
    # Expand a str formatted like this: 10-12, 32, 34
    # To a list: [10,11,12,32,34]
    modified_list = str(res).replace(' ', '').split(',')
    return_list = []
    for res in modified_list:
        res_list = str(res).split('-')
        if len(res_list) == 2:
            res_list = list(range(int(res_list[0]), int(res_list[1]) + 1))
        elif len(res_list) > 2:
            raise Exception('Has not been possible to change residues.')
        return_list.extend(map(int, res_list))
    return return_list


def renum_residues(res_list: List[int], mapping: Dict) -> List[int]:
    return [mapping[res] for res in res_list]


def rmsilent(file_path: str):
    # Remove file without error if it doesn't exist

    for file in glob.glob(file_path):
        try:
            os.remove(file)
        except OSError as e:
            if e.errno != errno.ENOENT:
                raise


def clean_files(input_dir: str):
    # Remove directory

    shutil.rmtree(input_dir)


def parse_aleph_annotate(file_path: str) -> Dict:
    # Read the output of aleph, return a dictionary containing:
    # {"ah": 59, "bs": 8, "number_total_residues": 1350}

    with open(file_path) as f_in:
        data = json.load(f_in)
    secondary_structure_dict = copy.deepcopy(data['annotation']['secondary_structure_content'])
    return secondary_structure_dict


def parse_aleph_ss(file_path: str) -> Dict:
    # Parse the aleph.txt file and get all the domains by chain

    chain_res_dict = {}
    with open(file_path) as f_in:
        lines = f_in.readlines()
        for line in lines:
            split_text = line.split()
            if len(split_text) == 12:
                chain = split_text[3]
                residues = expand_residues(f'{split_text[4]}-{split_text[6]}')
                try:
                    chain_res_dict[chain].append(residues)
                except KeyError:
                    chain_res_dict[chain] = [residues]
    return chain_res_dict

def parse_frobenius(file_path: str) -> Dict:
    # Parse the frobenius.txt file, get for each ranked, the coverage of angles and distances.
    
    results_dict = {'dist': {}, 'ang': {}}
    with open(file_path) as f_in:
        lines = f_in.readlines()
    index = [x for x in range(len(lines)) if 'Frobenius distances of CV' in lines[x]]
    distances = list(takewhile(lambda x: x!='\n', lines[index[0]+2:]))
    for distance in distances:
        distance_split = distance.split()
        results_dict['dist'][distance_split[0]] = distance_split[1]
    angles = list(takewhile(lambda x: x!='\n', lines[index[1]+2:]))
    for angle in angles:
        angle_split = angle.split()
        results_dict['ang'][angle_split[0]] = angle_split[1]
    return results_dict


def parse_pisa_general_multimer(pisa_output: str) -> List:
    # It parses the pisa output, in the following format:
    # List[Dict[
    #   area
    #   deltaG
    #   chain1
    #   chain2
    #   serial
    # ]]
    # It returns a list with all the interfaces found, each
    # interface contains the requiered information.

    return_list = []
    match1 = [m.start() for m in re.finditer(' LIST OF INTERFACES', pisa_output)][0]
    match2 = [m.start() for m in re.finditer(' ##: {2}serial number', pisa_output)][0]
    for line in pisa_output[match1:match2].split('\n')[4:-2]:
        line = line.split('|')
        area = line[3][:8].replace(' ', '')
        deltag = line[3][8:15].replace(' ', '')
        chain1 = line[1].replace(' ', '')
        chain2 = line[2].split()[0].replace(' ', '')
        serial = line[0][:4].replace(' ', '')
        return_list.append({'serial': serial, 'area': area, 'deltaG': deltag, 'chain1': chain1, 'chain2': chain2})

    return return_list


def parse_pisa_interfaces(pisa_output: str) -> Dict:
    # It parses the pisa output, in the following format:
    # List[Dict{
    #   solvation1
    #   solvation2
    #   se_gain1
    #   se_gain2
    #   chain1
    #   res_chain1
    #   chain2
    #   res_chain2
    # }]
    # It returns a list with the interface information, each
    # interface contains the required information.

    iter_list = iter(pisa_output.split('\n'))
    res_chain1 = []
    res_chain2 = []
    chain1 = chain2 = ''
    solvation1 = solvation2 = ''
    se_gain1 = se_gain2 = ''
    for line in iter_list:
        if 'Interfacing Residues: Structure' in line:
            next(iter_list)
            next(iter_list)
            next(iter_list)
            line = next(iter_list)
            res_list = []
            chain = ''
            while line != " -----'-'------------'--'----------------------":
                chain = line[10:11]
                energy = line[39:].replace(' ', '')
                if float(energy) != 0:
                    res_num = line[15:20].replace(' ', '')
                    res_list.append(int(res_num))
                line = next(iter_list)
            if not res_chain1:
                chain1 = chain
                res_chain1 = res_list
            else:
                chain2 = chain
                res_chain2 = res_list
        elif 'Solvation energy kcal/mol' in line:
            solvation1 = float(line.split('|')[1].replace(' ', ''))
            solvation2 = float(line.split('|')[2].replace(' ', ''))
        elif 'SE gain, kcal/mol' in line:
            se_gain1 = line.split('|')[1].replace(' ', '')
            se_gain2 = line.split('|')[2].replace(' ', '')

    return {'solvation1': solvation1, 'solvation2': solvation2, 'se_gain1': se_gain1,
            'se_gain2': se_gain2, 'chain1': chain1, 'res_chain1': res_chain1,
            'chain2': chain2, 'res_chain2': res_chain2}


def sort_by_digit(container: Any, item: int = 0):
    # Sort list or dictionary by a digit instead of str.
    # Dict can be like this:

    if isinstance(container, dict):
        return sorted(container.items(), key=lambda x: int("".join([i for i in x[item] if i.isdigit()])))
    elif isinstance(container, list):
        return sorted(container, key=lambda x: int("".join([i for i in x if i.isdigit()])))


def create_dir(dir_path: str, delete_if_exists: bool = False):
    # If directory not exists, create it
    # If directory exists, delete it and create it
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    elif delete_if_exists:
        shutil.rmtree(dir_path)
        os.makedirs(dir_path)


def remove_list_layer(input_list: List[List[str]]) -> List[str]:
    return [j for x in input_list for j in x]


def encode_data(input_data):
    return base64.b64encode(open(input_data, 'rb').read()).decode('utf-8')

def create_logger():
    # Create logger: The information will be stored in a buffer instead of a file. The buffer can be dumped to
    # a file later.

    logger = logging.getLogger()
    logger.handlers.clear()

    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

    logger.setLevel(logging.DEBUG)
    test = io.StringIO()
    stream_handler_ = logging.StreamHandler(test)
    stream_handler_.setLevel(logging.INFO)
    logger.addHandler(stream_handler_)

    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.INFO)
    logger.addHandler(stdout_handler)


def create_logger_dir(log_path: str):
    # Create logger in a working directory with a specific name:

    logger = logging.getLogger()
    logger_data = logger.handlers[0].stream.getvalue()
    logger.removeHandler(logger.handlers[0])
    with open(log_path, 'w+') as f_handle:
        f_handle.write(logger_data)
    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(logging.DEBUG)
    logger.addHandler(file_handler)
