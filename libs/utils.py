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
import pandas as pd
from pathlib import Path
from typing import Any, Dict, List

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

def get_mandatory_value(input_load: str, value: str) -> str:
    #Read value, Raise exception if value could not be found

    read_value = input_load.get(value)
    if read_value is None:
        raise Exception(f'{value} is mandatory')
    return read_value

def get_file_name(path: str) -> str:
    #Get the name of the file, without path or extension

    return os.path.splitext(os.path.basename(path))[0]

def get_readme() -> str:
    #Get README.md file

    return os.path.join(os.path.dirname(get_parent_folder(Path(__file__))), 'README.md')

def get_main_path() -> str:
    #Get the path of the main.py

    return Path(__file__).parent.parent.absolute()

def get_parent_folder(dir_path: str) -> str:

    return Path(dir_path).parent.absolute()

def get_working_dir() -> str:
    # Get working directory

    return os.getcwd()

def dict_values_to_list(input_dict: Dict):
    # Given a Dict, return all the values from the dict in a list

    return [value for value in input_dict.values()]

def get_key_for_value(value: str, search_dict: Dict) -> List:
    #Given a value, get the list of all keys that contains that value

    return list(search_dict.keys())[list(search_dict.values()).index(value)]

def get_positions_by_chain(path_list: List, chain: str) -> str:
    # Give a list of paths, return all the positions in the list
    # that contain the chain
    
    return [path_list.index(path) for path in get_paths_by_chain(path_list, chain)]

def get_paths_by_chain(path_list: List, search_chain: str) -> List:
    # Return all the paths that contain the chain

    return_list = []
    for path in path_list:
        if path is not None and get_chain_and_number(path)[0] == search_chain:
            return_list.append(path)
    return return_list

def get_chain_and_number(path_pdb: str) -> List:
    #Given a path: ../../template_A1.pdb return A and 1
    #Return CHAIN and NUMBER
    name = get_file_name(path_pdb)
    code = name.split('_')[-1]
    return code[0], int(code[1:])

def replace_last_number(text: str, value: str) -> str:
    #Replace the last number of text by the value

    return re.sub(r'\d+.pdb', str(value), str(text))+'.pdb'

def expand_residues(res: str) -> List:
    #Expand a str formatted like this: 10-12, 32, 34
    #To a list: [10,11,12,32,34]
    modified_list = str(res).replace(' ', '').split(',')
    return_list = []
    for res in modified_list:
        res_list = str(res).split('-')
        if len(res_list) == 2:
            res_list = list(range(int(res_list[0]), int(res_list[1])+1))
        elif len(res_list) > 2:
            raise Exception('Has not been possible to change residues.')
        return_list.extend(map(int,res_list))
    return return_list

def rmsilent(file_path: str):
    #Remove file without error if it doesn't exist

    for file in glob.glob(file_path):
        try:
            os.remove(file)
        except OSError as e:
            if e.errno != errno.ENOENT:
                raise

def clean_files(dir : str):
    #Remove directory

    shutil.rmtree(dir)

def parse_aleph_annotate(file_path: str) -> Dict:
    #Read the output of aleph, return a dictionary containing:
    #{"ah": 59, "bs": 8, "number_total_residues": 1350}

    secondary_structure_dict = {}
    with open(file_path) as f_in:
        data = json.load(f_in)
    secondary_structure_dict = copy.deepcopy(data['annotation']['secondary_structure_content'])
    return secondary_structure_dict

def parse_aleph_ss(file_path: str) -> Dict:
    #Parse the aleph.txt file and get all the domains by chain

    chain_res_dict = {}
    with open(file_path) as f_in:
        lines = f_in.readlines()
        for line in lines:
            splitted_text = line.split()
            if len(splitted_text) == 12:
                chain = splitted_text[3]
                residues = expand_residues(f'{splitted_text[4]}-{splitted_text[6]}')
                try:
                    chain_res_dict[chain].append(residues)
                except KeyError:
                    chain_res_dict[chain] = [residues]
    return chain_res_dict

def sort_by_digit(container: Any, item: int=0):
    #Sort list or dictionary by a digit instead of str.
    #Dict can be like this:

    if isinstance(container, dict):
        return sorted(container.items(), key=lambda x: int("".join([i for i in x[item] if i.isdigit()])))
    elif isinstance(container, list):
        return sorted(container, key=lambda x: int("".join([i for i in x if i.isdigit()])))

def create_dir(dir_path: str, delete_if_exists: bool = False):
    #If directory not exists, create it
    #If directory exists, delete it and create it
    # 
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    elif delete_if_exists:
        shutil.rmtree(dir_path)
        os.makedirs(dir_path)

def create_logger():
    #Create logger: The information will be stored in a buffer instead of a file. The buffer can be dumped to
    #a file later.

    logger = logging.getLogger()
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
    #Create logger in a working directory with a specific name:
    
    logger = logging.getLogger()
    logger_data = logger.handlers[0].stream.getvalue()
    logger.removeHandler(logger.handlers[0])
    with open(log_path, 'w+') as f_handle:
        f_handle.write(logger_data)
    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(logging.DEBUG)
    logger.addHandler(file_handler)


