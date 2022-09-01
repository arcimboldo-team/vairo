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

def get_mandatory_value(input_load: str, value: str) -> str:
    #Read value, Raise exception if value could not be found

    read_value = input_load.get(value)
    if read_value is None:
        raise Exception(f'{value} is mandatory')
    return read_value

def get_file_name(path: str) -> str:
    #Get the name of the file, without path or extension

    return os.path.splitext(os.path.basename(path))[0]

def get_main_path() -> str:
    #Get the path of the main.py

    return Path(__file__).parent.absolute()

def get_working_dir() -> str:
    # Get working directory

    return os.getcwd()

def get_key_for_value(value: str, search_dict: Dict) -> List:
    #Given a value, get the list of all keys that contains that value

    return list(search_dict.keys())[list(search_dict.values()).index(value)]

def get_chain_and_number(path_pdb: str) -> List:
    #Given a path: ../../template_A1.pdb return A and 1
    name = get_file_name(path_pdb)
    code = name.split("_", 1)[-1]
    return code[0], int(code[1:])

def replace_last_number(text: str, value: str) -> str:
    #Replace the last number of text by the value

    return re.sub(r'\d+.pdb', str(value), str(text))+'.pdb'


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
    #Create logger

    logstream = io.StringIO()
    formatter = logging.Formatter('%(message)s')

    logger = logging.getLogger()
    streamHandler = logging.StreamHandler(logstream)
    streamHandler.setLevel(logging.DEBUG)
    streamHandler.setFormatter(formatter)
    logger.addHandler(streamHandler)

    stdoutHandler = logging.StreamHandler(sys.stdout)
    stdoutHandler.setLevel(logging.DEBUG)
    stdoutHandler.setFormatter(formatter)
    logger.addHandler(stdoutHandler)

    logger.setLevel(logging.DEBUG)

