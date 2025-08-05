import base64, copy, errno, glob, io, json, logging, os, re, shutil, sys
import math
import random
import string
import subprocess
from itertools import groupby
from operator import itemgetter
from pathlib import Path
from typing import Any, Dict, List, Tuple
from sklearn import preprocessing
from libs import structures
from libs.global_variables import INPUT_PARAMETERS
from alphafold.common import residue_constants


def scale_values(input_list: List[int]) -> List[int]:
    max_value = max(input_list)
    new_list = []
    for value in input_list:
        if value <= 0:
            new_value = 0
        elif value >= max_value:
            new_value = 1
        else:
            new_value = round(math.log(value + 1) / math.log(max_value + 1), 2)
        new_list.append(new_value)
    return new_list


def check_external_programs():
    try:
        cmd = 'which gesamt'
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    except:
        raise Exception('The CCP4 programs cannot be found in the PATH. Please check if LSQKAB, PISA and PDBSET are '
                        'present in the PATH. If the CCP4 suite is not installed, please install it from the '
                        'following link: https://www.ccp4.ac.uk/')

    try:
        cmd = 'which maxit'
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    except:
        raise Exception('MAXIT cannot be found in the PATH. In order to continue, download and install maxit ('
                        'https://sw-tools.rcsb.org/apps/MAXIT/index.html).')


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
    logging.error('\n')
    logging.error(box)
    logging.error('\n')


def print_matrix(matrix: List):
    print('\n'.join('\t'.join(map(str, row)) for row in matrix))


def normalize_list(input_list: List):
    # Normalize list of values
    return preprocessing.normalize(input_list)[0]


def get_file_extension(path: str) -> str:
    # Get the extension of the file
    _, extension = os.path.splitext(path)
    return extension


def get_file_name(path: str) -> str:
    # Get the name of the file, without path or extension
    return os.path.splitext(os.path.basename(path))[0]


def get_readme() -> str:
    # Get README.md file
    return os.path.join(os.path.dirname(get_parent_folder(str(Path(__file__)))), 'helper.txt')


def get_main_path() -> Path:
    # Get the path of the main.py
    return Path(__file__).parent.parent.absolute()


def get_parent_folder(dir_path: str) -> Path:
    return Path(dir_path).parent.absolute()


def get_working_dir() -> str:
    # Get working directory
    return os.getcwd()


def dict_values_to_list(input_dict: Dict):
    # Given a Dict, return all the values from the dict in a list
    return [value for value in input_dict.values()]


def get_key_by_value(value: str, search_dict: Dict) -> List[str]:
    # Given a value, get the list of all keys that contains that value
    matching_keys = []
    for key, element in search_dict.items():
        if isinstance(element, str):
            if element == value:
                matching_keys.append(key)
        elif isinstance(element, int):
            if element == value:
                matching_keys.append(key)
        else:
            if value in element:
                matching_keys.append(key)

    return matching_keys


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


def print_consecutive_numbers(number_list: List[int]) -> str:
    # Given an integer list, return ranges of consecutive numbers
    if not number_list:
        return ""

    nums = sorted(set(number_list))  # Sort and remove duplicates
    ranges = []
    start = nums[0]

    for i in range(1, len(nums)):
        if nums[i] != nums[i - 1] + 1:
            ranges.append(f"{start}" if start == nums[i - 1] else f"{start}-{nums[i - 1]}")
            start = nums[i]

    ranges.append(f"{start}" if start == nums[-1] else f"{start}-{nums[-1]}")
    return ", ".join(ranges)


def get_chain_and_number(path_pdb: str) -> Tuple[str, int]:
    # Given a path: ../../template_A1.pdb return A and 1
    # Return CHAIN and NUMBER
    name = get_file_name(path_pdb)
    code = name.split('_')[-1]
    return code[0], int(code[1:])


def replace_last_number(text: str, value: int) -> str:
    # Replace the last amount of a text by the value
    return re.sub(r'\d+.pdb', str(value), str(text)) + '.pdb'


def expand_residues(res: str) -> List:
    # Expand a str formatted like this: 10-12, 32, 34
    # To a list: [10,11,12,32,34]
    if not res or res == '':
        return []
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


def expand_partition(res: str) -> List:
    # Expand a str formatted like this: 10-12, 32, 34
    # To a list of pairs [10-12, 32-32, 34-34]
    modified_list = str(res).replace(' ', '').split(',')
    return_list = []
    for partition in modified_list:
        partition_list = list(map(int, partition.split('-')))
        if len(partition_list) == 1:
            partition_list = [partition_list[0], partition_list[0]]
        return_list.append(partition_list)
    return return_list


def renum_residues(res_list: List[int], mapping: Dict) -> List[int]:
    return [mapping[res] for res in res_list]


def rmsilent(file_path: str):
    # Remove file without an error if it doesn't exist
    for file in glob.glob(file_path):
        try:
            os.remove(file)
        except OSError as e:
            if e.errno != errno.ENOENT:
                raise


def clean_files(input_dir: str):
    # Remove directory
    shutil.rmtree(input_dir)


def parse_cc_analysis(file_path: str) -> Dict:
    # Read the output of the cc analysis file
    return_dict = {}
    with open(file_path) as f_in:
        lines = f_in.readlines()
        for line in lines:
            split_text = line.split()
            coordinates = []
            if len(split_text) > 3:
                for pos in range(1, len(split_text) - 2):
                    coordinates.append(float(split_text[pos]))
                return_dict[split_text[0]] = structures.CCAnalysisOutput(coordinates, float(split_text[-2]),
                                                                         float(split_text[-1]))
            else:
                return_dict[split_text[0]] = structures.CCAnalysisOutput([split_text[1], 0], None, float(split_text[2]))
    return return_dict


def parse_hinges_chains(output: str) -> str:
    # Read the output of hines when using -p chains
    # Return the best chain combination
    lines = output.strip().split('\n')
    lowest_rmstot = float('inf')
    lowest_perm = ''
    found_hash = False
    for line in lines:
        if found_hash:
            split_line = line.split()
            rmstot = float(split_line[3])
            if rmstot < lowest_rmstot:
                lowest_rmstot = rmstot
                lowest_perm = split_line[1]
        elif line.startswith(' #'):
            found_hash = True

    return lowest_perm


def parse_hinges(output: str) -> structures.Hinges:
    # Read the output of hinges
    # Parse all the output of hinges with the following structure:
    # decreasing_rmsd: float
    # min_rmsd: str
    # groups: List[(int, int)]
    rmsd_list = []
    residues_list = []
    # split the data by ngroups
    groups = re.split(r'ngroups,\s*Ntot,\s*SUMSQ,\s*RMStot', output)[1:]
    # loop over each group and extract the data
    for i, group in enumerate(groups, start=1):
        ngroup_match = re.search(r"=\s*\d+\s+\d+\s+[\d.]+\s+([\d.]+)", group)
        if ngroup_match:
            rmsd_list.append(float(ngroup_match.group(1)))
            ngroup_data = []
            pairs = re.findall(r'A(\d+)\s+A(\d+)', group)
            for pair in pairs:
                ngroup_data.append((int(pair[0]), int(pair[1])))
            residues_list.append(ngroup_data)
    if not rmsd_list:
        return None

    # Extract counts from file 1 and file 2
    file_1_count_match = re.search(r'number of CA in file 1:\s+(\d+)', output)
    file_1_count = int(file_1_count_match.group(1)) if file_1_count_match else 0

    file_2_count_match = re.search(r'number of CA in file 2:\s+(\d+)', output)
    file_2_count = int(file_2_count_match.group(1)) if file_2_count_match else 0

    # Extract warnings for file 1 and file 2
    warnings_1_match = re.search(r'WARNING: (\d+) non-matching residue\(s\) in 1st sequence ignored', output)
    warnings_1 = int(warnings_1_match.group(1)) if warnings_1_match else 0

    warnings_2_match = re.search(r'WARNING: (\d+) non-matching residue\(s\) in 2nd sequence ignored', output)
    warnings_2 = int(warnings_2_match.group(1)) if warnings_2_match else 0

    # Calculate file ratios
    file1 = (file_1_count - warnings_1) / file_1_count if file_1_count != 0 else 0
    file2 = (file_2_count - warnings_2) / file_2_count if file_2_count != 0 else 0

    hinges_result = structures.Hinges(
        decreasing_rmsd_total=(rmsd_list[0] - rmsd_list[-1]) / rmsd_list[0] * 100 if rmsd_list[0] > 0 else 0,
        decreasing_rmsd_middle=(rmsd_list[0] - rmsd_list[len(rmsd_list) // 2]) / rmsd_list[0] * 100 if rmsd_list[
                                                                                                           0] > 0 else 0,
        one_rmsd=rmsd_list[0],
        middle_rmsd=rmsd_list[len(rmsd_list) // 2],
        min_rmsd=min(rmsd_list),
        overlap=file1 if file1 > file2 else file2,
        groups=residues_list)
    return hinges_result


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
    # interface contains the required information.

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
        nhb = line[3][15:22].replace(' ', '')
        return_list.append(
            {'serial': serial, 'area': area, 'deltaG': deltag, 'nhb': nhb, 'chain1': chain1, 'chain2': chain2})

    return return_list


def parse_pisa_interfaces(pisa_output: str) -> Dict:
    # It parses the pisa output, in the following format:
    # List[Dict{
    #   solvation1
    #   solvation2
    #   se_gain1
    #   se_gain2
    #   chain1
    #   chain2
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


def read_rankeds(input_path: str) -> List[str]:
    ranked_paths = [path for path in os.listdir(input_path) if check_ranked(path)]
    return [structures.Ranked(os.path.join(input_path, path)) for path in sort_by_digit(ranked_paths)]


def check_ranked(input_path: str) -> bool:
    return re.match('ranked_[0-9]+.pdb', input_path) or re.match('cluster_[0-9]+_ranked_[0-9]+.pdb', input_path)


def delete_old_rankeds(input_path: str):
    [os.remove(os.path.join(input_path, path)) for path in os.listdir(input_path) if check_ranked(path)]


def delete_old_html(input_path: str):
    [os.remove(os.path.join(input_path, path)) for path in os.listdir(input_path) if
     get_file_extension(path) == '.html']

def check_input(global_dict: Dict):
    all_keys = []
    [all_keys.extend(list(value.keys())) for key, value in INPUT_PARAMETERS.items()]

    def check_keys(data: Dict) -> str:
        if isinstance(data, dict):
            for key in data.keys():
                if key not in all_keys and key.lower() != 'all' and len(key) != 3:
                    raise Exception(f'Parameter {key} does not exist. Check the input file')
                if data[key] is None:
                    raise Exception(f'Paramter {key} does not have any value. Comment it or add input')
                if isinstance(data[key], list):
                    check_keys(data[key])

        if isinstance(data, list):
            for aux_dict in data:
                check_keys(aux_dict)

    check_keys(global_dict)


def get_input_value(name: str, section: str, input_dict: Dict, override_default=None):
    mapping = {
        'global': 'global_input',
        'sequence': 'sequence_input',
        'template': 'template_input',
        'modifications': 'modifications_input',
        'features': 'features_input',
        'mutations': 'mutations_input',
        'append_library': 'append_library_input',
    }
    chosen_dict = INPUT_PARAMETERS.get(mapping.get(section))
    value_dict = chosen_dict.get(name)
    value = input_dict.get(name)
    if value is None and value_dict['required']:
        raise Exception(f'{name} does not exist and it is a mandatory input parameter. Check the input file.')
    elif value is None:
        if override_default is not None:
            value = override_default
        else:
            value = value_dict['default']
    return value


def modification_list(query: List[int], target: List[int], length: int) -> List[int]:
    # Create a list of length LENGTH. Where each value, is the value that it should has in the target
    if query is None:
        query = '1'
    query = list(map(int, str(query).replace(' ', '').split(',')))
    if target is None:
        target = [(1, length)]
    else:
        target = target.replace(' ', '').split(',')
        target = [tuple(map(int, r.split('-'))) for r in target]
    if len(query) != len(target):
        raise ValueError('The number of query positions and library positions mismatch')
    return generate_modification_list(query=query, target=target, length=length)


def generate_modification_list(query: List[int], target: List[int], length: int) -> List[int]:
    result = [None] * length
    for query_value, target_range in zip(query, target):
        start, end = target_range
        for i in range(start, end + 1):
            if query_value - 1 < length:
                result[query_value - 1] = i
                query_value += 1
    return result


def print_dict(input_dict: Dict):
    for key, value in input_dict.items():
        if isinstance(value, list):
            logging.error(f'{key}: {" ".join(value)}')
        else:
            logging.error(f'{key}: {value}')


def create_logger():
    # Create logger: The information will be stored in a buffer instead of a file. The buffer can be dumped to
    # a file later.

    logger = logging.getLogger()
    logger.handlers.clear()

    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

    logger.setLevel(logging.DEBUG)
    test = io.StringIO()
    stream_handler_ = logging.StreamHandler(test)
    stream_handler_.setLevel(logging.NOTSET)
    logger.addHandler(stream_handler_)

    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.ERROR)
    logger.addHandler(stdout_handler)


def create_logger_dir(log_path: str, log_extended_path: str):
    # Create logger in a working directory with a specific name:
    logger = logging.getLogger()
    logger_data = logger.handlers[0].stream.getvalue()
    logger.removeHandler(logger.handlers[0])

    with open(log_path, 'w+') as f_handle:
        f_handle.write(logger_data)

    with open(log_extended_path, 'w+') as f_handle:
        f_handle.write(logger_data)

    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(logging.ERROR)
    logger.addHandler(file_handler)

    file_handler2 = logging.FileHandler(log_extended_path)
    file_handler2.setLevel(logging.DEBUG)
    logger.addHandler(file_handler2)


def generate_random_code(length: int):
    letters = string.ascii_letters
    random_code = ''.join(random.choice(letters) for i in range(length))
    return random_code


def read_mutations_dict(input_mutations: list):
    mutations_dict: Dict = {}
    for mutation in input_mutations:
        key = list(mutation.keys())[0]
        values = expand_residues(list(mutation.values())[0])
        if key not in list(residue_constants.restype_3to1.keys()):
            raise Exception(
                f'Mutation residues {"".join(values)} in {key} could not be possible. Residue {key} does not '
                f'exist')
        else:
            key = residue_constants.restype_3to1[key]
        mutations_dict.setdefault(key, []).extend(values)
    return mutations_dict








