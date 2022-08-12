import errno
import glob
import os
import logging
import io
import sys

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
    
    read_value = input_load.get(value)
    if read_value is None:
        raise Exception(f'{value} is mandatory')
    return read_value


def get_path_name(path: str) -> str:

    return os.path.splitext(os.path.basename(path))[0]

def rmsilent(file_path: str):
    for file in glob.glob(file_path):
        try:
            os.remove(file)
        except OSError as e:
            if e.errno != errno.ENOENT:
                raise

def clean_files(output_dir : str):
    rmsilent(f'{output_dir}/*.ffindex')
    rmsilent(f'{output_dir}/*.ffdata')
    rmsilent(f'{output_dir}/*.custom_output.hhr')
    rmsilent(f'{output_dir}/*.cif')
    rmsilent(f'{output_dir}/[A-Z].pdb')
    rmsilent(f'{output_dir}/*[A-Z]_template.pdb')
    rmsilent(f'{output_dir}/[0-9].pdb')

def create_logger():
    """
    Create logger: In case there is no working directory, the information
    will be stored in a buffer instead of a file. The buffer can be dumped to
    a file later.
    """
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

