import os
import shutil
import logging
import io
import sys
from Bio import SeqIO
from Bio.PDB import MMCIFIO, PDBIO, PDBList, PDBParser


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
    logger.info('\n')
    logger.info(box)
    logger.info('\n')

def get_mandatory_value(input_load: str, value: str) -> str:
    
    read_value = input_load.get(value)
    if read_value is None:
        raise Exception(f'{value} is mandatory')
    return read_value

def download_pdb(pdb_id: str, output_dir: str):

    logging.info(f'Downloading PDB {pdb_id}')

    pdbl = PDBList()
    result_ent = pdbl.retrieve_pdb_file(pdb_id, pdir=f'{output_dir}', file_format='pdb', obsolete=False)
    if not os.path.exists(result_ent):
        raise Exception(f'{pdb_id} could not be downloaded.')
    shutil.copy2(result_ent, f'{output_dir}/{pdb_id}.pdb')
    shutil.rmtree('obsolete')

def pdb2mmcif(pdb_in_path, cif_out_path):

    with open('/tmp/maxit.sh','w') as f:
        f.write(f'export RCSBROOT="/opt/maxit"\n')
        f.write(f'/opt/maxit/bin/maxit -input {pdb_in_path} -output {cif_out_path} -o 1\n')
    os.system('chmod +x /tmp/maxit.sh')
    os.system('bash /tmp/maxit.sh')


def pdb2cif(pdb_in_path: str, cif_out_path: str):
    p = PDBParser()
    structure = p.get_structure("pdb", pdb_in_path)
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(cif_out_path)

def cif2pdb(cif_in_path, pdb_out_path):

    p = MMCIFParser(QUIET=True)
    struc = p.get_structure('', f'{cif_in_path}')
    io = PDBIO()
    io.set_structure(struc)
    io.save(f'{pdb_out_path}')


def get_path_name(path: str) -> str:

    return os.path.splitext(os.path.basename(path))[0]


def extract_sequence(fasta_path: str) -> str:

    logging.info(f'Extracting sequence from {fasta_path}')

    try:
        record = SeqIO.read(fasta_path, "fasta")
    except:
        raise Exception(f'Not possible to extract the sequence from {fasta_path}')

    return str(record.seq)

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

