import logging
import os
import re
import shutil
import subprocess
import sys
from typing import List, Tuple
from Bio import SeqIO
from Bio.PDB import MMCIFIO, PDBIO, PDBList, PDBParser
from libs import utils


def download_pdb(pdb_id: str, output_dir: str):

    pdbl = PDBList()
    result_ent = pdbl.retrieve_pdb_file(pdb_id, pdir=f'{output_dir}', file_format='pdb', obsolete=False)
    if not os.path.exists(result_ent):
        raise Exception(f'{pdb_id} could not be downloaded.')
    shutil.copy2(result_ent, f'{output_dir}/{pdb_id}.pdb')
    shutil.rmtree('obsolete')

def pdb2mmcif(output_dir: str, pdb_in_path: str, cif_out_path: str):
    maxit_dir = f'{output_dir}/maxit'
    if not os.path.exists(maxit_dir):
        os.mkdir(maxit_dir)
        
    subprocess.Popen(['maxit', '-input', pdb_in_path, '-output', cif_out_path, '-o', '1'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    shutil.rmtree(maxit_dir)

def pdb2cif(pdb_id: str, pdb_in_path: str, cif_out_path: str):

    p = PDBParser()
    structure = p.get_structure(pdb_id, pdb_in_path)
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(cif_out_path)

def cif2pdb(cif_in_path: str, pdb_out_path: str):
    
    p = MMCIFParser(QUIET=True)
    struc = p.get_structure('', f'{cif_in_path}')
    io = PDBIO()
    io.set_structure(struc)
    io.save(f'{pdb_out_path}')

def extract_sequence(fasta_path: str) -> str:

    logging.info(f'Extracting sequence from {fasta_path}')

    try:
        record = SeqIO.read(fasta_path, "fasta")
    except:
        raise Exception(f'Not possible to extract the sequence from {fasta_path}')

    return str(record.seq)

def merge_pdbs(list_of_paths_of_pdbs_to_merge: str, merged_pdb_path: str):

    with open(merged_pdb_path, 'w') as f:
        counter = 0
        for pdb_path in list_of_paths_of_pdbs_to_merge:
            for line in open(pdb_path,'r').readlines():
                if line[:4] == 'ATOM':
                    counter += 1
                    f.write(line[:4] + str(counter).rjust(7) + line[11:])

def run_pisa(pdb_path: str):

    print(f'Generating REMARK 350 for {pdb_path} with PISA.')
    subprocess.Popen(['pisa', 'temp', '-analyse', pdb_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    pisa_output = subprocess.Popen(['pisa', 'temp', '-350'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0]
    return pisa_output.decode('utf-8')

def read_remark_350(pdb_path: str) -> Tuple[ List[str], List[float] ]:

    pdb_text = open(pdb_path, 'r').read()
    match_biomolecules = [m.start() for m in re.finditer(r'REMARK 350 BIOMOLECULE:', pdb_text)] # to know how many biomolecules there are.
    if len(match_biomolecules) == 0:
        pdb_text = run_pisa(pdb_path)
        match_biomolecules = [m.start() for m in re.finditer(r'REMARK 350 BIOMOLECULE:', pdb_text)] # to know how many biomolecules there are.

    if len(match_biomolecules) == 0:
        raise Exception(f'REMARK not found for template {pdb_path}.')
    elif len(match_biomolecules) == 1:
        match_last_350 = [m.start() for m in re.finditer(r'REMARK 350', pdb_text)][-1]
        match_end_in_last_350 = [m.end() for m in re.finditer(r'\n', pdb_text[match_last_350:])][-1]
        remark_350_text = pdb_text[match_biomolecules[0]:(match_last_350+match_end_in_last_350)]
    else:
        print('(It seem there is more than one biological assembly from REMARK 350. Only'
              ' "BIOMOLECULE 1" will be considered for the assembly generation)')
        remark_350_text = pdb_text[match_biomolecules[0]:match_biomolecules[1]-1]

    match_biomt1 = [m.start() for m in re.finditer(r'REMARK 350   BIOMT1', remark_350_text)]
    match_biomt3 = [m.end() for m in re.finditer(r'REMARK 350   BIOMT3', remark_350_text)]

    end_remark_350_block = [m.start() for m in re.finditer('\n', remark_350_text[match_biomt3[-1]:])]

    transformation_blocks_indices = match_biomt1 + [match_biomt3[-1] + end_remark_350_block[0] + 1]

    transformations_list = []
    for index in range(len(transformation_blocks_indices) - 1):
        block = remark_350_text[transformation_blocks_indices[index]:transformation_blocks_indices[index+1]]
        matrix = [item.split()[4:8] for item in block.split('\n')[:-1]]
        r11, r12, r13, t1 = matrix[0]
        r21, r22, r23, t2 = matrix[1]
        r31, r32, r33, t3 = matrix[2]
        transformation = [[r11, r12, r13], [r21, r22, r23], [r31, r32, r33], [t1, t2, t3]]
        transformations_list.append(transformation)

    chain_list = [line.split(':')[-1].replace(" ", "").split(',') for line in remark_350_text.split('\n')
                  if 'REMARK 350 APPLY THE FOLLOWING TO CHAINS:' in line][0]

    return chain_list, transformations_list

def rot_and_trans(pdb_path: str, out_pdb_path: str, rot_tra_matrix: str):

    r11, r12, r13 = rot_tra_matrix[0]
    r21, r22, r23 = rot_tra_matrix[1]
    r31, r32, r33 = rot_tra_matrix[2]
    t1, t2, t3 = rot_tra_matrix[3]

    with open('/tmp/pdbset.sh', 'w') as f:
        f.write(f'pdbset xyzin {pdb_path} xyzout {out_pdb_path} << eof\n')
        f.write(
            f'rotate {float(r11)} {float(r12)} {float(r13)} {float(r21)} {float(r22)} {float(r23)} {float(r31)} {float(r32)} {float(r33)}\n')
        f.write(f'shift {float(t1)} {float(t2)} {float(t3)}\n')
        f.write('end\n')
        f.write('eof')
    subprocess.Popen(['bash', '/tmp/pdbset.sh'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
    utils.rmsilent(f'/tmp/pdbset.sh')

def change_chain_and_apply_offset_in_single_chain(pdb_in_path, pdb_out_path, offset=None, chain=None):

    with open('/tmp/pdbset.sh', 'w') as f:
        f.write(f'pdbset xyzin {pdb_in_path} xyzout {pdb_out_path} << eof\n')
        if offset:
            f.write(f'renumber increment {offset}\n')
        if chain:
            f.write(f'chain {chain}\n')
        f.write('end\n')
        f.write('eof')
    subprocess.Popen(['bash', '/tmp/pdbset.sh'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
    utils.rmsilent(f'/tmp/pdbset.sh')