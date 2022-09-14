import logging
import os
import re
import shutil
import subprocess
from typing import Dict, List, Tuple
from Bio import SeqIO
from Bio.PDB import MMCIFIO, PDBIO, PDBList, PDBParser, Residue, Chain, Select, Selection, Structure
from libs import utils
from libs.alphafold_paths import AlphaFoldPaths
from scipy.spatial import distance
import numpy as np
import itertools

def download_pdb(pdb_id: str, output_dir: str):

    pdbl = PDBList()
    result_ent = pdbl.retrieve_pdb_file(pdb_id, pdir=f'{output_dir}', file_format='pdb', obsolete=False)
    if not os.path.exists(result_ent):
        raise Exception(f'{pdb_id} could not be downloaded.')
    shutil.copy2(result_ent, f'{output_dir}/{pdb_id}.pdb')
    os.remove(result_ent)
    shutil.rmtree('obsolete')

def pdb2mmcif(output_dir: str, pdb_in_path: str, cif_out_path: str):
    
    maxit_dir = f'{output_dir}/maxit'
    if not os.path.exists(maxit_dir):
        os.mkdir(maxit_dir)
    subprocess.Popen(['maxit', '-input', pdb_in_path, '-output', cif_out_path, '-o', '1'], cwd=maxit_dir,stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    shutil.rmtree(maxit_dir)

def pdb2cif(pdb_id: str, pdb_in_path: str, cif_out_path: str):

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_in_path)
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(cif_out_path)

def cif2pdb(cif_in_path: str, pdb_out_path: str):
    
    parser = MMCIFParser(QUIET=True)
    struc = parser.get_structure('', f'{cif_in_path}')
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

def run_pisa(pdb_path: str) -> str:

    logging.info(f'Generating REMARK 350 for {pdb_path} with PISA.')
    subprocess.Popen(['pisa', 'temp', '-analyse', pdb_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    pisa_output = subprocess.Popen(['pisa', 'temp', '-350'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0]
    return pisa_output.decode('utf-8')

def run_af2(output_dir:str, alphafold_paths:AlphaFoldPaths):
    
    logging.info('Running AF2')
    alphafold_paths.create_af2_script(output_dir=output_dir)
    af2_output = subprocess.Popen(['bash', alphafold_paths.run_alphafold_bash], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = af2_output.communicate()
    with open(alphafold_paths.run_alphafold_log, 'w') as f:
        f.write(stdout.decode('utf-8'))

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
        logging.info('It seem there is more than one biological assembly from REMARK 350. Only'
              ' "BIOMOLECULE 1" will be considered for the assembly generation')
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

    chain_list = [line.split(':')[-1].replace(' ', '').split(',') for line in remark_350_text.split('\n')
                  if 'REMARK 350 APPLY THE FOLLOWING TO CHAINS:' in line][0]

    return chain_list, transformations_list

def change_chain(pdb_in_path: str, pdb_out_path: str, rot_tra_matrix: str=None, offset: int=0, chain: str=None):

    with open('/tmp/pdbset.sh', 'w') as f:
        f.write(f'pdbset xyzin {pdb_in_path} xyzout {pdb_out_path} << eof\n')
        if rot_tra_matrix is not None:
            r11, r12, r13 = rot_tra_matrix[0]
            r21, r22, r23 = rot_tra_matrix[1]
            r31, r32, r33 = rot_tra_matrix[2]
            t1, t2, t3 = rot_tra_matrix[3]
            f.write(
                f'rotate {float(r11)} {float(r12)} {float(r13)} {float(r21)} {float(r22)} {float(r23)} {float(r31)} {float(r32)} {float(r33)}\n')
            f.write(f'shift {float(t1)} {float(t2)} {float(t3)}\n')
        f.write(f'renumber increment {offset}\n')
        if chain:
            f.write(f'chain {chain}\n')
        f.write('end\n')
        f.write('eof')
    subprocess.Popen(['bash', '/tmp/pdbset.sh'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
    utils.rmsilent(f'/tmp/pdbset.sh')

def remove_hydrogens(pdb_in_path: str, pdb_out_path:str):

    counter = 0
    with open(pdb_out_path, 'w') as f_out:
        with open(pdb_in_path, 'r') as f_in:
            for line in f_in.readlines():
                if line.split()[-1] in ['N', 'C', 'O', 'S']:
                    counter = counter + 1
                    f_out.write(line[:6] + str(counter).rjust(5) + line[11:])

def get_resseq(residue: Residue) -> int:
    #Return resseq number

    return residue.get_full_id()[3][1]

def get_chains(structure: Structure) -> List:
    #Return all chains from a PDB structure

    return [chain.get_id() for chain in structure.get_chains()]

def get_structure(pdb_path: str) -> Structure:
    #Get PDB structure

    pdb_id = utils.get_file_name(pdb_path)
    parser = PDBParser(QUIET=True)
    return parser.get_structure(pdb_id, pdb_path)

def extract_cryst_card_pdb(pdb_in_path: str) -> str:
    #Extract the crystal card from a pdb

    if os.path.isfile(pdb_in_path):
        with open(pdb_in_path, 'r') as f_in:
            pdb_lines = f_in.readlines()
        for line in pdb_lines:
            if line.startswith("CRYST1"):
                cryst_card = line
                return cryst_card
    return None


def get_atom_line(remark:str, num:int, name:str, res:int, chain:str, resseq, x:float, y:float, z:float,
                    occ:str, bfact:str, atype:str) -> str:
    #Given all elements of an atom, parse them in PDB format

    result = f'{remark:<6}{num:>5}  {name:<3}{res:>4} {chain}{resseq:>4}    {float(x):8.3f}{float(y):8.3f}{float(z):8.3f}{float(occ):6.2f}{float(bfact):6.2f}{atype:>12}\n'

    return result

def parse_pdb_line(line: str) -> Dict:
    #Parse all elements of an atom of a PDB line
    
    parsed_dict = {
            'remark':line[:6],
            'num': line[6:11],
            'name': line[12:16],
            'resname': line[17:20],
            'chain': line[21],
            'resseq': line[22:26],
            'x': line[30:38],
            'y': line[38:46],
            'z': line[46:54],
            'occ': line[54:60],
            'bfact': line[60:66],
            'atype': line[76:78]
    }
    for key, value in parsed_dict.items():
        parsed_dict[key] = value.replace(' ', '')

    return parsed_dict

def split_chains_assembly(pdb_in_path: str, pdb_out_path:str, a_air) -> bool:

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(utils.get_file_name(pdb_in_path), pdb_in_path)
    delete_residues = []
    total_length = len(a_air.query_sequence) + a_air.glycines
    chains = [chain.get_id() for chain in structure.get_chains()]
    chain = None
    chain_name = chains[0]

    if len(chains) > 1:
        logging.info(f'PDB: {pdb_in_path} is already splitted in several chains: {chains}')
    else:
        old_mod = -1
        for res in structure[0][chains[0]]:
            resseq = get_resseq(res)
            res_mod = resseq % total_length

            if res_mod < old_mod:
                chain_name = chr(ord(chain_name)+1)
                chain = Chain.Chain(chain_name)
                structure[0].add(chain)
 
            if res_mod == 0 or res_mod > len(a_air.query_sequence):
                delete_residues.append(res.id)
            elif resseq > total_length and chain is not None:
                delete_residues.append(res.id)
                res.parent = chain
                chain.add(res)
                chain[res.id].id = (' ', res_mod, ' ')
            old_mod = res_mod

        for id in delete_residues:
            structure[0][chains[0]].detach_child(id)

    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_out_path, preserve_atom_numbering = True)

def chain_splitter(pdb_in_path: str, pdb_out_path: str, chain: str):
    #Given a pdb_in and a chain, write a pdb_out containing only
    #the specified chain

    class ChainSelect(Select):
        def __init__(self, chain):
            self.chain = chain

        def accept_chain(self, chain):
            if chain.get_id() == self.chain:
                return 1
            else:          
                return 0

    structure = get_structure(pdb_path=pdb_in_path)
    io = PDBIO()               
    io.set_structure(structure)
    io.save(pdb_out_path, ChainSelect(chain))

def superpose_pdbs(query_pdb: str, target_pdb: str, output_pdb = None):

    superpose_input_list = ['superpose', f'{query_pdb}', '-s', '-all', f'{target_pdb}', '-s', '-all']
    if output_pdb is not None:
        superpose_input_list.extend(['-o', output_pdb])
    
    superpose_output = subprocess.Popen(superpose_input_list, stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
        
    for line in superpose_output.split('\n'):
        if 'r.m.s.d:' in line:
            rmsd = float(line.split()[1])
        if 'quality Q:' in line:
            quality_q = line.split()[2]
        if 'Nalign:' in line:
            nalign = line.split()[1]

    return rmsd, nalign, quality_q

def pdist(query_pdb: str, target_pdb: str) -> List[List]:

    structure_query = PDBParser(QUIET=1).get_structure('query', query_pdb)
    res_query_list = [res.id[1] for res in Selection.unfold_entities(structure_query, 'R')]

    structure_target = PDBParser(QUIET=1).get_structure('target', target_pdb)
    res_target_list = [res.id[1] for res in Selection.unfold_entities(structure_target, 'R')]

    common_res_list = list(set(res_query_list) & set(res_target_list))

    query_common_list = [res for res in Selection.unfold_entities(structure_query, 'R') if res.id[1] in common_res_list]
    query_matrix = calculate_distance_pdist(res_list=query_common_list)

    target_common_list = [res for res in Selection.unfold_entities(structure_target, 'R') if res.id[1] in common_res_list]
    target_matrix = calculate_distance_pdist(res_list=target_common_list)

    diff_pdist_matrix = np.abs(query_matrix - target_matrix)

    return diff_pdist_matrix.mean()

def calculate_distance_pdist(res_list: List):

    coords = [res['CA'].coord for res in res_list]
    pdist = distance.pdist(coords, "euclidean")
    return distance.squareform(pdist)    

def find_interface_from_pisa(pdb_in_path: str) -> Dict:

    interfaces_dict = {}

    subprocess.Popen(['pisa', 'temp', '-analyse', pdb_in_path],
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE).communicate()

    pisa_output = subprocess.Popen(['pisa', 'temp', '-list', 'interfaces'], stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE).communicate()[0].decode('utf-8')

    if 'NO INTERFACES FOUND' in pisa_output:
        return interfaces_dict

    match1 = [m.start() for m in re.finditer(' LIST OF INTERFACES', pisa_output)][0]
    match2 = [m.start() for m in re.finditer(' ##:  serial number', pisa_output)][0]

    for line in pisa_output[match1:match2].split('\n')[4:-2]:
        area = line.split('|')[3][:8].replace(' ', '')
        deltaG = line.split('|')[3][8:15].replace(' ', '')
        chain1 = line.split('|')[1].replace(' ', '')
        chain2 = line.split('|')[2].split()[0].replace(' ', '')
        interfaces_dict[f'{chain1}{chain2}'] =  (area, deltaG)

    return interfaces_dict

def calculate_auto_offset(input_list: List[List]) -> List:
    
    combinated_list = list(itertools.product(*input_list))
    trimmed_list = []
    for element in combinated_list:
        target_list = [target for _,target,_ in element]
        if len(target_list) == len(set(target_list)):
            trimmed_list.append(element)
    score_list = []
    for element in trimmed_list:
        score_list.append(sum(z for _,_,z in element))
    max_value = max(score_list)
    max_index = score_list.index(max_value)

    return trimmed_list[max_index]

def split_dimers_in_pdb(pdb_in_path: str, pdb_out_path: str, chain1: List, chain2: List):

    with open(pdb_in_path, 'r') as f_in:
        input_lines = f_in.readlines()

    with open(pdb_out_path, 'w') as f_out:
        for line in input_lines:
            if line[21:22] in [chain1, chain2]:
                f_out.write(line)