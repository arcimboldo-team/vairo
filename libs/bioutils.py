import logging
import os
import re
import shutil
import subprocess
from typing import Dict, List, Tuple
from Bio import SeqIO
from Bio.PDB import MMCIFIO, PDBIO, PDBList, PDBParser, Residue, Chain, Select
from libs import utils
from libs.alphafold_paths import AlphaFoldPaths

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

def create_af2_script(output_dir: str, script_path: str, alphafold_paths: AlphaFoldPaths):

    with open(script_path, 'w') as bash_file:
        previous_path_to_output_dir = '/'.join(output_dir.split('/')[:-1])
        name = output_dir.split('/')[-1]
        bash_file.write('#!/bin/bash\n')
        bash_file.write(f'python {os.path.dirname(os.path.abspath(__file__))}/ALPHAFOLD/run_alphafold.py \\\n')
        bash_file.write(f'--fasta_paths={name}.fasta \\\n')
        bash_file.write(f'--output_dir={previous_path_to_output_dir} \\\n')
        bash_file.write(f'--data_dir={alphafold_paths.af2_dbs_path} \\\n')
        bash_file.write(f'--uniref90_database_path={alphafold_paths.uniref90_db_path} \\\n')
        bash_file.write(f'--mgnify_database_path={alphafold_paths.mgnify_db_path} \\\n')
        bash_file.write(f'--template_mmcif_dir={alphafold_paths.mmcif_db_path} \\\n')
        bash_file.write('--max_template_date=2022-03-09 \\\n')
        bash_file.write(f'--obsolete_pdbs_path={alphafold_paths.obsolete_mmcif_db_path} \\\n')
        bash_file.write('--model_preset=monomer \\\n')
        bash_file.write(f'--bfd_database_path={alphafold_paths.bfd_db_path} \\\n')
        bash_file.write(f'--uniclust30_database_path={alphafold_paths.uniclust30_db_path} \\\n')
        bash_file.write(f'--pdb70_database_path={alphafold_paths.pdb70_db_path} \\\n')
        bash_file.write('--read_features_pkl=True\n')
        bash_file.close()

def run_af2(output_dir:str, alphafold_paths:AlphaFoldPaths):
    
    logging.info('Running AF2')

    script_path = f'{output_dir}/run_af2.sh'
    log_path = f'{output_dir}/af2_output.log'

    create_af2_script(output_dir=output_dir, script_path=script_path, alphafold_paths=alphafold_paths)
    af2_output = subprocess.Popen(['bash', script_path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = af2_output.communicate()
    with open(log_path, 'w') as f:
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

def rot_and_trans(pdb_in_path: str, pdb_out_path: str, rot_tra_matrix: str):

    r11, r12, r13 = rot_tra_matrix[0]
    r21, r22, r23 = rot_tra_matrix[1]
    r31, r32, r33 = rot_tra_matrix[2]
    t1, t2, t3 = rot_tra_matrix[3]

    with open('/tmp/pdbset.sh', 'w') as f:
        f.write(f'pdbset xyzin {pdb_in_path} xyzout {pdb_out_path} << eof\n')
        f.write(
            f'rotate {float(r11)} {float(r12)} {float(r13)} {float(r21)} {float(r22)} {float(r23)} {float(r31)} {float(r32)} {float(r33)}\n')
        f.write(f'shift {float(t1)} {float(t2)} {float(t3)}\n')
        f.write('end\n')
        f.write('eof')
    subprocess.Popen(['bash', '/tmp/pdbset.sh'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
    utils.rmsilent(f'/tmp/pdbset.sh')

def change_chain_and_apply_offset_in_single_chain(pdb_in_path: str , pdb_out_path: str , offset: str=None, chain: str=None):

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

def remove_hydrogens(pdb_in_path: str, pdb_out_path:str):

    counter = 0
    with open(pdb_out_path, 'w') as f_out:
        with open(pdb_in_path, 'r') as f_in:
            for line in f_in.readlines():
                if line.split()[-1] in ['N', 'C', 'O', 'S']:
                    counter = counter + 1
                    f_out.write(line[:6] + str(counter).rjust(5) + line[11:])

def convert_template_to_polyala(pdb_in_path: str, pdb_out_path:str, polyala_res):

    pdb_id = utils.get_file_name(pdb_in_path)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_in_path)
    chains = [chain.get_id() for chain in structure.get_chains()]

    if isinstance(polyala_res, dict):
        for key in polyala_res.keys():
            if key not in chains:
                raise Exception('Has not been possible to convert template to polyala. '
                                f'Chain: {key} does not exist. Available chains: {chains}.')
    elif isinstance(polyala_res, list):
        if len(chains) > 1:
            raise Exception('Has not been possible to convert template to polyala. '
                            'There is more than one chain available, select one chain in the configuration file. '
                            f'Available chains: {chains}.')           
        polyala_res = {chains[0]: polyala_res}
    else:
        raise Exception('Has not been possible to convert template to polyala.')
    
    polyala_res_dict = {}
    for key, value in polyala_res.items():
        polyala_res_list = []
        for res in value:
            res_list = str(res).replace(' ', '').split('-')
            if len(res_list) == 2:
                res_list = list(range(int(res_list[0]), int(res_list[1])+1))
            elif len(res_list) > 2:
                raise Exception('Has not been possible to convert template to polyala.')
            polyala_res_list.extend(map(int,res_list))

        polyala_res_dict[key] = list(set(polyala_res_list))

    ala_atoms_list = ['N', 'CA', 'C', 'CB', 'O']
    polyala_chains = polyala_res_dict.keys()
    atoms_del_list = []
    atoms_change_list = []

    logging.info(f'The following residues are going to be converted to polyala: {polyala_res_dict}')

    for chain in polyala_chains:
        for res in structure[0][chain]:
            if get_resseq(res) in polyala_res_dict[chain]:
                for atom in res:
                    res.resname = 'ALA'
                    if not atom.name in ala_atoms_list:
                        atoms_del_list.append(atom.get_serial_number())

    class Atom_select(Select):
            def accept_atom(self, atom):
                if atom.get_serial_number() in atoms_del_list:
                    return 0
                else:
                    return 1
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_out_path, select=Atom_select(), preserve_atom_numbering = True)

    cryst_card = extract_cryst_card_pdb(pdb_in_path=pdb_in_path)
    if cryst_card is not None:
        with open(pdb_out_path, "r+") as f_in: lines = f_in.read(); f_in.seek(0); f_in.write(f'{cryst_card}\n' + lines)

    return polyala_res_dict    

def get_resseq(residue: Residue) -> int:

    return residue.get_full_id()[3][1]

def extract_cryst_card_pdb(pdb_in_path: str) -> str:

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

    result = f'{remark:<6}{num:>5}  {name:<3}{res:>4} {chain}{resseq:>4}    {float(x):8.3f}{float(y):8.3f}{float(z):8.3f}{float(occ):6.2f}{float(bfact):6.2f}{atype:>12}\n'

    return result

def parse_pdb_line(line: str) -> Dict:
    parsed_dict = {
            'remark':line[:6],
            'num': line[6:11],
            'name': line[12:16],
            'resname': line[16:20],
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
        return

    for res in structure[0][chains[0]]:
        resseq = get_resseq(res)
        res_norm = resseq % total_length
        if res_norm == 0 or res_norm > len(a_air.query_sequence):
            if res_norm == len(a_air.query_sequence)+1:
                chain_name = chr(ord(chain_name)+1)
                chain = Chain.Chain(chain_name)
                structure[0].add(chain)
            delete_residues.append(res.id)
        elif resseq > total_length and chain is not None:
            delete_residues.append(res.id)
            res.parent = chain
            chain.add(res)
            chain[res.id].id = (' ', res_norm, ' ')

    for id in delete_residues:
        structure[0][chains[0]].detach_child(id)

    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_out_path, preserve_atom_numbering = True)

def superpose_pdbs(query_pdb: str, target_pdb: str, output_superposition: bool=True):

    # WARNING: this function is only for PDBs containing only one chain and has to be executed in the same
    # query_pdb and target_pdb path

    if output_superposition:
        superpose_output = subprocess.Popen(['superpose', f'{query_pdb}', '-s', '-all', f'{target_pdb}', '-s', '-all',
                                             '-o', f'{query_pdb[:-4]}_superposed.pdb'],
                                            stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
    else:
        superpose_output = subprocess.Popen(['superpose', f'{query_pdb}', '-s', '-all', f'{target_pdb}', '-s', '-all'],
                                            stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
    for line in superpose_output.split('\n'):
        if 'r.m.s.d:' in line:
            rmsd = float(line.split()[1])
        if 'quality Q:' in line:
            quality_q = line.split()[2]
        if 'Nalign:' in line:
            nalign = line.split()[1]

    try:
        match1 = [m.start() for m in re.finditer("TEXT:Residue alignment:", superpose_output)][0]
        match2 = [m.start() for m in re.finditer("`-------------'----------'-------------'", superpose_output)][0]
        alignment_output = superpose_output[match1:match2].split('\n')[5:]
        aligned_res_list = []
        for line in alignment_output:
            if line[23:25] == '**':
                aligned_res_list.append(int(line[36:39].replace(' ', '')))
    except:
        rmsd, nalign, quality_q, aligned_res_list = None, None, None, None

    return rmsd, nalign, quality_q, aligned_res_list

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

def split_dimers_in_pdb(pdb_in_path: str, pdb_out_path: str, chain1: List, chain2: List):

    with open(pdb_in_path, 'r') as f_in:
        input_lines = f_in.readlines()

    with open(pdb_out_path, 'w') as f_out:
        for line in input_lines:
            if line[21:22] in [chain1, chain2]:
                f_out.write(line)