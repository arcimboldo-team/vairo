import copy
import logging
import os
import re
import shutil
import subprocess
import numpy as np
import itertools
from typing import Any, Dict, List, Tuple, Union
from Bio import SeqIO
from Bio.PDB import PDBIO, PDBList, PDBParser, Residue, Chain, Select, Selection, Structure, Model
from libs import change_res, utils
from scipy.spatial import distance
from libs import sequence

def download_pdb(pdb_id: str, output_dir: str):

    pdbl = PDBList()
    result_ent = pdbl.retrieve_pdb_file(pdb_id, pdir=output_dir, file_format='pdb', obsolete=False)
    if not os.path.exists(result_ent):
        raise Exception(f'{pdb_id} could not be downloaded.')
    shutil.copy2(result_ent, os.path.join(output_dir,f'{pdb_id}.pdb'))
    os.remove(result_ent)
    shutil.rmtree('obsolete')

def pdb2mmcif(output_dir: str, pdb_in_path: str, cif_out_path: str):
    
    maxit_dir = os.path.join(output_dir, 'maxit')
    if not os.path.exists(maxit_dir):
        os.mkdir(maxit_dir)
    subprocess.Popen(['maxit', '-input', pdb_in_path, '-output', cif_out_path, '-o', '1'], cwd=maxit_dir,stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    shutil.rmtree(maxit_dir)
    return cif_out_path

def cif2pdb(cif_in_path: str, pdb_out_path: str):
    
    parser = MMCIFParser(QUIET=True)
    struc = parser.get_structure('', cif_in_path)
    io = PDBIO()
    io.set_structure(struc)
    io.save(pdb_out_path)

def check_pdb(pdb: str, output_dir: str) -> str:
    #Check if pdb is a path, and if it doesn't exist, download it.
    #If the pdb is a path, copy it to our input folder

    if not os.path.exists(pdb):
        download_pdb(pdb_id=pdb, output_dir=output_dir)
        pdb = os.path.join(output_dir, f'{pdb}.pdb')
    else:
        pdb_aux = os.path.join(output_dir, os.path.basename(pdb))
        if pdb != pdb_aux:
            shutil.copy2(pdb, pdb_aux)
            pdb = pdb_aux
    
    cryst_card = extract_cryst_card_pdb(pdb_in_path=pdb)
    remove_hetatm(pdb, pdb)
    if cryst_card is not None:
        add_cryst_card_pdb(pdb_in_path=pdb, cryst_card=cryst_card)

    return pdb

def add_cryst_card_pdb(pdb_in_path: str, cryst_card: str) -> bool:
    # Add a cryst1 record to a pdb file
    try:
        with open(pdb_in_path, 'r') as handle:
            pdb_dump = handle.read()
        with open(pdb_in_path, 'w') as handle:
            handle.write(cryst_card + "\n")
            handle.write(pdb_dump)
        return True
    except:
        logging.info(f'Something went wrong adding the CRYST1 record to the pdb at {pdb_in_path}')
        return False

def extract_sequence(fasta_path: str) -> str:

    logging.info(f'Extracting sequence from {fasta_path}')
    try:
        record = SeqIO.read(fasta_path, "fasta")
    except:
        raise Exception(f'Not possible to extract the sequence from {fasta_path}')
    return str(record.seq)

def extract_sequence_from_file(file_path: str) -> Union[str, None]:
    results = ''
    extension = utils.get_file_extension(file_path)
    if extension == '.pdb':
        extraction = 'pdb-atom'
    else:
        extraction = 'cif-atom'

    try:
        with open(file_path, 'r') as f_in:
            for record in SeqIO.parse(f_in, extraction):
                results += f'>{(record.id).replace("????", utils.get_file_name(file_path))}\n'
                results += f'{record.seq}\n'
        return results
    except:
        logging.info('Something went wrong extracting the fasta record from the pdb at', file_path)
        pass
    return None

def write_sequence(sequence: str, sequence_path: str) -> str:
    
    with open(sequence_path, 'w') as f_out:
        f_out.write(f'{sequence}\n')
    return sequence_path
    
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

def get_hetatm(residue: Residue) -> int:
    #Return hetatm 

    return residue.get_full_id()[3][0]

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

def split_chains_assembly(pdb_in_path: str, 
                        pdb_out_path:str, 
                        sequence_assembled: sequence.SequenceAssembled) -> List:
    #Split the assembly with serveral chains. The assembly is spitted 
    #by the query sequence length. Also, we have to take into account 
    #the glycines, So every query_sequence+glycines we can find a chain.
    #We return the list of chains.

    structure = get_structure(pdb_path=pdb_in_path)
    chains_return = []
    chains = [chain.get_id() for chain in structure.get_chains()]

    if len(chains) > 1:
        logging.info(f'PDB: {pdb_in_path} is already splitted in several chains: {chains}')
        shutil.copy2(pdb_in_path, pdb_out_path)
        chains_return = chains
    else:
        new_structure = Structure.Structure(structure.get_id)
        new_model = Model.Model(structure[0].id)
        new_structure.add(new_model)
        residues_list = list(structure[0][chains[0]].get_residues())
        idres_list = list([get_resseq(res) for res in residues_list])
        original_chain_name = chains[0]

        for i in range(sequence_assembled.total_copies):
            sequence_length = sequence_assembled.get_sequence_length(i)
            seq_with_glycines_length = sequence_length + sequence_assembled.glycines
            
            min = sequence_assembled.get_starting_length(i)
            max = min+sequence_length

            chain_name = chr(ord(original_chain_name)+i)
            chain = Chain.Chain(chain_name)
            new_structure[0].add(chain)
            for j in range(min+1, max+1):
                if j in idres_list:
                    res = residues_list[idres_list.index(j)]
                    new_res = copy.copy(res)
                    chain.add(new_res)
                    new_res.parent = chain
                    chain[new_res.id].id = (' ', j % seq_with_glycines_length, ' ')
            if list(new_structure[0][chain.id].get_residues()):
                chains_return.append(chain_name)
            else:
                chains_return.append(None)

        io = PDBIO()
        io.set_structure(new_structure)
        io.save(pdb_out_path, preserve_atom_numbering = True)
    return chains_return

def chain_splitter(pdb_path: str, chain: str = None) -> Dict:
    #Given a pdb_in and a optional chain, write one or serveral
    #pdbs containing each one a chain.
    #If chain is specified, only one file with the specific chain will be created
    #It will return a dictionary with the chain and the corresponding pdb

    return_chain_dict = {}

    structure = get_structure(pdb_path=pdb_path)
    chains = [chain.get_id() for chain in structure.get_chains()] if chain is None else [chain]

    for chain in chains:
        new_pdb = os.path.join(os.path.dirname(pdb_path), f'{utils.get_file_name(pdb_path)}_{chain}1.pdb')
        class ChainSelect(Select):
            def __init__(self, chain):
                self.chain = chain

            def accept_chain(self, chain):
                if chain.get_id() == self.chain:
                    return 1
                else:          
                    return 0

        io = PDBIO()               
        io.set_structure(structure)
        io.save(new_pdb, ChainSelect(chain))
        return_chain_dict[chain] = new_pdb
    
    return return_chain_dict

def generate_multimer_from_pdb(pdb_in_path: str, pdb_out_path: str):
    #Given a pdb_in, create the multimer and save it in pdb_out
    
    shutil.copy2(pdb_in_path, pdb_out_path)
    chain_dict = chain_splitter(pdb_out_path)
    multimer_chain_dict = dict(sorted(generate_multimer_chains(pdb_out_path, chain_dict).items()))
    chain_name = next(iter(multimer_chain_dict))
    result_chain_dict = {}
    for _, elements in multimer_chain_dict.items():
        for path in elements:
            result_chain_dict[chain_name] = path
            chain_name = chr(ord(chain_name)+1)
    change_chains(result_chain_dict)
    merge_pdbs(utils.dict_values_to_list(result_chain_dict), pdb_out_path)

def change_chains(chain_dict: Dict):
    #The Dict has to be: {A: path}
    #It will renaim the chains of the path to the
    #chain indicated in the key

    for key, value in chain_dict.items():
        change_chain(pdb_in_path=value,
                    pdb_out_path=value,
                    chain=key)        

def generate_multimer_chains(pdb_path:str, template_dict: Dict) -> Dict:
    #Read remark to get the transformations and the new chains
    #Apply transformations to generate the new ones
    #Rename chains with A1, A2...
    #Store a dict with the relation between old chains and new chains
    # Dict -> A: [path_to_A1, path_to_A2]

    chain_list, transformations_list = read_remark_350(pdb_path)
    multimer_dict = {}
    
    logging.info('Assembly can be build using chain(s) '+ str(chain_list) + ' by applying the following transformations:')
    for matrix in transformations_list:
        logging.info(str(matrix))

    for chain in chain_list:
        if isinstance(template_dict[chain], list):
            pdb_path = template_dict[chain][0]
        else:
            pdb_path = template_dict[chain]
        multimer_new_chains = []
        for i, transformation in enumerate(transformations_list):
            new_pdb_path = utils.replace_last_number(text=pdb_path, value=i+1)
            change_chain(pdb_in_path=pdb_path, 
                    pdb_out_path=new_pdb_path,
                    rot_tra_matrix=transformation)
            multimer_new_chains.append(new_pdb_path)
        multimer_dict[chain] = multimer_new_chains

    return multimer_dict

def remove_hetatm(pdb_in_path:str, pdb_out_path: str):
    #Tranform MSE HETATM to MSA ATOM
    #Remove HETATM from pdb

    class NonHetSelect(Select):
        def accept_residue(self, residue):
            return 1 if residue.id[0] == " " else 0
    structure = get_structure(pdb_path=pdb_in_path)
    for chain in structure[0]:
        for res in list(chain.get_residues()):
            if get_hetatm(res) == 'H_MSE':
                res.id = (' ', get_resseq(res), ' ')
                res.resname = 'MET'
                for atom in res:
                    if atom.element == 'SE':
                        atom.element = 'S'

    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_out_path, NonHetSelect())

def superpose_pdbs(pdb_list: List, output_pdb = None) -> List:
    
    superpose_input_list = ['superpose']
    for pdb in pdb_list:
        superpose_input_list.extend([pdb, '-s', '-all'])
    if output_pdb is not None:
        superpose_input_list.extend(['-o', output_pdb])
    
    superpose_output = subprocess.Popen(superpose_input_list, stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
    
    rmsd, quality_q, nalign = None,  None, None
    for line in superpose_output.split('\n'):
        if 'r.m.s.d:' in line:
            rmsd = float(line.split()[1])
        if 'quality Q:' in line:
            quality_q = line.split()[2]
        if 'Nalign:' in line:
            nalign = line.split()[1]
    return rmsd, nalign, quality_q

def pdist(query_pdb: str, target_pdb: str) -> List[List]:

    if query_pdb is None or target_pdb is None:
        return 1

    structure_query = get_structure(pdb_path=query_pdb)
    res_query_list = [res.id[1] for res in Selection.unfold_entities(structure_query, 'R')]

    structure_target = get_structure(pdb_path=target_pdb)
    res_target_list = [res.id[1] for res in Selection.unfold_entities(structure_target, 'R')]

    common_res_list = list(set(res_query_list) & set(res_target_list))
    if not common_res_list:
        return 0.9

    query_common_list = [res for res in Selection.unfold_entities(structure_query, 'R') if res.id[1] in common_res_list]
    query_matrix = calculate_distance_pdist(res_list=query_common_list)

    target_common_list = [res for res in Selection.unfold_entities(structure_target, 'R') if res.id[1] in common_res_list]
    target_matrix = calculate_distance_pdist(res_list=target_common_list)

    diff_pdist_matrix = np.abs(query_matrix - target_matrix)

    return diff_pdist_matrix.mean()

def calculate_distance_pdist(res_list: List) -> float:

    coords = [res['CA'].coord for res in res_list]
    pdist = distance.pdist(coords, "euclidean")
    return distance.squareform(pdist)

def find_interface_from_pisa(pdb_in_path: str, interfaces_path: str) -> List:

    subprocess.Popen(['pisa', 'temp', '-analyse', pdb_in_path],
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE).communicate()

    pisa_output = subprocess.Popen(['pisa', 'temp', '-list', 'interfaces'], stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE).communicate()[0].decode('utf-8')
    pisa_general_txt = os.path.join(interfaces_path, 'general_output.txt')
    with open (pisa_general_txt, 'w') as f_out:
        f_out.write(pisa_output)

    if 'NO INTERFACES FOUND' in pisa_output:
        logging.info('No interfaces found in pisa')
        return

    interfaces_list = utils.parse_pisa_general_multimer(pisa_output)
    interface_data_list = []
    for interface in interfaces_list:
        serial_output = subprocess.Popen(['pisa', 'temp', '-detail', 'interfaces', interface['serial']], stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE).communicate()[0].decode('utf-8')      
        
        interface_data = utils.parse_pisa_interfaces(serial_output)
        interface_data.update(interface)
        interface_data_list.append(interface_data)
        pisa_output_txt = os.path.join(interfaces_path, f'{utils.get_file_name(pdb_in_path)}_{interface_data["chain1"]}{interface_data["chain2"]}_interface.txt')
        with open (pisa_output_txt, 'w') as f_out:
            f_out.write(serial_output)

    return interface_data_list 

def create_interface_domain(pdb_in_path: str, interface: Dict, interfaces_path: str, domains_dict: Dict):
    add_domains_dict = {}
    bfactors_dict = {}
    for chain, residue in zip([interface['chain1'], interface['chain2']], [interface['res_chain1'], interface['res_chain2']]):
        added_res_list = []
        for domains in domains_dict[chain]:
            if bool(set(residue).intersection(domains)):
                added_res_list.extend(domains)
        added_res_list.extend(residue)
        add_domains_dict[chain] = list(set(added_res_list))
        bfactors_dict[chain] = [float(interface['bfactor'])] * len(add_domains_dict[chain])

    dimers_path = os.path.join(interfaces_path, f'{utils.get_file_name(pdb_in_path)}_{interface["chain1"]}{interface["chain2"]}.pdb')
    split_dimers_in_pdb(pdb_in_path=pdb_in_path,
                        pdb_out_path=dimers_path,
                        chain1=interface['chain1'],
                        chain2=interface['chain2'])

    change = change_res.ChangeResidues(chain_res_dict=add_domains_dict, chain_bfactors_dict=bfactors_dict)
    change.delete_residues_inverse(dimers_path, dimers_path)
    change.change_bfactors(dimers_path, dimers_path)

def calculate_auto_offset(input_list: List[List], length: int) -> List:
    
    combinated_list = list(itertools.product(*input_list))
    trimmed_list = []
    for element in combinated_list:
        sorted_list = sorted(element, key=lambda x:x[2])[:length]
        target_list = [target for _,target,_ in sorted_list]
        if len(target_list) == len(set(target_list)):
            trimmed_list.append(sorted_list)
    
    score_list = []
    for element in trimmed_list:
        score_list.append(sum(z for _,_,z in element))
    min_value = min(score_list)
    min_index = score_list.index(min_value)

    return trimmed_list[min_index]

def split_dimers_in_pdb(pdb_in_path: str, pdb_out_path: str, chain1: List, chain2: List):

    with open(pdb_in_path, 'r') as f_in:
        input_lines = f_in.readlines()

    with open(pdb_out_path, 'w') as f_out:
        for line in input_lines:
            if line[21:22] in [chain1, chain2]:
                f_out.write(line)