import copy
import itertools
import logging
import os
import re
import shutil
import subprocess
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
from Bio import SeqIO
from Bio.PDB import PDBIO, PDBList, PDBParser, Residue, Chain, Select, Selection, Structure, Model
from scipy.spatial import distance
from simtk import unit, openmm
from sklearn.cluster import KMeans

from libs import change_res, structures, utils, sequence


def download_pdb(pdb_id: str, output_dir: str):
    pdbl = PDBList(server='https://files.wwpdb.org')
    result_ent = pdbl.retrieve_pdb_file(pdb_code=pdb_id, pdir=output_dir, file_format='pdb', obsolete=False)
    if not os.path.exists(result_ent):
        raise Exception(f'{pdb_id} could not be downloaded.')
    shutil.copy2(result_ent, os.path.join(output_dir, f'{pdb_id}.pdb'))
    os.remove(result_ent)
    shutil.rmtree('obsolete')


def pdb2mmcif(output_dir: str, pdb_in_path: str, cif_out_path: str):
    maxit_dir = os.path.join(output_dir, 'maxit')
    if not os.path.exists(maxit_dir):
        os.mkdir(maxit_dir)
    subprocess.Popen(['maxit', '-input', pdb_in_path, '-output', cif_out_path, '-o', '1'], cwd=maxit_dir,
                     stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    shutil.rmtree(maxit_dir)
    return cif_out_path


def run_lsqkab(pdb_inf_path: str, pdb_inm_path: str, fit_ini: int, fit_end: int, match_ini: int, match_end: int,
               pdb_out: str, delta_out: str):
    # Run the program lsqkab. Write the superposed pdb in pdbout and the deltas in delta_out.
    # LSQKAB will match the CA atoms from the pdb_inf to fit in the pdb_inm.

    script_path = os.path.join(os.path.dirname(pdb_out), f'{utils.get_file_name(pdb_out)}_lsqkab.sh')
    with open(script_path, 'w') as f_in:
        f_in.write('lsqkab ')
        f_in.write(f'xyzinf {utils.get_file_name(pdb_inf_path)} ')
        f_in.write(f'xyzinm {utils.get_file_name(pdb_inm_path)} ')
        f_in.write(f'DELTAS {utils.get_file_name(delta_out)} ')
        f_in.write(f'xyzout {utils.get_file_name(pdb_out)} << END-lsqkab \n')
        f_in.write('title matching template and predictions \n')
        f_in.write('output deltas \n')
        f_in.write('output XYZ \n')
        f_in.write(f'fit RESIDU CA {match_ini} TO {match_end} CHAIN A \n')
        f_in.write(f'MATCH RESIDU {fit_ini} TO {fit_end} CHAIN A \n')
        f_in.write(f'end \n')
        f_in.write(f'END-lsqkab')
    subprocess.Popen(['bash', script_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                     cwd=os.path.dirname(pdb_out)).communicate()


def check_pdb(pdb: str, output_dir: str) -> str:
    # Check if pdb is a path, and if it doesn't exist, download it.
    # If the pdb is a path, copy it to our input folder

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


def check_sequence_path(path_in: str) -> str:
    if path_in is not None:
        if not os.path.exists(path_in):
            return path_in
        else:
            return extract_sequence(path_in)


def add_cryst_card_pdb(pdb_in_path: str, cryst_card: str) -> bool:
    # Add a cryst1 record to a pdb file
    try:
        with open(pdb_in_path, 'r') as handle:
            pdb_dump = handle.read()
        with open(pdb_in_path, 'w') as handle:
            handle.write(cryst_card + "\n")
            handle.write(pdb_dump)
        return True
    except Exception as e:
        logging.info(f'Something went wrong adding the CRYST1 record to the pdb at {pdb_in_path}')
        return False


def extract_sequence(fasta_path: str) -> str:
    logging.info(f'Extracting sequence from {fasta_path}')
    try:
        record = SeqIO.read(fasta_path, "fasta")
    except Exception as e:
        raise Exception(f'Not possible to extract the sequence from {fasta_path}')
    return str(record.seq)


def extract_sequences(fasta_path: str) -> Dict:
    logging.info(f'Extracting sequences from {fasta_path}')
    records = list(SeqIO.parse(fasta_path, 'fasta'))
    return dict([(rec.id, str(rec.seq)) for rec in records])


def extract_sequence_from_file(file_path: str) -> List[str]:
    results_list = []
    extension = utils.get_file_extension(file_path)
    if extension == '.pdb':
        extraction = 'pdb-atom'
    else:
        extraction = 'cif-atom'

    try:
        with open(file_path, 'r') as f_in:
            for record in SeqIO.parse(f_in, extraction):
                results = f'>{record.id.replace("????", utils.get_file_name(file_path)[:10])}\n'
                results += str(record.seq.replace("X", ""))
                results_list.append(results)
        return results_list
    except Exception as e:
        logging.info('Something went wrong extracting the fasta record from the pdb at', file_path)
        pass
    return results_list


def write_sequence(sequence_name: str, sequence_amino: str, sequence_path: str) -> str:
    with open(sequence_path, 'w') as f_out:
        f_out.write(f'>{sequence_name}\n')
        f_out.write(f'{sequence_amino}')
    return sequence_path


def split_pdb_in_chains(output_dir: str, pdb_in_path: str) -> Dict:
    aux_path = os.path.join(output_dir, os.path.basename(pdb_in_path))
    shutil.copy2(pdb_in_path, aux_path)
    chain_dict = chain_splitter(aux_path)
    extracted_chain_dict = {k: [v] for k, v in chain_dict.items()}
    return extracted_chain_dict


def merge_pdbs(list_of_paths_of_pdbs_to_merge: List[str], merged_pdb_path: str):
    with open(merged_pdb_path, 'w') as f:
        counter = 0
        for pdb_path in list_of_paths_of_pdbs_to_merge:
            for line in open(pdb_path, 'r').readlines():
                if line[:4] == 'ATOM':
                    counter += 1
                    f.write(line[:4] + str(counter).rjust(7) + line[11:])


def merge_pdbs_in_one_chain(list_of_paths_of_pdbs_to_merge: List[str], pdb_out_path: str):
    new_structure = Structure.Structure('struct')
    new_model = Model.Model('model')
    chain = Chain.Chain('A')
    new_structure.add(new_model)
    new_model.add(chain)

    count_res = 1
    for pdb_path in list_of_paths_of_pdbs_to_merge:
        structure = get_structure(pdb_path=pdb_path)
        residues_list = list(structure[0]['A'].get_residues())
        for residue in residues_list:
            new_res = copy.copy(residue)
            new_res.parent = None
            new_res.id = (' ', count_res, ' ')
            chain.add(new_res)
            new_res.parent = chain
            count_res += 1

    io = PDBIO()
    io.set_structure(new_structure)
    io.save(pdb_out_path)


def run_pisa(pdb_path: str) -> str:
    logging.info(f'Generating REMARK 350 for {pdb_path} with PISA.')
    subprocess.Popen(['pisa', 'temp', '-analyse', pdb_path], stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE).communicate()
    pisa_output = \
        subprocess.Popen(['pisa', 'temp', '-350'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0]
    erase_pisa(name='temp')
    return pisa_output.decode('utf-8')


def erase_pisa(name: str) -> str:
    subprocess.Popen(['pisa', name, '-erase'], stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE).communicate()


def read_remark_350(pdb_path: str) -> Tuple[List[str], List[List[List[Any]]]]:
    pdb_text = open(pdb_path, 'r').read()
    match_biomolecules = [m.start() for m in
                          re.finditer(r'REMARK 350 BIOMOLECULE:', pdb_text)]  # to know how many biomolecules there are.
    if len(match_biomolecules) == 0:
        pdb_text = run_pisa(pdb_path)
        match_biomolecules = [m.start() for m in re.finditer(r'REMARK 350 BIOMOLECULE:',
                                                             pdb_text)]  # to know how many biomolecules there are.

    if len(match_biomolecules) == 0:
        raise Exception(f'REMARK not found for template {pdb_path}.')
    elif len(match_biomolecules) == 1:
        match_last_350 = [m.start() for m in re.finditer(r'REMARK 350', pdb_text)][-1]
        match_end_in_last_350 = [m.end() for m in re.finditer(r'\n', pdb_text[match_last_350:])][-1]
        remark_350_text = pdb_text[match_biomolecules[0]:(match_last_350 + match_end_in_last_350)]
    else:
        logging.info('It seem there is more than one biological assembly from REMARK 350. Only'
                     ' "BIOMOLECULE 1" will be considered for the assembly generation')
        remark_350_text = pdb_text[match_biomolecules[0]:match_biomolecules[1] - 1]

    match_biomt1 = [m.start() for m in re.finditer(r'REMARK 350 {3}BIOMT1', remark_350_text)]
    match_biomt3 = [m.end() for m in re.finditer(r'REMARK 350 {3}BIOMT3', remark_350_text)]

    end_remark_350_block = [m.start() for m in re.finditer('\n', remark_350_text[match_biomt3[-1]:])]

    transformation_blocks_indices = match_biomt1 + [match_biomt3[-1] + end_remark_350_block[0] + 1]

    transformations_list = []
    for index in range(len(transformation_blocks_indices) - 1):
        block = remark_350_text[transformation_blocks_indices[index]:transformation_blocks_indices[index + 1]]
        matrix = [item.split()[4:8] for item in block.split('\n')[:-1]]
        r11, r12, r13, t1 = matrix[0]
        r21, r22, r23, t2 = matrix[1]
        r31, r32, r33, t3 = matrix[2]
        transformation = [[r11, r12, r13], [r21, r22, r23], [r31, r32, r33], [t1, t2, t3]]
        transformations_list.append(transformation)

    chain_list = [line.split(':')[-1].replace(' ', '').split(',') for line in remark_350_text.split('\n')
                  if 'REMARK 350 APPLY THE FOLLOWING TO CHAINS:' in line][0]

    return chain_list, transformations_list


def change_chain(pdb_in_path: str, pdb_out_path: str, rot_tra_matrix: List[List] = None, offset: Optional[int] = 0,
                 chain: Optional[str] = None):
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


def remove_hydrogens(pdb_in_path: str, pdb_out_path: str):
    counter = 0
    with open(pdb_out_path, 'w') as f_out:
        with open(pdb_in_path, 'r') as f_in:
            for line in f_in.readlines():
                if line.split()[-1] in ['N', 'C', 'O', 'S']:
                    counter = counter + 1
                    f_out.write(line[:6] + str(counter).rjust(5) + line[11:])


def get_resseq(residue: Residue) -> int:
    # Return resseq number
    return residue.get_full_id()[3][1]


def get_hetatm(residue: Residue) -> int:
    # Return hetatm
    return residue.get_full_id()[3][0]


def get_chains(pdb_path: str) -> List[str]:
    # Return all chains from a PDB structure
    structure = get_structure(pdb_path)
    return [chain.get_id() for chain in structure.get_chains()]


def get_structure(pdb_path: str) -> Structure:
    # Get PDB structure
    pdb_id = utils.get_file_name(pdb_path)
    parser = PDBParser(QUIET=True)
    return parser.get_structure(pdb_id, pdb_path)


def run_pdb2cc(templates_path: str, pdb2cc_path: str = None) -> str:
    cwd = os.getcwd()
    os.chdir(templates_path)
    output_path = 'cc_analysis.in'
    if pdb2cc_path is None:
        pdb2cc_path = 'pdb2cc'
    command_line = f'{pdb2cc_path} -m -i 10 -y 0.5 "orig.*.pdb" 0 {output_path}'
    p = subprocess.Popen(command_line, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    p.communicate()
    os.chdir(cwd)
    return os.path.join(templates_path, output_path)


def run_cc_analysis(input_path: str, n_clusters: int, cc_analysis_path: str = None) -> str:
    output_path = 'cc_analysis.out'
    if cc_analysis_path is None:
        cc_analysis_path = 'cc_analysis'
    cwd = os.getcwd()
    os.chdir(os.path.dirname(input_path))
    command_line = f'{cc_analysis_path} -dim {n_clusters} {os.path.basename(input_path)} {output_path}'
    p = subprocess.Popen(command_line, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    p.communicate()
    os.chdir(cwd)
    return f'{os.path.join(os.path.dirname(input_path), output_path)}'


def cc_analysis(paths_in: Dict, cc_analysis_paths: structures.CCAnalysis, cc_path: str, n_clusters: int = 2) -> List:
    utils.create_dir(cc_path, delete_if_exists=True)
    paths = [shutil.copy2(path, cc_path) for path in paths_in.values()]
    trans_dict = {}
    return_templates_cluster = [[] for _ in range(n_clusters)]
    clean_dict = {}

    for index, path in enumerate(paths):
        if utils.check_ranked(os.path.basename(path)):
            bfactors_dict = read_bfactors_from_residues(path)
            for chain, residues in bfactors_dict.items():
                for i in range(len(residues)):
                    if bfactors_dict[chain][i] is not None:
                        bfactors_dict[chain][i] = round(bfactors_dict[chain][i] - 70.0, 2)
            residues_dict = read_residues_from_pdb(path)
            change = change_res.ChangeResidues(chain_res_dict=residues_dict, chain_bfactors_dict=bfactors_dict)
            change.change_bfactors(path, path)
        new_path = os.path.join(cc_path, f'orig.{str(index)}.pdb')
        os.rename(os.path.join(cc_path, path), new_path)
        trans_dict[index] = utils.get_file_name(path)
    if trans_dict:
        with open(os.path.join(cc_path, 'labels.txt'), 'w+') as f:
            for key, value in trans_dict.items():
                f.write('%s:%s\n' % (key, value))
        output_pdb2cc = run_pdb2cc(templates_path=cc_path, pdb2cc_path=cc_analysis_paths.pd2cc_path)
        if os.path.exists(output_pdb2cc):
            output_cc = run_cc_analysis(input_path=output_pdb2cc,
                                        n_clusters=n_clusters,
                                        cc_analysis_path=cc_analysis_paths.cc_analysis_path)
            if os.path.exists(output_cc):
                cc_analysis_dict = utils.parse_cc_analysis(file_path=output_cc)
                for key, values in cc_analysis_dict.items():
                    if values.module > 0.1 or values.module < -0.1:
                        clean_dict[trans_dict[int(key) - 1]] = values
                points = np.array([[values.x, values.y] for values in clean_dict.values()])
                kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(points)
                lookup_table = {}
                counter = 0
                for i in kmeans.labels_:
                    if i not in lookup_table:
                        lookup_table[i] = counter
                        counter += 1
                conversion = [lookup_table[label] for label in kmeans.labels_]
                for i, label in enumerate(conversion):
                    return_templates_cluster[int(label)].append(paths_in[list(clean_dict.keys())[i]])

    return return_templates_cluster, clean_dict


def extract_cryst_card_pdb(pdb_in_path: str) -> Union[str, None]:
    # Extract the crystal card from a pdb

    if os.path.isfile(pdb_in_path):
        with open(pdb_in_path, 'r') as f_in:
            pdb_lines = f_in.readlines()
        for line in pdb_lines:
            if line.startswith("CRYST1"):
                cryst_card = line
                return cryst_card
    return None


def get_atom_line(remark: str, num: int, name: str, res: int, chain: str, resseq, x: float, y: float, z: float,
                  occ: str, bfact: str, atype: str) -> str:
    # Given all elements of an atom, parse them in PDB format

    result = f'{remark:<6}{num:>5}  {name:<3}{res:>4} {chain}{resseq:>4}    {float(x):8.3f}{float(y):8.3f}{float(z):8.3f}{float(occ):6.2f}{float(bfact):6.2f}{atype:>12}\n'
    return result


def parse_pdb_line(line: str) -> Dict:
    # Parse all elements of an atom of a PDB line

    parsed_dict = {
        'remark': line[:6],
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


def read_bfactors_from_residues(pdb_path: str) -> Dict:
    # Create a dictionary with each existing chain in the pdb.
    # In each chain, create a list of N length (corresponding to the number of residues)
    # Copy the bfactor in the corresponding residue number in the list.
    structure = get_structure(pdb_path=pdb_path)
    return_dict = {}
    for chain in structure[0]:
        return_dict[chain.get_id()] = []
        for res in list(chain.get_residues()):
            return_dict[chain.get_id()].append(res.get_unpacked_list()[0].bfactor)
    return return_dict


def read_residues_from_pdb(pdb_path: str) -> Dict:
    # Create a dictionary with each existing chain in the pdb.
    # In each chain, a list with the residue numbers
    structure = get_structure(pdb_path=pdb_path)
    return_dict = {}
    for chain in structure[0]:
        return_dict[chain.get_id()] = []
        for res in list(chain.get_residues()):
            return_dict[chain.get_id()].append(get_resseq(res))
    return return_dict


def split_chains_assembly(pdb_in_path: str,
                          pdb_out_path: str,
                          sequence_assembled: sequence.SequenceAssembled) -> Dict:
    # Split the assembly with several chains. The assembly is spitted
    # by the query sequence length. Also, we have to take into account
    # the glycines, So every query_sequence+glycines we can find a chain.
    # We return the list of chains.

    structure = get_structure(pdb_path=pdb_in_path)
    chains_return = {}
    chains = [chain.get_id() for chain in structure.get_chains()]

    if len(chains) > 1:
        logging.info(f'PDB: {pdb_in_path} is already split in several chains: {chains}')
        shutil.copy2(pdb_in_path, pdb_out_path)
    else:
        new_structure = Structure.Structure(structure.get_id)
        new_model = Model.Model(structure[0].id)
        new_structure.add(new_model)
        residues_list = list(structure[0][chains[0]].get_residues())
        idres_list = list([get_resseq(res) for res in residues_list])
        original_chain_name = chains[0]

        for i in range(sequence_assembled.total_copies):
            sequence_length = sequence_assembled.get_sequence_length(i)
            start_min = sequence_assembled.get_starting_length(i)
            start_max = start_min + sequence_length

            chain_name = chr(ord(original_chain_name) + i)
            chain = Chain.Chain(chain_name)
            new_structure[0].add(chain)
            mapping = {}
            for new_id, j in enumerate(range(start_min + 1, start_max + 1), start=1):
                if j in idres_list:
                    res = residues_list[idres_list.index(j)]
                    mapping[new_id] = j
                    new_res = copy.copy(res)
                    chain.add(new_res)
                    new_res.parent = chain
                    chain[new_res.id].id = (' ', new_id, ' ')
            chains_return[chain_name] = mapping

        io = PDBIO()
        io.set_structure(new_structure)
        io.save(pdb_out_path)
    return chains_return


def chain_splitter(pdb_path: str, chain: str = None) -> Dict:
    # Given a pdb_in and an optional chain, write one or several
    # pdbs containing each one a chain.
    # If chain is specified, only one file with the specific chain will be created
    # It will return a dictionary with the chain and the corresponding pdb

    return_chain_dict = {}
    structure = get_structure(pdb_path=pdb_path)
    chains = [chain.get_id() for chain in structure.get_chains()] if chain is None else [chain]

    for chain in chains:
        new_pdb = os.path.join(os.path.dirname(pdb_path), f'{utils.get_file_name(pdb_path)}_{chain}1.pdb')

        class ChainSelect(Select):
            def __init__(self, select_chain):
                self.chain = select_chain

            def accept_chain(self, select_chain):
                if select_chain.get_id() == self.chain:
                    return 1
                else:
                    return 0

        io = PDBIO()
        io.set_structure(structure)
        io.save(new_pdb, ChainSelect(chain))
        return_chain_dict[chain] = new_pdb

    return return_chain_dict


def generate_multimer_from_pdb(pdb_in_path: str, pdb_out_path: str):
    # Given a pdb_in, create the multimer and save it in pdb_out

    shutil.copy2(pdb_in_path, pdb_out_path)
    chain_dict = chain_splitter(pdb_out_path)
    multimer_chain_dict = dict(sorted(generate_multimer_chains(pdb_out_path, chain_dict).items()))
    chain_name = next(iter(multimer_chain_dict))
    result_chain_dict = {}
    for _, elements in multimer_chain_dict.items():
        for path in elements:
            result_chain_dict[chain_name] = path
            chain_name = chr(ord(chain_name) + 1)
    change_chains(result_chain_dict)
    merge_pdbs(utils.dict_values_to_list(result_chain_dict), pdb_out_path)


def change_chains(chain_dict: Dict):
    # The Dict has to be: {A: path}
    # It will rename the chains of the path to the
    # chain indicated in the key

    for key, value in chain_dict.items():
        change_chain(pdb_in_path=value,
                     pdb_out_path=value,
                     chain=key)


def generate_multimer_chains(pdb_path: str, template_dict: Dict) -> Dict:
    # Read remark to get the transformations and the new chains
    # Apply transformations to generate the new ones
    # Rename chains with A1, A2...
    # Store a dict with the relation between old chains and new chains
    # Dict -> A: [path_to_A1, path_to_A2]

    chain_list, transformations_list = read_remark_350(pdb_path)
    multimer_dict = {}

    logging.info(
        'Assembly can be build using chain(s) ' + str(chain_list) + ' by applying the following transformations:')
    for matrix in transformations_list:
        logging.info(str(matrix))

    for chain in chain_list:
        if isinstance(template_dict[chain], list):
            pdb_path = template_dict[chain][0]
        else:
            pdb_path = template_dict[chain]
        multimer_new_chains = []
        for i, transformation in enumerate(transformations_list):
            new_pdb_path = utils.replace_last_number(text=pdb_path, value=i + 1)
            change_chain(pdb_in_path=pdb_path,
                         pdb_out_path=new_pdb_path,
                         rot_tra_matrix=transformation)
            multimer_new_chains.append(new_pdb_path)
        multimer_dict[chain] = multimer_new_chains

    return multimer_dict


def remove_hetatm(pdb_in_path: str, pdb_out_path: str):
    # Transform MSE HETATM to MSA ATOM
    # Remove HETATM from pdb

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
                        atom.id = 'SD'
                        atom.fullname = 'SD'
                        atom.name = 'SD'

    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_out_path, NonHetSelect())


def run_pdbfixer(pdb_in_path: str, pdb_out_path: str):
    command_line = f'pdbfixer {os.path.abspath(pdb_in_path)} --output={pdb_out_path} --add-atoms=all ' \
                   f'--keep-heterogens=none --replace-nonstandard --add-residues --ph=7.0 '
    p = subprocess.Popen(command_line, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    p.communicate()


def run_arcimboldo_air(yml_path: str):
    arcimboldo_air_path = os.path.join(utils.get_main_path(), 'arcimboldo_air.py')
    command_line = f'{arcimboldo_air_path} {yml_path}'
    p = subprocess.Popen(command_line, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    logging.info('ARCIMBOLDO_AIR cluster run finished successfully.')


def run_openmm(pdb_in_path: str, pdb_out_path: str) -> structures.OpenmmEnergies:
    run_pdbfixer(pdb_in_path=pdb_in_path, pdb_out_path=pdb_out_path)
    protein_pdb = openmm.app.pdbfile.PDBFile(pdb_out_path)
    forcefield = openmm.app.ForceField('amber99sb.xml')
    system = forcefield.createSystem(protein_pdb.topology, constraints=openmm.app.HBonds)
    integrator = openmm.openmm.LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    simulation = openmm.app.simulation.Simulation(protein_pdb.topology, system, integrator)
    simulation.context.setPositions(protein_pdb.positions)
    simulation.minimizeEnergy()
    simulation.step(1000)
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    with open(pdb_out_path, 'w') as f_out:
        openmm.app.pdbfile.PDBFile.writeFile(protein_pdb.topology, state.getPositions(), file=f_out, keepIds=True)
    return structures.OpenmmEnergies(round(state.getKineticEnergy()._value, 2),
                                     round(state.getPotentialEnergy()._value, 2))


def superpose_pdbs(pdb_list: List, output_path: str = None) -> Tuple[Optional[float], Optional[str], Optional[str]]:
    superpose_input_list = ['superpose']
    for pdb in pdb_list:
        superpose_input_list.extend([pdb, '-s', '-all'])
    if output_path is not None:
        superpose_input_list.extend(['-o', output_path])

    superpose_output = subprocess.Popen(superpose_input_list, stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
    rmsd, quality_q, nalign = None, None, None
    for line in superpose_output.split('\n'):
        if 'r.m.s.d:' in line:
            rmsd = float(line.split()[1])
        if 'quality Q:' in line:
            quality_q = line.split()[2]
        if 'Nalign:' in line:
            nalign = line.split()[1]
    return rmsd, nalign, quality_q


def gesamt_pdbs(pdb_list: List, output_path: str = None) -> Tuple[Optional[float], Optional[str], Optional[str]]:
    superpose_input_list = ['gesamt']
    for pdb in pdb_list:
        superpose_input_list.extend([pdb, '-s', '-all'])
    if output_path is not None:
        superpose_input_list.extend(['-o', output_path])

    superpose_output = subprocess.Popen(superpose_input_list, stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
    rmsd, quality_q, nalign = None, None, None
    for line in superpose_output.split('\n'):
        if 'RMSD             :' in line:
            rmsd = float(line.split()[1])
        if 'Q-score          :' in line:
            quality_q = line.split()[2]
        if 'Aligned residues :' in line:
            nalign = line.split()[1]
    return rmsd, nalign, quality_q


def pdist(query_pdb: str, target_pdb: str) -> float:
    if query_pdb is None or target_pdb is None:
        return 1.0

    structure_query = get_structure(pdb_path=query_pdb)
    res_query_list = [res.id[1] for res in Selection.unfold_entities(structure_query, 'R')]

    structure_target = get_structure(pdb_path=target_pdb)
    res_target_list = [res.id[1] for res in Selection.unfold_entities(structure_target, 'R')]

    common_res_list = list(set(res_query_list) & set(res_target_list))
    if not common_res_list:
        return 0.9

    query_common_list = [res for res in Selection.unfold_entities(structure_query, 'R') if res.id[1] in common_res_list]
    query_matrix = calculate_distance_pdist(res_list=query_common_list)

    target_common_list = [res for res in Selection.unfold_entities(structure_target, 'R') if
                          res.id[1] in common_res_list]
    target_matrix = calculate_distance_pdist(res_list=target_common_list)

    diff_pdist_matrix = np.abs(query_matrix - target_matrix)

    return float(diff_pdist_matrix.mean())


def calculate_distance_pdist(res_list: List) -> List:
    coords = [res['CA'].coord for res in res_list]
    calculate_pdist = distance.pdist(coords, "euclidean")
    return distance.squareform(calculate_pdist)


def find_interface_from_pisa(pdb_in_path: str, interfaces_path: str) -> List[Union[Dict, None]]:
    interface_data_list = []
    pisa_text = subprocess.Popen(['pisa', 'temp', '-analyse', pdb_in_path],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE).communicate()[0].decode('utf-8')

    pisa_output = subprocess.Popen(['pisa', 'temp', '-list', 'interfaces'], stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE).communicate()[0].decode('utf-8')

    pisa_general_txt = os.path.join(interfaces_path, f'{utils.get_file_name(pdb_in_path)}_general_output.txt')
    with open(pisa_general_txt, 'w') as f_out:
        f_out.write(pisa_output)

    if 'NO INTERFACES FOUND' in pisa_output or 'no chains found in input file' in pisa_text:
        logging.info('No interfaces found in pisa')
    else:
        interfaces_list = utils.parse_pisa_general_multimer(pisa_output)
        for interface in interfaces_list:
            serial_output = \
                subprocess.Popen(['pisa', 'temp', '-detail', 'interfaces', interface['serial']], stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE).communicate()[0].decode('utf-8')

            interface_data = utils.parse_pisa_interfaces(serial_output)
            interface_data.update(interface)
            interface_data_list.append(interface_data)
            pisa_output_txt = os.path.join(interfaces_path,
                                           f'{utils.get_file_name(pdb_in_path)}_{interface_data["chain1"]}{interface_data["chain2"]}_interface.txt')
            with open(pisa_output_txt, 'w') as f_out:
                f_out.write(serial_output)

    erase_pisa(name='temp')

    return interface_data_list


def create_interface_domain(pdb_in_path: str, pdb_out_path: str, interface: Dict, domains_dict: Dict) \
        -> Dict[Any, List[Any]]:
    add_domains_dict = {}
    bfactors_dict = {}
    for chain, residue in zip([interface['chain1'], interface['chain2']],
                              [interface['res_chain1'], interface['res_chain2']]):
        added_res_list = []
        [added_res_list.extend(domains) for domains in domains_dict[chain] if bool(set(residue).intersection(domains))]
        added_res_list.extend(residue)
        add_domains_dict[chain] = list(set(added_res_list))
        bfactors_dict[chain] = [float(interface['bfactor'])] * len(add_domains_dict[chain])

    split_dimers_in_pdb(pdb_in_path=pdb_in_path,
                        pdb_out_path=pdb_out_path,
                        chain1=interface['chain1'],
                        chain2=interface['chain2'])

    change = change_res.ChangeResidues(chain_res_dict=add_domains_dict, chain_bfactors_dict=bfactors_dict)
    change.delete_residues_inverse(pdb_out_path, pdb_out_path)
    change.change_bfactors(pdb_out_path, pdb_out_path)

    return add_domains_dict


def calculate_auto_offset(input_list: List[List], length: int) -> List[int]:
    if length <= 0:
        return []

    combinated_list = list(itertools.product(*input_list))
    trimmed_list = []
    for element in combinated_list:
        aux_length = length
        element_aux = [x for x in element if x[3] is not False]
        if len(element_aux) <= 0:
            continue
        elif len(element_aux) < length:
            aux_length = len(element_aux)
        sorted_list = sorted(element_aux, key=lambda x: x[2])[:aux_length]
        target_list = [target for _, target, _, _ in sorted_list]
        if len(target_list) == len(set(target_list)):
            trimmed_list.append(sorted_list)

    score_list = []
    for element in trimmed_list:
        score_list.append(sum(z for _, _, z, _ in element))

    if score_list:
        min_value = min(score_list)
        min_index = score_list.index(min_value)
        return trimmed_list[min_index]
    else:
        return []


def split_dimers_in_pdb(pdb_in_path: str, pdb_out_path: str, chain1: List, chain2: List):
    with open(pdb_in_path, 'r') as f_in:
        input_lines = f_in.readlines()

    with open(pdb_out_path, 'w') as f_out:
        for line in input_lines:
            if line[21:22] in [chain1, chain2]:
                f_out.write(line)
