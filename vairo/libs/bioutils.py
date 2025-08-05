from collections import OrderedDict, defaultdict
from typing import Any, Dict, List, Optional, Tuple, Union
import copy, itertools, logging, os, re, shutil, subprocess, sys, tempfile
import numpy as np
import pandas as pd
import io
from Bio import SeqIO
from Bio.PDB import PDBIO, PDBList, PDBParser, Residue, Chain, Select, Selection, Structure, Model, PPBuilder, \
    Superimposer
from scipy.spatial import distance
from sklearn.cluster import KMeans
from ALEPH.aleph.core import ALEPH
from alphafold.common import residue_constants
from alphafold.relax import cleanup, amber_minimize
from simtk import unit
from libs import hhsearch, structures, utils, plots, global_variables, sequence, template_modifications


def download_pdb(pdb_id: str, pdb_path: str) -> str:
    pdbl = PDBList(server='https://files.wwpdb.org', verbose=False)
    result_ent = pdbl.retrieve_pdb_file(pdb_code=pdb_id, file_format='pdb', pdir='.', obsolete=False)
    if not os.path.exists(result_ent):
        raise Exception(f'{pdb_id} could not be downloaded.')
    shutil.copy2(result_ent, pdb_path)
    os.remove(result_ent)
    shutil.rmtree('obsolete')
    return pdb_path


def pdb2mmcif(pdb_in_path: str, cif_out_path: str) -> str:
    maxit_dir = os.path.join(os.path.dirname(cif_out_path), 'maxit')
    if not os.path.exists(maxit_dir):
        os.mkdir(maxit_dir)
    subprocess.Popen(['maxit', '-input', pdb_in_path, '-output', cif_out_path, '-o', '1'], cwd=maxit_dir,
                     stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    shutil.rmtree(maxit_dir)
    return cif_out_path


def superposition_by_chains(pdb1_in_path: str, pdb2_in_path: str) -> Dict:
    with tempfile.TemporaryDirectory() as tmpdirname:
        output_dict = defaultdict(lambda: defaultdict(dict))
        chains1_list = get_chains(pdb1_in_path)
        chains2_list = get_chains(pdb2_in_path)
        for chain1 in chains1_list:
            name1 = f'{utils.get_file_name(pdb1_in_path)}_{chain1}'
            for chain2 in chains2_list:
                name2 = f'{utils.get_file_name(pdb2_in_path)}_{chain2}'
                output_path = os.path.join(tmpdirname, 'aux.pdb')
                run_gesamt(pdb_reference=pdb1_in_path, pdb_superposed=pdb2_in_path, output_path=output_path,
                           reference_chains=[chain1], superposed_chains=[chain2])
                for chain3 in chains2_list:
                    name3 = f'{chain3}'
                    rmsd, _, _ = run_gesamt(pdb_reference=pdb1_in_path, pdb_superposed=output_path,
                                            reference_chains=[chain1], superposed_chains=[chain3])
                    output_dict[name1][name2][name3] = rmsd
    return output_dict


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
        f_in.write(f'fit RESIDUE CA {match_ini} TO {match_end} CHAIN A \n')
        f_in.write(f'MATCH RESIDUE {fit_ini} TO {fit_end} CHAIN A \n')
        f_in.write(f'end \n')
        f_in.write(f'END-lsqkab')
    subprocess.Popen(['bash', script_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                     cwd=os.path.dirname(pdb_out)).communicate()


def run_spong(pdb_in_path: str, spong_path: str) -> float:
    # Run Spong and return compactness

    store_old_dir = os.getcwd()
    os.chdir(os.path.dirname(pdb_in_path))
    command_line = f'{spong_path} {os.path.basename(pdb_in_path)}'
    output = subprocess.Popen(command_line, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE).communicate()[0].decode('utf-8')

    compactness = None
    pattern = r"COMPACTNESS=\s+(\d+\.?\d*)"
    match = re.search(pattern, output)
    if match:
        compactness = match.group(1)
    os.chdir(store_old_dir)
    try:
        file_created = f'{os.path.join(os.path.dirname(pdb_in_path), utils.get_file_name(pdb_in_path))}new.ent'
        os.remove(file_created)
    except:
        pass
    return float(compactness) >= 2, compactness


def check_pdb(pdb: str, pdb_out_path: str) -> str:
    # Check if pdb is a path, and if it doesn't exist, download it.
    # If the pdb is a path, copy it to our input folder
    if not os.path.exists(pdb):
        pdb_path = download_pdb(pdb_id=pdb, pdb_path=pdb_out_path)
    else:
        pdb_path = shutil.copy2(pdb, pdb_out_path)
    return pdb_path


def copy_positions_of_pdb(path_in: str, path_out: str, positions: List[str]) -> str:
    new_structure = Structure.Structure(utils.get_file_name(path_in))
    new_model = Model.Model('model')
    chain = Chain.Chain('A')
    new_structure.add(new_model)
    new_model.add(chain)
    structure = get_structure(path_in)
    chain_name = get_chains(path_in)[0]
    for m, pos in enumerate(positions):
        if pos != '-':
            residue = copy.copy(structure[0][chain_name][int(pos) + 1])
            residue.parent = None
            residue.id = (residue.id[0], m + 1, residue.id[2])
            chain.add(residue)
            residue.parent = chain

    class AtomSelect(Select):
        def accept_atom(self, atom):
            return atom.get_name() in global_variables.ATOM_TYPES

    pdb_io = PDBIO()
    pdb_io.set_structure(new_structure)
    pdb_io.save(path_out, select=AtomSelect())

    return path_out


def check_sequence_path(path_in: str) -> str:
    if path_in is not None:
        if not os.path.exists(path_in):
            return path_in
        else:
            return extract_sequence(path_in)
    return None


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


def extract_sequence_msa_from_pdb(pdb_path: str) -> dict:
    structure = get_structure(pdb_path)
    model = structure[0]
    sequences = {}
    for chain in model:
        sequence_with_gaps = ""
        prev_residue_number = 0
        for residue in chain:
            residue_number = residue.get_id()[1]
            if residue_number - prev_residue_number > 1:
                sequence_with_gaps += "-" * (residue_number - prev_residue_number - 1)
            try:
                # Convert MSE to M
                if residue.get_resname() == 'MSE':
                    sequence_with_gaps += 'M'
                else:
                    sequence_with_gaps += residue_constants.restype_3to1[residue.get_resname()]
            except KeyError:
                pass
            prev_residue_number = residue_number
        sequences[chain.id] = sequence_with_gaps
    return sequences


def extract_sequence(fasta_path: str) -> str:
    logging.error(f'Extracting sequence from {fasta_path}')
    try:
        record = SeqIO.read(fasta_path, "fasta")
    except Exception as e:
        raise Exception(f'Not possible to extract the sequence from {fasta_path}')
    return str(record.seq)


def extract_sequences(fasta_path: str) -> Dict:
    logging.info(f'Extracting sequences from {fasta_path}')
    records = list(SeqIO.parse(fasta_path, 'fasta'))
    return dict([(rec.id, str(rec.seq)) for rec in records])


def read_seqres(pdb_path: str) -> List[str]:
    sequences = {}
    results_list = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('SEQRES'):
                fields = line.split()
                chain_id = fields[2]
                sequence_ext = [residue_constants.restype_3to1[code] if code != 'MSE' else 'M' for code in fields[4:]]
                if chain_id in sequences:
                    sequences[chain_id] += ''.join(sequence_ext)
                else:
                    sequences[chain_id] = ''.join(sequence_ext)
    for chain, sequence_ext in sequences.items():
        results = f'>{utils.get_file_name(pdb_path)[:10]}:{chain}\n'
        results += sequence_ext
        results_list.append(results)
    return results_list


def extract_sequence_from_file(file_path: str) -> dict:
    results_dict = {}
    extension = utils.get_file_extension(file_path)
    if extension == '.cif':
        extraction = 'cif-atom'
    else:
        extraction = 'pdb-atom'

    try:
        with open(file_path, 'r') as f_in:
            for record in SeqIO.parse(f_in, extraction):
                key = f'>{record.id.replace("????", utils.get_file_name(file_path)[:10])}'
                value = str(record.seq.replace("X", "-"))
                results_dict[key] = value
    except Exception as e:
        logging.info('Something went wrong extracting the fasta record from the pdb at', file_path)
        pass
    return results_dict


def write_sequence(sequence_name: str, sequence_amino: str, sequence_path: str) -> str:
    with open(sequence_path, 'w') as f_out:
        f_out.write(f'>{sequence_name}\n')
        f_out.write(f'{sequence_amino}')
    return sequence_path


def merge_pdbs(list_of_paths_of_pdbs_to_merge: List[str], merged_pdb_path: str):
    with open(merged_pdb_path, 'w+') as f:
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
    tmp_name = utils.generate_random_code(6)
    logging.info(f'Generating REMARK 350 for {pdb_path} with PISA.')
    subprocess.Popen(['pisa', tmp_name, '-analyse', pdb_path], stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE).communicate()
    pisa_output = \
        subprocess.Popen(['pisa', tmp_name, '-350'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0]
    erase_pisa(tmp_name)
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
    try:
        tmp_file = tempfile.NamedTemporaryFile(delete=False)
        with open(tmp_file.name, 'w+') as f:
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
        subprocess.Popen(['bash', tmp_file.name], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
    finally:
        tmp_file.close()
        os.unlink(tmp_file.name)


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


def get_number_residues(pdb_path: str) -> int:
    return len([res for res in Selection.unfold_entities(get_structure(pdb_path), 'R')])


def run_pdb2cc(templates_dir: str, pdb2cc_path: str = None) -> str:
    try:
        cwd = os.getcwd()
        os.chdir(templates_dir)
        output_path = 'cc_analysis.in'
        if pdb2cc_path is None:
            pdb2cc_path = 'pdb2cc'
        command_line = f'{pdb2cc_path} -m -i 10 -y 0.5 "orig.*.pdb" 0 {output_path}'
        p = subprocess.Popen(command_line, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        p.communicate()
    finally:
        os.chdir(cwd)
    return os.path.join(templates_dir, output_path)


def run_cc_analysis(input_path: str, n_clusters: int, cc_analysis_path: str = None) -> str:
    # Run cc_analysis from cc_analysis_path and store the results in output_path
    # The number cluster to obtain is determined by n_clusters.
    # It will read all the pdb files from input_path and store the results inside the output_path
    # I change the directory so everything is stored inside the input_path

    output_path = 'cc_analysis.out'
    if cc_analysis_path is None:
        cc_analysis_path = 'cc_analysis'
    try:
        cwd = os.getcwd()
        os.chdir(os.path.dirname(input_path))
        command_line = f'{cc_analysis_path} -dim {n_clusters} {os.path.basename(input_path)} {output_path}'
        p = subprocess.Popen(command_line, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        p.communicate()
    finally:
        os.chdir(cwd)
    return f'{os.path.join(os.path.dirname(input_path), output_path)}'


def run_hinges(pdb1_path: str, pdb2_path: str, hinges_path: str = None, output_path: str = None) -> structures.Hinges:
    # Run hinges from hinges_path.
    # It needs two pdbs. Return the rmsd obtained.
    chains2_list = get_chains(pdb2_path)
    command_line = f'{hinges_path} {pdb1_path} {pdb2_path} -p {"".join(chains2_list)}'
    output = subprocess.Popen(command_line, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE).communicate()[0].decode('utf-8')
    best_chain_combination = utils.parse_hinges_chains(output)
    append = ''
    if best_chain_combination != '':
        append = ''
        for i, chain in enumerate(chains2_list):
            append += f'{chain}:{best_chain_combination[i]} '
        append = f'-r "{append}"'
    command_line = f'{hinges_path} {pdb1_path} {pdb2_path} {append} -p'
    output = subprocess.Popen(command_line, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE).communicate()[0].decode('utf-8')
    if output_path is not None:
        with open(output_path, 'w+') as f:
            f.write(output)
    return utils.parse_hinges(output)


def generate_ramachandran(pdb_path, output_dir: str = None) -> bool:
    # Ramachandran analysis, it generates the angles in degrees, and it generates the plot.
    valid_residues = ["MET", "SER", "ASN", "LEU", "GLU", "LYS", "GLN", "ILE", "ALA", "ARG",
                      "HIS", "CYS", "ASP", "THR", "GLY", "TRP", "PHE", "TYR", "PRO", "VAL"]

    # Create ramachandran plot and return analysis
    structure = get_structure(pdb_path)
    phi_angles = np.empty(0)
    psi_angles = np.empty(0)
    percentage_minimum = 5

    # Iterate over each polypeptide in the structure
    for pp in PPBuilder().build_peptides(structure):
        phi_psi = pp.get_phi_psi_list()
        for i, residue in enumerate(pp):
            if residue.resname not in valid_residues:
                continue
            # Get the phi and psi angles for the residue
            phi, psi = phi_psi[i]
            if phi and psi:
                phi_angles = np.append(phi_angles, np.degrees(phi))
                psi_angles = np.append(psi_angles, np.degrees(psi))

    phi_psi_angles = np.column_stack((phi_angles, psi_angles))

    if output_dir is not None:
        plots.plot_ramachandran(plot_path=os.path.join(output_dir, f'{utils.get_file_name(pdb_path)}.png'),
                                phi_psi_angles=phi_psi_angles)

    analysis = ramachandran_analysis(phi_psi_angles=phi_psi_angles)
    if len(phi_psi_angles) > 0:
        percentage = len(analysis) / len(phi_psi_angles) * 100
    else:
        percentage = 0
    logging.info(
        f'{round(percentage, 2)}% of outliers in the ramachandran analysis of {utils.get_file_name(pdb_path)}.')
    if percentage > percentage_minimum:
        return False, percentage
    return True, percentage


def ramachandran_analysis(phi_psi_angles: List[List[int]]) -> List[int]:
    # Do the ramachandran analysis. Given a matrix of PHI and PSI (X,Y), calculate the outliers given
    # a table from global_variables. If the value is different from 0, consider it not an outlier.
    outliers_list = []
    minimum_value = 1
    for phi, psi in phi_psi_angles:
        value = global_variables.RAMACHANDRAN_TABLE[int((psi + 180) / 10)][int((phi + 180) / 10)]
        if value < minimum_value:
            outliers_list.append(value)
    return outliers_list


def aleph_annotate(output_path: str, pdb_path: str) -> Union[None, Dict]:
    # Run aleph annotate. Given a path, it generates the annotation of the pdb (coil, bs and ah).
    # Also, it generates the domains.
    try:
        store_old_dir = os.getcwd()
        os.chdir(output_path)
        aleph_output_txt = os.path.join(output_path, f'aleph_{utils.get_file_name(pdb_path)}.txt')
        output_json = os.path.join(output_path, 'output.json')
        with open(aleph_output_txt, 'w') as sys.stdout:
            try:
                ALEPH.annotate_pdb_model(reference=pdb_path, strictness_ah=0.45, strictness_bs=0.2,
                                         peptide_length=3, width_pic=1, height_pic=1, write_graphml=False,
                                         write_pdb=True)
            except Exception as e:
                pass
        sys.stdout = sys.__stdout__
        if os.path.exists(output_json):
            return utils.parse_aleph_annotate(output_json), utils.parse_aleph_ss(aleph_output_txt)
        else:
            return None, None
    finally:
        os.chdir(store_old_dir)


def cc_and_hinges_analysis(pdbs: List[structures.Pdb], binaries_path: structures.BinariesPath, output_dir: str) -> List:
    templates_cluster2 = []
    templates_cluster = hinges(pdbs=pdbs,
                               binaries_path=binaries_path,
                               output_dir=os.path.join(output_dir, 'hinges'))

    pdbs_accepted_list = [template_in for template_list in templates_cluster for template_in in template_list]
    num_templates = len(pdbs_accepted_list)


    if num_templates >= 5:
        logging.info(
            f'Running ccanalysis with the following templates: {" ".join([pdb.name for pdb in pdbs_accepted_list])}')
        templates_cluster2, analysis_dict2 = cc_analysis(pdbs=pdbs_accepted_list,
                                                         cc_analysis_paths=binaries_path,
                                                         output_dir=os.path.join(output_dir, 'ccanalysis'))
    else:
        logging.info('Less than 5 templates recognised by hinges. Skipping ccanalysis.')

    # if len(templates_cluster) > 1 and templates_cluster2:
    if templates_cluster2:
        return templates_cluster2, analysis_dict2
    else:
        return templates_cluster, {}


def hinges(pdbs: List[structures.Pdb], binaries_path: structures.BinariesPath, output_dir: str) -> List:
    # Hinges algorithm does:
    # Check completeness and ramachandran of every template. If it is not at least 70% discard for hinges.
    # Do hinges 6 iterations in all for all the templates
    # If iter1 < 1 or iterMiddle < 4.5 or iter6 < 8 AND at least 70% of sequence length -> GROUP TEMPLATE
    # If it has generated more than one group with length > 1-> Return those groups that has generated
    # If there is no group generated -> Return the one more completed and the one more different to that template
    # Otherwise, return all the paths
    utils.create_dir(output_dir, delete_if_exists=True)
    threshold_completeness = 0.6
    threshold_completeness2 = 0.3
    threshold_rmsd_domains = 8
    threshold_rmsd_ss = 4.5
    threshold_rmsd_local = 1.5
    threshold_overlap = 0.7
    threshold_minimum = 3
    threshold_decrease = 40
    threshold_identity_upper = 90
    threshold_identity_down = 20

    logging.info('Starting hinges analysis')
    accepted_pdbs = []
    uncompleted_pdbs = []
    completed_pdbs = []

    # Do the analysis of the different templates. We are going to check:
    # Completeness respect the query size sequence
    # Ramachandran plot
    # And the compactness
    pdb_complete = ''
    pdb_complete_value = 0
    for pdb in pdbs:
        num_residues = sum(1 for _ in get_structure(pdb.split_path)[0].get_residues())
        # Validate using ramachandran, check the outliers
        validate_geometry, _ = generate_ramachandran(pdb_path=pdb.split_path, output_dir=output_dir)
        # Check the query sequence vs the number of residues of the pdb
        only_ca = check_not_only_CA(pdb_in_path=pdb.split_path)
        completeness = True
        identity = True
        if isinstance(pdb, structures.TemplateExtracted) and pdb.percentage_list:
            completeness = any(number > threshold_completeness for number in pdb.percentage_list)
            if pdb.identity > threshold_identity_upper or pdb.identity < threshold_identity_down:
                identity = False
        compactness_decision, _ = run_spong(pdb_in_path=pdb.split_path, spong_path=binaries_path.spong_path)
        if completeness and validate_geometry and compactness_decision and not only_ca and identity:
            accepted_pdbs.append(pdb)
            logging.info(f'PDB {pdb.name} has been accepted')
            if num_residues > pdb_complete_value:
                pdb_complete_value = num_residues
                pdb_complete = pdb
        else:
            uncompleted_pdbs.append(pdb)
            logging.info(f'PDB {pdb.name} has been filtered:')
            if not completeness:
                logging.info(f'    Not complete enough')
            if not validate_geometry:
                logging.info(f'    Ramachandran above limit')
            if not compactness_decision:
                logging.info(f'    Compactness below limit')
            if not identity:
                logging.info(f'    Too low/high identity with the query sequence ({pdb.identity})')

            if only_ca:
                logging.info(f'    Only CA')
        if isinstance(pdb, structures.TemplateExtracted) and pdb.percentage_list:
            if any(number > threshold_completeness2 for number in pdb.percentage_list):
                completed_pdbs.append(pdb)
        else:
            completed_pdbs.append(pdb)

    logging.info(f'There are {len(accepted_pdbs)} complete pdbs.')
    if len(accepted_pdbs) < 2:
        logging.info(f'Skipping hinges.')
        return [completed_pdbs]
    logging.info(f'Using hinges to create groups.')

    # Run hinges all-against-all, store the results in a dict.
    results_rmsd = {pdb.name: {} for pdb in accepted_pdbs}
    groups_names = {pdb.name: [] for pdb in accepted_pdbs}
    for pdb1 in accepted_pdbs:
        for pdb2 in accepted_pdbs:
            if pdb2.name not in results_rmsd[pdb1.name]:
                result_hinges = run_hinges(pdb1_path=pdb1.split_path, pdb2_path=pdb2.split_path,
                                           hinges_path=binaries_path.hinges_path,
                                           output_path=os.path.join(output_dir, f'{pdb1.name}_{pdb2.name}.txt'))
                results_rmsd[pdb1.name][pdb2.name] = result_hinges
                results_rmsd[pdb2.name][pdb1.name] = result_hinges
    results_rmsd = OrderedDict(sorted(results_rmsd.items(), key=lambda x: min(v.one_rmsd for v in x[1].values())))
    sorted_list = copy.deepcopy(results_rmsd)
    for key1, value in sorted_list.items():
        sorted_list[key1] = OrderedDict(
            sorted({k: v for k, v in value.items() if v is not None}.items(), key=lambda x: x[1].one_rmsd))
        selected_group = key1
        for key2, result in sorted_list[key1].items():
            if key1 != key2:
                group = utils.get_key_by_value(key2, groups_names)
                selected_for = None
                if group and result.overlap > threshold_overlap:
                    if result.one_rmsd <= threshold_rmsd_local:
                        selected_for = 'local changes'
                    elif result.middle_rmsd <= threshold_rmsd_ss and result.decreasing_rmsd_middle >= threshold_decrease:
                        selected_for = 'secondary structure changes'
                    elif result.min_rmsd <= threshold_rmsd_domains and result.decreasing_rmsd_total >= threshold_decrease:
                        selected_for = 'domain changes'
                    elif result.min_rmsd <= threshold_minimum:
                        selected_for = 'equivalence'
                    if selected_for is not None:
                        if selected_group not in groups_names or len(groups_names[group[0]]) > len(
                                groups_names[selected_group]):
                            selected_group = group[0]
                        logging.info(f'{key2} into group {selected_group} because of {selected_for} with {key1}')

        groups_names[selected_group].append(key1)
    groups_names = [values for values in groups_names.values() if len(values) > 1]

    tables_path = os.path.join(output_dir, 'tables.txt')

    with open(tables_path, 'w') as f_in:
        data1 = {'ranked': results_rmsd.keys()}
        data2 = {'ranked': results_rmsd.keys()}
        data3 = {'ranked': results_rmsd.keys()}
        for ranked in results_rmsd.values():
            for key, value in ranked.items():
                data1.setdefault(key, []).append(value.one_rmsd)
                data2.setdefault(key, []).append(value.middle_rmsd)
                data3.setdefault(key, []).append(value.min_rmsd)
        df = pd.DataFrame(data1)
        f_in.write('\n\n')
        f_in.write('Table first iteration\n')
        f_in.write(df.to_markdown())
        df = pd.DataFrame(data2)
        f_in.write('\n\n')
        f_in.write('Table middle iteration\n')
        f_in.write(df.to_markdown())
        df = pd.DataFrame(data3)
        f_in.write('\n\n')
        f_in.write('Table last iteration\n')
        f_in.write(df.to_markdown())

    if len(groups_names) > 1 or (len(groups_names) == 1 and len(groups_names[0]) > 1):
        # Return the groups that has generated
        logging.info(f'Hinges has created {len(groups_names)} group/s:')
        for i, values in enumerate(groups_names):
            logging.info(f'Group {i}: {",".join(values)}')
        return [[pdb for pdb in pdbs if pdb.name in group] for group in groups_names]
    elif len(list(results_rmsd.keys())) > 1:
        # Create two groups, more different and completes pdbs
        more_different = list(results_rmsd[pdb_complete.name].keys())[-1]
        pdb_diff = [pdb for pdb in pdbs if pdb.name == more_different][0]
        logging.info(f'Hinges could not create any groups')
        logging.info(f'Creating two groups: The more completed pdb: {pdb_complete.name} '
                     f'and the more different one: {pdb_diff.name}')

        return [[pdb_complete], [pdb_diff]]
    else:
        # Return the original list of pdbs
        logging.info('Not enough pdbs for hinges.')
        return [completed_pdbs]


def cc_analysis(pdbs: List[structures.Pdb], cc_analysis_paths: structures.BinariesPath, output_dir: str,
                n_clusters: int = 2) -> List:
    # CC_analysis. It is mandatory to have the paths of the programs to run pdb2cc and ccanalysis.
    # A dictionary with the different pdbs that are going to be analysed.

    utils.create_dir(output_dir, delete_if_exists=True)
    trans_dict = {}
    return_templates_cluster = [[] for _ in range(n_clusters)]
    clean_dict = {}

    for index, pdb in enumerate(pdbs):
        path = shutil.copy2(pdb.path, output_dir)
        # If it is ranked, it is mandatory to change the bfactors to VALUE-70.
        # We want to evaluate the residues that have a good PLDDT
        # PDB2CC ignore the residues with bfactors below 0
        if utils.check_ranked(os.path.basename(path)):
            bfactors_dict = read_bfactors_from_residues(path)
            for chain, residues in bfactors_dict.items():
                for i in range(len(residues)):
                    if bfactors_dict[chain][i] is not None:
                        bfactors_dict[chain][i] = round(bfactors_dict[chain][i] - 70.0, 2)

            modify_bfactors = template_modifications.TemplateModifications()
            modify_bfactors.append_modification(chains=list(bfactors_dict.keys()),
                                                bfactors=list(bfactors_dict.values()))
            modify_bfactors.modify_template(pdb_in_path=path, pdb_out_path=path, type_modify=['bfactors'])

        new_path = os.path.join(output_dir, f'orig.{str(index)}.pdb')
        os.rename(os.path.join(output_dir, path), new_path)
        # Create a translation dictionary, with and index and the pdb name
        trans_dict[index] = utils.get_file_name(path)
    if trans_dict:
        # Write the trans dict in order to be able to trace the pdbs in the output
        with open(os.path.join(output_dir, 'labels.txt'), 'w+') as f:
            for key, value in trans_dict.items():
                f.write('%s:%s\n' % (key, value))
        # run pdb2cc
        output_pdb2cc = run_pdb2cc(templates_dir=output_dir, pdb2cc_path=cc_analysis_paths.pd2cc_path)
        if os.path.exists(output_pdb2cc):
            # If pdb2cc has worked, launch cc analysis.
            output_cc = run_cc_analysis(input_path=output_pdb2cc,
                                        n_clusters=n_clusters,
                                        cc_analysis_path=cc_analysis_paths.cc_analysis_path)
            if os.path.exists(output_cc):
                # Parse the results of cc_analysis, we will have for each pdb, all the values given in ccanalysis
                cc_analysis_dict = utils.parse_cc_analysis(file_path=output_cc)
                for key, values in cc_analysis_dict.items():
                    # If the modules are higher than 0.1, keep the pdb, otherwise discard it
                    if values.module is None or values.module > 0.1 or values.module < -0.1:
                        clean_dict[trans_dict[int(key) - 1]] = values
                if clean_dict:
                    # Get the positions given by ccanalysis
                    points = np.array([values.coord for values in clean_dict.values()])
                    # Generate n clusters groups with KMEANS
                    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(points)
                    lookup_table = {}
                    counter = 0
                    # The kmeans results have a list with all the positions belonging to the corresponding pdb.
                    # Translate the labels into groups, so they appear in sequentially (group 0 first, 1...)
                    # Otherwise they are chosen randomly.
                    for i in kmeans.labels_:
                        if i not in lookup_table:
                            lookup_table[i] = counter
                            counter += 1
                    # Replace the values for the ordered ones.
                    conversion = [lookup_table[label] for label in kmeans.labels_]
                    # Translate the kmeans, which only has the position of the pdbs for the real pdbs.
                    for i, label in enumerate(conversion):
                        selected_pdb = [pdb for pdb in pdbs if pdb.name == list(clean_dict.keys())[i]][
                            0]
                        return_templates_cluster[int(label)].append(selected_pdb)

    # Return the clusters, a list, where each position has a group of pdbs.
    # Also, clean_dict has the cc_analysis vectors, so is useful to create the plots.
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


def convert_residues(residues_list: List[List], sequence_assembled):
    # Given a list of list, which each one corresponding a position in the query sequence
    # Return the real position of that residue before splitting, as the residues of the chain
    # are already separated by chains
    for i in range(0, len(residues_list)):
        if residues_list[i] is not None:
            for residue in residues_list[i]:
                result = sequence_assembled.get_real_residue_number(i, residue)
                if result is not None:
                    residues_list.append(result)
    return residues_list


def get_group(res: str) -> str:
    # Given a residue letter, return if that letter belongs to any group.
    groups = ['GAVLI', 'FYW', 'CM', 'ST', 'KRH', 'DENQ', 'P']
    group = [s for s in groups if res in s]
    if group:
        return group[0]
    return res


def compare_sequences(sequence1: str, sequence2: str, only_match=False) -> List[str]:
    # Given two sequences with same length, return a list showing
    # if there is a match, a group match, they are different, or
    # they are not aligned
    # Also, return the changes in the sequence2
    gap = 0
    if only_match:
        match = 1
        group_match = 0
        mismatch = 0
    else:
        match = 6
        group_match = 4
        mismatch = 2
    return_list = []
    changes_dict = {}
    for i, (res1, res2) in enumerate(itertools.zip_longest(sequence1, sequence2, fillvalue='-')):
        if res1 == '-' or res2 == '-':
            return_list.append(gap)
        elif res1 == res2:
            return_list.append(match)
        elif get_group(res1) == get_group(res2):
            return_list.append(group_match)
        else:
            return_list.append(mismatch)
        if res1 != res2:
            changes_dict[i] = res2

    return return_list, changes_dict


def sequence_identity_regions(seq1, seq2, regions_list: List):
    identity = 0
    num_amino = 0
    amino_region1 = ['X'] * len(seq1)
    amino_region2 = ['X'] * len(seq2)
    for region in regions_list:
        min_i = min(region[1], len(seq1), len(seq2))
        for i in range(region[0] - 1, min_i):
            identity += 1 if seq1[i] == seq2[i] and seq1[i] != '-' and seq2[i] != '-' else 0
            num_amino += 1
            amino_region1[i] = seq1[i]
            amino_region2[i] = seq2[i]
    return (identity / num_amino) * 100, ''.join(map(str, amino_region1)), ''.join(map(str, amino_region2))


def sequence_identity(seq1, seq2) -> float:
    # Compare the identity of two sequences
    identical_count = sum(1 for a, b in zip(seq1, seq2) if a != '-' and b != '-' and a == b)
    identity = (identical_count / len(seq1)) * 100
    return identity


def convert_msa_sequence(number_list: List[int]) -> str:
    return ''.join([residue_constants.ID_TO_HHBLITS_AA[res] for res in number_list])


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
    chains = list(set(get_chains(pdb_in_path)))

    if len(chains) > 1:
        logging.info(f'PDB: {pdb_in_path} is already split in several chains: {chains}')
        try:
            shutil.copy2(pdb_in_path, pdb_out_path)
        except shutil.SameFileError:
            pass
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


def split_pdb_in_chains(pdb_path: str, chain: str = None, output_dir: str = None) -> Dict:
    # Given a pdb_in and an optional chain, write one or several
    # pdbs containing each one a chain.
    # If chain is specified, only one file with the specific chain will be created
    # It will return a dictionary with the chain and the corresponding pdb

    return_chain_dict = {}
    structure = get_structure(pdb_path=pdb_path)
    chains = get_chains(pdb_path) if chain is None else [chain]

    if output_dir is None:
        output_dir = os.path.dirname(pdb_path)

    for chain in chains:
        new_pdb = os.path.join(output_dir, f'{utils.get_file_name(pdb_path)}_{chain}1.pdb')

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
    try:
        shutil.copy2(pdb_in_path, pdb_out_path)
    except:
        pass
    chain_dict = split_pdb_in_chains(pdb_path=pdb_out_path)
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
        if chain in template_dict.keys():
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
    for res in structure[0].get_residues():
        if get_hetatm(res) == 'H_MSE' or res.resname == 'MSE':
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


def remove_hydrogens(pdb_in_path: str, pdb_out_path: str):
    # Remove the atoms that don't belong to the list atom_types
    structure = get_structure(pdb_path=pdb_in_path)
    # Remove hydrogen atoms from the structure
    for residue in structure[0].get_residues():
        atoms = residue.get_unpacked_list()
        for atom in atoms:
            if atom.element not in ['N', 'C', 'O', 'S'] and atom in residue.get_unpacked_list():
                residue.detach_child(atom.get_id())

    # Save the edited structure to a new PDB file
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_out_path)


def run_pdbfixer(pdb_in_path: str, pdb_out_path: str):
    try:
        pdb_text = open(pdb_in_path, 'r').read()
        pdb_file = io.StringIO(pdb_text)
        pdb_output = cleanup.fix_pdb(pdb_file, {})
        with open(pdb_out_path, 'w') as f_out:
            f_out.write(pdb_output)
    except:
        logging.info(f'PDBFixer did not finish correctly for {utils.get_file_name(pdb_in_path)}. Skipping.')
        shutil.copy2(pdb_in_path, pdb_out_path)
        pass


def run_vairo(yml_path: str, input_path: str):
    vairo_path = os.path.join(utils.get_main_path(), 'run_vairo.py')
    command_line = f'{vairo_path} {yml_path}'
    p = subprocess.Popen(command_line, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if not utils.read_rankeds(input_path=input_path):
        raise Exception(
            f'VAIRO cluster did not finish correctly. Please check the log in the following directory: {input_path}')
    else:
        logging.error('VAIRO cluster run finished successfully.')


def run_openmm(pdb_in_path: str, pdb_out_path: str) -> float:
    energy = unit.kilocalories_per_mole
    length = unit.angstroms
    restraint_set = "non_hydrogen"
    max_iterations = 1
    tolerance = 2.39 * energy
    stiffness = 10.0 * energy / (length ** 2)
    run_pdbfixer(pdb_in_path=pdb_in_path, pdb_out_path=pdb_out_path)
    pdb_text = open(pdb_out_path, 'r').read()
    ret = amber_minimize._openmm_minimize(
        pdb_str=pdb_text,
        max_iterations=max_iterations,
        tolerance=tolerance,
        stiffness=stiffness,
        exclude_residues=[],
        restraint_set=restraint_set,
        use_gpu=True)
    with open(pdb_out_path, 'w+') as f:
        f.write(ret["min_pdb"])
    return round(ret["efinal"], 2)


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


def run_gesamt(pdb_reference: str, pdb_superposed: str, output_path: str = None, reference_chains: List[str] = [],
               superposed_chains: List[str] = []) -> Tuple[
    Optional[float], Optional[str], Optional[str]]:
    with tempfile.TemporaryDirectory() as tmpdirname:
        superpose_cmd = 'gesamt'
        logging_text = f'Superposing {pdb_reference} with {pdb_superposed}.'
        superpose_cmd += f' {pdb_reference}'
        if reference_chains:
            logging_text += f' Reference chains: {", ".join(reference_chains)}.'
            superpose_cmd += f' -s {",".join(reference_chains)}'
        superpose_cmd += f' {pdb_superposed}'
        if superposed_chains:
            f' Superposed chains: {", ".join(reference_chains)}.'
            superpose_cmd += f' -s {",".join(superposed_chains)}'
        if output_path is not None:
            superpose_cmd += f' -o {tmpdirname} -o-d'
        logging.info(logging_text)
        superpose_output = subprocess.Popen(superpose_cmd, stdout=subprocess.PIPE, shell=True).communicate()[0].decode(
            'utf-8')
        new_path = os.path.join(tmpdirname, f'{utils.get_file_name(pdb_superposed)}_2.pdb')
        if os.path.exists(new_path) and output_path:
            shutil.copy2(new_path, output_path)
        elif not os.path.exists(new_path) and output_path:
            logging.info(f'Not possible to superpose {pdb_reference} with {pdb_superposed}')
        rmsd, qscore, nalign = None, None, None
        for line in superpose_output.split('\n'):
            if 'RMSD             :' in line:
                rmsd = float(line.split()[2].strip())
            if 'Q-score          :' in line:
                qscore = float(line.split()[2].strip())
            if 'Aligned residues :' in line:
                nalign = int(line.split()[3].strip())
        return rmsd, nalign, qscore


def gesamt_pdbs(pdb_reference: str, pdb_superposed: str, output_path: str = None, check_chains: bool = True) -> Tuple[
    Optional[float], Optional[str], Optional[str]]:
    chains_superposed = []
    if check_chains:
        chains_r = get_chains(pdb_reference)
        chains_s = get_chains(pdb_superposed)
        if all(chain in chains_r for chain in chains_s):
            chains_superposed = chains_s

    rmsd, nalign, qscore = run_gesamt(pdb_reference=pdb_reference, pdb_superposed=pdb_superposed,
                                      output_path=output_path, reference_chains=chains_superposed, superposed_chains=[])
    return rmsd, nalign, qscore


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
    tmp_name = utils.generate_random_code(6)

    pisa_text = subprocess.Popen(['pisa', tmp_name, '-analyse', pdb_in_path],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE).communicate()[0].decode('utf-8')
    pisa_output = subprocess.Popen(['pisa', tmp_name, '-list', 'interfaces'], stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE).communicate()[0].decode('utf-8')

    pisa_general_txt = os.path.join(interfaces_path, f'{utils.get_file_name(pdb_in_path)}_general_output.txt')
    with open(pisa_general_txt, 'w') as f_out:
        f_out.write(pisa_output)
    if pisa_output == '' or 'NO INTERFACES FOUND' in pisa_output or 'no chains found in input file' in pisa_text:
        logging.info(f'No interfaces found in pisa for pdb {pdb_in_path}')
    else:
        interfaces_list = utils.parse_pisa_general_multimer(pisa_output)
        for interface in interfaces_list:
            serial_output = \
                subprocess.Popen(['pisa', tmp_name, '-detail', 'interfaces', interface['serial']],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE).communicate()[0].decode('utf-8')
            interface_data = utils.parse_pisa_interfaces(serial_output)
            new_interface = structures.Interface(name=f'{interface_data["chain1"]}-{interface_data["chain2"]}',
                                                 res_chain1=interface_data["res_chain1"],
                                                 res_chain2=interface_data["res_chain2"],
                                                 chain1=interface_data["chain1"],
                                                 chain2=interface_data["chain2"],
                                                 se_gain1=float(interface_data['se_gain1']),
                                                 se_gain2=float(interface_data['se_gain2']),
                                                 solvation1=float(interface_data['solvation1']),
                                                 solvation2=float(interface_data['solvation2']),
                                                 area=float(interface['area']),
                                                 deltaG=float(interface['deltaG']),
                                                 nhb=int(interface['nhb']),
                                                 )

            interface_data_list.append(new_interface)
            pisa_output_txt = os.path.join(interfaces_path,
                                           f'{utils.get_file_name(pdb_in_path)}_{interface_data["chain1"]}{interface_data["chain2"]}_interface.txt')
            with open(pisa_output_txt, 'w') as f_out:
                f_out.write(serial_output)

    erase_pisa(name=tmp_name)

    return interface_data_list


def parse_pdb_hits_hhr(hhr_text: str, pdb_name: str) -> Dict:
    pattern = rf'>{re.escape(pdb_name)}(.*)'
    match = re.search(pattern, hhr_text, re.MULTILINE | re.DOTALL)
    if match:
        protein_info = match.group(1).strip()
        e_value = re.search(r'E-value=(\S+)', protein_info).group(1)
        aligned_cols = re.search(r'Aligned_cols=(\d+)', protein_info).group(1)
        identity = re.search(r'Identities=(\d+)%', protein_info).group(1)
        total_residues = re.search(r'\((\d+)\)', protein_info).group(1)
        return e_value, aligned_cols, identity, total_residues
    else:
        return None, None, None, None


def check_not_only_CA(pdb_in_path: str) -> bool:
    structure = get_structure(pdb_in_path)
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_name() != 'CA':
                        return False
    return True


def create_interface_domain(pdb_in_path: str, pdb_out_path: str, interface: Dict, domains_dict: Dict) \
        -> Dict[Any, List[Any]]:
    # For a PDB and an interface (contains the chains and the residues involved in the interface).
    # Iterate the interface, selecting all those residues that belong to a domain (we have previously calculated
    # all the domains in the PDB).
    # We create a dictionary with all the residues of the interface, extending them to their whole domain.
    # Split the pdb into the chains of the interface and delete all the others residues.
    # Return the dictionary with the chains and the interface residues extended.
    add_domains_dict = {}
    for chain, residue in zip([interface.chain1, interface.chain2],
                              [interface.res_chain1, interface.res_chain2]):
        added_res_list = []
        [added_res_list.extend(domains) for domains in domains_dict[chain] if bool(set(residue).intersection(domains))]
        added_res_list.extend(residue)
        add_domains_dict[chain] = list(set(added_res_list))

    split_dimers_in_pdb(pdb_in_path=pdb_in_path,
                        pdb_out_path=pdb_out_path,
                        chain_list=[interface.chain1, interface.chain2])

    change = template_modifications.TemplateModifications()
    change.append_modification(chains=list(add_domains_dict.keys()), maintain_residues=list(add_domains_dict.values()))
    change.modify_template(pdb_in_path=pdb_out_path, pdb_out_path=pdb_out_path, type_modify=['delete'])

    return add_domains_dict



def calculate_auto_offset(input_list: List[List], length: int) -> List[int]:
    if length <= 0:
        return []
    trimmed_list = []
    for element in itertools.product(*input_list):
        if not element:
            continue

        sorted_list = sorted(element, key=lambda x: x[2])
        aux_list, x_list, y_list = [], set(), set()
        for tup in sorted_list:
            if tup[0] not in x_list and tup[1] not in y_list:
                aux_list.append(tup)
                x_list.add(tup[0])
                y_list.add(tup[1])
            if len(aux_list) == min(length, len(element)):
                break
        trimmed_list.append(aux_list)

    if not trimmed_list:
        return []

    max_length = max(len(lst) for lst in trimmed_list)
    trimmed_list = [lst for lst in trimmed_list if len(lst) == max_length]
    score_list = [sum(z for _, _, z, _, _ in element) for element in trimmed_list]
    if not score_list:
        return []
    min_score_index = min(range(len(score_list)), key=score_list.__getitem__)
    return trimmed_list[min_score_index]


def split_dimers_in_pdb(pdb_in_path: str, pdb_out_path: str, chain_list: List[str]):
    # Given a PDB, keep those chains that are in the list.
    # The other chains are going to be deleted.
    class ChainSelector(Select):
        def accept_chain(self, chain):
            if chain.get_id() in chain_list:
                return True
            else:
                return False

    structure = get_structure(pdb_in_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_out_path, ChainSelector())


def align_pdb(pdb_in_path: str, pdb_out_path: str, sequences_list: List[str],
              databases):
    with tempfile.TemporaryDirectory() as tmpdirname:
        chain_dict = split_pdb_in_chains(pdb_path=pdb_in_path, output_dir=tmpdirname)
        chains_aligned = []
        if len(chain_dict) != len(sequences_list):
            return None
        else:
            for i, path in enumerate(chain_dict.values()):
                aligned_chain, _ = hhsearch.run_hh(output_dir=tmpdirname, database_dir=tmpdirname,
                                                   query_sequence_path=sequences_list[i],
                                                   chain_in_path=path, databases=databases)
                shutil.copy2(aligned_chain, path)
                chains_aligned.append(path)
            merge_pdbs(list_of_paths_of_pdbs_to_merge=chains_aligned, merged_pdb_path=pdb_out_path)
            return pdb_out_path


def conservation_pdb(pdb_in_path: str, pdb_out_path: str, msa_list: List[str]):
    # Change the bfactors with the conservation found in the msa. Just one chain, as the msa
    # it is in a chain
    sequences_dict = extract_sequence_msa_from_pdb(pdb_path=pdb_in_path)
    chain = get_chains(pdb_in_path)[0]
    whole_seq = "".join([seq for seq in sequences_dict.values()])
    whole_seq += '-' * (len(msa_list[0]) - len(whole_seq))
    conservation_list = calculate_coverage(query_seq=whole_seq, sequences=msa_list, only_match=True)
    conservation_list = conservation_list * 100
    modify_bfactors = template_modifications.TemplateModifications()
    modify_bfactors.append_modification(chains=[chain], bfactors=conservation_list.tolist())
    modify_bfactors.modify_template(pdb_in_path=pdb_in_path, pdb_out_path=pdb_out_path, type_modify=['bfactors'])


def calculate_coverage(query_seq: str, sequences: List[str], only_match: bool) -> List[str]:
    # Coverage of the sequences. It is divided by the number of sequences.
    add_sequences = np.zeros(len(query_seq))
    for seq in sequences:
        aligned_sequence, _ = compare_sequences(query_seq, seq, only_match)
        add_sequences += np.array(aligned_sequence)
    add_sequences /= len(sequences)
    return add_sequences


def calculate_coverage_scaled(query_seq: str, sequences: List[str]):
    sequences_coverage = calculate_coverage(query_seq=query_seq, sequences=sequences, only_match=False)
    new_sequences = utils.scale_values(sequences_coverage)
    return new_sequences


def shift_pdb(pdb_in_path: str, sequence_predicted_assembled, sequence_assembled):
    # Given a list of shifts, one for each chain, apply them to each chain and return de pdb.
    shifts = sequence_predicted_assembled.get_region_starting_shifts()
    new_structure = Structure.Structure(utils.get_file_name(pdb_in_path))
    new_model = Model.Model('model')
    new_chain = Chain.Chain('A')
    new_structure.add(new_model)
    new_model.add(new_chain)
    structure = get_structure(pdb_in_path)
    chain = next(structure.get_chains())
    residues = chain.get_unpacked_list()
    for residue in residues:
        seq_position = sequence_predicted_assembled.get_position_by_residue_number(residue.id[1])
        if seq_position is not None:
            starting_position_original = sequence_assembled.get_starting_length(seq_position)
            starting_position_predicted = sequence_predicted_assembled.get_starting_length(seq_position)
            residue = copy.copy(residue)
            residue.parent = None
            new_id = starting_position_original + shifts[seq_position] - starting_position_predicted + residue.id[1]
            residue.id = (residue.id[0], new_id, residue.id[2])
            new_chain.add(residue)
            residue.parent = new_chain

    io = PDBIO() 
    io.set_structure(new_structure) 
    io.save(pdb_in_path)
