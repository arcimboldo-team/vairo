import copy, itertools, logging, os, re, shutil, subprocess, sys, tempfile
from collections import OrderedDict
from typing import Any, Dict, List, Optional, Tuple, Union
import numpy as np
from Bio import SeqIO
from Bio.PDB import PDBIO, PDBList, PDBParser, Residue, Chain, Select, Selection, Structure, Model, PPBuilder
from scipy.spatial import distance
from simtk import unit, openmm
from sklearn.cluster import KMeans

from ALEPH.aleph.core import ALEPH
from alphafold.common import residue_constants
from libs import change_res, structures, utils, plots, global_variables, sequence
from alphafold.relax import cleanup
from libs.structures import Hinges


def download_pdb(pdb_id: str, pdb_path: str) -> str:
    pdbl = PDBList(server='https://files.wwpdb.org')
    result_ent = pdbl.retrieve_pdb_file(pdb_code=pdb_id, file_format='pdb', obsolete=False)
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


def check_pdb(pdb: str, pdb_out_path: str) -> str:
    # Check if pdb is a path, and if it doesn't exist, download it.
    # If the pdb is a path, copy it to our input folder
    if not os.path.exists(pdb):
        pdb_path = download_pdb(pdb_id=pdb, pdb_path=pdb_out_path)
    else:
        pdb_path = shutil.copy2(pdb, pdb_out_path)
    return pdb_path


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


def extract_sequence_msa_from_pdb(pdb_path: str) -> str:
    structure = get_structure(pdb_path)
    model = structure[0]
    sequences = {}
    for chain in model:
        residue_numbers = set()
        for residue in chain:
            residue_numbers.add(residue.get_id()[1])
        sequence_ext = ""
        prev_residue_number = 0
        for residue in chain:
            residue_number = residue.get_id()[1]
            if residue_number - prev_residue_number > 1:
                for missing_number in range(prev_residue_number + 1, residue_number):
                    sequence_ext += "-"
            try:
                sequence_ext += residue_constants.restype_3to1[residue.get_resname()]
            except:
                pass
            prev_residue_number = residue_number
        sequences[chain.id] = sequence_ext
    return sequences


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


def read_seqres(pdb_path: str) -> str:
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


def extract_sequence_from_file(file_path: str) -> List[str]:
    results_dict = {}
    extension = utils.get_file_extension(file_path)
    if extension == '.pdb':
        extraction = 'pdb-atom'
    else:
        extraction = 'cif-atom'
    try:
        with open(file_path, 'r') as f_in:
            for record in SeqIO.parse(f_in, extraction):
                key = f'>{record.id.replace("????", utils.get_file_name(file_path)[:10])}'
                value = str(record.seq.replace("X", ""))
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


def split_pdb_in_chains(output_dir: str, pdb_in_path: str) -> Dict:
    aux_path = os.path.join(output_dir, os.path.basename(pdb_in_path))
    shutil.copy2(pdb_in_path, aux_path)
    chain_dict = chain_splitter(aux_path)
    return chain_dict


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


def run_pdb2cc(templates_path: str, pdb2cc_path: str = None) -> str:
    try:
        cwd = os.getcwd()
        os.chdir(templates_path)
        output_path = 'cc_analysis.in'
        if pdb2cc_path is None:
            pdb2cc_path = 'pdb2cc'
        command_line = f'{pdb2cc_path} -m -i 10 -y 0.15 "orig.*.pdb" 0 {output_path}'
        p = subprocess.Popen(command_line, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        p.communicate()
    finally:
        os.chdir(cwd)
    return os.path.join(templates_path, output_path)


def run_cc_analysis(input_path: str, n_clusters: int, cc_analysis_path: str = None) -> str:
    # Run cc_analysis from cc_analysis_path and store the results in output_path
    # The number cluster to obtain is determinated by n_clusters.
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


def run_hinges(pdb1_path: str, pdb2_path: str, hinges_path: str = None, output_path: str = None) -> Hinges:
    # Run hinges from hinges_path.
    # It needs two pdbs. Return the rmsd obtained.
    chains2_list = [chain.id for chain in get_structure(pdb2_path).get_chains()]
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


def generate_ramachandran(pdb_path, output_path: str = None) -> bool:
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

    if output_path is not None:
        plots.plot_ramachandran(plot_path=os.path.join(output_path, f'{utils.get_file_name(pdb_path)}.png'),
                                phi_psi_angles=phi_psi_angles)

    analysis = ramachandran_analysis(phi_psi_angles=phi_psi_angles)
    percentage = len(analysis) / len(phi_psi_angles) * 100
    logging.info(
        f'{round(percentage, 2)}% of outliers in the ramachandran analysis of {utils.get_file_name(pdb_path)}.')
    if percentage > percentage_minimum:
        return False
    return True


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


def cc_and_hinges_analysis(paths_in: Dict, binaries_path: structures.CCAnalysis, output_path: str,
                           length_sequences: Dict = None) -> List:
    analysis_dict = None
    templates_cluster2 = []
    templates_cluster = hinges(paths_in=paths_in,
                               hinges_path=binaries_path.hinges_path,
                               output_path=os.path.join(output_path, 'hinges'),
                               length_sequences=length_sequences)

    templates_path_list = [template_in for template_list in templates_cluster for template_in in template_list]
    num_templates = len(templates_path_list)
    if num_templates >= 5:
        templates_cluster2, analysis_dict = cc_analysis(paths_in=templates_path_list,
                                                        cc_analysis_paths=binaries_path,
                                                        cc_path=os.path.join(output_path, 'ccanalysis'))

    if len(templates_cluster) > 1 and templates_cluster2:
        return templates_cluster2, analysis_dict
    else:
        return templates_cluster, analysis_dict


def hinges(paths_in: Dict, hinges_path: str, output_path: str, length_sequences: Dict = None) -> List:
    # Hinges algorithm does:
    # Check completeness and ramachandran of every template. If it is not at least 70% discard for hinges.
    # Do hinges 6 iterations in all for all the templates
    # If iter1 < 1 or iterMiddle < 3 or iter6 < 10 AND at least 70% of sequence length -> GROUP TEMPLATE
    # If it has generated more than one group with length > 1-> Return those groups that has generated
    # If there is no group generated -> Return the one more completed and the one more different to that template
    # Otherwise, return all the paths
    utils.create_dir(output_path, delete_if_exists=True)
    threshold_completeness = 0.6
    threshold_rmsd_domains = 10
    threshold_rmsd_ss = 3
    threshold_rmsd_local = 1
    threshold_overlap = 0.7

    logging.info('Starting hinges analysis')
    accepted_pdbs = {}
    uncompleted_pdbs = {}

    # Do the analysis of the different templates. We are going to check:
    # Completeness respect the query size sequence
    # Ramachandran plot
    # And the compactness
    pdb_complete = ''
    pdb_complete_value = 0
    for key, value in paths_in.items():
        num_residues = sum(1 for _ in get_structure(value)[0].get_residues())
        # Validate using ramachandran, check the outliers
        validate_geometry = generate_ramachandran(pdb_path=value, output_path=output_path)
        # Check the query sequence vs the number of residues of the pdb
        completeness = True
        if length_sequences is not None and key in length_sequences:
            completeness = any(number > threshold_completeness for number in length_sequences[key])

        if completeness and validate_geometry:
            accepted_pdbs[key] = value
            if num_residues > pdb_complete_value:
                pdb_complete_value = num_residues
                pdb_complete = value
        else:
            uncompleted_pdbs[key] = value

    logging.info(f'There are {len(accepted_pdbs)} complete pdbs.')
    if len(accepted_pdbs) < 2:
        logging.info(f'Skipping hinges.')
        return [[values for values in paths_in.values()]]
    logging.info(f'Using hinges to create groups.')

    # Run hinges all-against-all, store the results in a dict.
    results_rmsd = {key: {} for key in accepted_pdbs}
    for key1, value1 in accepted_pdbs.items():
        for key2, value2 in accepted_pdbs.items():
            if key2 not in results_rmsd[key1] and key1 != key2:
                result_hinges = run_hinges(pdb1_path=value1, pdb2_path=value2, hinges_path=hinges_path,
                                           output_path=os.path.join(output_path, f'{key1}_{key2}.txt'))
                results_rmsd[key1][key2] = result_hinges
                results_rmsd[key2][key1] = result_hinges
    groups_names = {key: [] for key in accepted_pdbs}
    results_rmsd = OrderedDict(sorted(results_rmsd.items(), key=lambda x: min(v.one_rmsd for v in x[1].values())))
    for key1, value in results_rmsd.items():
        results_rmsd[key1] = OrderedDict(
            sorted({k: v for k, v in value.items() if v is not None}.items(), key=lambda x: x[1].one_rmsd))
        selected_group = key1
        for key2, result in results_rmsd[key1].items():
            group = utils.get_key_by_value(key2, groups_names)
            if group and (
                    (
                            result.min_rmsd < threshold_rmsd_domains or
                            result.one_rmsd < threshold_rmsd_local or
                            result.middle_rmsd < threshold_rmsd_ss
                    )
                    and result.overlap > threshold_overlap
            ):
                if selected_group not in groups_names or len(groups_names[group[0]]) > len(
                        groups_names[selected_group]):
                    selected_group = group[0]
        groups_names[selected_group].append(key1)
    groups_names = [values for values in groups_names.values() if len(values) > 1]

    if len(groups_names) > 1 or (len(groups_names) == 1 and len(groups_names[0]) > 1):
        # Return the groups that has generated
        logging.info(f'Hinges has created {len(groups_names)} group/s:')
        for i, values in enumerate(groups_names):
            logging.info(f'Group {i}: {",".join(values)}')
        return [[path for key, path in paths_in.items() if key in group] for group in groups_names]
    elif len(list(results_rmsd.keys())) > 1:
        # Create two groups, more different and completes pdbs
        more_different = list(results_rmsd[utils.get_file_name(pdb_complete)].keys())[-1]
        path_diff = paths_in[more_different]
        logging.info(f'Hinges could not create any groups')
        logging.info(f'Creating two groups: The more completed pdb: {utils.get_file_name(pdb_complete)} '
                     f'and the more different one: {more_different}')
        return [[pdb_complete], [path_diff]]
    else:
        # Return the original list of pdbs
        logging.info('Not enough pdbs for hinges.')
        return [[values for values in paths_in.values()]]


def cc_analysis(paths_in: List[str], cc_analysis_paths: structures.CCAnalysis, cc_path: str,
                n_clusters: int = 2) -> List:
    # CC_analysis. It is mandatory to have the paths of the programs in order to run pdb2cc and ccanalysis.
    # A dictionary with the different pdbs that are going to be analysed.

    utils.create_dir(cc_path, delete_if_exists=True)
    trans_dict = {}
    return_templates_cluster = [[] for _ in range(n_clusters)]
    clean_dict = {}

    for index, path in enumerate(paths_in):
        path = shutil.copy2(path, cc_path)
        # If it is a ranked, it is mandatory to change the bfactors to VALUE-70.
        # We want to evaluate the residues that have a good PLDDT
        # PDB2CC ignore the residues with bfactors below 0
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
        # Create a translation dictionary, with and index and the pdb name
        trans_dict[index] = utils.get_file_name(path)
    if trans_dict:
        # Write the trans dict in order to be able to trace the pdbs in the output
        with open(os.path.join(cc_path, 'labels.txt'), 'w+') as f:
            for key, value in trans_dict.items():
                f.write('%s:%s\n' % (key, value))
        # run pdb2cc
        output_pdb2cc = run_pdb2cc(templates_path=cc_path, pdb2cc_path=cc_analysis_paths.pd2cc_path)
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
                # Get the positions given by ccanalysis
                points = np.array([values.coord for values in clean_dict.values()])
                # Generate n clusters groups with KMEANS
                kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(points)
                lookup_table = {}
                counter = 0
                # The kmeans results has a list with all the positions belonging to the corresponding pdb.
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
                    real_path = [path for path in paths_in if utils.get_file_name(path) == list(clean_dict.keys())[i]][
                        0]
                    return_templates_cluster[int(label)].append(real_path)
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


def compare_sequences(sequence1: str, sequence2: str) -> List[str]:
    # Given two sequences with same length, return a list showing
    # if there is a match, a group match, they are different, or
    # they are not aligned
    return_list = []

    for i in range(0, len(sequence1)):
        if i < len(sequence2):
            res1 = sequence1[i]
            res2 = sequence2[i]
            if res1 == res2:
                return_list.append('0')
            elif res2 == '-':
                return_list.append('-')
            elif get_group(res1) == get_group(res2):
                return_list.append('0.4')
            else:
                return_list.append('0.7')
        else:
            return_list.append('-')

    return return_list


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

    chains = list(set([chain.get_id() for chain in structure.get_chains()]))

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
    for res in structure[0].get_residues():
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


def remove_hydrogens(pdb_in_path: str, pdb_out_path: str):
    # Remove the atoms that don't belong to the list atom_types
    structure = get_structure(pdb_path=pdb_in_path)
    # Remove hydrogen atoms from the structure
    for residue in structure[0].get_residues():
        atoms = residue.get_unpacked_list()
        for atom in atoms:
            if atom.element not in ['N', 'C', 'O', 'S']:
                residue.detach_child(atom.get_id())

    # Save the edited structure to a new PDB file
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_out_path)


def run_pdbfixer(pdb_in_path: str, pdb_out_path: str):
    try:
        pdb_text = open(pdb_in_path, 'r').read()
        pdb_output = cleanup.fix_pdb(pdb_text, {})
        with open(pdb_out_path, 'w') as f_out:
            f_out.write(pdb_output)
    except:
        logging.info(f'PDBFixer did not finish correctly for {utils.get_file_name(pdb_in_path)}. Skipping.')
        shutil.copy2(pdb_in_path, pdb_out_path)
        pass


def run_arcimboldo_air(yml_path: str):
    arcimboldo_air_path = os.path.join(utils.get_main_path(), 'arcimboldo_air.py')
    command_line = f'{arcimboldo_air_path} {yml_path}'
    p = subprocess.Popen(command_line, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    logging.info('ARCIMBOLDO_AIR cluster run finished successfully.')


def run_openmm(pdb_in_path: str, pdb_out_path: str) -> structures.OpenmmEnergies:
    try:
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
    except:
        logging.info(f'Not possible to calculate energies for {utils.get_file_name(pdb_in_path)}')
        return structures.OpenmmEnergies(None, None)


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


def gesamt_pdbs(pdb_list: List[str], output_path: str = None) -> Tuple[Optional[float], Optional[str], Optional[str]]:
    name_folder = 'tmp_gesamt'
    utils.create_dir(name_folder, delete_if_exists=True)
    superpose_input_list = ['gesamt'] + pdb_list
    if output_path is not None:
        superpose_input_list.extend(['-o', name_folder, '-o-d'])
    superpose_output = subprocess.Popen(superpose_input_list, stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
    new_path = [file for file in os.listdir(name_folder) if file.startswith(utils.get_file_name(pdb_list[-1]))]
    if new_path:
        shutil.copy2(os.path.join(name_folder, new_path[0]), output_path)
    shutil.rmtree(name_folder)
    rmsd, quality_q, nalign = None, None, None
    for line in superpose_output.split('\n'):
        if 'RMSD             :' in line:
            rmsd = float(line.split()[2].strip())
        if 'Q-score          :' in line:
            quality_q = line.split()[2].strip()
        if 'Aligned residues :' in line:
            nalign = line.split()[3].strip()
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
    # For a PDB and an interface (contains the chains and the residues involved in the interface).
    # Iterate the interface, selecting all those residues that belong to a domain (we have previously calculated
    # all the domains in the PDB).
    # We create a dictionary with all the residues of the interface, extending them to their whole domain.
    # Split the pdb into the chains of the interface and delete all the others residues.
    # Return the dictionary with the chains and the interface residues extended.
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
                        chain_list=[interface['chain1'], interface['chain2']])

    change = change_res.ChangeResidues(chain_res_dict=add_domains_dict, chain_bfactors_dict=bfactors_dict)
    change.delete_residues_inverse(pdb_out_path, pdb_out_path)
    # change.change_bfactors(pdb_out_path, pdb_out_path)

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

    if trimmed_list:
        max_length = max(len(lst) for lst in trimmed_list)
    else:
        return []

    trimmed_list = [lst for lst in trimmed_list if len(lst) == max_length]

    score_list = []
    for element in trimmed_list:
        score_list.append(sum(z for _, _, z, _ in element))

    if score_list:
        min_value = min(score_list)
        min_index = score_list.index(min_value)
        return trimmed_list[min_index]
    else:
        return []


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
