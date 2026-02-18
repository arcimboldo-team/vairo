#! /usr/bin/env python3
import json
import re
import shutil
import subprocess
import sys
import pickle
import os
import logging
import tempfile
import collections
from typing import List
import xml.etree.ElementTree as ET
import numpy as np
import requests
import csv
from alphafold.data import parsers, templates, mmcif_parsing
from Bio.Blast import NCBIWWW
from alphafold.common import residue_constants
from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain

current_directory = os.path.dirname(os.path.abspath(__file__))
target_directory = os.path.abspath(os.path.join(current_directory, '..', '..'))
sys.path.append(target_directory)

from vairo.libs import features, bioutils, plots, utils, structures, template_modifications, global_variables

def shift_pdb(template_path, output_path, where, shift_amount):
    parser = PDBParser(QUIET=True)
    try:
        original_structure = parser.get_structure('structure', template_path)
    except Exception as e:
        print(f"Error parsing file: {e}")
        return

    new_structure = Structure.Structure("renumbered")
    for i, model in enumerate(original_structure):
        new_model = Model.Model(model.id)
        new_structure.add(new_model)
        for chain in model:
            new_chain = Chain.Chain(chain.id)
            new_model.add(new_chain)
            residues = list(chain)
            if not residues:
                continue
            shift = shift_amount if where in [chain.id, 'all'] else 0
            print(f"  Chain {chain.id}: shifting by {shift}...")
            for res in residues:
                res.parent = None
                old_id = res.id
                new_number = old_id[1] + int(shift)
                res.id = (old_id[0], new_number, old_id[2])
                new_chain.add(res)
    io = PDBIO()
    io.set_structure(original_structure)
    try:
        io.save(output_path)
        print(f"Successfully saved renumbered PDB to: {output_path}")
    except Exception as e:
        print(f"Error saving file: {e}")


def renumber_template(template_path, output_path):
    parser = PDBParser(QUIET=True)
    try:
        original_structure = parser.get_structure("original", template_path)
    except FileNotFoundError:
        print(f"Error: Could not find file {template_path}")
        return
    new_structure = Structure.Structure("renumbered")
    print(f"Processing {template_path}...")
    for i, model in enumerate(original_structure):
        new_model = Model.Model(model.id)
        new_structure.add(new_model)
        for chain in model:
            new_chain = Chain.Chain(chain.id)
            new_model.add(new_chain)
            residues = list(chain)
            if not residues:
                continue
            first_residue_num = residues[0].id[1]
            offset = first_residue_num - 1
            print(f"  Chain {chain.id}: Starts at {first_residue_num}. shifting by -{offset}...")
            for res in residues:
                res.parent = None
                old_id = res.id
                new_number = old_id[1] - offset
                res.id = (old_id[0], new_number, old_id[2])
                new_chain.add(res)

    io = PDBIO()
    io.set_structure(new_structure)
    io.save(output_path)
    print(f"Success! Renumbered PDB saved to: {output_path}")


def write_features(features_path: str, output_dir: str = None):
    with open(os.path.abspath(features_path), 'rb') as f:
        data = pickle.load(f)
    if output_dir is None:
        output_dir = os.getcwd()
    features.write_templates_in_features(data, output_dir)


def print_features(features_path: str):
    logging.error = print
    features.print_features_from_file(features_path)


def print_sequence_info(seq_dict: dict, seq_type: str, ini: int = 0, end: int = 100):
    seq_sorted = sorted(seq_dict.items(), key=lambda x: x[1]['identity'], reverse=True)
    accepted_identity_elements = {key: value for key, value in seq_sorted if end >= value['identity'] >= ini}

    print(f'{seq_type} that have between a {ini}% and {end}% identity percentage ({len(accepted_identity_elements)}):')
    for key, values in accepted_identity_elements.items():
        print(F'SEQUENCE {key}')
        print(
            f'ID: {key} || Identity: {values["identity"]}% || Global Identity: {values["global_identity"]}% || Coverage: {values["coverage"]}%\n{values["seq"]}\n')
    return accepted_identity_elements


def mutate_features(features_path: str):
    feature = features.create_features_from_file(pkl_in_path=features_path)
    for i in range(feature.get_msa_length()):
        feature.msa_features['msa'][i][190] = residue_constants.HHBLITS_AA_TO_ID['N']
        feature.msa_features['msa'][i][161] = residue_constants.HHBLITS_AA_TO_ID['N']
        feature.msa_features['msa'][i][162] = residue_constants.HHBLITS_AA_TO_ID['T']
        feature.msa_features['msa'][i][244] = residue_constants.HHBLITS_AA_TO_ID['Y']
        feature.msa_features['msa'][i][40] = residue_constants.HHBLITS_AA_TO_ID['I']
    feature.write_pkl(features_path)


def extract_features_info(features_path: str, region: str = None, ini_identity: int = 0, end_identity: int = 100,
                          run_uniprot: bool = False):
    if region is None or region == "":
        region = '1-10000'
    region_list = region.replace(" ", "").split(',')
    region_result = []
    for r in region_list:
        start, end = map(int, r.split('-'))
        region_result.append((int(start), int(end)))

    
    ini_identity = int(ini_identity)
    end_identity = int(end_identity)
                              
    features_info_dict, region_query, query = features.extract_features_info(pkl_in_path=features_path,
                                                                             regions_list=region_result)
    files_info = []
    print('\n================================')
    print(f'REGION {region}')
    print(f'And we are looking for these specific regions:')
    print(f'{region_query}')

    if features_info_dict['msa']:
        print('\nMSA:')
        features_info_dict['msa'] = print_sequence_info(features_info_dict['msa'], 'Sequences', ini=ini_identity,
                                                        end=end_identity)
    if features_info_dict['templates']:
        print('\nTEMPLATES:')
        features_info_dict['templates'] = print_sequence_info(features_info_dict['templates'], 'Templates',
                                                              ini=ini_identity, end=end_identity)

    store_fasta_path = os.path.join(os.getcwd(), f'accepted_sequences_{region}.fasta')
    merged_dict = {**features_info_dict['msa'], **features_info_dict['templates']}
    print(f'Accepted sequences tanking into account the identity have been stored in: {store_fasta_path}')
    with open(store_fasta_path, 'w') as file:
        duplicate_list = []
        for key, values in merged_dict.items():
            if values["seq"] not in duplicate_list:
                file.write(f'\n>{key}\n')
                file.write(f'{values["seq"]}')
                duplicate_list.append(values["seq"])
    print('\n================================')
    files_info.append(store_fasta_path)

    residues_list = []
    for r in region_result:
        residues_list.extend(list(range(r[0], r[1] + 1)))
    if run_uniprot:
        results_uniprot = run_uniprot_blast(store_fasta_path, residues_list)
    else:
        results_uniprot = {}
    features_info_dict['templates'] = dict(
        sorted(features_info_dict['templates'].items(), key=lambda item: item[1]['identity'], reverse=True))
    features_info_dict['msa'] = dict(
        sorted(features_info_dict['msa'].items(), key=lambda item: item[1]['identity'], reverse=True))

    templates_keys = list(features_info_dict['templates'].keys())
    msa_keys = list(features_info_dict['msa'].keys())
    uniprot_description_statistics = collections.defaultdict(int)
    uniprot_organism_statistics = collections.defaultdict(int)
    uniprot_desidentity_statistics = collections.defaultdict(int)
    uniprot_orgidentity_statistics = collections.defaultdict(int)
    for key, value in results_uniprot.items():
        if key in templates_keys:
            features_info_dict['templates'][key]['uniprot'] = value
        elif key in msa_keys:
            features_info_dict['msa'][key]['uniprot'] = value

        for uni_element in value:
            protein_description = uni_element['uniprot_protein_description']
            uniprot_description_statistics[protein_description] += 1
            uniprot_desidentity_statistics[protein_description] = max(
                uniprot_desidentity_statistics[protein_description],
                int(uni_element['uniprot_identity']))

            organism_description = uni_element['uniprot_organism']
            uniprot_organism_statistics[organism_description] += 1
            uniprot_orgidentity_statistics[organism_description] = max(
                uniprot_orgidentity_statistics[organism_description],
                int(uni_element['uniprot_identity']))

    combined_uniprot_dict = {k: {'description': v, 'identity': uniprot_desidentity_statistics[k]}
                             for k, v in uniprot_description_statistics.items()}

    combined_org_uniprot_dict = {k: {'organism': v, 'identity': uniprot_orgidentity_statistics[k]}
                                 for k, v in uniprot_organism_statistics.items()}

    features_info_dict['uniprot_description_statistics'] = dict(
        sorted(combined_uniprot_dict.items(), key=lambda item: item[1]['description'], reverse=True))
    features_info_dict['uniprot_organism_statistics'] = dict(
        sorted(combined_org_uniprot_dict.items(), key=lambda item: item[1]['organism'], reverse=True))

    features_info_dict['general_information'] = {}
    features_info_dict['general_information']['query_sequence'] = query
    features_info_dict['general_information']['query_search'] = region_query

    templates_coverage = [0] * len(query)
    msa_coverage = [0] * len(query)
    if len(features_info_dict['templates']):
        templates_seq_list = [seq['seq'] for seq in features_info_dict['templates'].values()]
        templates_coverage = bioutils.calculate_coverage_scaled(query_seq=query, sequences=templates_seq_list)
    if len(features_info_dict['msa']):
        msa_seq_list = [seq['seq'] for seq in features_info_dict['msa'].values()]
        msa_coverage = bioutils.calculate_coverage_scaled(query_seq=query, sequences=msa_seq_list)

    features_info_dict['coverage'] = {
        'msa_coverage': msa_coverage,
        'num_msa': len(features_info_dict['msa']),
        'templates_coverage': templates_coverage,
        'num_templates': len(features_info_dict['templates'])
    }
    return features_info_dict


def generate_features(query_path: str, fasta_path: str):
    path = os.path.join(os.getcwd(), 'features.pkl')
    query = bioutils.extract_sequence(query_path)
    sequences = bioutils.extract_sequences(fasta_path)
    feature = features.Features(query)
    [feature.append_row_in_msa(sequence=seq, sequence_id=seq_id) for seq_id, seq in sequences.items()]
    write_features(path)


def hinges(template_path: str):
    output_path = os.path.join(template_path, 'hinges')
    os.listdir(template_path)
    templates_dict = {utils.get_file_name(path): os.path.join(template_path, path) for path in os.listdir(template_path)
                      if path.endswith('.pdb')}

    binaries_path = structures.CCAnalysis(os.path.join(utils.get_main_path(), 'binaries'))
    templates_cluster = bioutils.hinges(paths_in=templates_dict,
                                        binaries_path=binaries_path,
                                        output_path=output_path)

    for i, values in enumerate(templates_cluster):
        print(f'Group {i}: {",".join(values)}')


def ccanalysis(template_path: str):
    output_path = os.path.join(template_path, 'ccanalysis')
    os.listdir(template_path)
    templates_list = [structures.Pdb(os.path.join(template_path, path)) for path in os.listdir(template_path)
                      if path.endswith('.pdb')]
    binaries_path = structures.BinariesPath(os.path.join(utils.get_main_path(), 'binaries'))
    templates_cluster_list, analysis_dict = bioutils.cc_analysis(pdbs=templates_list, cc_analysis_paths=binaries_path,
                                                                 output_dir=output_path, n_clusters=2)
    if analysis_dict:
        plots.plot_cc_analysis(plot_path=os.path.join(output_path, 'plot.png'), analysis_dict=analysis_dict,
                               clusters=templates_cluster_list)


def superposition_chains(pdb1_path: str, pdb2_path: str):
    ret_dict = bioutils.superposition_by_chains(pdb1_in_path=pdb1_path, pdb2_in_path=pdb2_path)
    for key3, i3 in ret_dict.items():
        for key2, i2 in i3.items():
            for key1, i1 in i2.items():
                print(key3, key2, key1, i1)


def run_minimize(pdb1_path: str, pdb2_path: str):
    bioutils.remove_hetatm(pdb1_path, pdb2_path)
    print(bioutils.run_openmm(pdb2_path, pdb2_path))


def renumber():
    def check_consecutive(numbers):
        # Check if the difference between each pair of consecutive numbers is equal to 1
        for i in range(len(numbers) - 1):
            if numbers[i + 1] - numbers[i] != 1:
                return False
        return True

    # Specify the folder path containing the PDB files
    folder_path = "/Users/pep/work/transfers/clusters_lib"
    # Get a list of all PDB files in the folder
    pdb_files = [os.path.join(folder_path, file) for file in os.listdir(folder_path) if file.endswith(".pdb")]

    list_pdbs = []

    # Loop through each PDB file
    for pdb_file in pdb_files:

        structure = bioutils.get_structure(pdb_file)

        # Initialize counters for CYS residues and consecutive positions
        cys_count = 0
        save_residues = []
        save_pdb = False

        # Iterate over all residues in the structure
        for model in structure:
            for chain in model:
                residues = list(chain.get_residues())
                residues = sorted(residues, key=lambda x: bioutils.get_resseq(x))
                for j, residue in enumerate(residues):
                    # Check if the residue is CYS
                    if residue.get_resname() == 'CYS':
                        cys_count += 1
                        if cys_count == 1:
                            try:
                                list_cys = [residues[j + i] for i in range(-5, 2)]
                                list_cys = [bioutils.get_resseq(res) - 1 for res in list_cys]
                                if check_consecutive(list_cys):
                                    save_residues.extend(list_cys)
                                else:
                                    raise Exception
                            except:
                                cys_count = 3
                                pass
                        if cys_count == 2:
                            try:
                                list_cys = [residues[j + i] for i in range(-5, 3)]
                                list_cys = [bioutils.get_resseq(res) - 1 for res in list_cys]
                                if check_consecutive(list_cys):
                                    save_residues.extend(list_cys)
                                    if utils.get_file_name(pdb_file)[:4] not in list_pdbs:
                                        list_pdbs.append(utils.get_file_name(pdb_file)[:4])
                                        save_pdb = True
                                else:
                                    raise Exception
                            except:
                                pass

        if save_pdb:
            if len(save_residues) != 15:
                raise Exception
            bioutils.copy_positions_of_pdb(pdb_file, os.path.join("/Users/pep/work/transfers/library",
                                                                  utils.get_file_name(pdb_file)) + '.pdb',
                                           save_residues)
            print(f"Renumbering complete for {pdb_file}. Renumbered file saved as {utils.get_file_name(pdb_file)}.")


def merge_pdbs(pdb1_path: str, pdb2_path: str, inf_ini, inf_end, inm_ini, inm_end):
    MIN_RMSD_SPLIT = 5

    best_rankeds_dir = os.path.join(os.getcwd(), 'merged_pdb')
    utils.create_dir(best_rankeds_dir, delete_if_exists=True)
    aux_pdb1_path = os.path.join(best_rankeds_dir, 'pdb1_trimmed.pdb')
    merge_pdbs_list = [aux_pdb1_path]
    shutil.copy2(pdb1_path, best_rankeds_dir)
    shutil.copy2(pdb2_path, best_rankeds_dir)

    pdb_out = os.path.join(best_rankeds_dir, 'superposed.pdb')
    delta_out = os.path.join(best_rankeds_dir, 'deltas.dat')

    bioutils.run_lsqkab(pdb_inf_path=pdb1_path,
                        pdb_inm_path=pdb2_path,
                        fit_ranges=[(inf_ini, inf_end)],
                        match_ranges=[(inm_ini, inm_end)],
                        pdb_out=pdb_out,
                        delta_out=delta_out
                        )
    best_list = []
    best_min = MIN_RMSD_SPLIT
    with open(delta_out, 'r') as f_in:
        lines = f_in.readlines()
        lines = [line.replace('CA', '').split() for line in lines]
        for deltas in zip(lines, lines[1:], lines[2:], lines[3:]):
            deltas_sum = sum([float(delta[0]) for delta in deltas])
            if deltas_sum <= best_min:
                best_list = deltas
                best_min = deltas_sum

    if not best_list:
        raise Exception('RMSD minimum requirements not met in order to merge the results in mosaic mode.')

    inf_cut = int(best_list[1][3])
    inm_cut = int(best_list[2][1])

    delete_residues = template_modifications.TemplateModifications()
    delete_residues.append_modification(chains=['A'], delete_residues=[*range(inf_cut + 1, 10000 + 1, 1)])
    delete_residues.modify_template(pdb_in_path=pdb1_path, pdb_out_path=aux_pdb1_path, type_modify='delete')
    delete_residues = template_modifications.TemplateModifications()
    delete_residues.append_modification(chains=['A'], delete_residues=[*range(1, inm_cut, 1)])
    delete_residues.modify_template(pdb_in_path=pdb_out, pdb_out_path=pdb_out, type_modify='delete')

    merge_pdbs_list.append(pdb_out)
    bioutils.merge_pdbs_in_one_chain(list_of_paths_of_pdbs_to_merge=merge_pdbs_list,
                                     pdb_out_path=os.path.join(best_rankeds_dir, 'merged.pdb'))


def align_pdb(hhr_path: str, pdb_path: str, fasta_path: str):
    query_sequence = bioutils.extract_sequence(fasta_path=fasta_path)
    output_dir = os.path.join(os.getcwd(), 'templates')
    cif_path = os.path.join(output_dir, f'{utils.get_file_name(pdb_path)}.cif')
    pdb_path = os.path.abspath(pdb_path)

    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)

    chain_dict = bioutils.split_pdb_in_chains(pdb_path=pdb_path)
    results_list = []
    for chain, path in chain_dict.items():
        bioutils.pdb2mmcif(pdb_in_path=path, cif_out_path=cif_path)
        pdb_id = utils.get_file_name(pdb_path).upper()
        hhr_text = open(hhr_path, 'r').read()
        matches = re.finditer(r'No\s+\d+', hhr_text)
        matches_positions = [match.start() for match in matches] + [len(hhr_text)]

        detailed_lines_list = []
        for i in range(len(matches_positions) - 1):
            detailed_lines_list.append(hhr_text[matches_positions[i]:matches_positions[i + 1]].split('\n')[:-3])

        hits_list = [detailed_lines for detailed_lines in detailed_lines_list if
                     pdb_id in detailed_lines[1]]

        detailed_lines = hits_list[0]

        try:
            hit = parsers._parse_hhr_hit(detailed_lines)
        except:
            return None, None, None, 0, 0, 0

        template_sequence = hit.hit_sequence.replace('-', '')
        mapping = templates._build_query_to_hit_index_mapping(
            hit.query, hit.hit_sequence, hit.indices_hit, hit.indices_query,
            query_sequence)
        mmcif_string = open(cif_path).read()
        parsing_result = mmcif_parsing.parse(file_id=pdb_id, mmcif_string=mmcif_string)
        template_features, _ = templates._extract_template_features(
            mmcif_object=parsing_result.mmcif_object,
            pdb_id=pdb_id,
            mapping=mapping,
            template_sequence=template_sequence,
            query_sequence=query_sequence,
            template_chain_id=chain,
            kalign_binary_path='kalign')

        template_features['template_sum_probs'] = np.array([[hit.sum_probs]])
        template_features['template_aatype'] = np.array([template_features['template_aatype']])
        template_features['template_all_atom_masks'] = np.array([template_features['template_all_atom_masks']])
        template_features['template_all_atom_positions'] = np.array([template_features['template_all_atom_positions']])
        template_features['template_domain_names'] = np.array([template_features['template_domain_names']])
        template_features['template_sequence'] = np.array([template_features['template_sequence']])
        features.write_templates_in_features(template_features=template_features, output_dir=output_dir)
        result_pdb = os.path.join(output_dir, f'{pdb_id}_{chain}1.pdb')
        bioutils.change_chain(pdb_in_path=result_pdb, pdb_out_path=result_pdb, chain=chain)
        results_list.append(result_pdb)

    bioutils.merge_pdbs(list_of_paths_of_pdbs_to_merge=results_list, merged_pdb_path='result.pdb')


def delete_msas(pkl_in_path: str, pkl_out_path: str, delete_str: str):
    delete_list = list(map(int, delete_str.split(',')))
    features.delete_seq_from_msa(pkl_in_path=pkl_in_path, pkl_out_path=pkl_out_path, delete_list=delete_list)


def delete_templates(pkl_in_path: str, pkl_out_path: str, delete_str: str):
    delete_list = delete_str.split(',')
    feat = features.create_features_from_file(pkl_in_path=pkl_in_path)
    feat.delete_by_id(delete_list)
    feat.write_pkl(pkl_out_path)


def select_csv(pkl_in_path: str, csv_path: str, min_input: float, max_input: float):
    accepted_list = []
    deleted_list = []
    min_aux = min([min_input, max_input])
    max_aux = max([min_input, max_input])
    new_features_path = os.path.join(os.path.dirname(pkl_in_path), f'features_{min_aux}-{max_aux}.pkl')
    with open(csv_path) as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        next(csvreader)
        for row in csvreader:
            y_value = row[2]
            name = int(row[6])
            if min_aux <= y_value <= max_aux:
                accepted_list.append(name)
            else:
                deleted_list.append(name)

    print(f'The deleted list is the following one: {", ".join(list(map(str, deleted_list)))}\n')
    print(f'The accepted list is the following one: {", ".join(list(map(str, accepted_list)))}\n')
    print(f'The features file with just the accepted sequences can be found in:')
    features.delete_seq_from_msa(pkl_in_path=pkl_in_path, pkl_out_path=new_features_path, delete_list=deleted_list)


def run_uniprot_blast(fasta_path: str, residues_list: List[int], use_server: bool = False):
    sequences_dict = bioutils.extract_sequences(fasta_path)
    print(f'Running BLASTP with the sequences inside the file {fasta_path}')
    results_dict = {}
    for id, seq in sequences_dict.items():
        results_dict[id] = []
        print('================================')
        print(F'SEQUENCE {id}')
        print(f'Sequence id {id} with length {len(seq)}:')
        print(f'{seq}')
        modified_content = seq.replace('-', 'X')
        if use_server:
            database = 'swissprot'
            result_handle = NCBIWWW.qblast("blastp", database, modified_content)
            root = ET.fromstring(result_handle.read())
        else:
            with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp_file:
                temp_file.write(modified_content)
                temp_file.flush()
                blastp_path = shutil.which('blastp')
                blast_folder = os.path.dirname(blastp_path)
                blastp_cmd = f'{blastp_path} -db {os.path.join(blast_folder, "swissprot")} -query {temp_file.name} -outfmt 5'
                result = subprocess.Popen(blastp_cmd, stdout=subprocess.PIPE, shell=True)
                blastp_output = result.communicate()[0].decode('utf-8')
                root = ET.fromstring(blastp_output)

        for iteration_elem in root.findall('.//Hit'):
            ini_query = iteration_elem.find('.//Hsp_query-from').text
            end_query = iteration_elem.find('.//Hsp_query-to').text
            search_range = range(int(ini_query), int(end_query) + 1)
            residues = list(set(search_range) & set(residues_list))
            if len(residues) > 1:
                check = True
            else:
                check = False
            hit_accession = iteration_elem.find('.//Hit_accession').text
            evalue = iteration_elem.find('.//Hsp_evalue').text
            residues_identity = iteration_elem.find('.//Hsp_identity').text
            aligned_identity = iteration_elem.find('.//Hsp_align-len').text
            hsp_qseq = iteration_elem.find('.//Hsp_qseq').text
            hsp_hseq = iteration_elem.find('.//Hsp_hseq').text
            if float(evalue) < float(0.01):
                url = f'https://rest.uniprot.org/uniprotkb/search?query=accession_id:{hit_accession}&fields=annotation_score,protein_name,organism_name'
                response = requests.get(url)
                print('--------------')
                json_response = response.json()
                annotation_score = json_response['results'][0]['annotationScore']
                protein_description = json_response['results'][0]['proteinDescription']['recommendedName']['fullName'][
                    'value']
                organism = json_response['results'][0]['organism']['scientificName']
                print(f'Accession ID: {hit_accession}')
                print(f'E-value: {evalue}')
                print(f'Identity residues: {residues_identity}')
                print(f'Aligned residues: {aligned_identity}')
                print(f'Annotation Score: {annotation_score}')
                print(f'Protein description: {protein_description}')
                print(f'Organism: {organism}')
                if check:
                    print(f'It shares residues {", ".join(map(str, residues))} with the searched range')
                    print(f'The search query matching the range is the following one:')
                    chain = 'X' * (int(ini_query) - 1)
                    for i, res in enumerate(hsp_qseq, start=int(ini_query)):
                        if i in residues:
                            chain += res
                        else:
                            chain += 'X'
                    print("".join(chain))
                    print(f'The found sequence matching range is the following one:')
                    chain = 'X' * (int(ini_query) - 1)
                    for i, res in enumerate(hsp_hseq, start=int(ini_query)):
                        if i in residues:
                            chain += res
                        else:
                            chain += 'X'
                    print("".join(chain))
                else:
                    print(f'It does not have any matching residue with the specified range')

                results_dict[id].append({
                    'uniprot_protein_description': protein_description,
                    'uniprot_annotation_score': annotation_score,
                    'uniprot_organism': organism,
                    'uniprot_accession_id': hit_accession,
                    'uniprot_identity': residues_identity,
                    'uniprot_evalue': evalue
                })

        print('================================')
    return results_dict


def analyse_scwrl(fasta_path, template_path, num_steps):

    scwrl_path = '/opt/scwrl4/Scwrl4'
    sequences = bioutils.extract_sequences(fasta_path)
    sequence_template_dict = bioutils.extract_sequence_msa_from_pdb(template_path)
    template_dict = {}
    run_results = []

    print(f"Loaded sequences: {sequences}")

    for i, (chain, template_sequence) in enumerate(sequence_template_dict.items()):
        first_pos = None
        last_pos = None
        for index, item in enumerate(sequence_template_dict[chain]):
            if '-' not in item:
                first_pos = index
                break
        for index in range(len(sequence_template_dict[chain]) - 1, -1, -1):
            if '-' not in sequence_template_dict[chain][index]:
                last_pos = index
                break
        if first_pos is not None and last_pos is not None:
            template_dict[chain] = (list(sequences.values())[i], first_pos, last_pos)
    print(f"Template mapping: {template_dict}")

    interfaces_pdbs =  os.path.join(os.getcwd(), "interfaces_pdbs")
    os.makedirs(interfaces_pdbs, exist_ok=True)

    aleph_dir = os.path.join(os.getcwd(), "aleph")
    os.makedirs(aleph_dir, exist_ok=True)

    for i in range(-int(num_steps), int(num_steps) + 1, 1):
        general_sequence = ''
        valid_step = True

        for temp in template_dict.values():
            try:
                start_slice = temp[1] + i
                end_slice = temp[2] + i + 1
                if start_slice < 0: raise IndexError
                segment = temp[0][start_slice: end_slice]
                general_sequence += segment
            except IndexError:
                print(f"Warning: Index out of bounds at step {i}. Skipping.")
                valid_step = False
                break

        if not valid_step: continue

        # Directories
        os.makedirs("results", exist_ok=True)
        step_dir = os.path.join(os.getcwd(), "results", f"step_{i}")
        os.makedirs(step_dir, exist_ok=True)
        interfaces_dir = os.path.join(step_dir, "interfaces")
        os.makedirs(interfaces_dir, exist_ok=True)

        new_fasta_path = os.path.join(step_dir, "new.fasta")
        name = f"model_step_{i}"
        output_pdb_path = os.path.join(step_dir, f"{name}.pdb")

        with open(new_fasta_path, 'w') as f:
            f.write(general_sequence)

        cmd = [scwrl_path, '-i', template_path, '-o', output_pdb_path, '-s', new_fasta_path]
        print(f"Running SCWRL for step {i}...")
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        energy_match = re.search(r"Total minimal energy of the graph\s*=\s*([-\d\.]+)", result.stdout)


        results_dict, domains_dict = bioutils.aleph_annotate(output_path=aleph_dir, pdb_path=output_pdb_path)

        if energy_match:
            energy_score = float(energy_match.group(1))
            print(f"  -> SCWRL Success. Energy: {energy_score}")

            print(f"  -> Running PISA analysis for step {i}...")
            interfaces_data_list = bioutils.find_interface_from_pisa(output_pdb_path, interfaces_dir)
            print(f"  -> Found {len(interfaces_data_list)} interfaces.")

            for interface in interfaces_data_list:
                if not interface.chain1 in domains_dict or not interface.chain2 in domains_dict:
                    continue
                code = f'{interface.chain1}-{interface.chain2}'
                dimers_path = os.path.join(interfaces_pdbs, f'{name}_{code}.pdb')
                bioutils.create_interface_domain(pdb_in_path=output_pdb_path,
                                                 pdb_out_path=dimers_path,
                                                 interface=interface,
                                                 domains_dict=domains_dict)
                interface.set_structure(dimers_path)

            run_results.append({
                'step': i,
                'scwrl_energy': energy_score,
                'pdb_path': output_pdb_path,
                'sequence': general_sequence,
                'interfaces': interfaces_data_list
            })
        else:
            print("  -> Warning: SCWRL energy parsing failed.")


    # --- Analysis & Ranking ---
    print("\n" + "=" * 30)
    print("ANALYSIS OF SCWRL RUNS & INTERFACES")
    print("=" * 30)

    if not run_results:
        print("No successful runs to analyze.")
        return

    interfaces_by_type = collections.defaultdict(list)

    for run in run_results:
        for interface in run['interfaces']:
            key = interface.name
            interfaces_by_type[key].append({
                'step': run['step'],
                'scwrl_energy': run['scwrl_energy'],
                'pisa_deltaG': interface.deltaG,
                'interface_obj': interface
            })

    for key in interfaces_by_type:
        interfaces_by_type[key].sort(key=lambda x: x['pisa_deltaG'])

    # Print Ranking Table
    print(f"{'Interface':<10} | {'Rank':<5} | {'Step':<5} | {'Delta G':<10} | {'SCWRL E':<10}")
    print("-" * 60)

    for key, ranked_list in interfaces_by_type.items():
        print(f"--- Type: {key} ---")
        for rank, data in enumerate(ranked_list, 1):
            print(
                f"{key:<10} | {rank:<5} | {data['step']:<5} | {data['pisa_deltaG']:<10.4f} | {data['scwrl_energy']:<10.4f}")

    # --- 3. Generate PyMOL Script ---
    print("\nGenerating PyMOL visualization script...")
    generate_pymol_script(interfaces_by_type, "view_interfaces.pse")


def generate_pymol_script(ranked_interfaces_dict, output_path):
    """
    Generates a PyMOL script to visualize the ranked interfaces.
    """
    script = """
from pymol import cmd
cmd.set("bg_rgb", "0xffffff")
cmd.set("antialias", '2')
cmd.set("ribbon_sampling", '10')
cmd.set("hash_max", '220')
cmd.set("dash_length", '0.10000')
cmd.set("dash_gap", '0.30000')
cmd.set("cartoon_sampling", '14')
cmd.set("cartoon_loop_quality", '6.00000')
cmd.set("cartoon_rect_length", '1.10000')
cmd.set("cartoon_oval_length", '0.80000')
cmd.set("cartoon_oval_quality", '10.00000')
cmd.set("cartoon_tube_quality", '9.00000')
cmd.set("dash_width", '3.00000')
cmd.set("transparency", '0.60000')
cmd.set("two_sided_lighting", '0')
cmd.set("sculpt_vdw_weight", '0.45000')
cmd.set("sculpt_field_mask", '2047')
cmd.set("ray_shadow", 'off')
cmd.set("auto_color_next", '2')
cmd.set("button_mode_name", '3-Button Viewing')
cmd.set("mouse_selection_mode", '2')
cmd.set("cartoon_nucleic_acid_mode", '2')
cmd.set("cartoon_putty_quality", '11.00000')
cmd.set("cartoon_ring_mode", '1')
cmd.set("cartoon_ladder_color", 'cyan')
cmd.set("cartoon_nucleic_acid_color", 'cyan')
cmd.set("ray_trace_mode", '1')
cmd.set("sculpt_min_weight", '2.25000')
cmd.set("mesh_negative_color", 'grey30')
cmd.set("ray_transparency_oblique_power", '1.00000')
cmd.set("movie_quality", '60')
cmd.set("use_shaders", 'on')
cmd.set("volume_bit_depth", '8')
cmd.set("mesh_as_cylinders", 'on')
cmd.set("line_as_cylinders", 'on')
cmd.set("ribbon_as_cylinders", 'on')
cmd.set("nonbonded_as_cylinders", 'on')
cmd.set("nb_spheres_quality", '3')
cmd.set("alignment_as_cylinders", 'on')
cmd.set("dot_as_spheres", 'on')
cmd.set("valence", 'off')
"""

    added_objs = []

    for interface_type, interfaces in ranked_interfaces_dict.items():
        for rank, data in enumerate(interfaces):
            interface = data['interface_obj']
            step = data['step']
            obj_name = f"{interface_type}_step{step}_rank{rank}"
            if interface.path and os.path.exists(interface.path):
                script += f'cmd.load("{interface.path}", "{obj_name}")\n'
                script += f'cmd.show_as("sticks", "{obj_name}")\n'
                script += f'cmd.show("surface", "{obj_name}")\n'
                script += f'cmd.set("transparency", "0.50000", "{obj_name}")\n'

                script += f'cmd.color("lime", "resn ALA and {obj_name}")\n'
                script += f'cmd.color("density", "resn ARG and {obj_name}")\n'
                script += f'cmd.color("deepsalmon", "resn ASN and {obj_name}")\n'
                script += f'cmd.color("warmpink", "resn ASP and {obj_name}")\n'
                script += f'cmd.color("paleyellow", "resn CYS and {obj_name}")\n'
                script += f'cmd.color("tv_red", "resn GLN and {obj_name}")\n'
                script += f'cmd.color("ruby", "resn GLU and {obj_name}")\n'
                script += f'cmd.color("slate", "resn HIS and {obj_name}")\n'
                script += f'cmd.color("forest", "resn ILE and {obj_name}")\n'
                script += f'cmd.color("smudge", "resn LEU and {obj_name}")\n'
                script += f'cmd.color("deepblue", "resn LYS and {obj_name}")\n'
                script += f'cmd.color("sand", "resn MET and {obj_name}")\n'
                script += f'cmd.color("gray40", "resn PHE and {obj_name}")\n'
                script += f'cmd.color("gray20", "resn PRO and {obj_name}")\n'
                script += f'cmd.color("tv_orange", "resn SER and {obj_name}")\n'
                script += f'cmd.color("brown", "resn THR and {obj_name}")\n'
                script += f'cmd.color("palegreen", "resn TRP and {obj_name}")\n'
                script += f'cmd.color("wheat", "resn TYR and {obj_name}")\n'
                script += f'cmd.color("pink", "resn VAL and {obj_name}")\n'

                if rank > 0:
                    script += f'cmd.disable("{obj_name}")\n'

                added_objs.append(obj_name)

    script += f'cmd.save("{output_path}")\n'
    script += 'cmd.quit()\n'

    pymol_script_path = os.path.join("results", "temp.py")
    with open(pymol_script_path, 'w') as f:
        f.write(script)
    cmd = 'which pymol'
    subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    cmd = f'pymol -ckq {pymol_script_path}'
    out, err = subprocess.Popen(cmd, shell=True, env={}, stdin=subprocess.PIPE, stdout=subprocess.DEVNULL,
                                stderr=subprocess.STDOUT).communicate()


if __name__ == "__main__":
    print('Usage: utilities.py function input')
    print('Functions: write_features, print_features')
    logging.error = print
    args = sys.argv
    globals()[args[1]](*args[2:])
