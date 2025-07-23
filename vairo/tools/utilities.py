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
from Bio.PDB import PDBParser, PDBIO
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
import matplotlib.ticker as ticker
from collections import defaultdict


current_directory = os.path.dirname(os.path.abspath(__file__))
target_directory = os.path.abspath(os.path.join(current_directory, '..', '..'))
sys.path.append(target_directory)
from vairo.libs import features, bioutils, plots, utils, structures, template_modifications, global_variables


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
                        fit_ini=inf_ini,
                        fit_end=inf_end,
                        match_ini=inm_ini,
                        match_end=inm_end,
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


def lsqkab():

    def process_and_plot_results(results_dict, reference, plot_file, plot_title):
        colors = ['orange', 'blue', 'm', 'g']
        plt.figure(figsize=(12, 6))
        reference_results = results_dict.get(reference, {})
        for i, (model, pdbs) in enumerate(reference_results.items()):
            sorted_items = sorted(pdbs.items(), key=lambda x: os.path.basename(x[0]))
            y_vals = [float(rmsd) for _, rmsd in sorted_items]
            x_vals = list(range(1, len(y_vals) + 1))
            plt.plot(x_vals, y_vals, marker='o', linestyle='-', color=colors[i], label=utils.get_file_name(model))

        ax = plt.gca()
        ax.set_xlabel('MODELS')
        ax.set_ylabel('RMSD')
        ax.set_title(plot_title)
        ax.grid(True)
        ax.set_ylim(0, 10)

        ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
        ax.xaxis.set_minor_locator(ticker.NullLocator())

        plt.legend(title='Models')
        plt.tight_layout()
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')

    def split_models_in_pdb(pdb_path):
        name = utils.get_file_name(pdb_path)
        output_dir = os.path.join(os.getcwd(), name)
        utils.create_dir(name, True)
        structure = bioutils.get_structure(pdb_path)
        io = PDBIO()
        for model in structure:
            io.set_structure(model)
            out_file = os.path.join(output_dir, f"model_{model.id:03d}.pdb")
            io.save(out_file)
        return output_dir

    def extract_residue_ranges(res_list):
        results_list = []
        for seq, start_resnum, chain_id in res_list:
            end_resnum = start_resnum + len(seq) - 1
            results_list.append((chain_id, start_resnum, end_resnum))
        return results_list

    def prepare_lsqkab_input(reference_ranges, target_ranges):
        lines = []
        for reference, target in zip(reference_ranges, target_ranges):
            ref_chain, ref_start, ref_end = reference
            target_chain, target_start, target_end = target
            lines.append(f'FIT RESIDUE CA {ref_start} to {ref_end} CHAIN {ref_chain}')
            lines.append(f'MATCH {target_start} {target_end} CHAIN {target_chain}')
        return lines

    def write_lsqkab_input_file(lsqkab_input_lines, pdb1, pdb2):
        with open('lsqkab_input.cards', 'w') as f:
            f.write('lsqkab ')
            f.write(f'xyzinf {pdb1} ')
            f.write(f'xyzinm {pdb2} ')
            f.write('DELTAS deltas.out ')
            f.write('xyzout output.pdb << END-lsqkab\n')
            f.write('output deltas \n')
            f.write('output XYZ \n')
            for line in lsqkab_input_lines:
                f.write(line + '\n')
            f.write(f'end \n')
            f.write(f'END-lsqkab')
        stdout, stderr = subprocess.Popen(['bash', 'lsqkab_input.cards'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        text = stdout.decode('utf-8', errors='ignore')
        lines = text.splitlines()
        pattern = r"RMS\s+XYZ\s+DISPLACEMENT\s*=\s*([\d\.]+)"
        for line in lines:
            m = re.search(pattern, line)
            if m:
                rmsd = m.group(1)
                return rmsd
        return None

    def rmsd_ca(atoms1, atoms2):
        coords1 = np.array([a.get_coord() for a in atoms1])
        coords2 = np.array([a.get_coord() for a in atoms2])
        diff = coords1 - coords2
        return np.sqrt(np.mean(np.sum(diff * diff, axis=1)))

    def calculate_translation_list(reference_numbering, reference_pdb, target_pdb):
        reference_struct = bioutils.get_structure(reference_pdb)
        target_struct = bioutils.get_structure(target_pdb)
        reference_sequence = bioutils.extract_sequence_msa_from_pdb(reference_pdb)
        target_sequence = bioutils.extract_sequence_msa_from_pdb(target_pdb)
        translation_list = []

        for reference_chain, query_start, query_end in reference_numbering:
            query_ca = [reference_struct[0][reference_chain][i]['CA'] for i in range(query_start-1, query_end)]
            query_seq = reference_sequence[reference_chain][query_start-1:query_end]
            best_hit = None
            for target_chain, target_chain_seq in target_sequence.items():
                target_start = target_chain_seq.find(query_seq)
                while target_start != -1:
                    target_end = target_start + len(query_seq)
                    target_ca = [
                        target_struct[0][target_chain][res_id]['CA']
                        for res_id in range(target_start+1, target_end+1)
                    ]
                    score = rmsd_ca(query_ca, target_ca)
                    if best_hit is None or score < best_hit[3]:
                        best_hit = (target_chain, target_start+1, target_end, score)
                    target_start = target_chain_seq.find(query_seq, target_start + 1)
            if best_hit is not None:
                translation_list.append((best_hit[0], best_hit[1], best_hit[2]))
            else:
                sys.exit(1)
        return translation_list



    input_dict = {
        'greengreen': {
            'references':{
                'reference': 'GreenTetramerCrystal_clean.pdb',
                'vairo': 'greeboundDimer_super.pdb'
            },
            'frames': {
                'vairo': ['frames_GGxtal_200.pdb', 'frames_GGvairo_200.pdb']
            },
            'superpose_list': [('YNGKTYTANLKAD', 84, 'A'),
                              ('YNGKTYTANLKAD', 84, 'B'),
                              ('DVSFNFGSEN', 128, 'A'),
                              ('DVSFNFGSEN', 128, 'B'),
                              ('LDQNGVASLTN', 173, 'A'),
                              ('LDQNGVASLTN', 173, 'B'),
                              ]
        },
        'blueblue': {
            'references': {
                'reference': 'BlueTetramerCrystal.pdb',
                'vairo': 'blueboundDimer_clean_cut4md_renumbered_super.pdb'
            },
            'frames': {
                'vairo': ['frames_BBvairo_200.pdb', 'frames_BBxtal_200.pdb']
            },
            'superpose_list': [('NVNFYDVTSGATVTNG', 199,'B'),
                               ('NVNFYDVTSGATVTNG', 199,'A'),
                               ('AAQYADKKLNTRTANT', 241,'B'),
                               ('AAQYADKKLNTRTANT', 241,'A'),
                              ]
        },
        'greentetramer': {
            'references': {
                'reference': 'GreenTetramerCrystal_clean.pdb',
                'vairo': 'greentetramerTilefromPiecesgreenSeqMSAnomaskR0_clean_cut4md_renumbered_super.pdb'
            },
            'frames': {
                'vairo': ['frames_GTxtal_200.pdb', 'frames_GTvairo_200.pdb']
            },
            'superpose_list': [('SAVAANTANNTPAIAGNL', 59, 'A'),
                               ('SAVAANTANNTPAIAGNL', 59, 'E'),
                               ('YAINTTDNSN', 190, 'A'),
                               ('YAINTTDNSN', 190, 'E'),
                               ('SVNADNQGQVNVANVVAAINSKYF', 217, 'B'),
                               ('SVNADNQGQVNVANVVAAINSKYF', 217, 'D'),
                               ('LKDQKIDVNSVGYFKAPHTFTV', 264, 'B'),
                               ('LKDQKIDVNSVGYFKAPHTFTV', 264, 'D'),
                               ('YNGKTYTANLKAD', 84, 'A'),
                               ('YNGKTYTANLKAD', 84, 'B'),
                               ('YNGKTYTANLKAD', 84, 'D'),
                               ('YNGKTYTANLKAD', 84, 'E'),
                               ('DVSFNFGSEN', 128, 'A'),
                               ('DVSFNFGSEN', 128, 'B'),
                               ('DVSFNFGSEN', 128, 'D'),
                               ('DVSFNFGSEN', 128, 'E'),
                               ('LDQNGVASLTN', 173, 'A'),
                               ('LDQNGVASLTN', 173, 'B'),
                               ('LDQNGVASLTN', 173, 'D'),
                               ('LDQNGVASLTN', 173, 'E'),
                                ]
        },
        'bluetetramer': {
            'references': {
                'reference': 'BlueTetramerCrystal.pdb',
                'vairo': 'blueboundTetramerfromPiecesSeqgbTR1_clean_cut4md_renumbered_super.pdb'
            },
            'frames': {
                'vairo': ['frames_BTxtal_200.pdb', 'frames_BTvairo_200.pdb'],
            },
            'superpose_list': [('SAVAANTANNTPAIAGNL', 59, 'E'),
                       ('SAVAANTANNTPAIAGNL', 59, 'A'),
                       ('YAINTTDNSN', 190, 'E'),
                       ('YAINTTDNSN', 190, 'A'),
                       ('SVNADNQGQVNVANVVAAINSKYF', 217,'B'),
                       ('SVNADNQGQVNVANVVAAINSKYF', 217,'D'),
                       ('LKDQKIDVNSVGYFKAPHTFTV', 264,'B'),
                       ('LKDQKIDVNSVGYFKAPHTFTV', 264,'D'),
                       ('NVNFYDVTSGATVTNG', 199,'B'),
                       ('NVNFYDVTSGATVTNG', 199,'D'),
                       ('NVNFYDVTSGATVTNG', 199,'A'),
                       ('NVNFYDVTSGATVTNG', 199,'E'),
                       ('AAQYADKKLNTRTANT', 241,'B'),
                       ('AAQYADKKLNTRTANT', 241,'D'),
                       ('AAQYADKKLNTRTANT', 241,'A'),
                       ('AAQYADKKLNTRTANT', 241 ,'E'),
                       ]
        },
    }
    generate_list = ['greengreen', 'blueblue', 'bluetetramer', 'greemtetramer']
    #generate_list = ['greengreen']
    for generate in generate_list:
        old_path = os.getcwd()
        path = os.path.join(os.getcwd(), generate)
        os.chdir(path)
        superpose_list = input_dict[generate]['superpose_list']
        references = input_dict[generate]['references']
        frames = input_dict[generate]['frames']

        reference_superpose_list = extract_residue_ranges(superpose_list)
        results_dict = {ref: {} for ref in references.values()}
        superpose_translation_dict = {'reference': reference_superpose_list}
        superpose_translation_dict.update({
            key: calculate_translation_list(reference_superpose_list, references['reference'], ref)
            for key, ref in references.items()
            if key != 'reference'
        })

        for frame_type, frame_list in frames.items():
            for frame_value in frame_list:
                dir_path = split_models_in_pdb(frame_value)
                for reference_type, reference_values in references.items():
                    results_dict[reference_values][frame_value] = {}
                    lsqkab_input = prepare_lsqkab_input(superpose_translation_dict[frame_type], superpose_translation_dict[reference_type])
                    for pdb in os.listdir(dir_path):
                        pdb_path = os.path.join(dir_path, pdb)
                        rmsd = write_lsqkab_input_file(lsqkab_input, reference_values, pdb_path)
                        results_dict[reference_values][frame_value][pdb_path] = rmsd
        process_and_plot_results(results_dict, references['reference'], 'reference_crystal','Reference Crystal')
        process_and_plot_results(results_dict, references['vairo'], 'reference_vairo', 'Reference Vairo')
        os.chdir(old_path)

if __name__ == "__main__":
    print('Usage: utilities.py function input')
    print('Functions: write_features, print_features')
    logging.error = print
    args = sys.argv
    globals()[args[1]](*args[2:])
