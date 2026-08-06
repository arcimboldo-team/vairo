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
import itertools
import ast
from alphafold.data import parsers, templates, mmcif_parsing
from Bio.Blast import NCBIWWW
from alphafold.common import residue_constants
from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain
import yaml

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


def renumber_template(template_path, output_path, chain_offsets=None):
    parser = PDBParser(QUIET=True)
    original_structure = parser.get_structure("original", template_path)
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
            offset = int((chain_offsets or {}).get(chain.id, 0))
            print(f"  Chain {chain.id}: {residues[0].id[1]}..{residues[-1].id[1]}, shift {offset} "
                  f"-> native {offset:+d} = {residues[0].id[1] + offset}..{residues[-1].id[1] + offset}")
            for res in residues:
                res.parent = None
                old_id = res.id
                new_number = old_id[1] + offset
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


def _write_shift_summary(results_dir, run_results):
    summary_path = os.path.join(results_dir, "shift_summary.csv")
    with open(summary_path, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(["combo_name", "sequence"])
        for run in run_results:
            writer.writerow([run['combo_name'], run.get('sequence', '')])
    print(f"  -> Wrote shift summary: {summary_path}")
    return summary_path


def _grouped_combinations(step_ranges_list, chain_order, chain_groups):
    range_by_chain = dict(zip(chain_order, step_ranges_list))
    groups, seen = [], set()
    for g in chain_groups:
        members = [c for c in g if c in range_by_chain and c not in seen]
        if members:
            groups.append(members)
            seen.update(members)
    for c in chain_order:
        if c not in seen:
            groups.append([c])
            seen.add(c)
    per_group = [] 
    for g in groups:
        ranges = [range_by_chain[c] for c in g]
        max_len = max((len(r) for r in ranges), default=0)
        padded = [r + [r[-1]] * (max_len - len(r)) if r else [0] * max_len for r in ranges]
        per_group.append([dict(zip(g, vals)) for vals in zip(*padded)])
    combos = []
    for pick in itertools.product(*per_group):
        merged = {}
        for d in pick:
            merged.update(d)
        combos.append(tuple(merged.get(c, 0) for c in chain_order))
    return combos


def generate_shift_models(fasta_path, template_path, chain_steps_dict, bor_file='', use_product=True,
                          one_chain=False, output_dir=None, chain_query_offsets=None, chain_groups=None,
                          use_scwrl=False):

    sequences_dict = dict(bioutils.extract_sequences(fasta_path).items())
    sequence_template_dict = dict(bioutils.extract_sequence_msa_from_pdb(template_path).items())

    results_dir = output_dir if output_dir else os.path.join(os.getcwd(), "results")
    os.makedirs(results_dir, exist_ok=True)

    scwrl_candidates = ('Scwrl4', 'SCWRL4', 'scwrl4', 'scwrl', 'SCWRL')
    scwrl_path = next((shutil.which(name) for name in scwrl_candidates if shutil.which(name)), None)
    if use_scwrl and not scwrl_path:
        raise FileNotFoundError("use_scwrl is set but SCWRL was not found on PATH "
                                f"(looked for {scwrl_candidates})")

    output_path = os.path.join(results_dir, 'fixed.pdb')
    shutil.copy2(template_path, output_path)
    bioutils.run_pdbfixer(output_path, output_path)
    template_path = output_path

    template_struct = bioutils.get_structure(template_path)
    chain_resnums = {ch.id: [res.id[1] for res in ch if res.id[0] == ' '] for ch in template_struct[0]}

    template_dict = {}
    run_results = []
    generated_configs = []
    for i, (chain, template_sequence) in enumerate(sequence_template_dict.items()):
        first_pos = None
        last_pos = None

        for index, item in enumerate(template_sequence):
            if '-' not in item:
                first_pos = index
                break
        for index in range(len(template_sequence) - 1, -1, -1):
            if '-' not in template_sequence[index]:
                last_pos = index
                break

        if first_pos is not None and last_pos is not None:
            if chain in chain_steps_dict:
                template_dict[chain] = {
                    'is_target': True,
                    'seq': sequences_dict[chain],
                    'first': first_pos,
                    'last': last_pos
                }
            else:
                wt_sequence = "".join([res for res in template_sequence if '-' not in res])
                template_dict[chain] = {
                    'is_target': False,
                    'seq': wt_sequence
                }

    print(f"Template mapping: {template_dict}")

    chain_order = list(template_dict.keys())
    step_ranges_list = []

    for chain in chain_order:
        if chain not in chain_steps_dict:
            step_ranges_list.append([0])
            continue
        spec = chain_steps_dict[chain]
        if isinstance(spec, (list, tuple)):
            lo, hi = int(spec[0]), int(spec[1])
        else:
            n = int(spec)
            lo, hi = (0, n) if n >= 0 else (n, 0)
        if hi < lo:
            lo, hi = hi, lo
        temp = template_dict[chain]
        expected_len = temp['last'] - temp['first'] + 1
        lo_feas = -temp['first']
        hi_feas = len(temp['seq']) - expected_len - temp['first']
        clo, chi = max(lo, lo_feas), min(hi, hi_feas)
        if (clo, chi) != (lo, hi):
            if chi >= clo:
                used = f"[{clo},{chi}]"
            else:
                used = "[0,0] (requested range is entirely outside the feasible window; using shift 0 only)"
            print(f"Warning: shifts for chain {chain} clamped from [{lo},{hi}] to {used}; "
                  f"feasible window [{lo_feas},{hi_feas}] (fasta len {len(temp['seq'])}, region {expected_len}).")
        step_ranges_list.append(list(range(clo, chi + 1)) if chi >= clo else [0])

    if chain_groups:
        all_combinations = _grouped_combinations(step_ranges_list, chain_order, chain_groups)
    elif use_product:
        all_combinations = list(itertools.product(*step_ranges_list))
    else:
        max_len = max((len(r) for r in step_ranges_list), default=0)
        padded_ranges = [
            r + [r[-1]] * (max_len - len(r)) if r else '' * max_len
            for r in step_ranges_list
        ]
        all_combinations = list(zip(*padded_ranges))

    print(f"Total combinatorial steps to evaluate: {len(all_combinations)}")

    for combo in all_combinations:
        current_steps = dict(zip(chain_order, combo))
        general_sequence = ''
        valid_step = True
        step_name_parts = []

        for chain in chain_order:
            temp = template_dict[chain]
            shift = current_steps[chain]
            if chain in chain_steps_dict:
                step_name_parts.append(f"{chain}{shift}")
                expected_len = temp['last'] - temp['first'] + 1
                start_slice = temp['first'] + shift
                end_slice = start_slice + expected_len
                if start_slice < 0 or end_slice > len(temp['seq']):
                    print(f"Warning: Index out of bounds for chain {chain} at shift {shift}. Skipping combination.")
                    valid_step = False
                    break
                general_sequence += temp['seq'][start_slice: end_slice]
            else:
                general_sequence += temp['seq']

        if not valid_step:
            continue

        if one_chain and chain_query_offsets:
            occupied = set()
            overlap = False
            for chain in chain_order:
                if chain not in chain_steps_dict:
                    continue
                off = chain_query_offsets.get(chain)
                if off is None or off < 0:
                    continue
                temp = template_dict[chain]
                shift = current_steps[chain]
                length = temp['last'] - temp['first'] + 1
                start_q = off + temp['first'] + shift
                numbers = set(range(start_q, start_q + length))
                if numbers & occupied:
                    overlap = True
                    break
                occupied |= numbers
            if overlap:
                print(f"Warning: segments overlap in the query sequence for combination "
                      f"{'_'.join(step_name_parts)}. Skipping combination.")
                continue

        combo_name = "_".join(step_name_parts) if step_name_parts else "NM"
        name = f"model_{combo_name}"

        new_fasta_path = os.path.join(results_dir, f"{name}.fasta")
        output_pdb_path = os.path.join(results_dir, f"{name}.pdb")

        with open(new_fasta_path, 'w') as f:
            f.write(general_sequence)

        if bor_file != '':
            try:
                with open(bor_file, 'r') as file:
                    config = yaml.safe_load(file)
                    if 'output_dir' in config:
                        config['output_dir'] = os.path.join(results_dir, f"{name}_output")
                    if 'templates' in config:
                        if 'pdb' in config['templates'][0]:
                            config['templates'][0]['pdb'] = os.path.join(results_dir, f"{name}.pdb")
                output_bor_path = os.path.join(results_dir, f"{name}.bor")
                with open(output_bor_path, 'w') as file:
                    yaml.dump(config, file, default_flow_style=False, sort_keys=False)
                generated_configs.append(output_bor_path)
            except Exception as e:
                print('Not possible to read bor file.')


        if use_scwrl:
            cmd = [scwrl_path, '-i', template_path, '-o', output_pdb_path, '-s', new_fasta_path]
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            print(f"  -> SCWRL threaded {combo_name}")
        else:
            modifications = template_modifications.TemplateModifications()
            for chain in chain_order:
                if chain not in chain_steps_dict:
                    continue
                temp = template_dict[chain]
                expected_len = temp['last'] - temp['first'] + 1
                start_slice = temp['first'] + current_steps[chain]
                chain_seq = temp['seq'][start_slice: start_slice + expected_len]
                mutations = []
                for k, resseq in enumerate(chain_resnums.get(chain, [])):
                    if k >= len(chain_seq):
                        break
                    three = residue_constants.restype_1to3.get(chain_seq[k])
                    if three:
                        mutations.append(template_modifications.ResidueMutate(
                            mutate_residues_number=[resseq], mutate_with=three))
                if mutations:
                    modifications.append_chain_modification(
                        template_modifications.ChainModifications(chain=chain, mutations=mutations))

            modifications.modify_template(pdb_in_path=template_path, pdb_out_path=output_pdb_path,
                                          type_modify=['mutate'])
            print(f"  -> threaded {combo_name}: mutated {len(modifications.modifications_list)} chain(s)")

        run_results.append({'combo_name': combo_name, 'sequence': general_sequence})

        if one_chain:
            base_offsets = chain_query_offsets or {}
            merge_offsets = {chain: int(base_offsets.get(chain, 0) or 0) + int(shift)
                             for chain, shift in current_steps.items()}
            bioutils.merge_all_chains_into_one(output_pdb_path, output_pdb_path,
                                               chain_offsets=merge_offsets)
        else:
            renumber_template(output_pdb_path, output_pdb_path, chain_offsets=current_steps)

    print(f"\nGenerated {len(run_results)} shifted models under {results_dir}")
    _write_shift_summary(results_dir, run_results)
    return generated_configs


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
    print('Usage: utilities.py function [positional_args] [key=value_args]')
    import sys
    import ast
    args = sys.argv
    func_name = args
    pos_args = []
    kwargs = {}
    for arg in args[2:]:
        if '=' in arg:
            key, value = arg.split('=', 1)
            try:
                kwargs[key] = ast.literal_eval(value)
                print(f"Parsed kwarg: {key} -> {kwargs[key]}")
            except (ValueError, SyntaxError):
                kwargs[key] = value
        else:
            try:
                parsed_arg = ast.literal_eval(arg)
                pos_args.append(parsed_arg)
            except (ValueError, SyntaxError):
                pos_args.append(arg)

    print(f"Executing: {func_name}")
    globals()[func_name[1]](*pos_args, **kwargs)