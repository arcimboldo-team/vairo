import logging
import os
import re
import shutil
import sys
import statistics
from typing import Dict, List
import pandas as pd
from ALEPH.aleph.core import ALEPH
from libs import bioutils, utils
from Bio.PDB import PDBParser, Selection
import matplotlib.pyplot as plt

PERCENTAGE_FILTER = 0.8

def plot_plddt(plot_path: str, ranked_models_dict: Dict) -> Dict:

    return_plddt_dict = {}

    plt.clf()
    for ranked, ranked_path in ranked_models_dict.items():
        plddt_list = []
        with open(ranked_path) as f:
            for line in f.readlines():
                if line[:4] == 'ATOM' and line[13:16] == 'CA ':
                    plddt_list.append(float(line[60:66].replace(" ", "")))
        res_list = [int(item) for item in range(1, len(plddt_list) + 1)]
        return_plddt_dict[ranked] = statistics.median(map(float, plddt_list))
        plt.plot(res_list, plddt_list, label=ranked)
    plt.legend()
    plt.xlabel('residue number')
    plt.ylabel('pLDDT')
    plt.savefig(plot_path)
    
    return return_plddt_dict

def analyse_output(a_air):
    
    plots_path = f'{a_air.output_dir}/plots'
    templates_path = f'{a_air.output_dir}/templates'
    interfaces_path = f'{a_air.output_dir}/interfaces'
    analysis_path = f'{plots_path}/analysis.txt'
    aleph_results_path = f'{a_air.run_dir}/output.json'
    plddt_plot_path = f'{plots_path}/plddt.png'

    utils.create_dir(dir_path=plots_path,delete_if_exists=True)
    utils.create_dir(dir_path=templates_path,delete_if_exists=True)
    utils.create_dir(dir_path=interfaces_path,delete_if_exists=True)

    ##Read all templates and rankeds, if there are no ranked, raise an error
    template_dict = a_air.features.write_all_templates_in_features(output_dir=templates_path)
    ranked_models_dict = {utils.get_file_name(ranked): os.path.join(a_air.run_dir, ranked) for ranked in os.listdir(a_air.run_dir) if re.match('ranked_[0-9]+.pdb', ranked)}
    if not bool(ranked_models_dict):
        logging.info('No ranked PDBs found')
        return

    ##Create a plot with the ranked plddts, also, calculate the maximum plddt
    plddt_dict = plot_plddt(plot_path=plddt_plot_path, ranked_models_dict=ranked_models_dict)
    max_plddt = max(plddt_dict.values())

    ##Split the templates with chains
    for template, template_path in template_dict.items():
        bioutils.split_chains_assembly(pdb_in_path=template_path, pdb_out_path=template_path, a_air=a_air)

    ##Filter rankeds, split them in chains.
    ranked_filtered = []
    for ranked, ranked_path in ranked_models_dict.items():
        new_pdb_path = os.path.join(a_air.run_dir, f'splitted_{os.path.basename(ranked_path)}')
        bioutils.split_chains_assembly(pdb_in_path=ranked_path, pdb_out_path=new_pdb_path, a_air=a_air)
        if plddt_dict[ranked] >= (PERCENTAGE_FILTER*max_plddt):
            ranked_filtered.append(ranked)
        ranked_models_dict[ranked] = new_pdb_path

    ##Superpose each template with all the rankeds.
    rmsd_dict = {}
    if template_dict:
        for ranked, ranked_path in utils.sort_by_digit(ranked_models_dict):
            res_list_length = len([res for res in Selection.unfold_entities(PDBParser().get_structure(template, template_path), 'R')])
            for template, template_path in template_dict.items():
                rmsd, nalign, quality_q = bioutils.superpose_pdbs([ranked_path, template_path])
                try:
                    rmsd_dict[f'{template}'].append(f'{rmsd}, {nalign} ({res_list_length})')
                except KeyError:
                    rmsd_dict[f'{template}'] = [f'{rmsd}, {nalign} ({res_list_length})']

            if ranked in ranked_filtered:
                output_pdb = os.path.join(a_air.output_dir, f'{ranked}_superposed.pdb')
                bioutils.superpose_pdbs([ranked_path] + list(template_dict.values()), output_pdb=output_pdb)
    
    ##Use aleph to generate domains and calculate secondary structure percentage
    secondary_dict = {}
    for ranked, ranked_path in utils.sort_by_digit(ranked_models_dict):
        aleph_file = os.path.join(a_air.run_dir, f'aleph_{ranked}.txt')
        with open(aleph_file, 'w') as sys.stdout:
            ALEPH.annotate_pdb_model(reference=ranked_path, strictness_ah=0.45, strictness_bs=0.2, peptide_length=3, 
                                    width_pic=1, height_pic=1, write_graphml=False, write_pdb=True)
        sys.stdout = sys.__stdout__
        if os.path.exists(aleph_results_path):
            secondary_dict[ranked] = utils.parse_aleph_annotate(file_path=aleph_results_path)
            aleph_txt_path = f'{a_air.run_dir}/aleph_{ranked}.txt'
            domains_dict = utils.parse_aleph_ss(aleph_txt_path)
        else:
            break

        if ranked in ranked_filtered and os.path.exists(aleph_txt_path):
            interfaces_data_list = bioutils.find_interface_from_pisa(ranked_path, interfaces_path)
            if interfaces_data_list:
                deltas_list = [interface['deltaG'] for interface in interfaces_data_list]
                deltas_list = utils.normalize_list([deltas_list])
                for i, interface in enumerate(interfaces_data_list):
                    interface['bfactor'] = deltas_list[i]
                    if not ((float(interface['se_gain1']) >= 0) == (float(interface['se_gain2']) >= 0)):
                        interface['bfactor'] = abs(max(deltas_list)) * 2
                    bioutils.create_interface_domain(ranked_path, interface, interfaces_path, domains_dict)

    ##Superpose the experimental pdb with all the rankeds and templates
    experimental_dict = {}
    if a_air.experimental_pdb is not None:
        for ranked, ranked_path in ranked_models_dict.items():
            rmsd, nalign, quality_q = bioutils.superpose_pdbs([ranked_path, a_air.experimental_pdb])
            experimental_dict[ranked] = rmsd
        for template, template_path in template_dict.items():
            rmsd, nalign, quality_q = bioutils.superpose_pdbs([template_path, a_air.experimental_pdb])
            experimental_dict[template] = rmsd       

    with open(analysis_path, 'w') as f_in:

        if bool(rmsd_dict):
            f_in.write('\n')
            f_in.write('Superpositions of rankeds and templates\n')
            rows = []
            for key in rmsd_dict.keys():
                rows.append([key] + rmsd_dict[key])
            df = pd.DataFrame(rows, columns=['template'] + utils.sort_by_digit(list(ranked_models_dict)))
            f_in.write(df.to_markdown())

        if bool(secondary_dict):
            f_in.write('\n\n')
            f_in.write('Secondary structure percentages calculated with ALEPH\n')
            rows = []
            for key in secondary_dict.keys():
                rows.append([key] + list(secondary_dict[key].values()))
            df = pd.DataFrame(rows, columns=['ranked'] + list(secondary_dict[key]))
            f_in.write(df.to_markdown())

        if bool(plddt_dict):
            f_in.write('\n\n')
            f_in.write('rankeds PLDDT\n')
            rows = []
            for key in plddt_dict.keys():
                rows.append([key, plddt_dict[key]])
            df = pd.DataFrame(rows, columns=['ranked', 'plddt'])
            f_in.write(df.to_markdown())

        if bool(experimental_dict):
            f_in.write('\n\n')
            f_in.write(f'Superposition with experimental structure {a_air.experimental_pdb}\n')
            rows = []
            for key in experimental_dict.keys():
                rows.append([key, experimental_dict[key]])
            df = pd.DataFrame(rows, columns=['pdb', 'rmsd'])
            f_in.write(df.to_markdown())

        f_in.write('\n\n')




