import os
import re
import shutil
import subprocess
import sys
from typing import Dict
import pandas as pd
from ALEPH.aleph.core import ALEPH
from libs import bioutils, utils
from Bio.PDB import PDBParser, Selection
import matplotlib.pyplot as plt


def plot_plddt(plot_path: str, ranked_models_dict: Dict):
    
    plt.clf()
    for ranked, ranked_path in utils.sort_by_digit(ranked_models_dict):
        plddt_list = []
        with open(ranked_path) as f:
            for line in f.readlines():
                if line[:4] == 'ATOM' and line[13:16] == 'CA ':
                    plddt_list.append(float(line[60:66].replace(" ", "")))
        res_list = [int(item) for item in range(1, len(plddt_list) + 1)]
        plt.plot(res_list, plddt_list, label=ranked)
    plt.legend()
    plt.xlabel('residue number')
    plt.ylabel('pLLDT')
    plt.savefig(plot_path)

def find_interface_from_pisa(pdb_in_path: str) -> Dict:

    subprocess.Popen(['pisa', 'temp', '-analyse', pdb_in_path],
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE).communicate()

    pisa_output = subprocess.Popen(['pisa', 'temp', '-list', 'interfaces'], stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE).communicate()[0].decode('utf-8')

    match1 = [m.start() for m in re.finditer(" LIST OF INTERFACES", pisa_output)][0]
    match2 = [m.start() for m in re.finditer(" ##:  serial number", pisa_output)][0]

    interfaces_dict = {}
    for line in pisa_output[match1:match2].split('\n')[4:-2]:
        area = line.split('|')[3][:8].replace(' ', '')
        deltaG = line.split('|')[3][8:15].replace(' ', '')
        chain1 = line.split('|')[1].replace(' ', '')
        chain2 = line.split('|')[2].split()[0].replace(' ', '')
        interfaces_dict[f'{chain1}{chain2}'] =  (area, deltaG)

    return interfaces_dict

def split_dimers_in_pdb(pdb_in_path, pdb_out_path, chain1, chain2):

    with open(pdb_in_path, 'r') as f_in:
        input_lines = f_in.readlines()

    with open(pdb_out_path, 'w') as f_out:
        for line in input_lines:
            if line[21:22] in [chain1, chain2]:
                f_out.write(line)

def analyse_output(output_dir: str, run_dir: str, a_air):
    
    plots_path = f'{output_dir}/plots'
    templates_path = f'{output_dir}/templates'
    interfaces_path = f'{output_dir}/interfaces'
    analysis_path = f'{plots_path}/analysis.txt'
    aleph_results_path = f'{run_dir}/output.json'
    plddt_plot_path = f'{plots_path}/plddt.png'

    utils.create_dir(dir_path=plots_path,delete_if_exists=True)
    utils.create_dir(dir_path=plots_path,delete_if_exists=True)

    if os.path.exists(interfaces_path):
        shutil.rmtree(interfaces_path)
    os.makedirs(interfaces_path)

    if os.path.exists(templates_path):
        shutil.rmtree(templates_path)
    os.makedirs(templates_path)

    ranked_models_dict = {utils.get_file_name(ranked): os.path.join(run_dir, ranked) for ranked in os.listdir(run_dir) if re.match('ranked_[0-9]+.pdb', ranked)}
    if not bool(ranked_models_dict):
        raise Exception('No ranked PDBs found')
    
    template_dict = a_air.features.write_all_templates_in_features(templates_path)

    for ranked, ranked_path in ranked_models_dict.items():
        shutil.copy(ranked_path, output_dir)
        ranked_models_dict[ranked] = os.path.join(output_dir, os.path.basename(ranked_path))
    
    plot_plddt(plot_path=plddt_plot_path, ranked_models_dict=ranked_models_dict)

    rmsd_dict = {}
    num = 0
    for template, template_path in template_dict.items():
        res_list_length = len([res for res in Selection.unfold_entities(PDBParser().get_structure(template, template_path), 'R')])
        results_list = []
        for ranked, ranked_path in utils.sort_by_digit(ranked_models_dict):
            rmsd, nalign, quality_q, aligned_res_list = bioutils.superpose_pdbs(query_pdb=ranked_path,
                                                                       target_pdb=template_path,
                                                                       output_superposition=False)
            results_list.append(f'{rmsd}, {nalign} ({res_list_length})')
        if template in rmsd_dict:
            num = num + 1
            rmsd_dict[f'{template}_({num})'] = results_list
        else:
            rmsd_dict[f'{template}'] = results_list

    secondary_dict = {}
    for ranked, ranked_path in utils.sort_by_digit(ranked_models_dict):
        aleph_file = os.path.join(run_dir, f'aleph_{ranked}.txt')
        with open(aleph_file, 'w') as sys.stdout:
            ALEPH.annotate_pdb_model_with_aleph(ranked_path)
        sys.stdout = sys.__stdout__
        if os.path.exists(aleph_results_path):
            secondary_dict[ranked] = utils.parse_aleph_annotate(file_path=aleph_results_path)

    with open(analysis_path, 'w') as f_in:

        ##rmsd
        if bool(rmsd_dict):
            rows = []
            for key in rmsd_dict.keys():
                rows.append([key] + rmsd_dict[key])
            df = pd.DataFrame(rows, columns=['template'] + utils.sort_by_digit(list(ranked_models_dict)))
            f_in.write(df.to_markdown())

        f_in.write('\n\n')

        if bool(secondary_dict):
            rows = []
            for key in secondary_dict.keys():
                rows.append([key] + list(secondary_dict[key].values()))
            df = pd.DataFrame(rows, columns=['ranked'] + list(secondary_dict[key]))
            f_in.write(df.to_markdown())

    for ranked, ranked_path in utils.sort_by_digit(ranked_models_dict):
        bioutils.split_chains_assembly(pdb_in_path=ranked_path, pdb_out_path=ranked_path, a_air=a_air)
        interfaces_dict = find_interface_from_pisa(ranked_path)
        for interfaces in interfaces_dict:
            chain1, chain2 = list(interfaces)
            dimers_path = os.path.join(interfaces_path, f'{ranked}_{chain1}{chain2}.pdb')
            split_dimers_in_pdb(pdb_in_path=ranked_path,
                                pdb_out_path=dimers_path,
                                chain1=chain1,
                                chain2=chain2)


    #ranked_0, ranked_0_path = template_dict.pop('ranked_0')
    #for ranked, ranked_path in template_dict:
    #    bioutils.superpose_pdbs(query_pdb=ranked_path,
    #                   target_pdb=ranked_0_path,
    #                   output_superposition=True)




