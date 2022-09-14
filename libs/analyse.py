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
    plt.ylabel('pLLDT')
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

    template_dict = a_air.features.write_all_templates_in_features(output_dir=templates_path)
    ranked_models_dict = {utils.get_file_name(ranked): os.path.join(a_air.run_dir, ranked) for ranked in os.listdir(a_air.run_dir) if re.match('ranked_[0-9]+.pdb', ranked)}
    if not bool(ranked_models_dict):
        raise Exception('No ranked PDBs found')

    pllddt_dict = plot_plddt(plot_path=plddt_plot_path, ranked_models_dict=ranked_models_dict)
    max_plldt = max(pllddt_dict.values())


    ranked_filtered = []
    for ranked, ranked_path in ranked_models_dict.items():
        new_ranked_path = os.path.join(a_air.run_dir, f'trimmed_{os.path.basename(ranked_path)}')
        ranked_models_dict[ranked] = new_ranked_path
        bioutils.split_chains_assembly(pdb_in_path=ranked_path, pdb_out_path=new_ranked_path, a_air=a_air)

        if pllddt_dict[ranked] >= (PERCENTAGE_FILTER*max_plldt):
            ranked_filtered.append(ranked)
            interfaces_dict = bioutils.find_interface_from_pisa(new_ranked_path)
            for interfaces in interfaces_dict:
                chain1, chain2 = list(interfaces)
                dimers_path = os.path.join(interfaces_path, f'{ranked}_{chain1}{chain2}.pdb')
                bioutils.split_dimers_in_pdb(pdb_in_path=new_ranked_path,
                                    pdb_out_path=dimers_path,
                                    chain1=chain1,
                                    chain2=chain2)

    for template, template_path in template_dict.items():
        bioutils.split_chains_assembly(pdb_in_path=template_path, pdb_out_path=template_path, a_air=a_air)

    rmsd_dict = {}
    num = 0
    for template, template_path in template_dict.items():
        res_list_length = len([res for res in Selection.unfold_entities(PDBParser().get_structure(template, template_path), 'R')])
        results_list = []
        for ranked, ranked_path in utils.sort_by_digit(ranked_models_dict):
            output_superposition = True if ranked in ranked_filtered else False
            if ranked in ranked_filtered:
                output_pdb = os.path.join(a_air.output_dir, f'{ranked}_{template}.pdb')
            else:
                output_pdb = None
            rmsd, nalign, quality_q = bioutils.superpose_pdbs(query_pdb=ranked_path,
                                                            target_pdb=template_path,
                                                            output_pdb=output_pdb)
                                                            
            results_list.append(f'{rmsd}, {nalign} ({res_list_length})')
        if template in rmsd_dict:
            num = num + 1
            rmsd_dict[f'{template}_({num})'] = results_list
        else:
            rmsd_dict[f'{template}'] = results_list

    secondary_dict = {}
    for ranked, ranked_path in utils.sort_by_digit(ranked_models_dict):
        aleph_file = os.path.join(a_air.run_dir, f'aleph_{ranked}.txt')
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




