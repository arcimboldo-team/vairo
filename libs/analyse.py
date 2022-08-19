import os
import shutil
import sys
import pandas as pd
from ALEPH.aleph.core import ALEPH
from libs import bioutils, plots, utils
from libs.features import Features
from Bio.PDB import PDBParser, Selection

def analyse_output(output_dir: str, run_dir: str, features: Features):
    
    plots_path = f'{output_dir}/plots'
    templates_path = f'{output_dir}/templates'
    analysis_path = f'{plots_path}/analysis.txt'

    utils.create_dir(dir_path=plots_path,delete_if_exists=True)
    utils.create_dir(dir_path=plots_path,delete_if_exists=True)

    if os.path.exists(templates_path):
        shutil.rmtree(templates_path)
    os.makedirs(templates_path)

    ranked_models_dict = {utils.get_file_name(ranked): os.path.join(run_dir, ranked) for ranked in os.listdir(run_dir) if ranked.startswith('ranked_') and ranked.endswith('.pdb')}
    template_dict = features.write_all_templates_in_features(templates_path)

    for ranked, ranked_path in ranked_models_dict.items():
        shutil.copy2(ranked_path, output_dir)
    
    plots.create_plots(plots_path=plots_path, ranked_models_dict=ranked_models_dict)

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
        secondary_dict[ranked] = utils.parse_aleph_annotate(file_path=aleph_file)

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


    #ranked_0, ranked_0_path = template_dict.pop('ranked_0')
    #for ranked, ranked_path in template_dict:
    #    bioutils.superpose_pdbs(query_pdb=ranked_path,
    #                   target_pdb=ranked_0_path,
    #                   output_superposition=True)




