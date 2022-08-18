import os
import shutil
import pandas as pd
from libs import bioutils, plots
from libs.features import Features
from Bio.PDB import PDBParser, Selection

def analyse_output(output_dir: str, run_dir: str, features: Features):
    
    plots_path = f'{output_dir}/plots'
    templates_path = f'{output_dir}/templates'
    rmsd_path = f'{plots_path}/rmsd.log'

    if os.path.exists(plots_path):
        shutil.rmtree(plots_path)
    os.makedirs(plots_path)

    if os.path.exists(templates_path):
        shutil.rmtree(templates_path)
    os.makedirs(templates_path)

    ranked_models_dict = {os.path.basename(ranked): os.path.join(run_dir, ranked) for ranked in os.listdir(run_dir) if ranked.startswith('ranked_')}
    template_dict = features.write_all_templates_in_features(templates_path)

    for ranked, ranked_path in ranked_models_dict.items():
        shutil.copy2(ranked_path, output_dir)
    
    plots.create_plots(plots_path=plots_path, ranked_models_dict=ranked_models_dict)

    results_dict = {}
    num = 0
    for template, template_path in template_dict.items():
        res_list_length = len([res for res in Selection.unfold_entities(PDBParser().get_structure(template, template_path), 'R')])
        results_list = []
        for ranked, ranked_path in sorted(ranked_models_dict.items(), key=lambda x: int("".join([i for i in x[0] if i.isdigit()]))):
            rmsd, nalign, quality_q, aligned_res_list = bioutils.superpose_pdbs(query_pdb=ranked_path,
                                                                       target_pdb=template_path,
                                                                       output_superposition=False)
            results_list.append(f'{rmsd}, {nalign} ({res_list_length})')
        if template in results_dict:
            num = num + 1
            results_dict[f'{template}_({num})'] = results_list
        else:
            results_dict[f'{template}'] = results_list

    with open(rmsd_path, 'w') as f:
        rows = []
        for key in results_dict.keys():
            rows.append([key] + results_dict[key])
        df = pd.DataFrame(rows, columns=['template'] + list((sorted(ranked_models_dict.keys(), key=lambda x: int("".join([i for i in x if i.isdigit()]))))))
        f.write(df.to_markdown())
    f.close()

    #ranked_0, ranked_0_path = template_dict.pop('ranked_0')
    #for ranked, ranked_path in template_dict:
    #    bioutils.superpose_pdbs(query_pdb=ranked_path,
    #                   target_pdb=ranked_0_path,
    #                   output_superposition=True)




