#! /usr/bin/env python3

import pickle
import subprocess
from libs import alphafold_paths, analyse, arcimboldo_air, utils
import os
import sys
import toml
import logging


def create_af2(output_dir: str, alphafold_paths: alphafold_paths.AlphaFoldPaths):

    with open(f'{output_dir}/run_af2.sh', 'w') as bash_file:
        previous_path_to_output_dir = '/'.join(output_dir.split('/')[:-1])
        name = output_dir.split('/')[-1]
        bash_file.write('#!/bin/bash\n')
        bash_file.write(f'python {os.path.dirname(os.path.abspath(__file__))}/ALPHAFOLD/run_alphafold.py \\\n')
        bash_file.write(f'--fasta_paths={name}.fasta \\\n')
        bash_file.write(f'--output_dir={previous_path_to_output_dir} \\\n')
        bash_file.write(f'--data_dir={alphafold_paths.af2_dbs_path} \\\n')
        bash_file.write(f'--uniref90_database_path={alphafold_paths.uniref90_db_path} \\\n')
        bash_file.write(f'--mgnify_database_path={alphafold_paths.mgnify_db_path} \\\n')
        bash_file.write(f'--template_mmcif_dir={alphafold_paths.mmcif_db_path} \\\n')
        bash_file.write('--max_template_date=2022-03-09 \\\n')
        bash_file.write(f'--obsolete_pdbs_path={alphafold_paths.obsolete_mmcif_db_path} \\\n')
        bash_file.write('--model_preset=monomer \\\n')
        bash_file.write(f'--bfd_database_path={alphafold_paths.bfd_db_path} \\\n')
        bash_file.write(f'--uniclust30_database_path={alphafold_paths.uniclust30_db_path} \\\n')
        bash_file.write(f'--pdb70_database_path={alphafold_paths.pdb70_db_path} \\\n')
        bash_file.write('--read_features_pkl=True\n')
        bash_file.close()

def main():

    utils.create_logger()

    if len(sys.argv) != 2:
        raise Exception('Wrong command-line arguments.')

    input_path = os.path.abspath(sys.argv[1])

    logging.info(f'Reading the configuration file for ARCIMBOLDO_AIR at {input_path}')

    if not os.path.exists(input_path):
        raise Exception('The given path for the configuration file either does not exist or you do not have the permissions to read it')
    try:
        input_load = toml.load(input_path)
    except:
        raise Exception('It has not been possible to read the input file')

    a_air = arcimboldo_air.ArcimboldoAir(parameters_dict=input_load)

    for template in a_air.templates:
        template.generate_features(a_air=a_air)
        if template.add_to_msa:
            sequence_from_template = template.template_features['template_sequence'][0].decode('utf-8')
            a_air.features.append_row_in_msa(sequence=sequence_from_template, msa_uniprot_accession_identifiers=template.pdb_id)
            logging.info(f'Sequence from template \"{template.pdb_id}\" was added to msa.')
        if template.add_to_templates:
            a_air.features.append_new_template_features(new_template_features=template.template_features, custom_sum_prob=template.sum_prob)
            logging.info(f'Template \"{template.pdb_id}\" was added to templates.')

    a_air.features.write_pkl(output_dir=f'{a_air.output_dir}/features.pkl')

    if a_air.run_af2:
        create_af2(a_air.output_dir, a_air.alphafold_paths)
        logging.info('Running AF2')
        af2_output = subprocess.Popen(['bash', f'{a_air.output_dir}/run_af2.sh'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout, stderr = af2_output.communicate()
        with open(f'{a_air.output_dir}/af2_output.log', 'w') as f:
            f.write(stdout.decode('utf-8'))
        analyse.af_summary(af_job_path=f'{a_air.output_dir}')

    features_file = open(f'{a_air.output_dir}/features.pkl','rb')
    features = pickle.load(features_file)
    features_file.close()

    for key in features.keys():
        print(key, features[key].shape)
    
    utils.clean_files(a_air.output_dir)

if __name__ == "__main__":
    main()







