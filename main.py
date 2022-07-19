#! /usr/bin/env python3

from libs import arcimboldo_air, utils
import os
import sys
import toml
import logging


def run_alphafold():

    delattr(ALPHAFOLD.run_alphafold.flags.FLAGS,'data_dirs')
    ALPHAFOLD.run_alphafold.flags.DEFINE_string('data_dir', "data", 'Path to directory of supporting data.')

def af2_sh(output_dir):
    #try:
    #    print('\n')

    #except:
    #    raise ValueError('Please, set all AF2 databases manually')

    with open(f'{output_dir}/run_af2.sh', 'w') as bash_file:
        previous_path_to_output_dir = '/'.join(output_dir.split('/')[:-1])
        name = output_dir.split('/')[-1]
        bash_file.write('#!/bin/bash\n')
        bash_file.write(f'python ./alphafold/run_alphafold.py \\\n')
        bash_file.write(f'--fasta_paths={name}.fasta \\\n')
        bash_file.write(f'--output_dir={previous_path_to_output_dir} \\\n')
        bash_file.write(f'--data_dir={AF2_DBS_PATH} \\\n')
        bash_file.write(f'--uniref90_database_path={UNIREF90_DB_PATH} \\\n')
        bash_file.write(f'--mgnify_database_path={MGNIFY_DB_PATH} \\\n')
        bash_file.write(f'--template_mmcif_dir={MMCIF_DB_PATH} \\\n')
        bash_file.write('--max_template_date=2022-03-09 \\\n')
        bash_file.write(f'--obsolete_pdbs_path={OBSOLETE_MMCIF_DB_PATH} \\\n')
        bash_file.write('--model_preset=monomer \\\n')
        bash_file.write(f'--bfd_database_path={BFD_DB_PATH} \\\n')
        bash_file.write(f'--uniclust30_database_path={UNICLUST30_DB_PATH} \\\n')
        bash_file.write(f'--pdb70_database_path={PDB70_DB_PATH} \\\n')
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

    a_air = arcimboldo_air.ArcimboldoAir(input_load)

    for template in a_air.templates:
        logging.info(f'Looking for alignment between {template.pdb_id} and query sequence using hhsearch:')

        template.generate_features(a_air=a_air)

        if template.sequence_from_template:
            a_air.features.append_row_in_msa(sequence=template.sequence_from_template, msa_uniprot_accession_identifiers=template.pdb_name)
            print(f'Sequence from template \"{template.pdb_id}\" was already added to features.')
        if template.template_features:
            a_air.features.append_new_template_features(new_template_features=template.template_features, custom_sum_prob=template.sum_prob)
            print(f'Template \"{template.pdb_id}\" was already added to features.')

    a_air.features.write_pkl(output_dir=f'{a_air.output_dir}/features.pkl')
    
    sys.exit(1)
    
    if RUN_AF2:
        af2_sh(output_dir=OUTPUT_DIR)
        print('Running AF2')
        af2_output = subprocess.Popen(['bash', f'{OUTPUT_DIR}/run_af2.sh'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout, stderr = af2_output.communicate()
        with open(f'{OUTPUT_DIR}/af2_output.log', 'w') as f:
            f.write(stdout.decode('utf-8'))
        os.system(f'rm -r {OUTPUT_DIR}/msas')
        af_summary(af_job_path=f'{OUTPUT_DIR}')

    features_file = open(f'{OUTPUT_DIR}/features.pkl','rb')
    features = pickle.load(features_file)
    features_file.close()

    print('\n')
    for key in features.keys():
        print(key, features[key].shape)


if __name__ == "__main__":
    main()







