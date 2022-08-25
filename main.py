#! /usr/bin/env python3

import shutil
from libs import analyse, arcimboldo_air, bioutils, features, utils
import os
import sys
import toml
import logging

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
    os.chdir(a_air.run_dir)
    shutil.copy2(input_path, a_air.input_dir)

    for template in a_air.templates:
        template.generate_features(a_air=a_air)
        if template.add_to_msa:
            sequence_from_template = template.template_features['template_sequence'][0].decode('utf-8')
            a_air.features.append_row_in_msa(sequence=sequence_from_template, msa_uniprot_accession_identifiers=template.pdb_id)
            logging.info(f'Sequence from template \"{template.pdb_id}\" was added to msa.')
        if template.add_to_templates:
            a_air.features.append_new_template_features(new_template_features=template.template_features, custom_sum_prob=template.sum_prob)
            logging.info(f'Template \"{template.pdb_id}\" was added to templates.')

    a_air.features.write_pkl(output_dir=f'{a_air.run_dir}/features.pkl')

    if a_air.run_af2:
        #bioutils.run_af2(output_dir=a_air.run_dir, alphafold_paths=a_air.alphafold_paths)
        a_air.check_if_assembly()
        analyse.analyse_output(output_dir=a_air.output_dir, run_dir=a_air.run_dir, a_air=a_air)
    
    if not a_air.verbose:
        utils.clean_files(dir=a_air.run_dir)
    else:
        features.print_features(f'{a_air.run_dir}/features.pkl')

if __name__ == "__main__":
    main()







