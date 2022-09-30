#! /usr/bin/env python3

from libs import analyse, bioutils, features, structure_air, utils
import os
import sys
import logging
import yaml
import shutil

def main():

    utils.create_logger()
    logging.info('')
    logging.info('ARCIBOLDO_AIR')
    logging.info('--------------')
    logging.info('')

    try:
        input_path = os.path.abspath(sys.argv[1])
        if not os.path.isfile(input_path):
            raise Exception
    except:
        logging.info('USAGE')
        logging.info('------')
        logging.info(open(utils.get_readme()).read())
        raise SystemExit

    logging.info('Starting ARCIMBOLDO_AIR...')
    if not os.path.exists(input_path):
        raise Exception('The given path for the configuration file either does not exist or you do not have the permissions to read it')
    logging.info(f'Reading the configuration file for ARCIMBOLDO_AIR at {input_path}')

    try:
        with open(input_path) as f:
            input_load = yaml.load(f, Loader=yaml.SafeLoader)
    except:
        raise Exception('It has not been possible to read the input file')

    a_air = structure_air.StructureAir(parameters_dict=input_load)
    
    utils.create_logger_dir(a_air.log_path)
    os.chdir(a_air.run_dir)
    shutil.copy2(input_path, a_air.input_dir)

    if a_air.use_features:
        logging.info('Generating features.pkl for AlphaFold2')
        for template in a_air.templates_list:
            template.generate_features(a_air=a_air)
            if template.add_to_msa:
                sequence_from_template = template.template_features_dict['template_sequence'][0].decode('utf-8')
                a_air.features.append_row_in_msa(sequence=sequence_from_template, msa_uniprot_accession_identifiers=template.pdb_id)
                logging.info(f'Sequence from template \"{template.pdb_id}\" was added to msa')
            if template.add_to_templates:
                a_air.features.append_new_template_features(new_template_features=template.template_features_dict, custom_sum_prob=template.sum_prob)
                logging.info(f'Template \"{template.pdb_id}\" was added to templates')
        a_air.features.write_pkl(output_dir=f'{a_air.run_dir}/features.pkl')
    else:
        logging.info('No features.pkl added, default AlphaFold2 run')
    
    if a_air.run_af2:
        logging.info('Running AlphaFold2')
        a_air.run_alphafold()
        logging.info('AlphaFold2 has finshed succesfully. Proceeding to analyse the results')
        a_air.check_if_assembly()
        analyse.analyse_output(a_air=a_air)

    if not a_air.verbose:
        utils.clean_files(dir=a_air.run_dir)
    else:
        features.print_features(f'{a_air.run_dir}/features.pkl')

    logging.info('ARCIMBOLDO_AIR has finished succesfully')

if __name__ == "__main__":
    try:
        main()
    except SystemExit as e:
        sys.exit(e)
    except:
        logging.error('ERROR:', exc_info=True)








