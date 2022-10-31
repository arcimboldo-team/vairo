#! /usr/bin/env python3

import os
import sys
import logging
import yaml
import shutil
from libs import analyse, features, structure_air, utils


def main():
    utils.create_logger()
    logging.info('')
    logging.info('ARCIMBOLDO_AIR')
    logging.info('--------------')
    logging.info('')

    try:
        input_path = os.path.abspath(sys.argv[1])
        if not os.path.isfile(input_path):
            raise Exception
    except Exception as e:
        logging.info('USAGE')
        logging.info('------')
        logging.info(open(utils.get_readme()).read())
        raise SystemExit

    logging.info('Starting ARCIMBOLDO_AIR...')
    if not os.path.exists(input_path):
        raise Exception(
            'The given path for the configuration file either does not exist or you do not have the permissions to '
            'read it')
    logging.info(f'Reading the configuration file for ARCIMBOLDO_AIR at {input_path}')

    try:
        with open(input_path) as f:
            input_load = yaml.load(f, Loader=yaml.SafeLoader)
    except Exception as e:
        raise Exception('It has not been possible to read the input file')

    a_air = structure_air.StructureAir(parameters_dict=input_load)

    utils.create_logger_dir(a_air.log_path)
    os.chdir(a_air.run_dir)
    shutil.copy2(input_path, a_air.input_dir)

    features_list = []
    if a_air.custom_features:
        logging.info('Generating features.pkl for AlphaFold2')
        feature = features.Features(query_sequence=a_air.sequence_assembled.sequence_assembled)
        for template in a_air.templates_list:
            alignment_dict = {}
            database_dir = os.path.join(a_air.run_dir, template.pdb_id)
            utils.create_dir(database_dir)
            template.generate_database(output_dir=database_dir, database_path=a_air.alphafold_paths.bfd_db_path)
            for sequence in a_air.sequence_assembled.sequence_list:
                alignment_dir = os.path.join(a_air.run_dir, sequence.name)
                utils.create_dir(alignment_dir)
                alignment_dict[sequence.name] = template.align(output_dir=alignment_dir,
                                                               fasta_path=sequence.fasta_path)
            results_path_position = template.generate_features(
                output_dir=a_air.run_dir,
                alignment_dict=alignment_dict,
                global_reference=a_air.reference,
                sequence_assembled=a_air.sequence_assembled)
            a_air.append_line_in_templates(results_path_position)
            if template.add_to_msa:
                sequence_from_template = template.template_features['template_sequence'][0].decode('utf-8')
                feature.append_row_in_msa(sequence=sequence_from_template,
                                          msa_uniprot_accession_identifiers=template.pdb_id)
                logging.info(f'Sequence from template \"{template.pdb_id}\" was added to msa')
            if template.add_to_templates:
                feature.append_new_template_features(new_template_features=template.template_features,
                                                     custom_sum_prob=template.sum_prob)
                logging.info(f'Template \"{template.pdb_id}\" was added to templates')
        feature.write_pkl(os.path.join(a_air.run_dir, 'features.pkl'))
        if a_air.mosaic is not None:
            features_list = feature.slicing_features(mosaic=a_air.mosaic)
        else:
            features_list.append(feature)
    else:
        features_list.append(None)
        logging.info('No features.pkl added, default AlphaFold2 run')

    if a_air.run_af2:
        logging.info('Start running AlphaFold2')
        a_air.run_alphafold(features_list=features_list)
        a_air.merge_results()
        os.chdir(a_air.run_dir)
        analyse.analyse_output(a_air=a_air)

    if not a_air.verbose:
        utils.clean_files(input_dir=a_air.run_dir)

    logging.info('ARCIMBOLDO_AIR has finished succesfully')


if __name__ == "__main__":
    try:
        main()
    except SystemExit as e:
        sys.exit(e)
    except Exception as e:
        logging.error('ERROR:', exc_info=True)
