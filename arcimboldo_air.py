#! /usr/bin/env python3

import os
#os.environ['XLA_PYTHON_CLIENT_PREALLOCATE'] = 'false'
import sys
import logging
import yaml
from libs import features, structure_air, utils

def main():
    try:
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
        a_air.write_input_file()
        a_air.generate_output()

        features_list = []
        if a_air.custom_features:
            logging.info('Generating features.pkl for AlphaFold2')
            a_air.set_feature(feature = features.Features(query_sequence=a_air.sequence_assembled.sequence_assembled))

            if a_air.features_input:
                feat_aux = features.create_features_from_file(pkl_in_path=a_air.features_input.path)
                if a_air.features_input.keep_msa:
                    a_air.feature.set_msa_features(new_msa=feat_aux.msa_features, starting=0)
                if a_air.features_input.keep_templates:
                    a_air.feature.set_template_features(new_templates=feat_aux.template_features)
                        
            a_air.change_state(state=1)
            a_air.generate_output()
            for template in a_air.templates_list:
                if not template.aligned:
                    database_dir = os.path.join(a_air.run_dir, template.pdb_id)
                    utils.create_dir(database_dir)
                    template.generate_database(output_dir=database_dir, database_path=a_air.alphafold_paths.bfd_db_path)
                for sequence in a_air.sequence_assembled.sequence_list:
                    alignment_dir = os.path.join(a_air.run_dir, sequence.name)
                    utils.create_dir(alignment_dir)
                    template.align(output_dir=alignment_dir, sequence=sequence)
                template.generate_features(
                    output_dir=a_air.run_dir,
                    global_reference=a_air.reference,
                    sequence_assembled=a_air.sequence_assembled)
                a_air.append_line_in_templates(template.results_path_position)
                if template.add_to_msa:
                    sequence_from_template = template.template_features['template_sequence'][0].decode('utf-8')
                    a_air.feature.append_row_in_msa(sequence=sequence_from_template,
                                            sequence_id=template.pdb_id)
                    logging.info(f'Sequence from template \"{template.pdb_id}\" was added to msa')
                if template.add_to_templates:
                    a_air.feature.append_new_template_features(new_template_features=template.template_features,
                                                        custom_sum_prob=template.sum_prob)
                    logging.info(f'Template {template.pdb_id} was added to templates')
            a_air.feature.write_pkl(os.path.join(a_air.run_dir, 'features.pkl'))
            features_list = a_air.feature.slicing_features(mosaic=a_air.mosaic, overlap=a_air.mosaic_overlap)
        
        else:
            features_list = [None] * a_air.mosaic
            logging.info('No features.pkl added, default AlphaFold2 run')

        a_air.change_state(state=2)
        a_air.generate_output()
        if a_air.run_af2:
            logging.info('Start running AlphaFold2')
            a_air.run_alphafold(features_list=features_list)
            a_air.merge_results()
            features_path = os.path.join(a_air.run_dir, 'features.pkl')
            if a_air.feature is None and os.path.exists(features_path):
                new_features = features.create_features_from_file(features_path)
                a_air.set_feature(new_features)
            os.chdir(a_air.run_dir)
            a_air.output.set_run_dir(run_dir=a_air.run_dir)
            a_air.output.analyse_output(sequence_assembled=a_air.sequence_assembled, feature=a_air.feature, experimental_pdb=a_air.experimental_pdb, custom_features=a_air.custom_features)
            if a_air.cluster_templates:
                a_air.dendogram_clustering()

        a_air.change_state(state=3)
        a_air.generate_output()
        logging.info('ARCIMBOLDO_AIR has finished successfully')
        
    except SystemExit as e:
        sys.exit(e)
    except Exception as e:
        a_air.change_state(-1)
        a_air.generate_output()
        logging.error('ERROR:', exc_info=True)

if __name__ == "__main__":
    main()

