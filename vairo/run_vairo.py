#! /usr/bin/env python3

import os
os.environ['XLA_PYTHON_CLIENT_PREALLOCATE'] = 'false'
import sys
import logging
import yaml
import argparse

from datetime import datetime
from libs import features, main_structure, utils, bioutils, pymol_script


def main():
    try:
        utils.create_logger()
        logging.error('')
        logging.error('VAIRO')
        logging.error('--------------')
        logging.error('')
        parser = argparse.ArgumentParser(description='Process configuration file and optional check flag.')
        parser.add_argument('input_path', nargs='?', help='Path to the input file (.yml)')
        parser.add_argument('-check', action='store_true', help='Check the input file and exit')
        args = parser.parse_args()
        # If config_file is not provided, print the README and exit
        if not args.input_path:
            logging.error('USAGE')
            logging.error('------')
            logging.error(open(utils.get_readme()).read())
            raise SystemExit

        # Retrieve the values
        input_path = args.input_path
        check_mode = args.check

        logging.error('Starting VAIRO...')
        logging.error(f'Timestamp: {datetime.now()}')
        utils.check_external_programs()
        if not os.path.exists(input_path):
            raise Exception(
                'The given path for the configuration file either does not exist or you do not have the permissions to '
                'read it')
        logging.error(f'Reading the configuration file for VAIRO at {input_path}')

        try:
            with open(input_path) as f:
                input_load = yaml.load(f, Loader=yaml.SafeLoader)
        except Exception as e:
            raise Exception('It has not been possible to read the input file')
        if check_mode:
            logging.error(f'Input file is CORRECT')
            raise SystemExit

        utils.check_input(input_load)
        a_air = main_structure.MainStructure(parameters_dict=input_load)
        os.chdir(a_air.run_dir)
        a_air.write_input_file()
        if a_air.custom_features:
            logging.error('Guided mode selected: Generating custom features.pkl for AlphaFold2')
            a_air.set_feature(
                feature=features.Features(query_sequence=a_air.sequence_assembled.sequence_mutated_assembled))
            a_air.change_state(state=1)
            a_air.generate_output()
            for template in a_air.templates_list:
                logging.error(f'Reading template {template.pdb_id}')
                if template.add_to_templates or template.add_to_msa:
                    template.generate_chains(sequence_assembled=a_air.sequence_assembled)
                    if not template.aligned:
                        template.align(databases=a_air.alphafold_paths)
                    template.generate_features(
                        global_reference=a_air.reference,
                        sequence_assembled=a_air.sequence_assembled)

                    a_air.append_line_in_templates(template.results_path_position)
                    if template.add_to_msa:
                        sequence_from_template = template.get_old_sequence(
                            sequence_list=a_air.sequence_assembled.sequence_list_expanded,
                            glycines=a_air.glycines)
                        a_air.feature.append_row_in_msa(sequence_in=sequence_from_template,
                                                        sequence_id=template.pdb_id)
                        logging.error(f'     Adding the template sequence to the MSA')
                    if template.add_to_templates:
                        a_air.feature.append_new_template_features(new_template_features=template.template_features,
                                                                   custom_sum_prob=template.sum_prob)
                        logging.error(f'     Adding template to templates')
            for feat in a_air.features_input:
                logging.error(f'Reading features {feat.path}')
                feat_aux = features.create_features_from_file(pkl_in_path=feat.path)
                num_msa = 0
                num_templates = 0
                modifications_list = utils.modification_list(query=feat.numbering_query, target=feat.numbering_features,
                                                             length=a_air.sequence_assembled.length)
                # Delete the residues before expanding, so we avoid shifting them
                feat_aux.delete_residues_msa(delete_positions=feat.msa_mask)
                feat_aux.replace_sequence_template(sequence_in=feat.replace_sequence)
                feat_aux.mutate_residues(mutation_dict=feat.mutate_residues)
                # Cut and expand the features, in order to fit the general features.pkl
                feat_aux = feat_aux.cut_expand_features(query_sequence=a_air.sequence_assembled.sequence_assembled,
                                                        modifications_list=modifications_list)
                if feat.keep_msa != 0:
                    # Send without masking features, as we have deleted them
                    num_msa = a_air.feature.set_msa_features(new_msa=feat_aux.msa_features, start=1,
                                                             finish=feat.keep_msa,
                                                             delete_positions=[])
                    logging.error(f'     Adding {num_msa} sequence/s to the MSA')
                if feat.keep_templates != 0:
                    num_templates = a_air.feature.set_template_features(new_templates=feat_aux.template_features,
                                                                        finish=feat.keep_templates)
                    logging.error(f'     Adding {num_templates} template/s to templates')
                feat.add_information(num_msa=num_msa, num_templates=num_templates)

            for i, library in enumerate(a_air.library_list):
                logging.error(f'Reading library {library.path}')
                aux_list = [os.path.join(library.path, file) for file in os.listdir(library.path)] if os.path.isdir(
                    library.path) else [library.path]
                paths = [path for path in aux_list if utils.get_file_extension(path) in ['.pdb', '.fasta']]
                modifications_list = utils.modification_list(query=library.numbering_query,
                                                             target=library.numbering_library,
                                                             length=a_air.sequence_assembled.length)
                lib_feat = features.Features(query_sequence=a_air.sequence_assembled.sequence_mutated_assembled)
                for aux_path in paths:
                    if library.add_to_templates:
                        if utils.get_file_extension(aux_path) == '.fasta' and library.add_to_templates:
                            logging.error(f'Ignoring add_to_templates to True for fasta file {aux_path}')
                        else:
                            template_path = f'{os.path.join(a_air.input_dir, utils.get_file_name(aux_path))}.pdb'
                            bioutils.remove_hetatm(aux_path, template_path)
                            bioutils.remove_hydrogens(template_path, template_path)
                            template_features = features.extract_template_features_from_aligned_pdb_and_sequence(
                                query_sequence=a_air.sequence_assembled.sequence_assembled,
                                pdb_path=template_path,
                                pdb_id=utils.get_file_name(aux_path),
                                chain_id='A')
                            lib_feat.append_new_template_features(new_template_features=template_features)

                    if library.add_to_msa:
                        extension = utils.get_file_extension(aux_path)
                        sequence_list =  []
                        if extension in ['.pdb']:
                            sequence_list = bioutils.extract_sequence_msa_from_pdb(aux_path)
                            sequence_list = list(sequence_list.values())
                        if extension == '.fasta':
                            sequence_list = list(bioutils.extract_sequences(aux_path).values())
                        for num, sequence in enumerate(sequence_list):
                            lib_feat.append_row_in_msa(sequence, f'lib_{i}-{num}_{utils.get_file_name(aux_path)}', 1)

                num_msa = 0
                num_templates = 0
                lib_feat = lib_feat.cut_expand_features(query_sequence=a_air.sequence_assembled.sequence_assembled,
                                                        modifications_list=modifications_list)
                if library.add_to_msa:
                    num_msa = a_air.feature.set_msa_features(new_msa=lib_feat.msa_features)
                if library.add_to_templates:
                    num_templates = a_air.feature.set_template_features(new_templates=lib_feat.template_features)

                if num_templates > 0:
                    logging.error(f'     Adding {num_templates} template/s to templates')
                if num_msa > 0:
                    logging.error(f'     Adding {num_msa} sequence/s to the MSA')
                library.add_information(num_msa=num_msa, num_templates=num_templates)
        
            features_list = a_air.partition_mosaic()
        else:
            a_air.generate_output()
            features_list = [None] * a_air.mosaic
            a_air.partition_mosaic()
            logging.error('Naive mode selected: No custom features.pkl generated')

        a_air.change_state(state=2)
        a_air.generate_output()
        logging.error('All input information has been processed correctly')
        a_air.run_alphafold(features_list=features_list)        
        if len(features_list) > 1:
            a_air.merge_results()

        features_path = os.path.join(a_air.results_dir, 'features.pkl')
        if a_air.feature is None and os.path.exists(features_path):
            new_features = features.create_features_from_file(features_path)
            a_air.set_feature(new_features)

        # a_air.align_experimental_pdbs()
        if a_air.mode == 'naive' and a_air.run_af2 and a_air.cluster_templates:
            # store results features before trimming
            old_features_path = os.path.join(a_air.results_dir, 'alphafold_features.pkl')
            a_air.feature.write_pkl(pkl_path=old_features_path)
            a_air.feature.select_msa_templates(sequence_assembled=a_air.sequence_predicted_assembled)
            a_air.extract_results(region_predicted=True)
            a_air.templates_clustering()
        else:
            a_air.extract_results(region_predicted=True)
            if a_air.sequence_predicted_assembled.mutated:
                a_air.delete_mutations()

        a_air.expand_features_predicted_sequence()
        a_air.extract_results(region_predicted=False)
        a_air.analyse_output()
        a_air.change_state(state=3)
        pymol_script.create_pymol_session(a_air)
        logging.error(f'Timestamp: {datetime.now()}')
        logging.error('VAIRO has finished successfully')
        if a_air.feature is not None:
            a_air.feature.set_extra_info()
            a_air.feature.write_pkl(pkl_path=features_path)
        a_air.generate_output()

    except SystemExit as e:
        sys.exit()
    except Exception as e:
        logging.error(f'Timestamp: {datetime.now()}')
        logging.error('ERROR:', exc_info=True)
        try:
            a_air.change_state(-1)
            a_air.generate_output()
        except Exception as e2:
            pass


if __name__ == "__main__":
    main()
