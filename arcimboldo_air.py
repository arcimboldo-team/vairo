#! /usr/bin/env python3
import copy
import os
import shutil

os.environ['XLA_PYTHON_CLIENT_PREALLOCATE'] = 'false'
import sys
import logging
import yaml
from libs import features, structure_air, utils, bioutils


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

        utils.check_input(input_load)
        a_air = structure_air.StructureAir(parameters_dict=input_load)

        utils.create_logger_dir(a_air.log_path)
        os.chdir(a_air.run_dir)
        a_air.write_input_file()
        a_air.generate_output()
        if a_air.custom_features:
            logging.info('Generating features.pkl for AlphaFold2')
            a_air.set_feature(
                feature=features.Features(query_sequence=a_air.sequence_assembled.sequence_mutated_assembled))
            a_air.change_state(state=1)
            a_air.generate_output()
            for template in a_air.templates_list:
                if template.add_to_templates or template.add_to_msa:
                    if not template.aligned:
                        database_dir = os.path.join(a_air.run_dir, template.pdb_id)
                        utils.create_dir(database_dir)
                        template.generate_database(output_dir=database_dir, databases=a_air.alphafold_paths)
                    for sequence in a_air.sequence_assembled.sequence_list:
                        alignment_dir = os.path.join(a_air.run_dir, sequence.name)
                        utils.create_dir(alignment_dir)
                        template.align(output_dir=alignment_dir, sequence_in=sequence, databases=a_air.alphafold_paths)
                    template.generate_features(
                        output_dir=a_air.run_dir,
                        global_reference=a_air.reference,
                        sequence_assembled=a_air.sequence_assembled)
                    a_air.append_line_in_templates(template.results_path_position)
                    if template.add_to_msa:
                        sequence_from_template = template.get_old_sequence(
                            sequence_list=a_air.sequence_assembled.sequence_list_expanded,
                            glycines=a_air.glycines)
                        a_air.feature.append_row_in_msa(sequence_in=sequence_from_template,
                                                        sequence_id=template.pdb_id)
                        logging.info(f'Sequence from template \"{template.pdb_id}\" was added to msa')
                    if template.add_to_templates:
                        a_air.feature.append_new_template_features(new_template_features=template.template_features,
                                                                   custom_sum_prob=template.sum_prob)
                    logging.info(f'Template {template.pdb_id} was added to templates')
            for feat in a_air.features_input:
                feat_aux = features.create_features_from_file(pkl_in_path=feat.path)
                positions = a_air.sequence_assembled.get_range_residues(position_ini=feat.positions[0] - 1,
                                                                        position_end=feat.positions[-1] - 1)
                num_msa = 0
                num_templates = 0
                if feat.keep_msa != 0:
                    num_msa = a_air.feature.set_msa_features(new_msa=feat_aux.msa_features, start=1,
                                                   finish=feat.keep_msa,
                                                   delete_positions=feat.msa_delete,
                                                   positions=positions)
                if feat.keep_templates != 0:
                    num_templates = a_air.feature.set_template_features(new_templates=feat_aux.template_features,
                                                        finish=feat.keep_templates,
                                                        positions=positions,
                                                        sequence_in=feat.sequence)
                feat.add_information(extracted_sequences=num_msa, extracted_templates=num_templates)
            for i, library in enumerate(a_air.libraries):
                num_msa = 0
                num_templates = 0
                if library.add_to_templates:
                    template_path = f'{os.path.join(a_air.input_dir, utils.get_file_name(library.path))}.pdb'
                    if library.positions:
                        template_path = bioutils.copy_positions_of_pdb(path_in=library.path, path_out=template_path, positions=library.positions_list)
                    else:
                        shutil.copy2(library.path, template_path)
                    bioutils.remove_hetatm(template_path, template_path)
                    bioutils.remove_hydrogens(template_path, template_path)
                    template_features = features.extract_template_features_from_aligned_pdb_and_sequence(
                        query_sequence=a_air.sequence_assembled.sequence_assembled,
                        pdb_path=template_path,
                        pdb_id=utils.get_file_name(library.path),
                        chain_id='A')
                    a_air.feature.append_new_template_features(new_template_features=template_features)
                    num_templates = 1

                if library.add_to_msa:
                    extension = utils.get_file_extension(library.path)
                    if extension in ['.pdb']:
                        sequence_list = bioutils.extract_sequence_msa_from_pdb(library.path)
                        sequence_list = list(sequence_list.values())
                    if extension == '.fasta':
                        sequence_list = list(bioutils.extract_sequences(library.path).values())

                    for sequence in sequence_list:
                        seq_aux = sequence
                        if library.positions_list:
                            aux_library_list = copy.copy(library.positions_list)
                            for m, pos in enumerate(aux_library_list):
                                if pos != '-':
                                    aux_library_list[m] = seq_aux[pos]
                            seq_aux = aux_library_list
                        seq_aux = ''.join(seq_aux)
                        a_air.feature.append_row_in_msa(seq_aux, f'lib_{i}_{utils.get_file_name(library.path)}', 1)
                        num_msa += 1
                library.add_information(extracted_sequences=num_msa, extracted_templates=num_templates)

            features_list = a_air.partition_mosaic()
        else:
            features_list = [None] * a_air.mosaic
            a_air.partition_mosaic()
            logging.info('No features.pkl added, default AlphaFold2 run')

        a_air.change_state(state=2)
        a_air.generate_output()

        logging.info('Start running AlphaFold2')
        a_air.run_alphafold(features_list=features_list)
        if len(features_list) > 1:
            a_air.merge_results()
        features_path = os.path.join(a_air.results_dir, 'features.pkl')
        if a_air.feature is None and os.path.exists(features_path):
            new_features = features.create_features_from_file(features_path)
            a_air.set_feature(new_features)
        a_air.output.analyse_output(results_dir=a_air.results_dir,
                                    sequence_assembled=a_air.sequence_assembled,
                                    feature=a_air.feature,
                                    experimental_pdbs=a_air.experimental_pdbs,
                                    binaries_paths=a_air.binaries_paths,
                                    cluster_templates=a_air.cluster_templates)
        if a_air.cluster_templates:
            if a_air.run_af2:
                a_air.templates_clustering()
            a_air.output.analyse_output(results_dir=a_air.results_dir,
                                        sequence_assembled=a_air.sequence_assembled,
                                        feature=a_air.feature,
                                        experimental_pdbs=a_air.experimental_pdbs,
                                        binaries_paths=a_air.binaries_paths)
        a_air.change_state(state=3)
        a_air.generate_output()
        logging.info('ARCIMBOLDO_AIR has finished successfully')

    except SystemExit as e:
        sys.exit(e)
    except Exception as e:
        try:
            a_air.change_state(-1)
            a_air.generate_output()
        except Exception as e2:
            pass
        logging.error('ERROR:', exc_info=True)


if __name__ == "__main__":
    main()
