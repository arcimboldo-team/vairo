import features_pipeline
import analyse_pipeline
from analyse_pipeline import af_summary
from energy_analysis import run_pisa
from absl import app
from absl import flags
from absl import logging
import os
import subprocess


FLAGS = flags.FLAGS

flags.DEFINE_enum('feature_mode', None, ['from_pkl', 'from_chain_in_pdb', 'from_assembly_in_pdb', 'from_custom_pdb'],
                  'Could be: read features pkl file (from_pkl), generate features from specific chain in deposited PDB '
                  '(from_pdb_and_chain), generate features from assembly in deposited PDB (from_assembly_in_pdb)')
flags.DEFINE_string('output_dir', None, 'Path to the output directory.')
flags.DEFINE_string('fasta_path', None, 'Fasta path.')
flags.DEFINE_string('features_pkl_path', None, 'Path to features.pkl.')
flags.DEFINE_string('template_pdb_id', None, 'Template PDB id.')
flags.DEFINE_string('template_pdb_chain_id', None, 'Template PDB chain id.')
flags.DEFINE_string('hhr_path', None, 'Path to HHR file provided by the user. HHR file has to contain the target'
                                      ' template in its headers.')
flags.DEFINE_string('custom_pdb_path', None, 'Path to custom PDB. Coordinates must be aligned with the query sequence!')
flags.DEFINE_string('poly_ala_list_of_res_ranges', None, 'Convert template as polyalanine model. Can be "all" o defined'
                                                         ' as 1-40,45-50,400-500')
# ALPHAFOLD DBS FLAGS:
flags.DEFINE_string('mmcif_db', None, 'Path to mmcif db from AF databases.')
flags.DEFINE_string('pdb70_db', None, 'Path to PDB70 db from AF.')



def run_af2(features_pkl_directory):

    out_dir = '/'.join(features_pkl_directory.split('/')[:-1])
    name = features_pkl_directory.split('/')[-1]

    fasta_name, template_pdb_id = name.split('-')[0], name.split('-')[1][:4]
    if FLAGS.feature_mode == 'from_chain_in_pdb':
        template_chain_id = name.split('-')[1][-1]
        bash_file = open(f'{out_dir}/{fasta_name}-{template_pdb_id}{template_chain_id}/{fasta_name}-{template_pdb_id}{template_chain_id}.sh','w')
    if FLAGS.feature_mode == 'from_assembly_in_pdb':
        bash_file = open(f'{out_dir}/{fasta_name}-{template_pdb_id}/{fasta_name}-{template_pdb_id}.sh','w')
    if FLAGS.feature_mode == 'from_custom_pdb':
        bash_file = open(f'{out_dir}/{fasta_name}-custom/{fasta_name}-custom.sh', 'w')
    else:
        logging.fatal('You must set a feature mode (--feature_mode)')


    bash_file.write('#!/bin/bash\n')
    bash_file.write('python ./alphafold/run_alphafold.py \\\n')
    bash_file.write('--use_precomputed_msas=False \\\n')
    bash_file.write(f'--fasta_paths=fake.fasta \\\n')
    bash_file.write(f'--output_dir={out_dir} \\\n')
    bash_file.write('--data_dir=/alpha/alphauser/af2_dbs \\\n')
    bash_file.write('--uniref90_database_path=/alpha/alphauser/af2_dbs/uniref90/uniref90.fasta \\\n')
    bash_file.write('--mgnify_database_path=/alpha/alphauser/af2_dbs/mgnify/mgy_clusters_2018_12.fa \\\n')
    bash_file.write(f'--template_mmcif_dir={FLAGS.mmcif_db} \\\n')
    bash_file.write('--max_template_date=2021-10-22 \\\n')
    bash_file.write(f'--obsolete_pdbs_path=/alpha/alphauser/af2_dbs/pdb_mmcif/obsolete.dat \\\n')
    bash_file.write('--bfd_database_path=/alpha/alphauser/af2_dbs/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \\\n')
    bash_file.write('--uniclust30_database_path=/alpha/alphauser/af2_dbs/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \\\n')
    bash_file.write('--model_preset=monomer \\\n')
    bash_file.write('--pdb70_database_path=/alpha/alphauser/af2_dbs/pdb70/pdb70 \\\n')
    bash_file.write('--read_features_pkl=True\n')
    bash_file.close()


def main(argv):

    f = features_pipeline.Features()

    if FLAGS.feature_mode == 'from_pkl':
        if FLAGS.features_pkl_path:
            f.features_from_features_pkl(features_pkl_path=FLAGS.features_pkl_path)
        else:
            logging.fatal('Please, provide path to features.pkl!')

    if FLAGS.feature_mode == 'from_chain_in_pdb':

        logging.info(f'MODE: features from {FLAGS.template_pdb_id}{FLAGS.template_pdb_chain_id} will be used as template.')
        fasta_name, query_sequence = features_pipeline.extract_query_sequence(FLAGS.fasta_path)
        isExist = os.path.exists(f'{FLAGS.output_dir}/{fasta_name}-{FLAGS.template_pdb_id}{FLAGS.template_pdb_chain_id}')
        if not isExist:
            os.mkdir(f'{FLAGS.output_dir}/{fasta_name}-{FLAGS.template_pdb_id}{FLAGS.template_pdb_chain_id}')

        f.features_from_pdb_id(query_sequence=query_sequence,
                               pdb_id=FLAGS.template_pdb_id,
                               chain_id=FLAGS.template_pdb_chain_id,
                               out_dir=FLAGS.output_dir,
                               mmcif_db=FLAGS.mmcif_db,
                               pdb70_db=FLAGS.pdb70_db,
                               hhr_path=FLAGS.hhr_path)

        features_pipeline.write_pkl_from_features(features=f.features,
                                                  out_path=f'{FLAGS.output_dir}/{fasta_name}-{FLAGS.template_pdb_id}'
                                                  f'{FLAGS.template_pdb_chain_id}/features.pkl')

        run_af2(features_pkl_directory=f'{FLAGS.output_dir}/{fasta_name}-{FLAGS.template_pdb_id}{FLAGS.template_pdb_chain_id}')

    if FLAGS.feature_mode == 'from_assembly_in_pdb':

        logging.info(f'MODE: features from the {FLAGS.template_pdb_id} assembly will be used as template.')
        fasta_name, query_sequence = features_pipeline.extract_query_sequence(FLAGS.fasta_path)
        isExist = os.path.exists(
            f'{FLAGS.output_dir}/{fasta_name}-{FLAGS.template_pdb_id}')
        if not isExist:
            os.mkdir(f'{FLAGS.output_dir}/{fasta_name}-{FLAGS.template_pdb_id}')

        f.features_for_query_sequence_and_experimental_assembly_in_pdb(query_subunit_sequence=query_sequence,
                                                                       pdb_id=FLAGS.template_pdb_id,
                                                                       mmcif_db=FLAGS.mmcif_db,
                                                                       pdb70_db=FLAGS.pdb70_db,
                                                                       out_dir=f'{FLAGS.output_dir}/{fasta_name}-{FLAGS.template_pdb_id}',
                                                                       hhr_path=FLAGS.hhr_path)

        features_pipeline.write_pkl_from_features(features=f.features,
                                                  out_path=f'{FLAGS.output_dir}/{fasta_name}-{FLAGS.template_pdb_id}/features.pkl')

        run_af2(features_pkl_directory=f'{FLAGS.output_dir}/{fasta_name}-{FLAGS.template_pdb_id}')

    if FLAGS.feature_mode == 'from_custom_pdb':

        logging.info(f'MODE: features from custom pdb will be used as template.')
        fasta_name, query_sequence = features_pipeline.extract_query_sequence(FLAGS.fasta_path)
        isExist = os.path.exists(f'{FLAGS.output_dir}/{fasta_name}-custom')
        if not isExist:
            os.mkdir(f'{FLAGS.output_dir}/{fasta_name}-custom')

        if FLAGS.poly_ala_list_of_res_ranges:

            if FLAGS.poly_ala_list_of_res_ranges == 'all':
                poly_ala_list_of_res_ranges = [[-10000, 10000]]
            else:
                poly_ala_list_of_res_ranges = [item.split('-') for item in FLAGS.poly_ala_list_of_res_ranges.split(',')]

            features_pipeline.convert_template_to_polyala(pdb_in_path=FLAGS.custom_pdb_path,
                                                          pdb_out_path=f'{FLAGS.output_dir}/{fasta_name}-custom/{FLAGS.custom_pdb_path.split("/")[-1][:-4]}_polyala.pdb',
                                                          list_of_res_ranges=poly_ala_list_of_res_ranges)
            f.features_for_custom_pdb(query_sequence=query_sequence,
                                      pdb_path=f'{FLAGS.output_dir}/{fasta_name}-custom/{FLAGS.custom_pdb_path.split("/")[-1][:-4]}_polyala.pdb')

            if FLAGS.poly_ala_list_of_res_ranges == 'all':
                f.delete_rows_in_msa(list_to_remove=[-1])
            else:
                list_to_remove = []
                if len(poly_ala_list_of_res_ranges) > 1:
                    for item in poly_ala_list_of_res_ranges:
                        list_to_remove.extend(list(range(int(item[0])-1, int(item[1])-1)))

                else:
                    list_to_remove = list(range(int(poly_ala_list_of_res_ranges[0][0])-1, int(poly_ala_list_of_res_ranges[0][1])-1))

                f.replace_columns_by_gaps_in_msa(list_to_remove=list_to_remove, specific_row=-1)

        else:

            f.features_for_custom_pdb(query_sequence=query_sequence,
                                      pdb_path=f'{FLAGS.custom_pdb_path}')

        features_pipeline.write_pkl_from_features(features=f.features,
                                                  out_path=f'{FLAGS.output_dir}/{fasta_name}-custom/features.pkl')

        run_af2(features_pkl_directory=f'{FLAGS.output_dir}/{fasta_name}-custom')


###################################

    out_dir = '/cri4/albert/Desktop/1ixc-3fxq'
    query_fasta_path = '/cri4/albert/Desktop/rcsb_pdb_1IXC.fasta'
    template_list = ['3fxq']
    msa_list = ['3fxq']

    isExist = os.path.exists(out_dir)
    if not isExist:
        os.mkdir(out_dir)

    non_redundant_msa_list = []
    for pdb in msa_list:
        if pdb not in template_list:
            non_redundant_msa_list.append(pdb)

    query_subunit_sequence = features_pipeline.extract_query_sequence(query_fasta_path)[1]

    if len(template_list) != 0:
        f.features_for_query_sequence_and_experimental_assembly_in_pdb(query_subunit_sequence=query_subunit_sequence,
                                                                       pdb_id=template_list[0],
                                                                       mmcif_db='/alpha/alphauser/af2_dbs/pdb_mmcif/mmcif_files',
                                                                       pdb70_db='/alpha/alphauser/af2_dbs/pdb70/pdb70',
                                                                       out_dir=out_dir,
                                                                       hhr_path=False)
        if template_list[0] not in msa_list:
            f.delete_rows_in_msa(list_to_remove=[-1])
        for pdb_id in template_list[1:]:
            g = features_pipeline.Features()
            g.features_for_query_sequence_and_experimental_assembly_in_pdb(query_subunit_sequence=query_subunit_sequence,
                                                                           pdb_id=pdb_id,
                                                                           mmcif_db='/alpha/alphauser/af2_dbs/pdb_mmcif/mmcif_files',
                                                                           pdb70_db='/alpha/alphauser/af2_dbs/pdb70/pdb70',
                                                                           out_dir=out_dir,
                                                                           hhr_path=f'{out_dir}/output.hhr')
            f.features = features_pipeline.append_template_features_from_pdb(features=f.features,
                                                                             features_for_append=g.features)
            if pdb_id in msa_list:
                sequence_to_append = g.features['template_sequence'][0].decode('utf-8')
                f.append_row_in_msa(sequence=sequence_to_append)
    else:
        logging.info('Template list is empty.')
        f.features_for_query_sequence_and_experimental_assembly_in_pdb(query_subunit_sequence=query_subunit_sequence,
                                                                       pdb_id=non_redundant_msa_list[0],
                                                                       mmcif_db='/alpha/alphauser/af2_dbs/pdb_mmcif/mmcif_files',
                                                                       pdb70_db='/alpha/alphauser/af2_dbs/pdb70/pdb70',
                                                                       out_dir=out_dir,
                                                                       hhr_path=f'{out_dir}/output.hhr')
        # TODO make remove_row_in_templates in features_pipeline

    if len(template_list) == 0:
        non_redundant_msa_list = non_redundant_msa_list[1:]

    for pdb_id in non_redundant_msa_list:

        if '6xtu' in pdb_id: # TODO exception for 6xtu
            hhr_path = '/cri4/albert/Desktop/atzr_vs_6xtu.hhr'
        else:
            hhr_path = f'{out_dir}/output.hhr'

        h = features_pipeline.Features()
        h.features_for_query_sequence_and_experimental_assembly_in_pdb(query_subunit_sequence=query_subunit_sequence,
                                                                       pdb_id=pdb_id,
                                                                       mmcif_db='/alpha/alphauser/af2_dbs/pdb_mmcif/mmcif_files',
                                                                       pdb70_db='/alpha/alphauser/af2_dbs/pdb70/pdb70',
                                                                       out_dir=out_dir,
                                                                       hhr_path=hhr_path)
        sequence_to_append = h.features['template_sequence'][0].decode('utf-8')
        f.append_row_in_msa(sequence=sequence_to_append)

    ID_TO_HHBLITS_AA = {0: 'A', 1: 'C', 2: 'D', 3: 'E', 4: 'F', 5: 'G', 6: 'H',
                        7: 'I', 8: 'K', 9: 'L', 10: 'M', 11: 'N', 12: 'P', 13: 'Q',
                        14: 'R', 15: 'S', 16: 'T', 17: 'V', 18: 'W', 19: 'Y',
                        20: 'X', 21: 'o'}

    for key in f.features.keys():
        print(key, f.features[key].shape)
    for seq in f.features['msa']:
        print(''.join([ID_TO_HHBLITS_AA[res] for res in list(seq)]))
    print(f.features['template_domain_names'])
    os.system(f'rm {out_dir}/*_template.pdb')
    f.write_all_templates_in_features(output_path=out_dir)
    features_pipeline.write_pkl_from_features(features=f.features,
                                              out_path=f'{out_dir}/features.pkl')


###################################

if __name__ == '__main__':

    app.run(main)


