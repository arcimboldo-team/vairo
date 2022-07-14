#! /usr/bin/env python3

from libs import utils, analyse, arcimboldo_air, features, hhsearch
import copy
import glob, os
import sys
import toml
import shutil
import logging

def generate_features(output_dir, num_of_seqs, query_sequence, mode, template, add_to_msa, add_to_templates,
                      pdb70_db_path, sum_prob):

    assembly_sequence_with_linkers = ''
    for num in range(num_of_seqs):
        assembly_sequence_with_linkers = assembly_sequence_with_linkers + query_sequence + 50 * 'G'
    assembly_sequence_with_linkers = assembly_sequence_with_linkers[:-50]


    if num_of_seqs == 1:

        if mode == 'monomer_from_pdb':

            pdb_id, chain = template[:4], template[-1]

            print_msg_box(msg=f'Monomer from PDB: {pdb_id.upper()}_{chain}',
                          indent=2, title=None)

            if os.path.exists(f'{output_dir}/{template}.pdb'):
                print(f'{template}.pdb already exists in {output_dir}.')
            else:
                download_pdb(pdb_id=pdb_id, out_dir=output_dir)
                os.system(f'mv {output_dir}/pdb{pdb_id}.ent {output_dir}/{pdb_id}.pdb')

            print(f'Looking for alignment between {pdb_id}_{chain} and query sequence using hhsearch:')
            if os.path.exists(f'{output_dir}/output.hhr'):
                print(f'hhr file already present in {output_dir}')
            else:
                run_hhsearch(query_fasta_path=QUERY_FASTA_PATH, pdb70_db=f'{pdb70_db_path}/pdb70',
                             output_path=f'{output_dir}/output.hhr')

            g = Features(query_sequence=query_sequence)
            try:
                template_features = g.extract_template_features_from_pdb(hhr_path=f'{output_dir}/output.hhr',
                                                                         pdb_id=pdb_id.upper(),
                                                                         chain_id=chain,
                                                                         mmcif_db=MMCIF_DB_PATH)
            except:
                generate_hhsearch_db(template_cif_path=f'{MMCIF_DB_PATH}/{pdb_id.lower()}.cif', output_dir=output_dir)

                run_hhsearch(query_fasta_path=QUERY_FASTA_PATH, pdb70_db=f'{output_dir}/pdb70',
                             output_path=f'{output_dir}/{template}_output.hhr')
                template_features = g.extract_template_features_from_pdb(hhr_path=f'{output_dir}/{template}_output.hhr',
                                                                         pdb_id=pdb_id,
                                                                         chain_id=chain,
                                                                         mmcif_db=MMCIF_DB_PATH)
                os.system(f'rm {output_dir}/*.ffindex')
                os.system(f'rm {output_dir}/*.ffdata')

            g.append_new_template_features(new_template_features=template_features, custom_sum_prob=sum_prob)
            g.write_all_templates_in_features(output_path=output_dir)
            _template_features = copy.deepcopy(template_features)

        if mode == 'monomer_from_custom_template':

            name = template.split("/")[-1][:-4]

            print_msg_box(msg=f'Monomer from custom template {name}', indent=2, title=None)

            pdb2mmcif(pdb_in_path=template, cif_out_path=f'{output_dir}/{name}.cif')
            generate_hhsearch_db(template_cif_path=f'{output_dir}/{name}.cif', output_dir=output_dir)

            run_hhsearch(query_fasta_path=QUERY_FASTA_PATH, pdb70_db=f'{output_dir}/pdb70',
                         output_path=f'{output_dir}/custom_output.hhr')

            g = Features(query_sequence=query_sequence)

            template_features = g.extract_template_features_from_pdb(
                hhr_path=f'{output_dir}/custom_output.hhr',
                pdb_id=name,
                chain_id='A',
                mmcif_db=output_dir)

            g.append_new_template_features(new_template_features=template_features, custom_sum_prob=sum_prob)
            # g.write_all_templates_in_features(output_path=output_dir)
            _template_features = copy.deepcopy(template_features)

            os.system(f'rm {output_dir}/*ffdata')
            os.system(f'rm {output_dir}/*ffindex')
            os.system(f'rm {output_dir}/custom_output.hhr')
            os.system(f'rm {output_dir}/{name}.cif')

        if mode == 'monomer_from_aligned_custom_template':

            print_msg_box(msg=f'Monomer from aligned custom template: {template.split("/")[-1]} ({num_of_seqs} subunits)',
                          indent=2, title=None)

            name = template.split("/")[-1][:-4]
            os.system(f'cp {template} {output_dir}/{name}_template.pdb')

            g = Features(query_sequence=query_sequence)

            template_features = g.extract_template_features_from_aligned_pdb_and_sequence(
                pdb_path=f'{output_dir}/{name}_template.pdb',
                chain_ID='A')

            g.append_new_template_features(new_template_features=template_features, custom_sum_prob=sum_prob)
            # g.write_all_templates_in_features(output_path=output_dir)
            _template_features = copy.deepcopy(template_features)

        if add_to_templates:
            pass
        else:
            template_features = None
        if add_to_msa:
            sequence_from_template = _template_features['template_sequence'][0].decode('utf-8')
        else:
            sequence_from_template = None

        return sequence_from_template, template_features

    else:
        if mode == 'assembly_from_pdb':

            print_msg_box(msg=f'Assembly from PDB: {template.upper()} ({num_of_seqs} subunits)',
                          indent=2, title=None)

            if os.path.exists(f'{output_dir}/{template}.pdb'):
                print(f'{template}.pdb already exists in {output_dir}.')
            else:
                download_pdb(pdb_id=template, out_dir=output_dir)
                os.system(f'mv {output_dir}/pdb{template}.ent {output_dir}/{template}.pdb')



            chain_list, transformations_list = read_remark_350(f'{output_dir}/{template}.pdb', use_pisa=False)
            print('Assembly can be build using chain(s)', *chain_list, 'by applying the following transformations:')
            for matrix in transformations_list:
                print(*matrix)
            new_chain_list = list(string.ascii_uppercase)[:len(transformations_list) * len(chain_list)]
            if len(new_chain_list) != num_of_seqs:
                print(f'Assembly description from REMARK 350 contains {len(new_chain_list)} subunits. Please, try to'
                      f'generate a new REMARK 350 (manually or using e.g. PISA) for considering a new assembly with '
                      f'{num_of_seqs} subunits')
                sys.exit()

            query_seq_length = len(query_sequence)
            g = Features(query_sequence=query_sequence)
            counter = 0
            for chain in chain_list:
                try:
                    template_features = g.extract_template_features_from_pdb(hhr_path=f'{output_dir}/output.hhr',
                                                                             pdb_id=template.upper(),
                                                                             chain_id=chain,
                                                                             mmcif_db=MMCIF_DB_PATH)
                except:
                    generate_hhsearch_db(template_cif_path=f'{MMCIF_DB_PATH}/{template.lower()}.cif', output_dir=output_dir)

                    run_hhsearch(query_fasta_path=QUERY_FASTA_PATH, pdb70_db=f'{output_dir}/pdb70',
                                 output_path=f'{output_dir}/{template}_output.hhr')
                    template_features = g.extract_template_features_from_pdb(hhr_path=f'{output_dir}/{template}_output.hhr',
                                                                             pdb_id=template,
                                                                             chain_id=chain,
                                                                             mmcif_db=MMCIF_DB_PATH)
                    os.system(f'rm {output_dir}/*.ffindex')
                    os.system(f'rm {output_dir}/*.ffdata')


                g.append_new_template_features(new_template_features=template_features, custom_sum_prob=sum_prob)
                g.write_all_templates_in_features(output_path=output_dir)
                for num, transformation in enumerate(transformations_list):
                    counter += 1
                    rot_and_trans(pdb_path=f'{output_dir}/{template}{chain}_template.pdb',
                                  out_pdb_path=f'{output_dir}/{counter}.pdb',
                                  rot_tra_matrix=transformation)
                    if counter == 1:
                        change_chain_and_apply_offset_in_single_chain(pdb_in_path=f'{output_dir}/{counter}.pdb',
                                                                      pdb_out_path=f'{output_dir}/{new_chain_list[counter - 1]}.pdb',
                                                                      offset=None,
                                                                      chain='A')
                    else:
                        offset = query_seq_length * (counter - 1) + 50 * (counter - 1)  # 50 glycines as offset !!!
                        change_chain_and_apply_offset_in_single_chain(pdb_in_path=f'{output_dir}/{counter}.pdb',
                                                                      pdb_out_path=f'{output_dir}/{new_chain_list[counter - 1]}.pdb',
                                                                      offset=offset,
                                                                      chain='A')
                    os.system(f'rm {output_dir}/{counter}.pdb')

            list_of_paths_of_pdbs_to_merge = [f'{output_dir}/{ch}.pdb' for ch in new_chain_list]
            merge_pdbs(list_of_paths_of_pdbs_to_merge=list_of_paths_of_pdbs_to_merge,
                       merged_pdb_path=f'{output_dir}/{template}_template.pdb')
            #os.system(f'rm {output_dir}/[A-Z].pdb')
            #os.system(f'rm {output_dir}/{template}[A-Z]_template.pdb')

            name = template

        if mode == 'assembly_from_custom_template':

            print_msg_box(msg=f'Assembly from custom template: {template.split("/")[-1]} ({num_of_seqs} subunits)',
                          indent=2, title=None)

            name = template.split("/")[-1][:-4]

            pdb2mmcif(pdb_in_path=template, cif_out_path=f'{output_dir}/{name}.cif')
            generate_hhsearch_db(template_cif_path=f'{output_dir}/{name}.cif', output_dir=output_dir)

            run_hhsearch(fasta_path=QUERY_FASTA_PATH, pdb70_db=f'{output_dir}/pdb70',
                         output_path=f'{output_dir}/custom_output.hhr')

            chain_list, transformations_list = read_remark_350(pdb_path=f'{output_dir}/{name}.cif', use_pisa=True)
            new_chain_list = list(string.ascii_uppercase)[:len(transformations_list) * len(chain_list)]
            print('Assembly can be build using chain(s)', *chain_list, 'by applying the following transformations:')
            for matrix in transformations_list:
                print(*matrix)

            query_seq_length = len(query_sequence)
            g = Features(query_sequence=query_sequence)
            counter = 0 
            for chain in chain_list:
                template_features = g.extract_template_features_from_pdb(
                    hhr_path=f'{output_dir}/custom_output.hhr',
                    pdb_id=name,
                    chain_id=chain,
                    mmcif_db=output_dir)
                g.append_new_template_features(new_template_features=template_features, custom_sum_prob=sum_prob)
                g.write_all_templates_in_features(output_path=output_dir)
                for num, transformation in enumerate(transformations_list):
                    counter += 1
                    rot_and_trans(pdb_path=f'{output_dir}/{name}{chain}_template.pdb',
                                  out_pdb_path=f'{output_dir}/{counter}.pdb',
                                  rot_tra_matrix=transformation)

                    if counter == 1:
                        change_chain_and_apply_offset_in_single_chain(pdb_in_path=f'{output_dir}/{counter}.pdb',
                                                                      pdb_out_path=f'{output_dir}/{new_chain_list[counter - 1]}.pdb',
                                                                      offset=None,
                                                                      chain='A')
                    else:
                        offset = query_seq_length * (counter - 1) + 50 * (counter - 1)  # 50 glycines as offset !!!
                        change_chain_and_apply_offset_in_single_chain(pdb_in_path=f'{output_dir}/{counter}.pdb',
                                                                      pdb_out_path=f'{output_dir}/{new_chain_list[counter - 1]}.pdb',
                                                                      offset=offset,
                                                                      chain='A')
            list_of_paths_of_pdbs_to_merge = [f'{output_dir}/{ch}.pdb' for ch in new_chain_list]
            merge_pdbs(list_of_paths_of_pdbs_to_merge=list_of_paths_of_pdbs_to_merge,
                       merged_pdb_path=f'{output_dir}/{name}_template.pdb')
            os.system(f'rm {output_dir}/[A-Z].pdb')
            os.system(f'rm {output_dir}/[0-9].pdb')
            os.system(f'rm {output_dir}/{name}[A-Z]_template.pdb')
            os.system(f'rm {output_dir}/*ffdata')
            os.system(f'rm {output_dir}/*ffindex')
            os.system(f'rm {output_dir}/custom_output.hhr')
            os.system(f'rm {output_dir}/{name}.cif')

        if mode == 'assembly_from_aligned_custom_template':

            print_msg_box(msg=f'Assembly from aligned custom template: {template.split("/")[-1]} ({num_of_seqs} subunits)',
                          indent=2, title=None)

            name = template.split("/")[-1][:-4]
            os.system(f'cp {template} {output_dir}/{name}_template.pdb')

        h = Features(query_sequence=assembly_sequence_with_linkers)
        template_features = h.extract_template_features_from_aligned_pdb_and_sequence(
            pdb_path=f'{output_dir}/{name}_template.pdb',
            chain_ID='A')
        _template_features = copy.deepcopy(template_features)
        if add_to_templates:
            pass
        else:
            template_features = None
        if add_to_msa:
            sequence_from_template = _template_features['template_sequence'][0].decode('utf-8')
            h.append_row_in_msa(sequence=sequence_from_template, msa_uniprot_accession_identifiers=template)
            if add_to_templates:
                pass
            else:
                template_features = None
        else:
            sequence_from_template = None

        return sequence_from_template, template_features


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




########################################################################################################################

'''
num_of_seqs=1, query_fasta_path=fasta_path, mode=monomer_from_pdb, template=pdb_chain_id, add_to_msa=Bool, add_to_templates=Bool, polyala_res_list=[], sum_prob = None
num_of_seqs=1, query_fasta_path=fasta_path, mode=monomer_from_custom_template, template=path, add_to_msa=Bool, add_to_templates=Bool, polyala_res_list=[], sum_prob = None
num_of_seqs=1, query_fasta_path=fasta_path, mode=monomer_from_aligned_custom_template, template=path, add_to_msa=Bool, add_to_templates=Bool, polyala_res_list=[], sum_prob = None

num_of_seqs=n, query_fasta_path=fasta_path, mode=assembly_from_pdb, template=pdb_id, add_to_msa=Bool, add_to_templates=Bool, polyala_res_list=[], sum_prob = None]
num_of_seqs=n, query_fasta_path=fasta_path, mode=assembly_from_custom_template, template=path, add_to_msa=Bool, add_to_templates=Bool, polyala_res_list=[], sum_prob = None
num_of_seqs=n, query_fasta_path=fasta_path, mode=assembly_from_aligned_custom_template, template=path, add_to_msa=Bool, add_to_templates=Bool, polyala_res_list=[], sum_prob = None
'''


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

    arcimboldo = arcimboldo_air.ArcimboldoAir(input_load)

    for template in filter(lambda x: x.aligned, arcimboldo.templates):
        logging.info(f'Looking for alignment between {template.pdb} and query sequence using hhsearch:')

        name_pdb = utils.get_path_name(template.pdb)
        cif_path = f'{arcimboldo.output_dir}/{name_pdb}.cif'
        utils.pdb2cif(pdb_in_path=template.pdb, cif_out_path=cif_path) 
        hhsearch.generate_hhsearch_db(template_cif_path=cif_path, output_dir=arcimboldo.output_dir)
        hhsearch.run_hhsearch(fasta_path=arcimboldo.fasta_path, pdb70_db=f'{arcimboldo.output_dir}/pdb70',
                         output_path=f'{arcimboldo.output_dir}/custom_output.hhr')




        
    sys.exit(1)
    sequence_from_template, template_features = generate_features(output_dir=OUTPUT_DIR,
                                                                    num_of_seqs=NUM_OF_SEQS,
                                                                    query_sequence=query_sequence,
                                                                    mode=mode,
                                                                    template=template,
                                                                    add_to_msa=add_to_msa,
                                                                    add_to_templates=add_to_templates,
                                                                    pdb70_db_path=PDB70_DB_PATH,
                                                                    sum_prob=sum_prob)
    if 'assembly' not in mode:
        name = template.split("/")[-1][:-4]
    if 'monomer_from_pdb' in mode:
        name = f'{template[:-1]}_{template[-1]}'
    else:
        name = template
    if sequence_from_template:
        f.append_row_in_msa(sequence=sequence_from_template, msa_uniprot_accession_identifiers=name)
        print(f'Sequence from template \"{name}\" was already added to features.')
    if template_features:
        f.append_new_template_features(new_template_features=template_features, custom_sum_prob=sum_prob)
        print(f'Template \"{name}\" was already added to features.')

    features = merge_features(sequence_features=f.sequence_features,
                            msa_features=f.msa_features,
                            template_features=f.template_features)
    write_pkl_from_features(features=features, out_path=f'{OUTPUT_DIR}/features.pkl')









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







