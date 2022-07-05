#! /usr/bin/env python3

from libs.features import *
from libs.analyse import af_summary
import copy
import glob, os
import sys

def print_msg_box(msg, indent=1, title=None):

    lines = msg.split('\n')
    space = " " * indent

    width = max(map(len, lines))
    box = f'╔{"═" * (width + indent * 2)}╗\n'  # upper_border
    if title:
        box += f'║{space}{title:<{width}}{space}║\n'  # title
        box += f'║{space}{"-" * len(title):<{width}}{space}║\n'  # underscore
    box += ''.join([f'║{space}{line:<{width}}{space}║\n' for line in lines])
    box += f'╚{"═" * (width + indent * 2)}╝'  # lower_border
    print('\n')
    print(box)
    print('\n')


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

            print(f'Looking for alignment between {template} and query sequence using hhsearch:')
            if os.path.exists(f'{output_dir}/output.hhr'):
                print(f'hhr file already present in {output_dir}')
            else:
                run_hhsearch(query_fasta_path=QUERY_FASTA_PATH, pdb70_db=f'{pdb70_db_path}/pdb70',
                             output_path=f'{output_dir}/output.hhr')

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

            run_hhsearch(query_fasta_path=QUERY_FASTA_PATH, pdb70_db=f'{output_dir}/pdb70',
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

def af2_sh(output_dir):

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

# OUTPUT_DIR = '/home/albert/Desktop/atzr_DNA/pep'
# QUERY_FASTA_PATH = '/home/albert/Desktop/atzr_DNA/atzr.fasta'
# NUM_OF_SEQS = 4

# BOR_PATH ='/home/albert/repos/arcimboldo-air/config.bor'
# AF2_DBS_PATH = '/storage2/af2/af2_dbs'

# RUN_AF2 = False


BOR_PATH = sys.argv[1] 

bor_text = open(BOR_PATH, 'r').read()

print('\n')
for line in bor_text.split('\n'):
    if line.startswith('OUTPUT_DIR'):
        OUTPUT_DIR = line.split('=')[1].replace(' ', '')
        print('Output directory:', OUTPUT_DIR)
    if line.startswith('QUERY_FASTA_PATH'):
        QUERY_FASTA_PATH = line.split('=')[1].replace(' ', '')
        print('Query fasta path:', QUERY_FASTA_PATH)
    if line.startswith('NUM_OF_SEQS'):
        NUM_OF_SEQS = int(line.split('=')[1].replace(' ', ''))
        print('Number of sequences:', NUM_OF_SEQS)
    if line.startswith('AF2_DBS_PATH'):
        AF2_DBS_PATH = line.split('=')[1].replace(' ', '')
    if line.startswith('RUN_AF2'):
        RUN_AF2 = line.split('=')[1].replace(' ', '')
        if RUN_AF2 == 'True':
            RUN_AF2 = True
            print('Run AF2:', RUN_AF2)
        else:
            RUN_AF2 = False
            print('Run AF2:', RUN_AF2)
         
try:
    print('\n')
    for db in os.listdir(f'{AF2_DBS_PATH}'):
        if 'mgnify' in db:
            MGNIFY_DB_PATH = glob.glob(f'{AF2_DBS_PATH}/{db}/*.fa', recursive=True)[0]
            print('Mgnify DB path:', MGNIFY_DB_PATH)
        elif 'uniref90' in db:
            UNIREF90_DB_PATH = glob.glob(f'{AF2_DBS_PATH}/{db}/*.fasta', recursive=True)[0]
            print('Uniref90 DB path', UNIREF90_DB_PATH)
        elif 'pdb_mmcif' in db:
            MMCIF_DB_PATH = f'{AF2_DBS_PATH}/{db}/mmcif_files'
            OBSOLETE_MMCIF_DB_PATH = f'{AF2_DBS_PATH}/{db}/obsolete.dat'
            print('mmCIF DB path:', MMCIF_DB_PATH)
            print('Obsolte mmCIF path:', OBSOLETE_MMCIF_DB_PATH)
        elif 'bfd' in db:
            BFD_DB_PATH = '_'.join(glob.glob(f'{AF2_DBS_PATH}/{db}/*', recursive=True)[0].split('_')[:-1])
            print('BFD DB path:', BFD_DB_PATH)
        elif 'uniclust30' in db:
            for file in glob.glob(f'{AF2_DBS_PATH}/{db}/**/*', recursive=True):
                if '.cs219' in file[-6:]:
                    UNICLUST30_DB_PATH = file.split('.')[:-1][0]
                    print('Uniclust30 DB path:', UNICLUST30_DB_PATH)
        elif 'pdb70' in db:
            PDB70_DB_PATH = f'{AF2_DBS_PATH}/{db}/pdb70'
            print('PDB70 DB path:', PDB70_DB_PATH)
except:
    print('Please, set all AF2 databases manually')
    sys.exit()


isExist = os.path.exists(OUTPUT_DIR)
if not isExist:
    os.makedirs(OUTPUT_DIR)
os.system(f'cp {BOR_PATH} {OUTPUT_DIR}')

fasta_name, query_sequence = extract_query_sequence(query_fasta_path=QUERY_FASTA_PATH)
assembly_sequence_with_linkers = ''
for num in range(NUM_OF_SEQS):
    assembly_sequence_with_linkers = assembly_sequence_with_linkers + query_sequence + 50 * 'G'
assembly_sequence_with_linkers = assembly_sequence_with_linkers[:-50]

f = Features(query_sequence=assembly_sequence_with_linkers)

job_parameters = bor_text[[m.end() for m in re.finditer(r'\[ARCIMBOLDO_AIR\]', bor_text)][0] + 1:]
job_lines = job_parameters.split('\n')[:-1]

for num, line in enumerate(job_lines):
    mode = line[:-1].split(',')[0].split('=')[1].replace(' ', '')
    template = line[:-1].split(',')[1].split('=')[1].replace(' ', '')
    add_to_msa = line[:-1].split(',')[2].split('=')[1].replace(' ', '')
    if add_to_msa == 'True':
        add_to_msa = True
    else:
        add_to_msa = False
    add_to_templates = line[:-1].split(',')[3].split('=')[1].replace(' ', '')
    if add_to_templates == 'True':
        add_to_templates = True
    else:
        add_to_templates = False
    polyala_res_list = line[:-1].split(',')[4].split('=')[1].replace(' ', '')
    sum_prob = line[:-1].split(',')[5].split('=')[1].replace(' ', '')
    if sum_prob == 'None':
        sum_prob = None

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








