import logging
import os
import shutil
import subprocess
from typing import List
import alphafold
from libs import utils, alphafold_classes, bioutils, features, structures


def create_a3m(fasta_path, databases: alphafold_classes.AlphaFoldPaths, output_dir: str) -> str:
    path = os.path.join(output_dir, f'{utils.get_file_name(fasta_path)}.a3m')
    hhblits = alphafold.data.tools.hhblits.HHBlits(binary_path='hhblits',
                                                   databases=[databases.bfd_db_path, databases.uniclust30_db_path])
    result = hhblits.query(fasta_path)
    with open(path, 'w+') as f_in:
        f_in.write(result[0]['a3m'])
    return path


def create_database_from_pdb(fasta_path: str, databases: alphafold_classes.AlphaFoldPaths, output_dir: str) -> str:
    name = utils.get_file_name(fasta_path)
    database_dir = os.path.join(output_dir, f'{name}_database')
    data_name = os.path.join(database_dir, name)
    utils.create_dir(database_dir, delete_if_exists=True)
    a3m_path = create_a3m(fasta_path, databases, database_dir)
    try:
        store_old_dir = os.getcwd()
        os.chdir(database_dir)

        subprocess.call(
            ['ffindex_build', '-as', f'{name}_a3m.ffdata', f'{name}_a3m.ffindex', os.path.basename(a3m_path)],
            stdout=subprocess.PIPE)
        subprocess.call(['ffindex_apply', f'{name}_a3m.ffdata', f'{name}_a3m.ffindex', '-i',
                         f'{name}_hhm.ffindex', '-d', f'{name}_hhm.ffdata', '--', 'hhmake',
                         '-i', 'stdin', '-o', 'stdout', '-v', '0'], stdout=subprocess.PIPE)
        subprocess.call(['cstranslate', '-f', '-x', '0.3', '-c', '4', '-I', 'a3m', '-i', f'{name}_a3m', '-o',
                         f'{name}_cs219'], stdout=subprocess.PIPE)
    finally:
        os.chdir(store_old_dir)
    if not os.path.exists(f'{data_name}_cs219.ffindex'):
        raise Exception(f'Could not create alignment for chain {utils.get_file_name(fasta_path)}.')
    return data_name


def run_hhsearch(a3m_path: str, database_path: str, output_path: str) -> str:
    out = subprocess.Popen(['hhsearch', '-i', a3m_path, '-o', output_path, '-maxseq',
                            '1000000', '-d', database_path, '-p', '20', '-Z', '250', '-loc', '-z', '1',
                            '-b', '1', '-B', '250', '-ssm', '2', '-sc', '1', '-seq', '1', '-dbstrlen', '10000',
                            '-norealign', '-maxres', '32000'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = out.communicate()
    hhr = stdout.decode('utf-8')

    return hhr


def run_hhalign(fasta_ref_path: str, fasta_aligned_path: str, output_path: str) -> str:
    out = subprocess.Popen(['hhalign', '-i', fasta_ref_path, '-t', fasta_aligned_path, '-o', output_path],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = out.communicate()
    hhr = stdout.decode('utf-8')

    return hhr


def run_hh(output_dir: str, database_dir: str, query_sequence_path: str, chain_in_path: str,
           databases: alphafold_classes.AlphaFoldPaths, name: str = None):
    query_sequence = bioutils.extract_sequence(query_sequence_path)
    template_sequence = bioutils.extract_sequence_from_file(file_path=chain_in_path)
    aux_key = list(template_sequence.keys())[0].replace('>', '')
    sequence_name = aux_key.split(':')[0]
    sequence_chain = aux_key.split(':')[1]
    database_chain_dir = os.path.join(database_dir, sequence_chain)

    if name is None:
        name = utils.get_file_name(chain_in_path)

    template_fasta_path = os.path.join(database_chain_dir, f'{name}.fasta')
    cif_path = os.path.join(output_dir, f'{name}.cif')
    hhr_path = os.path.join(output_dir, f'{name}.hhr')
    bioutils.pdb2mmcif(pdb_in_path=chain_in_path, cif_out_path=cif_path)

    if not os.path.exists(database_chain_dir):
        utils.create_dir(database_chain_dir)
        bioutils.write_sequence(sequence_name=f'{sequence_name}:{sequence_chain}',
                                sequence_amino=list(template_sequence.values())[0],
                                sequence_path=template_fasta_path)
        create_database_from_pdb(fasta_path=template_fasta_path, databases=databases, output_dir=database_chain_dir)

    run_hhalign(fasta_ref_path=query_sequence_path, fasta_aligned_path=template_fasta_path, output_path=hhr_path)

    template_features, mapping, identities, aligned_columns, total_columns, evalue = \
        features.extract_template_features_from_pdb(
            query_sequence=query_sequence,
            hhr_path=hhr_path,
            cif_path=cif_path,
            chain_id=sequence_chain
        )

    if int(aligned_columns) < int(total_columns * 0.95):
        hhr_path2 = os.path.join(output_dir, f'{utils.get_file_name(template_fasta_path)}2.hhr')
        a3m_path = os.path.join(output_dir, f'{utils.get_file_name(template_fasta_path)}.a3m')
        if not os.path.exists(a3m_path):
            a3m_path = create_a3m(fasta_path=query_sequence_path,
                                  databases=databases,
                                  output_dir=output_dir)
        run_hhsearch(a3m_path=a3m_path, database_path=database_chain_dir, output_path=hhr_path2)

        try:
            template_features2, mapping2, identities2, aligned_columns2, total_columns2, evalue2 = \
                features.extract_template_features_from_pdb(
                    query_sequence=query_sequence,
                    hhr_path=hhr_path2,
                    cif_path=cif_path,
                    chain_id=sequence_chain)
        except:
            aligned_columns2 = 0
            pass

        if int(aligned_columns) < int(aligned_columns2):
            os.remove(hhr_path)
            shutil.move(hhr_path2, hhr_path)
            template_features, mapping, identities, aligned_columns, total_columns, evalue = \
                template_features2, mapping2, identities2, aligned_columns2, total_columns2, evalue2

    extracted_chain_path = None

    logging.info(
        f'Alignment results for pdb {utils.get_file_name(chain_in_path)} and chain {sequence_chain}, with sequence {utils.get_file_name(query_sequence_path)}:')
    logging.info(
        f'Aligned columns: {aligned_columns} ({total_columns}), Evalue: {evalue}, Identities: {identities}')
    if template_features is not None:
        g = features.Features(query_sequence=query_sequence)
        g.append_new_template_features(new_template_features=template_features)
        aux_dict = g.write_all_templates_in_features(output_dir=output_dir, chain=sequence_chain)
        extracted_chain_path = list(aux_dict.values())[0]

    alignment_chain_struct = structures.Alignment(hhr_path=hhr_path, identities=identities,
                                                  aligned_columns=aligned_columns,
                                                  total_columns=total_columns, evalue=evalue,
                                                  mapping=mapping, chain=sequence_chain)
    return extracted_chain_path, alignment_chain_struct
