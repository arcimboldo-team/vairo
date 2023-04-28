import os
import subprocess
import alphafold
from libs import utils, alphafold_classes


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
    subprocess.call(['ffindex_build', '-as', f'{data_name}_a3m.ffdata', f'{data_name}_a3m.index', a3m_path])
    subprocess.call(['ffindex_apply', f'{data_name}_a3m.ffdata', f'{data_name}_a3m.ffindex', '-i',
                     f'{data_name}_hhm.ffindex', '-d', f'{data_name}_hhm.ffdata', '--', 'hhmake',
                     '-i', 'stdin', '-o', 'stdout', '-v', '0'])
    subprocess.call(['cstranslate', '-f', '-x', '0.3', '-c', '4', '-I', 'a3m', '-i', f'{data_name}_a3m', '-o',
                     f'{data_name}_cs219'])
    return data_name


def run_hhsearch(a3m_path: str, database_path: str, output_path: str) -> str:
    out = subprocess.Popen(['hhsearch', '-i', a3m_path, '-o', output_path, '-maxseq',
                            '1000000', '-d', database_path, '-p', '20', '-Z', '250', '-loc', '-z', '1',
                            '-b', '1', '-B', '250', '-ssm', '2', '-sc', '1', '-seq', '1', '-dbstrlen', '10000',
                            '-norealign', '-maxres', '32000'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = out.communicate()
    hhr = stdout.decode('utf-8')

    return hhr
