import os
import shutil
import subprocess

from libs import utils


def create_database_from_pdb(fasta_path: str, database_path: str, output_dir: str) -> str:
    name = utils.get_file_name(fasta_path)
    database_dir = os.path.join(output_dir, f'{name}_database')
    data_name = os.path.join(database_dir, name)
    if not os.path.exists(f'{data_name}_cs219.ffindex'):
        utils.create_dir(database_dir, delete_if_exists=True)
        subprocess.call(
            ['ffindex_from_fasta', '-s', f'{data_name}.ffdata', f'{data_name}.ffindex', fasta_path])
        subprocess.call(['hhblits_omp', '-i', f'{data_name}', '-d', database_path, '-oa3m',
                         f'{data_name}_a3m_wo_ss', '-n', '2', '-cpu', '4', '-v', '0'])
        shutil.copy2(os.path.join(database_dir, f'{name}_a3m_wo_ss.ffdata'),
                     os.path.join(database_dir, f'{name}_a3m.ffdata'))
        shutil.copy2(os.path.join(database_dir, f'{name}_a3m_wo_ss.ffindex'),
                     os.path.join(database_dir, f'{name}_a3m.ffindex'))
        subprocess.call(['ffindex_apply', f'{data_name}_a3m.ffdata', f'{data_name}_a3m.ffindex', '-i',
                         f'{data_name}_hhm.ffindex', '-d', f'{data_name}_hhm.ffdata', '--', 'hhmake',
                         '-i', 'stdin', '-o', 'stdout', '-v', '0'])
        subprocess.call(['cstranslate', '-f', '-x', '0.3', '-c', '4', '-I', 'a3m', '-i', f'{data_name}_a3m', '-o',
                         f'{data_name}_cs219'])
    return os.path.join(database_dir, name)


def run_hhsearch(fasta_path: str, database_path: str, output_path: str) -> str:
    out = subprocess.Popen(['hhsearch', '-i', fasta_path, '-o', output_path, '-maxseq',
                            '1000000', '-d', database_path, '-p', '20', '-Z', '250', '-loc', '-z', '1',
                            '-b', '1', '-B', '250', '-ssm', '2', '-sc', '1', '-seq', '1', '-dbstrlen', '10000',
                            '-norealign', '-maxres', '32000'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = out.communicate()
    hhr = stdout.decode('utf-8')

    return hhr
