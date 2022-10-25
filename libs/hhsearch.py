import logging
import os
import shutil
import subprocess
from libs import utils

def create_database_from_pdb(fasta_path: str, database_path: str, output_dir: str) -> str:

    name = utils.get_file_name(fasta_path)
    database_dir = os.path.join(output_dir, f'{name}_database')
    utils.create_dir(database_dir, delete_if_exists=True)
    subprocess.Popen(['ffindex_from_fasta','-s', f'{database_dir}/template.ffdata', f'{database_dir}/template.ffindex', fasta_path])
    subprocess.Popen(['hhblits_omp','-i', f'{database_dir}/template', '-d', database_path, '-oa3m', f'{database_dir}/template_a3m_wo_ss','-n','2','-v','0'])
    shutil.copy2(os.path.join(database_dir, 'template_a3m_wo_ss.ffdata'), os.path.join(database_dir, 'template_a3m.ffdata'))
    shutil.copy2(os.path.join(database_dir, 'template_a3m_wo_ss.ffindex'), os.path.join(database_dir, 'template_a3m.ffindex'))
    subprocess.Popen(['ffindex_apply',f'{database_dir}/template_a3m.ff{{data,index}}','-i', f'{database_dir}/template_hhm.ffindex','-d', f'{database_dir}/template_hhm.ffdata','--','hhmake','-i','stdin','-o','stdout','-v','0'])
    subprocess.Popen(['cstranslate','-f','-x','0.3','-c','4','-I','a3m','-i',f'{database_dir}/template_a3m','-o',f'{database_dir}/template_cs219'])
    return database_dir

def run_hhsearch(fasta_path: str, database_path, output_path: str) -> str:

    out = subprocess.Popen(['hhsearch','-i','fasta_path','-o',output_path,'-maxseq',
                            '1000000','-d',database_path,'-p','20','-Z','250','-loc','-z','1',
                            '-b','1','-B','250','-ssm','2','-sc','1','-seq','1','-dbstrlen','10000',
                            '-norealign','-maxres','32000'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = out.communicate()
    hhr = stdout.decode('utf-8')

    return hhr