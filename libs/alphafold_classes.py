import glob
import os
import logging
import shutil
import subprocess
from libs import utils, features

class AlphaFoldRun:
    def __init__ (self, output_dir:str, fasta_path: str, use_features: bool, feature: features.Features = None):
        self.run_alphafold_bash: str
        self.results_dir: str
        self.fasta_path: str
        self.use_Features: bool
        self.feature: features.Features = None

        self.feature = feature
        self.use_features = use_features
        self.results_dir = output_dir
        utils.create_dir(self.results_dir,delete_if_exists=True)
        self.fasta_path = os.path.join(self.results_dir, f'{os.path.basename(output_dir)}.fasta')
        self.run_alphafold_bash = os.path.join(self.results_dir, 'run_af2.sh')
        shutil.copy2(fasta_path, self.fasta_path)


    def run_af2(self):

        af2_output = subprocess.Popen(['bash', self.run_alphafold_bash], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        while True:
            line = af2_output.stdout.readline()
            if not line:
                break
            logging.debug(line.decode('utf-8'))
        af2_output.wait()
        return_code = af2_output.poll()
        logging.debug(f'AlphaFold2 return code is {return_code}')
        if af2_output.returncode != 0:
            raise Exception('AlphaFold2 stopped abruptly. Check the logfile')
        logging.info('AlphaFold2 has finshed succesfully. Proceeding to analyse the results')


    def create_af2_script(self, alphafold_paths):
        #Create the script to launch alphafold. It contins all the databases,
        #paths to the outputdir and fasta.

        previous_path = utils.get_parent_folder(dir_path=self.results_dir)

        with open(self.run_alphafold_bash, 'w') as bash_file:
            bash_file.write('#!/bin/bash\n')
            bash_file.write(f'python {alphafold_paths.run_alphafold_script} \\\n')
            bash_file.write(f'--fasta_paths={self.fasta_path} \\\n')
            bash_file.write(f'--output_dir={previous_path} \\\n')
            bash_file.write(f'--data_dir={alphafold_paths.af2_dbs_path} \\\n')
            bash_file.write(f'--uniref90_database_path={alphafold_paths.uniref90_db_path} \\\n')
            bash_file.write(f'--mgnify_database_path={alphafold_paths.mgnify_db_path} \\\n')
            bash_file.write(f'--template_mmcif_dir={alphafold_paths.mmcif_db_path} \\\n')
            bash_file.write('--max_template_date=2022-03-09 \\\n')
            bash_file.write(f'--obsolete_pdbs_path={alphafold_paths.obsolete_mmcif_db_path} \\\n')
            bash_file.write('--model_preset=monomer \\\n')
            bash_file.write(f'--bfd_database_path={alphafold_paths.bfd_db_path} \\\n')
            bash_file.write(f'--uniclust30_database_path={alphafold_paths.uniclust30_db_path} \\\n')
            bash_file.write(f'--pdb70_database_path={alphafold_paths.pdb70_db_path} \\\n')
            bash_file.write(f'--read_features_pkl={self.use_features}\n')
            bash_file.close()

class AlphaFoldPaths:

    def __init__ (self, af2_dbs_path: str):
        self.run_alphafold_script: str
        self.run_alphafold_bash: str
        self.af2_dbs_path: str
        self.mgnify_db_path: str
        self.uniref90_db_path: str
        self.mmcif_db_path: str
        self.obsolete_mmcif_db_path: str
        self.bfd_db_path: str
        self.uniclust30_db_path: str
        self.pdb70_db_path: str

        self.af2_dbs_path = af2_dbs_path
        self.run_alphafold_script = f'{utils.get_main_path()}/ALPHAFOLD/run_alphafold.py'

        for db in os.listdir(f'{self.af2_dbs_path}'):
            if 'mgnify' in db:
                self.mgnify_db_path = glob.glob(f'{self.af2_dbs_path}/{db}/*.fa', recursive=True)[0]
                logging.info(f'Mgnify DB path: {self.mgnify_db_path}')
            elif 'uniref90' in db:
                self.uniref90_db_path = glob.glob(f'{self.af2_dbs_path}/{db}/*.fasta', recursive=True)[0]
                logging.info(f'Uniref90 DB path {self.uniref90_db_path}')
            elif 'pdb_mmcif' in db:
                self.mmcif_db_path = f'{self.af2_dbs_path}/{db}/mmcif_files'
                self.obsolete_mmcif_db_path = f'{self.af2_dbs_path}/{db}/obsolete.dat'
                logging.info(f'mmCIF DB path: {self.mmcif_db_path}')
                logging.info(f'Obsolte mmCIF path: {self.obsolete_mmcif_db_path}')
            elif 'bfd' in db:
                self.bfd_db_path = '_'.join(glob.glob(f'{self.af2_dbs_path}/{db}/*', recursive=True)[0].split('_')[:-1])
                logging.info(f'BFD DB path: {self.bfd_db_path}')
            elif 'uniclust30' in db:
                for file in glob.glob(f'{self.af2_dbs_path}/{db}/**/*', recursive=True):
                    if '.cs219' in file[-6:]:
                        self.uniclust30_db_path = file.split('.')[:-1][0]
                        logging.info(f'Uniclust30 DB path: {self.uniclust30_db_path}')
            elif 'pdb70' in db:
                self.pdb70_db_path = f'{self.af2_dbs_path}/{db}/pdb70'
                logging.info(f'PDB70 DB path: {self.pdb70_db_path}')
    
    def __repr__(self):
        return f' \
        run_alphafold_script: {self.run_alphafold_script} \n \
        af2_dbs_path: {self.af2_dbs_path} \n \
        mgnify_db_path: {self.mgnify_db_path} \n \
        uniref90_db_path: {self.uniref90_db_path} \n \
        mmcif_db_path: {self.mmcif_db_path} \n \
        obsolete_mmcif_db_path: {self.obsolete_mmcif_db_path} \n \
        bfd_db_path: {self.bfd_db_path} \n \
        uniclust30_db_path: {self.uniclust30_db_path} \n \
        pdb70_db_path: {self.pdb70_db_path}'