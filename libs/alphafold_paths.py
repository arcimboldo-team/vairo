import glob
import os
import logging
from libs import utils

class AlphaFoldPaths:

    def __init__ (self, af2_dbs_path: str, output_dir:str):
        self.run_alphafold_script: str
        self.run_alphafold_bash: str
        self.run_alphafold_log: str
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

        self.run_alphafold_bash = f'{output_dir}/run_af2.sh'
        self.run_alphafold_log = f'{output_dir}/af2_output.log'

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

    def create_af2_script(self, output_dir: str):

        with open(self.run_alphafold_bash, 'w') as bash_file:
            output_name = utils.get_file_name(output_dir)
            bash_file.write('#!/bin/bash\n')
            bash_file.write(f'python {self.run_alphafold_script} \\\n')
            bash_file.write(f'--fasta_paths={output_name}.fasta \\\n')
            bash_file.write(f'--output_dir={output_dir} \\\n')
            bash_file.write(f'--data_dir={self.af2_dbs_path} \\\n')
            bash_file.write(f'--uniref90_database_path={self.uniref90_db_path} \\\n')
            bash_file.write(f'--mgnify_database_path={self.mgnify_db_path} \\\n')
            bash_file.write(f'--template_mmcif_dir={self.mmcif_db_path} \\\n')
            bash_file.write('--max_template_date=2022-03-09 \\\n')
            bash_file.write(f'--obsolete_pdbs_path={self.obsolete_mmcif_db_path} \\\n')
            bash_file.write('--model_preset=monomer \\\n')
            bash_file.write(f'--bfd_database_path={self.bfd_db_path} \\\n')
            bash_file.write(f'--uniclust30_database_path={self.uniclust30_db_path} \\\n')
            bash_file.write(f'--pdb70_database_path={self.pdb70_db_path} \\\n')
            bash_file.write('--read_features_pkl=True\n')
            bash_file.close()
    
    def __repr__(self):
        return f' \
        run_alphafold_script: {self.run_alphafold_script} \n \
        run_alphafold_bash: {self.run_alphafold_bash} \n \
        run_alphafold_log: {self.run_alphafold_log} \n \
        af2_dbs_path: {self.af2_dbs_path} \n \
        mgnify_db_path: {self.mgnify_db_path} \n \
        uniref90_db_path: {self.uniref90_db_path} \n \
        mmcif_db_path: {self.mmcif_db_path} \n \
        obsolete_mmcif_db_path: {self.obsolete_mmcif_db_path} \n \
        bfd_db_path: {self.bfd_db_path} \n \
        uniclust30_db_path: {self.uniclust30_db_path} \n \
        pdb70_db_path: {self.pdb70_db_path}'