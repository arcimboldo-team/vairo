import glob
import os
import logging
import run_alphafold
from typing import Any, Union
from libs import bioutils, features, utils


class AlphaFoldRun:
    def __init__(self, results_dir: str, sequence: str, custom_features: bool, cluster_templates: bool, small_bfd: bool, start_chunk: int,
                 end_chunk: int, feature: features.Features = None):
        self.run_alphafold_bash: str
        self.results_dir: str
        self.fasta_path: str
        self.custom_features: bool
        self.cluster_templates: bool
        self.small_bfd: bool
        self.feature: Union[features.Features, None] = None
        self.start_chunk: int
        self.end_chunk: int

        self.feature = feature
        self.custom_features = custom_features
        self.cluster_templates = cluster_templates
        self.small_bfd = small_bfd
        self.results_dir = results_dir
        self.start_chunk = start_chunk
        self.end_chunk = end_chunk
        utils.create_dir(self.results_dir, delete_if_exists=False)
        self.fasta_path = os.path.join(self.results_dir, f'{os.path.basename(results_dir)}.fasta')
        bioutils.write_sequence(sequence_name=utils.get_file_name(self.fasta_path), sequence_amino=sequence,
                                sequence_path=self.fasta_path)
        self.run_alphafold_bash = os.path.join(self.results_dir, 'run_af2.sh')

    def run_af2(self, alphafold_paths):

        logging.info(f'Running AlphaFold2 in directory {self.results_dir}')

        previous_path = utils.get_parent_folder(dir_path=self.results_dir)
        if self.custom_features:
            self.feature.write_pkl(os.path.join(self.results_dir, 'features.pkl'))
        try:
            run_alphafold.launch_alphafold2(
                fasta_path=[self.fasta_path],
                output_dir=previous_path,
                data_dir=alphafold_paths.af2_dbs_path,
                max_template_date='2022-10-10',
                model_preset='monomer',
                uniref90_database_path=alphafold_paths.uniref90_db_path,
                mgnify_database_path=alphafold_paths.mgnify_db_path,
                template_mmcif_dir=alphafold_paths.mmcif_db_path,
                obsolete_pdbs_path=alphafold_paths.obsolete_mmcif_db_path,
                bfd_database_path=alphafold_paths.bfd_db_path,
                uniclust30_database_path=alphafold_paths.uniclust30_db_path,
                pdb70_database_path=alphafold_paths.pdb70_db_path,
                small_bfd_database_path=alphafold_paths.small_bfd_path,
                small_bfd=self.small_bfd,
                read_features_pkl=self.custom_features,
                stop_after_msa=self.cluster_templates)

        except SystemExit as e:
            pass
        except:
            raise Exception('AlphaFold2 stopped abruptly. Check the logfile')

        logging.info('AlphaFold2 has finished successfully. Proceeding to analyse the results')


class AlphaFoldPaths:

    def __init__(self, af2_dbs_path: str):
        self.af2_dbs_path: str
        self.mgnify_db_path: str
        self.uniref90_db_path: str
        self.mmcif_db_path: str
        self.obsolete_mmcif_db_path: str
        self.bfd_db_path: str = ''
        self.uniclust30_db_path: str = ''
        self.pdb70_db_path: str
        self.small_bfd_path: str = ''

        self.af2_dbs_path = af2_dbs_path

        for db in os.listdir(f'{self.af2_dbs_path}'):
            if 'mgnify' == db:
                self.mgnify_db_path = glob.glob(f'{self.af2_dbs_path}/{db}/*.fa', recursive=True)[0]
                logging.info(f'Mgnify DB path: {self.mgnify_db_path}')
            elif 'uniref90' == db:
                self.uniref90_db_path = glob.glob(f'{self.af2_dbs_path}/{db}/*.fasta', recursive=True)[0]
                logging.info(f'Uniref90 DB path {self.uniref90_db_path}')
            elif 'pdb_mmcif' == db:
                self.mmcif_db_path = f'{self.af2_dbs_path}/{db}/mmcif_files'
                self.obsolete_mmcif_db_path = f'{self.af2_dbs_path}/{db}/obsolete.dat'
                logging.info(f'mmCIF DB path: {self.mmcif_db_path}')
                logging.info(f'Obsolte mmCIF path: {self.obsolete_mmcif_db_path}')
            elif 'bfd' == db:
                self.bfd_db_path = '_'.join(glob.glob(f'{self.af2_dbs_path}/{db}/*', recursive=True)[0].split('_')[:-1])
                logging.info(f'BFD DB path: {self.bfd_db_path}')
            elif 'uniclust30' == db:
                for file in glob.glob(f'{self.af2_dbs_path}/{db}/**/*', recursive=True):
                    if '.cs219' in file[-6:]:
                        self.uniclust30_db_path = file.split('.')[:-1][0]
                        logging.info(f'Uniclust30 DB path: {self.uniclust30_db_path}')
            elif 'pdb70' == db:
                self.pdb70_db_path = f'{self.af2_dbs_path}/{db}/pdb70'
                logging.info(f'PDB70 DB path: {self.pdb70_db_path}')
            elif 'small_bfd' == db:
                self.small_bfd_path = glob.glob(f'{self.af2_dbs_path}/{db}/*.fasta', recursive=True)[0]
                logging.info(f'Small BFD path {self.small_bfd_path}')

    def __repr__(self):
        return f' \
        af2_dbs_path: {self.af2_dbs_path} \n \
        mgnify_db_path: {self.mgnify_db_path} \n \
        uniref90_db_path: {self.uniref90_db_path} \n \
        mmcif_db_path: {self.mmcif_db_path} \n \
        obsolete_mmcif_db_path: {self.obsolete_mmcif_db_path} \n \
        bfd_db_path: {self.bfd_db_path} \n \
        uniclust30_db_path: {self.uniclust30_db_path} \n \
        pdb70_db_path: {self.pdb70_db_path}  \n \
        small_bfd_path: {self.small_bfd_path}'
