import glob
import os
import logging

class AlphaFoldPaths:

    def __init__ (self, af2_dbs_path: str):
        self.af2_dbs_path: str
        self.mgnify_db_path: str
        self.uniref90_db_path: str
        self.mmcif_db_path: str
        self.obsolete_mmcif_db_path: str
        self.bfd_db_path: str
        self.uniclust30_db_path: str
        self.pdb70_db_path: str

        self.af2_dbs_path = af2_dbs_path

        for db in os.listdir(f'{self.af2_dbs_path}'):
            if 'mgnify' in db:
                self.mgnify_db_path = glob.glob(f'{self.af2_dbs_path}/{db}/*.fa', recursive=True)[0]
                logging.info('Mgnify DB path:', self.mgnify_db_path)
            elif 'uniref90' in db:
                self.uniref90_db_path = glob.glob(f'{self.af2_dbs_path}/{db}/*.fasta', recursive=True)[0]
                logging.info('Uniref90 DB path', self.uniref90_db_path)
            elif 'pdb_mmcif' in db:
                self.mmcif_db_path = f'{self.af2_dbs_path}/{db}/mmcif_files'
                self.obsolete_mmcif_db_path = f'{self.af2_dbs_path}/{db}/obsolete.dat'
                logging.info('mmCIF DB path:', self.mmcif_db_path)
                logging.info('Obsolte mmCIF path:', self.obsolete_mmcif_db_path)
            elif 'bfd' in db:
                self.bfd_db_path = '_'.join(glob.glob(f'{self.af2_dbs_path}/{db}/*', recursive=True)[0].split('_')[:-1])
                logging.info('BFD DB path:', self.bfd_db_path)
            elif 'uniclust30' in db:
                for file in glob.glob(f'{self.af2_dbs_path}/{db}/**/*', recursive=True):
                    if '.cs219' in file[-6:]:
                        self.uniclust30_db_path = file.split('.')[:-1][0]
                        logging.info('Uniclust30 DB path:', self.uniclust30_db_path)
            elif 'pdb70' in db:
                self.pdb70_db_path = f'{self.af2_dbs_path}/{db}/pdb70'
                logging.info('PDB70 DB path:', self.pdb70_db_path)
    
    def __repr__(self):
        logging.info('af2_dbs_path: ', self.af2_dbs_path)
        logging.info('mgnify_db_path: ', self.mgnify_db_path)
        logging.info('uniref90_db_path: ', self.uniref90_db_path)
        logging.info('mmcif_db_path: ', self.mmcif_db_path)
        logging.info('obsolete_mmcif_db_path: ', self.obsolete_mmcif_db_path)
        logging.info('bfd_db_path: ', self.bfd_db_path)
        logging.info('uniclust30_db_path: ', self.uniclust30_db_path)
        logging.info('pdb70_db_path: ', self.pdb70_db_path)