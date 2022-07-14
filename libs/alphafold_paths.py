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

        for db in os.listdir(f'{af2_dbs_path}'):
            if 'mgnify' in db:
                self.mgnify_db_path = glob.glob(f'{af2_dbs_path}/{db}/*.fa', recursive=True)[0]
                logger.info('Mgnify DB path:', self.mgnify_db_path)
            elif 'uniref90' in db:
                self.uniref90_db_path = glob.glob(f'{af2_dbs_path}/{db}/*.fasta', recursive=True)[0]
                logger.info('Uniref90 DB path', self.uniref90_db_path)
            elif 'pdb_mmcif' in db:
                self.mmcif_db_path = f'{af2_dbs_path}/{db}/mmcif_files'
                self.obsolete_mmcif_db_path = f'{af2_dbs_path}/{db}/obsolete.dat'
                logger.info('mmCIF DB path:', self.mmcif_db_path)
                logger.info('Obsolte mmCIF path:', self.obsolete_mmcif_db_path)
            elif 'bfd' in db:
                self.bfd_db_path = '_'.join(glob.glob(f'{af2_dbs_path}/{db}/*', recursive=True)[0].split('_')[:-1])
                logger.info('BFD DB path:', self.bfd_db_path)
            elif 'uniclust30' in db:
                for file in glob.glob(f'{af2_dbs_path}/{db}/**/*', recursive=True):
                    if '.cs219' in file[-6:]:
                        self.uniclust30_db_path = file.split('.')[:-1][0]
                        logger.info('Uniclust30 DB path:', self.uniclust30_db_path)
            elif 'pdb70' in db:
                self.pdb70_db_path = f'{af2_dbs_path}/{db}/pdb70'
                logger.info('PDB70 DB path:', self.pdb70_db_path)