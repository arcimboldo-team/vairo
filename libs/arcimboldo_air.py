import os
from typing import List
from libs import template, utils, features, alphafold_paths, template

class ArcimboldoAir:

    def __init__ (self, parameters_list: List):
        self.output_dir: str
        self.fasta_path: str
        self.num_of_copies: int
        self.features: features.Features
        self.alphafold_paths: afpaths.AlphaFoldPaths
        self.templates: List[template.Template] = []
        self.run_af2: bool = False

        self.output_dir = utils.get_mandatory_value(input_load = parameters_list, value = 'output_directory')
        self.fasta_path = utils.get_mandatory_value(input_load = parameters_list, value = 'fasta_path')
        self.af2_dbs_path = utils.get_mandatory_value(input_load = parameters_list, value = 'af2_dbs_path')
        self.num_of_copies = utils.get_mandatory_value(input_load = parameters_list, value = 'num_of_copies')
        self.run_af2 = parameters_list.get('run_alphafold', self.run_af2)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        if not os.path.exists(self.fasta_path):
            raise Exception('fasta_path does not exist')
        if not os.path.exists(self.af2_dbs_path):
            raise Exception('af2_dbs_path does not exist')

        for parameters_template in parameters_list['template']:
            self.templates.append(template.Template(parameters_template, self.output_dir))

        self.features = features.Features(fasta_path=self.fasta_path, num_of_copies=self.num_of_copies)
        self.alphafold_paths = alphafold_paths.AlphaFoldPaths(self.af2_dbs_path)