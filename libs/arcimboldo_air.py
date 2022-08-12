import os
from typing import List, Dict
from libs import bioutils, template, utils, features, alphafold_paths, template

class ArcimboldoAir:

    def __init__ (self, parameters_dict: Dict):
        self.output_dir: str
        self.fasta_path: str
        self.query_sequence: str
        self.query_sequence_assembled: str
        self.num_of_copies: int
        self.features: features.Features
        self.alphafold_paths: alphafold_paths.AlphaFoldPaths
        self.templates: List[template.Template] = []
        self.run_af2: bool = False
        
        self.output_dir = utils.get_mandatory_value(input_load = parameters_dict, value = 'output_directory')
        self.fasta_path = utils.get_mandatory_value(input_load = parameters_dict, value = 'fasta_path')
        self.num_of_copies = utils.get_mandatory_value(input_load = parameters_dict, value = 'num_of_copies')
        af2_dbs_path = utils.get_mandatory_value(input_load = parameters_dict, value = 'af2_dbs_path')
        self.run_af2 = parameters_dict.get('run_alphafold', self.run_af2)
        self.query_sequence = bioutils.extract_sequence(fasta_path=self.fasta_path)
        self.query_sequence_assembled = (self.query_sequence + 50 * 'G') * (int(self.num_of_copies)-1) + self.query_sequence

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        if not os.path.exists(self.fasta_path):
            raise Exception('fasta_path does not exist')
        if not os.path.exists(af2_dbs_path):
            raise Exception('af2_dbs_path does not exist')
        if not 'template' in parameters_dict:
            raise Exception('No templates detected. Check if the [[template]] tag exists.')

        for parameters_template in parameters_dict['template']:
            self.templates.append(template.Template(parameters_template, self.output_dir, self.num_of_copies))

        self.features = features.Features(query_sequence=self.query_sequence_assembled)
        self.alphafold_paths = alphafold_paths.AlphaFoldPaths(af2_dbs_path)

    def __repr__(self):
        return f' \
        output_dir: {self.output_dir} \n \
        fasta_path: {self.fasta_path} \n \
        query_sequence: {self.query_sequence} \n \
        query_sequence_assembled: {self.query_sequence_assembled} \n \
        num_of_copies: {self.num_of_copies} \n \
        run_af2: {self.run_af2}'