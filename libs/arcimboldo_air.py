import os
import shutil
from typing import List, Dict
from libs import bioutils, template, utils, features, alphafold_paths, template

class ArcimboldoAir:

    def __init__ (self, parameters_dict: Dict):
        self.output_dir: str
        self.run_dir: str
        self.input_dir: str
        self.fasta_path: str
        self.query_sequence: str
        self.query_sequence_assembled: str
        self.num_of_copies: int
        self.features: features.Features
        self.alphafold_paths: alphafold_paths.AlphaFoldPaths
        self.templates: List[template.Template] = []
        self.run_af2: bool = False
        self.verbose: bool = True
        self.glycines: int = 50
        
        self.output_dir = utils.get_mandatory_value(input_load = parameters_dict, value = 'output_dir')
        self.run_dir = parameters_dict.get('run_dir', os.path.join(self.output_dir, "run"))
        self.input_dir = os.path.join(self.run_dir, "input")

        utils.create_dir(self.output_dir)
        utils.create_dir(self.run_dir)
        utils.create_dir(self.input_dir)

        fasta_path = utils.get_mandatory_value(input_load = parameters_dict, value = 'fasta_path')
        
        if not os.path.exists(fasta_path):
            raise Exception(f'{fasta_path} does not exist')
        else:
            self.fasta_path = os.path.join(self.input_dir, os.path.basename(fasta_path))
            shutil.copy2(fasta_path, self.fasta_path)
        
        self.num_of_copies = utils.get_mandatory_value(input_load = parameters_dict, value = 'num_of_copies')
        af2_dbs_path = utils.get_mandatory_value(input_load = parameters_dict, value = 'af2_dbs_path')
        self.run_af2 = parameters_dict.get('run_alphafold', self.run_af2)
        self.verbose = parameters_dict.get('verbose', self.verbose)
        self.query_sequence = bioutils.extract_sequence(fasta_path=self.fasta_path)
        self.query_sequence_assembled = self.generate_query_assembled()

        if not os.path.exists(af2_dbs_path):
            raise Exception('af2_dbs_path does not exist')
        if not 'template' in parameters_dict:
            raise Exception('No templates detected. Check if the [[template]] tag exists.')

        for parameters_template in parameters_dict['template']:
            self.templates.append(template.Template(parameters_template, self.run_dir, self.input_dir, self.num_of_copies))

        self.features = features.Features(query_sequence=self.query_sequence_assembled)
        self.alphafold_paths = alphafold_paths.AlphaFoldPaths(af2_dbs_path)

    def check_if_assembly(self):

        if self.num_of_copies == 1:
            split_result = self.query_sequence.split('G'*self.glycines)
            if len(split_result)>1 and len(set(split_result)) == 1:
                self.query_sequence = split_result[0]
                self.num_of_copies = len(split_result)
                self.query_sequence_assembled = self.generate_query_assembled()

    def generate_query_assembled(self) -> str:

        return (self.query_sequence + self.glycines * 'G') * (int(self.num_of_copies)-1) + self.query_sequence

    def __repr__(self):
        return f' \
        output_dir: {self.output_dir} \n \
        fasta_path: {self.fasta_path} \n \
        query_sequence: {self.query_sequence} \n \
        query_sequence_assembled: {self.query_sequence_assembled} \n \
        num_of_copies: {self.num_of_copies} \n \
        run_af2: {self.run_af2}'