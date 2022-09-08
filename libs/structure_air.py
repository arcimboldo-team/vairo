import os
import shutil
from typing import List, Dict
from urllib.parse import non_hierarchical
from libs import bioutils, template, utils, features, alphafold_paths

class StructureAir:

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
        self.templates_list: List[template.Template] = []
        self.run_af2: bool = False
        self.verbose: bool = True
        self.glycines: int = 50
        self.template_positions_list: List = [List]
        self.reference: template.Template = None
        
        self.output_dir = utils.get_mandatory_value(input_load=parameters_dict, value='output_dir')
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
        if not 'templates' in parameters_dict:
            raise Exception('No templates detected. Check if the [[template]] tag exists.')

        reference = parameters_dict.get('reference')
        for parameters_template in parameters_dict.get('templates'):
            new_template = template.Template(parameters_template, self.run_dir, self.input_dir, self.num_of_copies)
            self.templates_list.append(new_template)
            if new_template.pdb_id == reference:
                self.reference = new_template

        for element in self.templates_list:
            element.set_reference_templates(self)

        self.order_templates_with_restrictions()

        if self.reference is None:
            self.reference = self.templates_list[0]

        self.features = features.Features(query_sequence=self.query_sequence_assembled)
        self.alphafold_paths = alphafold_paths.AlphaFoldPaths(af2_dbs_path, self.run_dir)

    def get_template_by_id(self, pdb_id: str) -> template.Template:
        #Return the template matching the pdb_id

        for template in self.templates_list:
            if template.pdb_id == pdb_id:
                return template
        return None


    def order_templates_with_restrictions(self):
        # Order the templates list in order to meet the requiered dependencies
        # All the templates are goingt to be in order, so the references will be calculated
        # before needed

        new_templates_list = []
        old_templates_list = self.templates_list

        if self.reference is not None:
            new_templates_list.append(self.reference)
            old_templates_list.remove(self.reference)
            
        while old_templates_list:
            deleted_items = []
            for template in old_templates_list:
                reference_list = template.get_reference_list()
                if set(reference_list).issubset(new_templates_list):
                    new_templates_list.append(template)
                    deleted_items.append(template)
            old_templates_list = [x for x in old_templates_list if (x not in deleted_items)]
            if not deleted_items:
                raise Exception('The match conditions could not be applied, there is an endless loop')
        self.templates_list = new_templates_list
        
    def append_line_in_templates(self, new_list: List):

        self.template_positions_list.append(new_list)

    def check_if_assembly(self):
        #Check if it was assembled before the execution of the program
        #We check if our 'G' pattern matches the input pdb
        #If so, we can change the parameters to convert it to a multimer

        if self.num_of_copies == 1:
            split_result = self.query_sequence.split('G'*self.glycines)
            if len(split_result)>1 and len(set(split_result)) == 1:
                self.query_sequence = split_result[0]
                self.num_of_copies = len(split_result)
                self.query_sequence_assembled = self.generate_query_assembled()

    def generate_query_assembled(self) -> str:
        #Transform the query_sequence to an assembled one
        #Repeat the sequence num_of_copies times, and put 'G' between the copies

        return (self.query_sequence + self.glycines * 'G') * (int(self.num_of_copies)-1) + self.query_sequence

    def __repr__(self):
        return f' \
        output_dir: {self.output_dir} \n \
        fasta_path: {self.fasta_path} \n \
        query_sequence: {self.query_sequence} \n \
        query_sequence_assembled: {self.query_sequence_assembled} \n \
        num_of_copies: {self.num_of_copies} \n \
        run_af2: {self.run_af2}'