import logging
import os
from typing import List, Dict, Union
from libs import alphafold_classes, bioutils, template, utils, features, sequence


class StructureAir:

    def __init__(self, parameters_dict: Dict):

        self.output_dir: str
        self.run_dir: str
        self.input_dir: str
        self.log_path: str
        self.sequence_assembled = sequence.SequenceAssembled
        self.afrun_list: List[alphafold_classes.AlphaFoldRun] = []
        self.alphafold_paths: alphafold_classes.AlphaFoldPaths
        self.templates_list: List[template.Template] = []
        self.run_af2: bool = True
        self.verbose: bool = True
        self.glycines: int = 50
        self.template_positions_list: List = [List]
        self.reference: Union[template.Template, None] = None
        self.custom_features: bool = True
        self.experimental_pdb: Union[str, None] = None
        self.mosaic: Union[int, None] = None

        self.output_dir = utils.get_mandatory_value(input_load=parameters_dict, value='output_dir')
        self.run_dir = parameters_dict.get('run_dir', os.path.join(self.output_dir, 'run'))
        self.input_dir = os.path.join(self.run_dir, 'input')
        self.log_path = os.path.join(self.output_dir, 'output.log')

        utils.create_dir(self.output_dir)
        utils.create_dir(self.run_dir)
        utils.create_dir(self.input_dir)

        af2_dbs_path = utils.get_mandatory_value(input_load=parameters_dict, value='af2_dbs_path')
        self.run_af2 = parameters_dict.get('run_alphafold', self.run_af2)
        self.verbose = parameters_dict.get('verbose', self.verbose)
        self.custom_features = parameters_dict.get('custom_features', self.custom_features)
        self.mosaic = parameters_dict.get('mosaic', self.mosaic)

        experimental_pdb = parameters_dict.get('experimental_pdb', self.experimental_pdb)
        if experimental_pdb is not None:
            experimental_pdb = bioutils.check_pdb(experimental_pdb, self.input_dir)
            self.experimental_pdb = os.path.join(self.run_dir, os.path.basename(experimental_pdb))
            bioutils.generate_multimer_from_pdb(experimental_pdb, self.experimental_pdb)

        sequence_list = []
        if 'sequences' not in parameters_dict:
            raise Exception('No sequences detected. Mandatory input')
        else:
            for parameters_sequence in parameters_dict.get('sequences'):
                new_sequence = sequence.Sequence(parameters_sequence, self.input_dir)
                sequence_list.append(new_sequence)
        self.sequence_assembled = sequence.SequenceAssembled(sequence_list, self.glycines)

        if not os.path.exists(af2_dbs_path):
            raise Exception('af2_dbs_path does not exist')
        if 'templates' not in parameters_dict:
            logging.info('No templates detected')
        else:
            reference = parameters_dict.get('reference')
            for parameters_template in parameters_dict.get('templates'):
                new_template = template.Template(parameters_dict=parameters_template, output_dir=self.run_dir,
                                                 input_dir=self.input_dir,
                                                 num_of_copies=self.sequence_assembled.total_copies)
                self.templates_list.append(new_template)
                if new_template.pdb_id == reference:
                    self.reference = new_template

            for element in self.templates_list:
                element.set_reference_templates(self)

            self.order_templates_with_restrictions()

            if self.reference is None:
                self.reference = self.templates_list[0]

        self.alphafold_paths = alphafold_classes.AlphaFoldPaths(af2_dbs_path=af2_dbs_path)

    def get_template_by_id(self, pdb_id: str) -> Union[template.Template, None]:
        # Return the template matching the pdb_id

        for temp in self.templates_list:
            if temp.pdb_id == pdb_id:
                return temp
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
            for temp in old_templates_list:
                reference_list = temp.get_reference_list()
                if set(reference_list).issubset(new_templates_list):
                    new_templates_list.append(temp)
                    deleted_items.append(temp)
            old_templates_list = [x for x in old_templates_list if (x not in deleted_items)]
            if not deleted_items:
                raise Exception('The match conditions could not be applied, there is an endless loop')
        self.templates_list = new_templates_list

    def append_line_in_templates(self, new_list: List):
        # Add line to the template's matrix.
        # The list contains the position of the chains

        self.template_positions_list.append(new_list)

    def run_alphafold(self, features_list: List[features.Features]):
        # Create the script and run alphafold
        partitions = utils.chunk_string(len(self.sequence_assembled.sequence_assembled), self.mosaic)
        for i, feature in enumerate(features_list):
            name = f'results_{i}'
            path = os.path.join(self.run_dir, name)
            sequence_chunk = self.sequence_assembled.sequence_assembled[partitions[i][0]:partitions[i][1]]
            afrun = alphafold_classes.AlphaFoldRun(output_dir=path,
                                                sequence=sequence_chunk,
                                                custom_features=self.custom_features,
                                                feature=feature)
             
            self.afrun_list.append(afrun)
            afrun.run_af2(alphafold_paths=self.alphafold_paths)

    def merge_results(self):
        if len(self.afrun_list) == 1:
            self.run_dir = self.afrun_list[0].results_dir
        else:
            for self.afrun in self.afrun_list:
                print('change things')

    def __repr__(self) -> str:
        return f' \
        output_dir: {self.output_dir} \n \
        run_af2: {self.run_af2}'