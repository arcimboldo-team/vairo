import copy
import logging
import os
import shutil
from typing import List, Dict, Union
from libs import alphafold_classes, bioutils, change_res, output_air, template, utils, features, sequence, structures, \
    plots
from jinja2 import Environment, FileSystemLoader

MIN_RMSD_SPLIT = 5


class StructureAir:

    def __init__(self, parameters_dict: Dict):

        self.mode: str
        self.output_dir: str
        self.run_dir: str
        self.results_dir: str
        self.name_results_dir = 'results'
        self.input_dir: str
        self.input_path: str
        self.log_path: str
        self.log_extended_path: str
        self.binaries_path: str
        self.cluster_path: str
        self.cluster_list: List[structures.Cluster] = []
        self.af2_dbs_path: str
        self.binaries_paths: structures.BinariesPath
        self.sequence_assembled = sequence.SequenceAssembled
        self.afrun_list: List[alphafold_classes.AlphaFoldRun] = []
        self.alphafold_paths: alphafold_classes.AlphaFoldPaths
        self.templates_list: List[template.Template] = []
        self.run_af2: bool
        self.small_bfd: bool
        self.cluster_templates: bool
        self.cluster_templates_msa: int
        self.cluster_templates_msa_delete: List[int]
        self.cluster_templates_sequence: str
        self.glycines: int
        self.template_positions_list: List[List] = []
        self.reference: Union[template.Template, None]
        self.custom_features: bool
        self.experimental_pdbs: List[str] = []
        self.mosaic: Union[int, None]
        self.mosaic_overlap: int = 150
        self.mosaic_partition: List[int]
        self.mosaic_seq_partition: List[int]
        self.feature: Union[features.Features, None] = None
        self.output: output_air.OutputAir
        self.state: int = 0
        self.features_input: List[structures.FeaturesInput] = []
        self.features_list: List[features.Features] = []
        self.library_list: List[structures.Library] = []
        self.chunk_list: List[int] = []
        self.pymol_show_list: List[str] = []
        self.template_html_path: str
        self.pymol_template_path: str
        self.pymol_keys_dict: dict = {}


        self.output_dir = utils.get_input_value(name='output_dir', section='global', input_dict=parameters_dict)
        utils.create_dir(self.output_dir)
        self.log_path = os.path.join(self.output_dir, 'output.log')
        self.log_extended_path = os.path.join(self.output_dir, 'output_extended.log')
        utils.create_logger_dir(self.log_path, self.log_extended_path)
        self.mode = utils.get_input_value(name='mode', section='global', input_dict=parameters_dict)
        self.custom_features = False if self.mode == 'naive' else True
        self.run_dir = utils.get_input_value(name='run_dir', section='global', input_dict=parameters_dict)
        if self.run_dir is None:
            self.run_dir = os.path.join(self.output_dir, 'run')
        self.input_dir = os.path.join(self.output_dir, 'input')
        self.experimental_dir: str = os.path.join(self.output_dir, 'experimental_pdbs')
        self.results_dir = os.path.join(self.run_dir, self.name_results_dir)
        self.cluster_path = os.path.join(self.output_dir, 'clustering')
        self.input_path = os.path.join(self.input_dir, 'config.yml')
        self.binaries_path = os.path.join(utils.get_main_path(), 'binaries')
        self.binaries_paths = structures.BinariesPath(self.binaries_path)
        self.output = output_air.OutputAir(output_dir=self.output_dir)
        self.dir_templates_path = f'{utils.get_main_path()}/templates'
        self.template_html_path = os.path.join(self.dir_templates_path, 'output.html')
        self.pymol_template_path = os.path.join(self.dir_templates_path, 'pymol_script.py')

        utils.create_dir(self.run_dir)
        utils.create_dir(self.input_dir)
        utils.create_dir(self.experimental_dir, delete_if_exists=True)
        utils.delete_old_rankeds(self.output_dir)
        utils.delete_old_html(self.output_dir)
        if os.path.exists(self.output.pymol_script_path):
            os.remove(self.output.pymol_script_path)

        self.af2_dbs_path = utils.get_input_value(name='af2_dbs_path', section='global', input_dict=parameters_dict)
        if not os.path.exists(self.af2_dbs_path):
            raise Exception(f'Path {self.af2_dbs_path} does not exist. Check the af2_dbs_path input parameter')
        self.run_af2 = utils.get_input_value(name='run_af2', section='global', input_dict=parameters_dict)
        self.glycines = utils.get_input_value(name='glycines', section='global', input_dict=parameters_dict)
        self.mosaic = utils.get_input_value(name='mosaic', section='global', input_dict=parameters_dict)
        self.small_bfd = utils.get_input_value(name='small_bfd', section='global', input_dict=parameters_dict)
        pyoml_show_str = utils.get_input_value(name='pymol_show', section='global', input_dict=parameters_dict)
        if pyoml_show_str:
            self.pymol_show_list = pyoml_show_str.replace(' ', '').split(',')

        if self.mode == 'naive':
            self.cluster_templates = utils.get_input_value(name='cluster_templates', section='global', input_dict=parameters_dict, override_default=True)
        else:
            self.cluster_templates = utils.get_input_value(name='cluster_templates', section='global', input_dict=parameters_dict)
        
        self.cluster_templates_msa = utils.get_input_value(name='cluster_templates_msa', section='global',
                                                           input_dict=parameters_dict)
        self.cluster_templates_msa_delete = utils.expand_residues(
            utils.get_input_value(name='cluster_templates_msa_delete', section='global', input_dict=parameters_dict))

        self.cluster_templates_sequence = bioutils.check_sequence_path(
            utils.get_input_value(name='cluster_templates_sequence', section='global', input_dict=parameters_dict))
        self.mosaic_partition = utils.get_input_value(name='mosaic_partition', section='global',
                                                      input_dict=parameters_dict)
        self.mosaic_seq_partition = utils.get_input_value(name='mosaic_seq_partition', section='global',
                                                          input_dict=parameters_dict)

        experimental_string = utils.get_input_value(name='experimental_pdbs', section='global',
                                                    input_dict=parameters_dict)
        if experimental_string:
            experimental_list = experimental_string.replace(' ', '').split(',')
            for pdb in experimental_list:
                pdb_path = bioutils.check_pdb(pdb,
                                              f'{os.path.join(self.experimental_dir, utils.get_file_name(pdb))}.pdb')
                self.experimental_pdbs.append(os.path.join(self.experimental_dir, os.path.basename(pdb_path)))
                try:
                    bioutils.generate_multimer_from_pdb(self.experimental_pdbs[-1], self.experimental_pdbs[-1])
                except Exception as e:
                    logging.debug(
                        f'Not possible to generate the multimer for {utils.get_file_name(self.experimental_pdbs[-1])}')

        sequence_list = []
        logging.error('Building query sequence')
        for parameters_sequence in utils.get_input_value(name='sequences', section='global',
                                                         input_dict=parameters_dict):
            new_sequence = sequence.Sequence(parameters_sequence, self.input_dir)
            sequence_list.append(new_sequence)
        self.sequence_assembled = sequence.SequenceAssembled(sequence_list, self.glycines)

        for library in utils.get_input_value(name='append_library', section='global', input_dict=parameters_dict):
            path = utils.get_input_value(name='path', section='append_library', input_dict=library)
            aligned = utils.get_input_value(name='aligned', section='append_library', input_dict=library)
            add_to_msa = utils.get_input_value(name='add_to_msa', section='append_library', input_dict=library)
            add_to_templates = utils.get_input_value(name='add_to_templates', section='append_library', input_dict=library)
            positions = utils.get_input_value(name='positions', section='append_library', input_dict=library)
            positions_list = None
            if positions:
                positions_list = ['-'] * self.sequence_assembled.length
                for pos in positions:
                    for pos_lib, pos_query_list in pos.items():
                        pos_lib_exp = utils.expand_residues(pos_lib)
                        for pos_query in pos_query_list.split(','):
                            pos_query_exp = utils.expand_residues(pos_query)
                            if len(pos_lib_exp) != len(pos_query_exp):
                                raise Exception(
                                    'Wrong format in the positions in append_library. Residue range numbers '
                                    'mismatch')
                            for i, pos_aux in enumerate(pos_query_exp):
                                positions_list[pos_aux - 1] = pos_lib_exp[i] - 1

            if os.path.exists(path):
                self.library_list.append(structures.Library(path=path, aligned=aligned,
                                                        add_to_msa=add_to_msa,
                                                        add_to_templates=add_to_templates,
                                                        positions=positions,
                                                        positions_list=positions_list))
            else:
                raise Exception(f'Path {path} does not exist. Check the input append_library parameter.')

        for parameters_features in utils.get_input_value(name='features', section='global', input_dict=parameters_dict):
            positions = utils.get_input_value(name='positions', section='features', input_dict=parameters_features)
            if positions is None:
                positions = f'1-{self.sequence_assembled.total_copies}'
            self.features_input.append(structures.FeaturesInput(
                path=utils.get_input_value(name='path', section='features', input_dict=parameters_features),
                keep_msa=utils.get_input_value(name='keep_msa', section='features', input_dict=parameters_features),
                keep_templates=utils.get_input_value(name='keep_templates', section='features',
                                                     input_dict=parameters_features),
                msa_delete=utils.expand_residues(
                    utils.get_input_value(name='msa_delete', section='features', input_dict=parameters_features)),
                positions=utils.expand_residues(positions),
                sequence=bioutils.check_sequence_path(
                    utils.get_input_value(name='sequence', section='features', input_dict=parameters_features))
            ))

        if self.mosaic_partition:
            self.mosaic_partition = utils.expand_partition(self.mosaic_partition)
            self.mosaic = len(self.mosaic_partition)
        elif self.mosaic_seq_partition:
            expanded = utils.expand_partition(self.mosaic_seq_partition)
            for exp in expanded:
                self.mosaic_partition.append([self.sequence_assembled.get_starting_length(exp[0] - 1) + 1,
                                              self.sequence_assembled.get_starting_length(
                                                  exp[1] - 1) + self.sequence_assembled.get_sequence_length(
                                                  exp[1] - 1)])
            self.mosaic = len(self.mosaic_partition)

        self.reference = utils.get_input_value(name='reference', section='global', input_dict=parameters_dict)
        templates = utils.get_input_value(name='templates', section='global', input_dict=parameters_dict)
        if templates:
            counter = 0
            for parameters_template in templates:
                new_name = None
                pdb = utils.get_input_value(name='pdb', section='template', input_dict=parameters_template)
                result = self.get_template_by_id(utils.get_file_name(pdb))
                if result is not None:
                    counter += 1
                    new_name = f'{result.pdb_id}_{counter}'
                new_template = template.Template(parameters_dict=parameters_template, output_dir=self.run_dir,
                                                 num_of_copies=self.sequence_assembled.total_copies, new_name=new_name)
                self.templates_list.append(new_template)
                self.reference = new_template if new_template.pdb_id == self.reference else self.reference
            [element.set_reference_templates(self) for element in self.templates_list]
            self.order_templates_with_restrictions()
            self.reference = self.templates_list[0] if self.reference is None else self.reference
        self.alphafold_paths = alphafold_classes.AlphaFoldPaths(af2_dbs_path=self.af2_dbs_path)

    def partition_mosaic(self) -> List[features.Features]:
        if not self.mosaic_partition:
            self.chunk_list = self.sequence_assembled.partition(number_partitions=self.mosaic,
                                                                overlap=self.mosaic_overlap)
        else:
            [self.chunk_list.append((partition[0] - 1, partition[1])) for partition in self.mosaic_partition]
        if self.feature is not None:
            self.features_list = self.feature.slicing_features(chunk_list=self.chunk_list)
        return self.features_list


    def render_output(self, reduced: bool):
        render_dict = {}

        template_str = open(self.template_html_path, 'r').read()
        jinja_template = Environment(loader=FileSystemLoader(self.dir_templates_path)).from_string(
            template_str)
        
        if reduced:
            accepted_templates = self.output.templates_selected
            if os.path.exists(self.output.html_complete_path):
                render_dict['complete_html'] = self.output.html_complete_path
        else:
            accepted_templates = [temp for temp in self.output.templates_dict.keys()]

        render_dict['frobenius_equation'] = utils.encode_data(
            input_data=f'{utils.get_main_path()}/templates/frobenius_equation.png')
        render_dict['frobenius_equation2'] = utils.encode_data(
            input_data=f'{utils.get_main_path()}/templates/frobenius_equation2.png')
        render_dict['custom_features'] = self.custom_features
        render_dict['mosaic'] = self.mosaic
        render_dict['total_copies'] = self.sequence_assembled.total_copies
        render_dict['number_alignments'] = len(
            [template_path for template_list in self.template_positions_list for template_path in template_list if
             template_path is not None])

        with open(self.input_path, 'r') as f_in:
            render_dict['bor_text'] = f_in.read()

        with open(self.log_path, 'r') as f_in:
            render_dict['log_text'] = f_in.read()

        if os.path.exists(self.output.pymol_script_path):
            render_dict['pymol'] = {'script': self.output.pymol_script_path, 'keys_dict': self.pymol_keys_dict}

        if self.feature is not None:
            self.create_plot_gantt(reduced=reduced)
            if reduced:
                if self.output.gantt_plots is not None:
                    render_dict['gantt'] = self.output.gantt_plots
            else:
                if self.output.gantt_complete_plots is not None:
                    render_dict['gantt'] = self.output.gantt_complete_plots

        if os.path.exists(self.output.plddt_plot_path):
            render_dict['plddt'] = utils.encode_data(self.output.plddt_plot_path)

        if os.path.exists(self.output.sequence_plot_path):
            render_dict['sequence_plot'] = utils.encode_data(input_data=self.output.sequence_plot_path)

        if self.cluster_templates:
            if os.path.exists(self.output.analysis_plot_path):
                render_dict['clustering_plot'] = utils.encode_data(input_data=self.output.analysis_plot_path)

            if os.path.exists(self.output.analysis_ranked_plot_path):
                render_dict['clustering_ranked_plot'] = utils.encode_data(input_data=self.output.analysis_ranked_plot_path)

            if self.cluster_list:
                render_dict['cluster_list'] = self.cluster_list
        
        if self.feature:
            if self.mode != 'naive':
                info_input_list = []
                sum_msa = 1
                sum_templates = 0
                normal_input = {'num_templates': 0, 'num_msa': 0, 'type': 'user input'}
                for template in self.templates_list:
                    if template.add_to_msa:
                        normal_input['num_msa'] += 1
                    if template.add_to_templates:
                        normal_input['num_templates'] += 1
                sum_msa += normal_input['num_msa']
                sum_templates += normal_input['num_templates']
                info_input_list.append(normal_input)
                for library in self.library_list:
                    info_input_list.append({'num_templates': library.num_templates, 'num_msa': library.num_msa, 'type': 'library', 'path': library.path})
                    sum_msa += library.num_msa
                    sum_templates += library.num_templates
                for feature in self.features_input:
                    info_input_list.append({'num_templates': feature.num_templates, 'num_msa': feature.num_msa, 'type': 'features', 'path': feature.path})
                    sum_msa += feature.num_msa
                    sum_templates += feature.num_templates                
                render_dict['info_input'] = info_input_list
                render_dict['num_msa'] = sum_msa
                render_dict['num_templates'] = sum_templates
            else:
                render_dict['num_msa'] = self.feature.get_msa_length()
                render_dict['num_templates'] = self.feature.get_templates_length()

        if self.output.ranked_list:
            render_dict['table'] = {}
            plddt_dict = {}
            secondary_dict = {}
            rmsd_dict = {}
            ranked_rmsd_dict = {}
            energies_dict = {}
            interfaces_dict = {}
            frobenius_dict = {}
            conclusion_dict = {}

            for ranked in self.output.ranked_list:
                ranked_rmsd_dict[ranked.name] = {}
                for ranked2 in self.output.ranked_list:
                    if ranked.name == ranked2.name:
                        ranked_rmsd_dict[ranked.name][ranked.name] = 0
                    else:
                        ranked_rmsd_dict[ranked.name][ranked2.name] = ranked.rmsd_dict[ranked2.name]

                plddt_dict[ranked.name] = {'plddt': ranked.plddt, 'compactness': ranked.compactness, 'ramachandran': ranked.ramachandran}
                try:
                    secondary_dict[ranked.name] = {'ah': ranked.ah, 'bs': ranked.bs,
                                                   'number_total_residues': ranked.total_residues}
                except:
                    pass
                if ranked.potential_energy is not None:
                    energies_dict[ranked.name] = ranked.potential_energy


                if ranked.superposition_experimental:
                    conclusion_dict.setdefault(ranked.superposition_experimental[0].pdb, []).append(ranked.name)
                    conclusion_type = 'experimental'
                elif ranked.superposition_templates:
                    conclusion_dict.setdefault(ranked.superposition_templates[0].pdb, []).append(ranked.name)
                    conclusion_type = 'template'

                if ranked.superposition_templates and any(ranked_template.pdb in accepted_templates for ranked_template in ranked.superposition_templates):
                    rmsd_dict[ranked.name] = {}
                    for ranked_template in ranked.superposition_templates:
                        if ranked_template.pdb in accepted_templates:
                            rmsd_dict[ranked.name][ranked_template.pdb] = {'rmsd': ranked_template.rmsd,
                                                                                'aligned_residues': ranked_template.aligned_residues,
                                                                                'total_residues': ranked_template.total_residues
                                                                                }                

                if ranked.filtered and ranked.interfaces:
                    interfaces_list = [copy.deepcopy(interface) for interface in ranked.interfaces if interface.interface_template and any([inter.template in accepted_templates for inter in interface.interface_template])]
                    if interfaces_list:
                        for inter in interfaces_list:
                            inter.interface_template = [int_temp for int_temp in inter.interface_template if int_temp.template in accepted_templates]
                        interfaces_dict[ranked.name] = interfaces_list
                if ranked.frobenius_plots:
                    new_frobenius_plots = [plts for plts in ranked.frobenius_plots if plts.template in accepted_templates]
                    if new_frobenius_plots:
                        ordered_list = sorted(new_frobenius_plots, key=lambda x: x.core, reverse=True)
                        frobenius_plots_list = [ordered_list.pop(0)]
                        if ordered_list:
                            frobenius_plots_list.append(ordered_list.pop())
                        frobenius_dict[ranked.name] = frobenius_plots_list + ordered_list

            render_dict['bests_dict'] = {ranked.name: ranked for ranked in self.output.ranked_list if ranked.best}
            render_dict['filtered_dict'] = {ranked.name: ranked for ranked in self.output.ranked_list if
                                            ranked.filtered}

            if self.templates_list:
                render_dict['templates_list'] = self.templates_list
            if self.output.ranked_list:
                render_dict['ranked_list'] = self.output.ranked_list
            if self.output.group_ranked_by_rmsd_dict:
                render_dict['ranked_by_rmsd'] = self.output.group_ranked_by_rmsd_dict
            if self.output.template_interfaces:
                render_dict['template_interfaces'] = self.output.template_interfaces
            if conclusion_dict:
                render_dict['conclusion_dict'] = conclusion_dict
                render_dict['conclusion_type'] = conclusion_type
            if ranked_rmsd_dict:
                render_dict['table']['ranked_rmsd_dict'] = ranked_rmsd_dict
            if secondary_dict:
                render_dict['table']['secondary_dict'] = secondary_dict
            if rmsd_dict:
                render_dict['table']['rmsd_dict'] = rmsd_dict
            if energies_dict:
                render_dict['table']['energies_dict'] = energies_dict
            if interfaces_dict:
                render_dict['interfaces_dict'] = interfaces_dict
            if frobenius_dict:
                render_dict['frobenius_dict'] = frobenius_dict
            
            if self.output.experimental_dict:
                new_dict = copy.deepcopy(self.output.experimental_dict)
                for key, inner_dict in new_dict.items():
                    new_dict[key] = {k: v for k, v in inner_dict.items() if k in accepted_templates or k in ranked_rmsd_dict.keys()}
                render_dict['table']['experimental_dict'] = new_dict

            self.output.write_tables(rmsd_dict=rmsd_dict, ranked_rmsd_dict=ranked_rmsd_dict,
                                     secondary_dict=secondary_dict, plddt_dict=plddt_dict,
                                     energies_dict=energies_dict)

        render_dict['state'] = self.get_state_text()

        if reduced:
            write_output = self.output.html_path
        else:
            write_output = self.output.html_complete_path
        
        with open(write_output, 'w') as f_out:
            f_out.write(jinja_template.render(data=render_dict))


    def generate_output(self):
        self.generate_pymol_script()
        if self.feature and self.feature.get_templates_length() > 20:
            self.render_output(reduced=True)
            self.render_output(reduced=False)
        else:
            self.render_output(reduced=True)


    def generate_pymol_script(self):
        # Read the pymol script found in the templates directory. It has the configuration
        # Add information regarding the tempaltes and experimental pdbs
        # Set if user has selected a specific zone to zoom
        if self.output.ranked_list:
            pdb_list = [ranked.split_path for ranked in self.output.ranked_list]
            pdb_list.extend(self.experimental_pdbs)
            template_str = open(self.pymol_template_path, 'r').read()
            for i, pdb_path in enumerate(pdb_list):
                template_str += f'\ncmd.load("{pdb_path}", "{utils.get_file_name(pdb_path)}")'
                if i != 0:
                    template_str += f'\ncmd.disable("{utils.get_file_name(pdb_path)}")'
            
            for i, zoom in enumerate(self.pymol_show_list, 1):
                key = f'F{i}'
                self.pymol_keys_dict[key] = f'Zoom into residues {zoom}'
                template_str += f'\ncmd.set_key("{key}", lambda: cmd.zoom("resi {zoom}"))'
            
            with open(self.output.pymol_script_path, 'w') as f_out:
                f_out.write(template_str)        


    def get_template_by_id(self, pdb_id: str) -> Union[template.Template, None]:
        # Return the template matching the pdb_id
        for temp in self.templates_list:
            if temp.pdb_id == pdb_id:
                return temp
        return None

    def order_templates_with_restrictions(self):
        # Order the templates list in order to meet the required dependencies
        # All the templates are going to be in order, so the references will be calculated
        # before needed
        new_templates_list = []
        old_templates_list = self.templates_list

        if self.reference is not None:
            new_templates_list.append(self.reference)
            old_templates_list.remove(self.reference)

        while old_templates_list:
            deleted_items = []
            for temp in old_templates_list:
                reference_list = temp.get_templates_references()
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
        for i, feature in enumerate(features_list):
            if len(features_list) == 1:
                name = self.name_results_dir
            else:
                name = f'{self.name_results_dir}{i}'
            path = os.path.join(self.run_dir, name)
            if self.cluster_templates and self.mode == 'naive':
                sequence_chunk = self.sequence_assembled.sequence_assembled[
                                 self.chunk_list[i][0]:self.chunk_list[i][1]]
            else:
                sequence_chunk = self.sequence_assembled.sequence_mutated_assembled[
                                 self.chunk_list[i][0]:self.chunk_list[i][1]]
            run_af2 = False if self.mode == 'guided' and self.cluster_templates else self.run_af2
            afrun = alphafold_classes.AlphaFoldRun(results_dir=path,
                                                   sequence=sequence_chunk,
                                                   custom_features=self.custom_features,
                                                   cluster_templates=self.cluster_templates,
                                                   small_bfd=self.small_bfd,
                                                   start_chunk=self.chunk_list[i][0],
                                                   end_chunk=self.chunk_list[i][1],
                                                   run=run_af2,
                                                   feature=feature
                                                   )
            self.afrun_list.append(afrun)
            afrun.run_af2(alphafold_paths=self.alphafold_paths)

    def merge_results(self):
        best_rankeds_dir = os.path.join(self.results_dir, 'best_rankeds')

        utils.create_dir(self.results_dir, delete_if_exists=True)
        utils.create_dir(best_rankeds_dir, delete_if_exists=True)

        best_ranked_list = []
        for afrun in self.afrun_list:
            ranked_list = utils.read_rankeds(input_path=afrun.results_dir)
            if not ranked_list:
                logging.error('No predictions found')
                return
            plot_path = os.path.join(afrun.results_dir, 'plddt.png')
            plots.plot_plddt(plot_path=plot_path, ranked_list=ranked_list)
            ranked_list.sort(key=lambda x: x.plddt, reverse=True)

            new_ranked_path = os.path.join(best_rankeds_dir,
                                           f'ranked_{afrun.start_chunk + 1}-{afrun.end_chunk}.pdb')
            shutil.copy2(ranked_list[0].path, new_ranked_path)
            ranked_list[0].set_path(path=new_ranked_path)
            best_ranked_list.append(ranked_list[0])

        inf_path = best_ranked_list[0].path
        merge_pdbs_list = [inf_path]
        for i, ranked in enumerate(best_ranked_list[1:]):
            len_sequence = len(bioutils.extract_sequence(self.afrun_list[i].fasta_path))

            if self.mosaic_partition:
                self.mosaic_overlap = self.mosaic_partition[i][1] - self.mosaic_partition[i + 1][0] + 1

            inf_ini = len_sequence - self.mosaic_overlap + 1
            inf_end = len_sequence
            inm_ini = 1
            inm_end = self.mosaic_overlap
            pdb_out = os.path.join(best_rankeds_dir, f'{utils.get_file_name(ranked.path)}_superposed.pdb')
            delta_out = os.path.join(best_rankeds_dir, f'{utils.get_file_name(ranked.path)}_deltas.dat')
            bioutils.run_lsqkab(pdb_inf_path=inf_path,
                                pdb_inm_path=ranked.path,
                                fit_ini=inf_ini,
                                fit_end=inf_end,
                                match_ini=inm_ini,
                                match_end=inm_end,
                                pdb_out=pdb_out,
                                delta_out=delta_out
                                )
            best_list = []
            best_min = MIN_RMSD_SPLIT
            with open(delta_out, 'r') as f_in:
                lines = f_in.readlines()
                lines = [line.replace('CA', '').split() for line in lines]
                for deltas in zip(lines, lines[1:], lines[2:], lines[3:]):
                    deltas_sum = sum([float(delta[0]) for delta in deltas])
                    if deltas_sum <= best_min:
                        best_list = deltas
                        best_min = deltas_sum

            if not best_list:
                raise Exception('RMSD minimum requirements not met in order to merge the results in mosaic mode.')

            inf_cut = int(best_list[1][3])
            inm_cut = int(best_list[2][1])
            delete_residues = change_res.ChangeResidues(
                chain_res_dict={'A': [*range(inf_cut + 1, len_sequence + 1, 1)]})
            delete_residues.delete_residues(pdb_in_path=inf_path, pdb_out_path=inf_path)
            delete_residues = change_res.ChangeResidues(chain_res_dict={'A': [*range(1, inm_cut, 1)]})
            delete_residues.delete_residues(pdb_in_path=pdb_out, pdb_out_path=pdb_out)
            merge_pdbs_list.append(pdb_out)
            inf_path = pdb_out

        bioutils.merge_pdbs_in_one_chain(list_of_paths_of_pdbs_to_merge=merge_pdbs_list,
                                         pdb_out_path=os.path.join(self.results_dir, 'ranked_0.pdb'))

    def set_feature(self, feature: features.Features):
        self.feature = feature

    def change_state(self, state: int):
        self.state = state

    def get_state_text(self):
        return {
            '-1': 'Finished with errors, not completed',
            '0': 'Starting',
            '1': 'Template alignment',
            '2': 'Running AlphaFold2',
            '3': 'Finished'
        }[str(self.state)]

    def templates_clustering(self):
        counter = 0
        utils.create_dir(self.cluster_path, delete_if_exists=False)
        templates_cluster, _ = bioutils.cc_and_hinges_analysis(paths_split_in=self.output.templates_dict,
                                                               paths_nonsplit_in=self.output.templates_nonsplit_dict,
                                                               binaries_path=self.binaries_paths,
                                                               output_path=self.results_dir,
                                                               length_sequences=self.output.percentage_sequences)
        if templates_cluster:
            logging.error(
                f'The templates can be grouped in {len(templates_cluster)} clusters')
            for templates in templates_cluster:
                name_job = f'cluster_{counter}'
                label_job = f'Cluster {counter}'
                new_path = os.path.join(self.cluster_path, name_job)
                logging.error(f'Launching an ARCIMBOLDO_AIR job in {new_path} with the following templates:')
                logging.error(', '.join([utils.get_file_name(template_in) for template_in in templates]))
                counter += 1
                yml_path = self.create_cluster(job_path=new_path, templates=templates)
                bioutils.run_arcimboldo_air(yml_path=yml_path)
                rankeds = utils.read_rankeds(input_path=new_path)
                results_path = os.path.join(new_path, os.path.basename(self.run_dir),
                                            os.path.basename(self.results_dir))
                rankeds_path_list = []
                for ranked in rankeds:
                    rankeds_path_list.append(ranked.path)
                    nonsplit_path = os.path.join(results_path, f'{ranked.name}.pdb')
                    new_name = f'{name_job}_{ranked.name}.pdb'
                    shutil.copy2(nonsplit_path, os.path.join(self.results_dir, new_name))
                self.cluster_list.append(structures.Cluster(
                    name=name_job,
                    label=label_job,
                    path=new_path,
                    relative_path=os.path.join(os.path.basename(self.output_dir),
                                               os.path.relpath(new_path, self.output_dir),
                                               os.path.basename(self.output.html_path)),
                    rankeds={utils.get_file_name(ranked_path): ranked_path for ranked_path in rankeds_path_list},
                    templates={utils.get_file_name(template_in): template_in for template_in in templates}
                ))

    def create_plot_gantt(self, reduced: bool):
        gantt_plots_both, legend_both = plots.plot_gantt(plot_type='both', plot_path=self.output.plots_path,
                                                        a_air=self, reduced=reduced)
        gantt_plots_template, legend_template = plots.plot_gantt(plot_type='templates', plot_path=self.output.plots_path,
                                                                a_air=self, reduced=reduced)
        gantt_plots_msa, legend_msa = plots.plot_gantt(plot_type='msa', plot_path=self.output.plots_path, a_air=self)

        struct = structures.GanttPlot(plot_both=utils.encode_data(gantt_plots_both),
                                                    legend_both=legend_both,
                                                    plot_template=utils.encode_data(gantt_plots_template),
                                                    legend_template=legend_template,
                                                    plot_msa=utils.encode_data(gantt_plots_msa),
                                                    legend_msa=legend_msa)

        if reduced:
            self.output.gantt_plots = struct
        else:
            self.output.gantt_complete_plots = struct

        if self.sequence_assembled.total_copies > 1:
            plots.plot_sequence(plot_path=self.output.sequence_plot_path, a_air=self)


    def create_cluster(self, job_path: str, templates: List[str]) -> str:
        yml_path = os.path.join(job_path, 'config.yml')
        features_path = os.path.join(job_path, 'features.pkl')
        utils.create_dir(dir_path=job_path, delete_if_exists=False)
        new_features = features.Features(self.sequence_assembled.sequence_assembled)
        for template_in in templates:
            index = self.feature.get_index_by_name(utils.get_file_name(template_in))
            template_dict = self.feature.get_template_by_index(index)
            new_features.set_template_features(new_templates=template_dict,
                                               sequence_in=self.cluster_templates_sequence)
        total_msa = self.feature.get_msa_length() if self.cluster_templates_msa == -1 else self.cluster_templates_msa + 1
        if self.cluster_templates_msa != 0:
            new_features.set_msa_features(new_msa=self.feature.msa_features, start=1, finish=total_msa,
                                          delete_positions=self.cluster_templates_msa_delete)
        new_features.write_pkl(features_path)

        with open(yml_path, 'w') as f_out:
            f_out.write(f'mode: guided\n')
            f_out.write(f'output_dir: {job_path}\n')
            f_out.write(f'run_dir: {os.path.join(job_path, os.path.basename(self.run_dir))}\n')
            f_out.write(f'af2_dbs_path: {self.af2_dbs_path}\n')
            f_out.write(f'glycines: {self.glycines}\n')
            f_out.write(f'mosaic: {self.mosaic}\n')
            if self.mosaic_partition:
                txt_aux = []
                for partition in self.mosaic_partition:
                    txt_aux.append("-".join(map(str, partition)))
                f_out.write(f'mosaic_partition: {",".join(map(str, txt_aux))}\n')
            f_out.write(f'\nsequences:\n')
            for sequence_in in self.sequence_assembled.sequence_list:
                f_out.write('-')
                f_out.write(f' fasta_path: {sequence_in.fasta_mutated_path}\n')
                f_out.write(f'  num_of_copies: {sequence_in.num_of_copies}\n')
                new_positions = [position + 1 if position != -1 else position for position in sequence_in.positions]
                f_out.write(f'  positions: {",".join(map(str, new_positions))}\n')
            f_out.write(f'\nfeatures:\n')
            f_out.write('-')
            f_out.write(f' path: {features_path}\n')
            f_out.write(f'  keep_msa: -1\n')
            f_out.write(f'  keep_templates: -1\n')
        return yml_path

    def write_input_file(self):
        with open(self.input_path, 'w') as f_out:
            f_out.write(f'mode: {self.mode}\n')
            f_out.write(f'output_dir: {self.output_dir}\n')
            f_out.write(f'run_dir: {self.run_dir}\n')
            f_out.write(f'af2_dbs_path: {self.af2_dbs_path}\n')
            f_out.write(f'glycines: {self.glycines}\n')
            f_out.write(f'run_af2: {self.run_af2}\n')
            if self.reference is not None:
                f_out.write(f'reference: {self.reference.pdb_path}\n')
            if self.experimental_pdbs:
                f_out.write(f'experimental_pdbs: {",".join(map(str, self.experimental_pdbs))}\n')
            f_out.write(f'small_bfd: {self.small_bfd}\n')
            f_out.write(f'mosaic: {self.mosaic}\n')
            if self.mosaic_partition:
                txt_aux = []
                for partition in self.mosaic_partition:
                    txt_aux.append("-".join(map(str, partition)))
                f_out.write(f'mosaic_partition: {",".join(map(str, txt_aux))}\n')
            f_out.write(f'cluster_templates: {self.cluster_templates}\n')
            if self.cluster_templates:
                f_out.write(f'cluster_templates_msa: {self.cluster_templates_msa}\n')
                if self.cluster_templates_msa_delete:
                    f_out.write(
                        f'cluster_templates_msa_delete: {",".join(map(str, self.cluster_templates_msa_delete))}\n')
                if self.cluster_templates_sequence is not None:
                    f_out.write(f'cluster_templates_sequence: {self.cluster_templates_sequence}\n')
            if self.library_list:
                f_out.write(f'\nappend_library:\n')
                for library in self.library_list:
                    f_out.write('-')
                    f_out.write(f' path: {library.path}\n')
                    f_out.write(f'  aligned: {library.aligned}\n')
                    if library.add_to_msa:
                        f_out.write(f'  add_to_msa: {library.add_to_msa}\n')
                    if library.add_to_templates:
                        f_out.write(f'  add_to_templates: {library.add_to_templates}\n')
                    if library.positions:
                        f_out.write(f'  positions:\n')
                        for positions in library.positions:
                            for lib, query in positions.items():
                                f_out.write(f'  -{lib}: {query}\n')
            if self.features_input:
                f_out.write(f'\nfeatures:\n')
                for feat in self.features_input:
                    f_out.write('-')
                    f_out.write(f' path: {feat.path}\n')
                    f_out.write(f'  keep_msa: {feat.keep_msa}\n')
                    f_out.write(f'  keep_templates: {feat.keep_templates}\n')
                    if feat.msa_delete:
                        f_out.write(f'  msa_delete: {",".join(map(str, feat.msa_delete))}\n')
                    if feat.positions:
                        f_out.write(f'  positions: {",".join(map(str, feat.positions))}\n')
                    if feat.sequence is not None:
                        f_out.write(f'  sequence: {feat.sequence}\n')
            f_out.write(f'\nsequences:\n')
            for sequence_in in self.sequence_assembled.sequence_list:
                f_out.write('-')
                f_out.write(f' fasta_path: {sequence_in.fasta_path}\n')
                f_out.write(f'  num_of_copies: {sequence_in.num_of_copies}\n')
                new_positions = [position + 1 if position != -1 else position for position in sequence_in.positions]
                f_out.write(f'  positions: {",".join(map(str, new_positions))}\n')
                if sequence_in.mutations_dict.items():
                    f_out.write(f'  mutations:\n')
                    for residue, values in sequence_in.mutations_dict.items():
                        f_out.write(f'  -{residue}: {",".join(map(str, values))}\n')
            if self.templates_list:
                f_out.write(f'\ntemplates:\n')
                for template_in in self.templates_list:
                    f_out.write('-')
                    f_out.write(f' pdb: {template_in.pdb_path}\n')
                    f_out.write(f'  add_to_msa: {template_in.add_to_msa}\n')
                    f_out.write(f'  add_to_templates: {template_in.add_to_templates}\n')
                    f_out.write(f'  generate_multimer: {template_in.generate_multimer}\n')
                    f_out.write(f'  aligned: {template_in.aligned}\n')
                    f_out.write(f'  legacy: {template_in.legacy}\n')
                    f_out.write(f'  strict: {template_in.strict}\n')
                    if template_in.reference is not None:
                        f_out.write(f'  reference: {self.reference}\n')
                    if template_in.change_res_struct.change_residues_list:
                        f_out.write(f'  change_res:\n')
                        for change in template_in.change_res_struct.change_residues_list:
                            f_out.write('  -')
                            if change.resname is not None:
                                f_out.write(f' resname: {change.resname}\n')
                            elif change.sequence is not None:
                                if change.fasta_path is not None:
                                    f_out.write(f' fasta_path: {change.fasta_path}\n')
                                else:
                                    f_out.write(f' fasta_path: {change.sequence}\n')
                            f_out.write(f'    when: {change.when}\n')
                            for key, value in change.chain_group_res_dict.items():
                                f_out.write(f'    {key}: {", ".join(map(str, value))}\n')

                    if template_in.match_restrict_struct.match_restrict_list:
                        f_out.write(f'  match:\n')
                        for match in template_in.match_restrict_struct.match_restrict_list:
                            f_out.write('  -')
                            f_out.write(f' chain: {match.chain}\n')
                            f_out.write(
                                f'    position: {match.position + 1 if match.position != -1 else match.position}\n')
                            if match.residues is not None:
                                f_out.write(
                                    f'    residues: {",".join(map(str, list(match.residues.chain_res_dict.values())[0]))}\n')
                            if match.reference is not None:
                                f_out.write(f'    reference: {match.reference}\n')
                            if match.reference_chain is not None:
                                f_out.write(f'    reference_chain: {match.reference_chain}\n')

    def __repr__(self) -> str:
        return f' \
        output_dir: {self.output_dir} \n \
        run_af2: {self.run_af2}'
