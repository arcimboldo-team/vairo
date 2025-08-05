import copy
import itertools
import logging
import os
import shutil
from typing import List, Dict, Union
from libs import alphafold_classes, bioutils, output, template, utils, features, sequence, structures, plots, \
    template_modifications
from jinja2 import Environment, FileSystemLoader

MIN_RMSD_SPLIT = 5


class MainStructure:

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
        self.sequence_predicted_assembled = sequence.SequenceAssembled
        self.afrun_list: List[alphafold_classes.AlphaFoldRun] = []
        self.alphafold_paths: alphafold_classes.AlphaFoldPaths
        self.templates_list: List[template.Template] = []
        self.run_af2: bool
        self.small_bfd: bool
        self.cluster_templates: bool
        self.cluster_templates_msa: int
        self.cluster_templates_msa_mask: List[int]
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
        self.output: output.OutputStructure
        self.state: int = 0
        self.features_input: List[structures.FeaturesInput] = []
        self.features_list: List[features.Features] = []
        self.library_list: List[structures.Library] = []
        self.chunk_list: List[int] = []
        self.pymol_show_list: List[str] = []
        self.template_html_path: str

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
        self.output = output.OutputStructure(output_dir=self.output_dir)
        self.dir_templates_path = f'{utils.get_main_path()}/templates'
        self.template_html_path = os.path.join(self.dir_templates_path, 'output.html')

        utils.create_dir(self.run_dir)
        utils.create_dir(self.input_dir)
        utils.create_dir(self.experimental_dir, delete_if_exists=True)
        utils.delete_old_rankeds(self.output_dir)
        utils.delete_old_html(self.output_dir)

        if os.path.exists(self.output.pymol_session_path):
            os.remove(self.output.pymol_session_path)

        self.af2_dbs_path = utils.get_input_value(name='af2_dbs_path', section='global', input_dict=parameters_dict)
        if not os.path.exists(self.af2_dbs_path):
            raise Exception(f'Path {self.af2_dbs_path} does not exist. Check the af2_dbs_path input parameter')
        self.run_af2 = utils.get_input_value(name='run_af2', section='global', input_dict=parameters_dict)
        self.glycines = utils.get_input_value(name='glycines', section='global', input_dict=parameters_dict)
        self.mosaic = utils.get_input_value(name='mosaic', section='global', input_dict=parameters_dict)
        self.small_bfd = utils.get_input_value(name='small_bfd', section='global', input_dict=parameters_dict)
        pyoml_show_str = utils.get_input_value(name='show_pymol', section='global', input_dict=parameters_dict)

        if pyoml_show_str:
            self.pymol_show_list = str(pyoml_show_str).split(',')

        if self.mode == 'naive':
            self.cluster_templates = utils.get_input_value(name='cluster_templates', section='global',
                                                           input_dict=parameters_dict, override_default=True)
        else:
            self.cluster_templates = utils.get_input_value(name='cluster_templates', section='global',
                                                           input_dict=parameters_dict)

        self.cluster_templates_msa = utils.get_input_value(name='cluster_templates_msa', section='global',
                                                           input_dict=parameters_dict)
        self.cluster_templates_msa_mask = utils.expand_residues(
            utils.get_input_value(name='cluster_templates_msa_mask', section='global', input_dict=parameters_dict))

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
                    logging.info(
                        f'Not possible to generate the multimer for {utils.get_file_name(self.experimental_pdbs[-1])}')

        sequence_list = []
        sequence_prediced_list = []
        logging.error('Building query sequence')
        for parameters_sequence in utils.get_input_value(name='sequences', section='global',
                                                         input_dict=parameters_dict):
            new_sequence = sequence.Sequence(parameters_sequence, self.input_dir, self.run_dir, predicted=False)
            sequence_list.append(new_sequence)
            new_sequence = sequence.Sequence(parameters_sequence, self.input_dir, self.run_dir, predicted=True)
            sequence_prediced_list.append(new_sequence)

        self.sequence_assembled = sequence.SequenceAssembled(sequence_list, self.glycines)
        self.sequence_predicted_assembled = sequence.SequenceAssembled(sequence_prediced_list, self.glycines)

        for library in utils.get_input_value(name='append_library', section='global', input_dict=parameters_dict):
            path = utils.get_input_value(name='path', section='append_library', input_dict=library)
            aligned = utils.get_input_value(name='aligned', section='append_library', input_dict=library)
            add_to_msa = utils.get_input_value(name='add_to_msa', section='append_library', input_dict=library)
            add_to_templates = utils.get_input_value(name='add_to_templates', section='append_library',
                                                     input_dict=library)
            numbering_query = utils.get_input_value(name='numbering_query', section='append_library',
                                                    input_dict=library)
            numbering_library = utils.get_input_value(name='numbering_library', section='append_library',
                                                      input_dict=library)


            if os.path.exists(path):
                self.library_list.append(structures.Library(path=path, aligned=aligned,
                                                            add_to_msa=add_to_msa,
                                                            add_to_templates=add_to_templates,
                                                            numbering_query=numbering_query,
                                                            numbering_library=numbering_library))
            else:
                raise Exception(f'Path {path} does not exist. Check the input append_library parameter.')

        for parameters_features in utils.get_input_value(name='features', section='global', input_dict=parameters_dict):
            numbering_features = utils.get_input_value(name='numbering_features', section='features',
                                                       input_dict=parameters_features)

            positions = utils.get_input_value(name='positions', section='features', input_dict=parameters_features)
            numbering_query = utils.get_input_value(name='numbering_query', section='features',
                                                    input_dict=parameters_features)
            if numbering_query is None and positions is not None:
                numbering_query = f'{self.sequence_assembled.get_starting_length(positions - 1) + 1}'

            mutations = utils.get_input_value(name='mutations', section='features', input_dict=parameters_features)
            if mutations:
                mutations = utils.read_mutations_dict(mutations)
            self.features_input.append(structures.FeaturesInput(
                path=utils.get_input_value(name='path', section='features', input_dict=parameters_features),
                keep_msa=utils.get_input_value(name='keep_msa', section='features', input_dict=parameters_features),
                keep_templates=utils.get_input_value(name='keep_templates', section='features',
                                                     input_dict=parameters_features),
                msa_mask=utils.expand_residues(
                    utils.get_input_value(name='msa_mask', section='features', input_dict=parameters_features)),
                numbering_features=numbering_features,
                numbering_query=numbering_query,
                mutate_residues=mutations,
                replace_sequence=bioutils.check_sequence_path(
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
            translation_dict = {}
            for parameters_template in templates:
                pdb = utils.get_input_value(name='pdb', section='template', input_dict=parameters_template)
                pdb_name = utils.get_file_name(pdb)
                if pdb_name in translation_dict:
                    translation_dict[pdb_name].append(translation_dict[pdb_name][-1] + 1)
                else:
                    translation_dict[pdb_name] = [1]
            for parameters_template in templates:
                pdb = utils.get_input_value(name='pdb', section='template', input_dict=parameters_template)
                new_name = None
                pdb_name = utils.get_file_name(pdb)
                if pdb_name in translation_dict and sum(translation_dict[pdb_name]) != 1:
                    value = translation_dict[pdb_name].pop(0)
                    new_name = f'{pdb_name}_{value}'

                new_template = template.Template(parameters_dict=parameters_template, output_dir=self.run_dir,
                                                 num_of_copies=self.sequence_assembled.total_copies, new_name=new_name)
                self.templates_list.append(new_template)
                self.reference = new_template if new_template.pdb_id == self.reference else self.reference

            if self.reference is not None:
                self.templates_list.insert(0, self.templates_list.pop(self.reference))
            else:
                self.reference = self.templates_list[0]
        self.alphafold_paths = alphafold_classes.AlphaFoldPaths(af2_dbs_path=self.af2_dbs_path)

    def resize_features_predicted_sequence(self):
        numbering_query = []
        numbering_target = []
        for i, sequence_in in enumerate(self.sequence_predicted_assembled.sequence_list_expanded):
            if sequence_in.predict_region:
                ini = self.sequence_assembled.get_real_residue_number(i=i, residue=sequence_in.predict_region[0])
                end = self.sequence_assembled.get_real_residue_number(i=i, residue=sequence_in.predict_region[1])
                starting_seq = self.sequence_predicted_assembled.get_starting_length(i)
                numbering_query.append(starting_seq + 1)
                numbering_target.append(tuple([ini, end]))

        if numbering_query and numbering_target:
            modifications_list = utils.generate_modification_list(query=numbering_query, target=numbering_target,
                                                                length=self.sequence_predicted_assembled.length)

            self.feature = self.feature.cut_expand_features(self.sequence_predicted_assembled.sequence_assembled,
                                                            modifications_list)

    def expand_features_predicted_sequence(self):
        numbering_query = []
        numbering_target = []
        for i, sequence_in in enumerate(self.sequence_predicted_assembled.sequence_list_expanded):
            if sequence_in.predict_region:
                ini = self.sequence_predicted_assembled.get_starting_length(i) + 1
                end = self.sequence_predicted_assembled.get_finishing_length(i) + 1
                starting_seq = self.sequence_assembled.get_starting_length(i) + sequence_in.predict_region[0] - 1
                numbering_query.append(starting_seq + 1)
                numbering_target.append(tuple([ini, end]))

        if numbering_query and numbering_target:
            modifications_list = utils.generate_modification_list(query=numbering_query, target=numbering_target,
                                                                length=self.sequence_assembled.length)
            self.feature = self.feature.cut_expand_features(self.sequence_assembled.sequence_assembled, modifications_list)

    def partition_mosaic(self) -> List[features.Features]:
        if not self.mosaic_partition:
            self.chunk_list = self.sequence_predicted_assembled.partition(number_partitions=self.mosaic,
                                                                          overlap=self.mosaic_overlap)
        else:
            [self.chunk_list.append((partition[0] - 1, partition[1])) for partition in self.mosaic_partition]
        if self.feature is not None:
            self.resize_features_predicted_sequence()
            self.feature.select_msa_templates(sequence_assembled=self.sequence_predicted_assembled)
            self.features_list = self.feature.slicing_features(chunk_list=self.chunk_list)
        return self.features_list

    def render_output(self, reduced: bool):
        render_dict = {}

        template_str = open(self.template_html_path, 'r').read()
        jinja_template = Environment(loader=FileSystemLoader(self.dir_templates_path)).from_string(
            template_str)

        accepted_templates = self.output.templates_selected
        if reduced and os.path.exists(self.output.html_complete_path):
            render_dict['complete_html'] = self.output.html_complete_path

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

        if self.output.dendogram_struct:
            render_dict['dendogram_struct'] = self.output.dendogram_struct

        if self.cluster_templates:
            if os.path.exists(self.output.analysis_plot_path):
                render_dict['clustering_plot'] = utils.encode_data(input_data=self.output.analysis_plot_path)

            if os.path.exists(self.output.analysis_ranked_plot_path):
                render_dict['clustering_ranked_plot'] = utils.encode_data(
                    input_data=self.output.analysis_ranked_plot_path)

            if self.cluster_list:
                render_dict['cluster_list'] = self.cluster_list

        if self.templates_list:
            render_dict['templates_list'] = self.templates_list

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
                    info_input_list.append(
                        {'num_templates': library.num_templates, 'num_msa': library.num_msa, 'type': 'library',
                         'path': library.path})
                    sum_msa += library.num_msa
                    sum_templates += library.num_templates
                for feature in self.features_input:
                    info_input_list.append(
                        {'num_templates': feature.num_templates, 'num_msa': feature.num_msa, 'type': 'features',
                         'path': feature.path})
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
            ranked_qscore_dict = {}
            energies_dict = {}
            frobenius_dict = {}
            conclusion_dict = {}
            interfaces_dict = {'interfaces': {}, 'pdbs': {}}

            for pdb_in in self.output.ranked_filtered_list + self.output.experimental_list + self.output.templates_list:
                interfaces_dict['pdbs'][pdb_in.name] = pdb_in.interfaces
                for interface in pdb_in.interfaces:
                    if interface.name not in interfaces_dict['interfaces']:
                        interfaces_dict['interfaces'][interface.name] = None
                    if self.output.best_experimental is not None and pdb_in.name == self.output.best_experimental:
                        interfaces_dict['interfaces'][interface.name] = interface.deltaG

            render_dict['num_interfaces'] = self.output.num_interfaces

            for ranked in self.output.ranked_list:
                ranked_qscore_dict[ranked.name] = {}
                for ranked2 in self.output.ranked_list:
                    if ranked.name == ranked2.name:
                        ranked_qscore_dict[ranked.name][ranked.name] = 0
                    else:
                        ranked_qscore_dict[ranked.name][ranked2.name] = ranked.qscore_dict[ranked2.name]

                plddt_dict[ranked.name] = {'plddt': ranked.plddt, 'compactness': ranked.compactness,
                                           'ramachandran': ranked.ramachandran}
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

                if ranked.superposition_templates and any(
                        ranked_template.pdb in accepted_templates for ranked_template in
                        ranked.superposition_templates):
                    rmsd_dict[ranked.name] = {}
                    for ranked_template in ranked.superposition_templates:
                        if ranked_template.pdb in accepted_templates:
                            rmsd_dict[ranked.name][ranked_template.pdb] = {'qscore': ranked_template.qscore,
                                                                           'rmsd': ranked_template.rmsd,
                                                                           'aligned_residues': ranked_template.aligned_residues,
                                                                           'total_residues': ranked_template.total_residues
                                                                           }

                if ranked.frobenius_plots:
                    new_frobenius_plots = [plts for plts in ranked.frobenius_plots if
                                           plts.template in accepted_templates]
                    if new_frobenius_plots:
                        ordered_list = sorted(new_frobenius_plots, key=lambda x: x.core, reverse=True)
                        frobenius_plots_list = [ordered_list.pop(0)]
                        if ordered_list:
                            frobenius_plots_list.append(ordered_list.pop())
                        frobenius_dict[ranked.name] = frobenius_plots_list + ordered_list

            render_dict['bests_dict'] = {ranked.name: ranked for ranked in self.output.ranked_list if ranked.best}
            render_dict['filtered_dict'] = {ranked.name: ranked for ranked in self.output.ranked_filtered_list}

            if self.output.ranked_list:
                render_dict['ranked_list'] = self.output.ranked_list
            if self.output.group_ranked_by_qscore_dict:
                render_dict['ranked_by_qscore'] = self.output.group_ranked_by_qscore_dict
            if self.output.filtered_ranked_reason_dict:
                render_dict['filtered_ranked_reason'] = self.output.filtered_ranked_reason_dict
            if conclusion_dict:
                render_dict['conclusion_dict'] = conclusion_dict
                render_dict['conclusion_type'] = conclusion_type
            if ranked_qscore_dict:
                render_dict['table']['ranked_qscore_dict'] = ranked_qscore_dict
            if secondary_dict:
                render_dict['table']['secondary_dict'] = secondary_dict
            if rmsd_dict:
                render_dict['table']['rmsd_dict'] = rmsd_dict
            if energies_dict:
                render_dict['table']['energies_dict'] = energies_dict
            if interfaces_dict['interfaces']:
                render_dict['interfaces_dict'] = interfaces_dict
            if frobenius_dict:
                render_dict['frobenius_dict'] = frobenius_dict

            if self.output.experimental_dict:
                new_dict = copy.deepcopy(self.output.experimental_dict)
                for key, inner_dict in new_dict.items():
                    new_dict[key] = {k: v for k, v in inner_dict.items() if
                                     k in accepted_templates or k in ranked_qscore_dict.keys()}
                render_dict['table']['experimental_dict'] = new_dict

            self.output.write_tables(rmsd_dict=rmsd_dict, ranked_qscore_dict=ranked_qscore_dict,
                                     secondary_dict=secondary_dict, plddt_dict=plddt_dict,
                                     energies_dict=energies_dict)

        render_dict['state'] = self.get_state_text()
        render_dict['mode'] = self.mode.capitalize()

        if reduced:
            write_output = self.output.html_path
        else:
            write_output = self.output.html_complete_path

        if os.path.exists(self.output.pymol_session_path):
            render_dict['pymol'] = self.output.pymol_session_path
        elif self.output.ranked_list:
            render_dict['pymol'] = None

        jinja_template.globals['print_consecutive_numbers'] = utils.print_consecutive_numbers
        with open(write_output, 'w') as f_out:
            f_out.write(jinja_template.render(data=render_dict))

    def generate_output(self):
        if self.feature and self.feature.get_templates_length() > 20:
            self.render_output(reduced=True)
            self.render_output(reduced=False)
        else:
            self.render_output(reduced=True)

    def get_template_by_id(self, pdb_id: str) -> Union[template.Template, None]:
        # Return the template matching the pdb_id
        for temp in self.templates_list:
            if temp.pdb_id == pdb_id:
                return temp
        return None

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
                sequence_chunk = self.sequence_predicted_assembled.sequence_assembled[
                                 self.chunk_list[i][0]:self.chunk_list[i][1]]
            else:
                sequence_chunk = self.sequence_predicted_assembled.sequence_mutated_assembled[
                                 self.chunk_list[i][0]:self.chunk_list[i][1]]
            run_af2 = False if self.mode == 'guided' and self.cluster_templates else self.run_af2
            stop_after_msa = True if self.mode == 'naive' and self.cluster_templates else False
            afrun = alphafold_classes.AlphaFoldRun(results_dir=path,
                                                   sequence=sequence_chunk,
                                                   custom_features=self.custom_features,
                                                   stop_after_msa=stop_after_msa,
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
        for j, afrun in enumerate(self.afrun_list):
            best_ranked_list.append([])
            ranked_list = utils.read_rankeds(input_path=afrun.results_dir)
            if not ranked_list:
                logging.error('No predictions found')
                return
            [ranked.set_plddt() for ranked in ranked_list]
            ranked_list.sort(key=lambda x: x.plddt, reverse=True)

            # Select the two best rankeds; I don't think a superposition would help here.
            for i in range(2):
                new_ranked_path = os.path.join(best_rankeds_dir,
                                               f'ranked_{afrun.start_chunk + 1}-{afrun.end_chunk}_v{i}.pdb')
                shutil.copy2(ranked_list[i].path, new_ranked_path)
                ranked_list[i].set_path(path=new_ranked_path)
                best_ranked_list[j].append(ranked_list[i])

        if best_ranked_list:
            combinations = list(itertools.product(*best_ranked_list))
            for num, combination in enumerate(combinations):
                inf_path = combination[0].path
                merge_pdbs_list = [inf_path]
                for i, ranked in enumerate(combination[1:]):
                    len_sequence = len(bioutils.extract_sequence(self.afrun_list[i].fasta_path))

                    if self.mosaic_partition:
                        self.mosaic_overlap = self.mosaic_partition[i][1] - self.mosaic_partition[i + 1][0] + 1

                    inf_ini = len_sequence - self.mosaic_overlap + 1
                    inf_end = len_sequence
                    inm_ini = 1
                    inm_end = self.mosaic_overlap
                    pdb_out = os.path.join(best_rankeds_dir,
                                           f'{utils.get_file_name(ranked.path)}_{num}-{i}_superposed.pdb')
                    delta_out = os.path.join(best_rankeds_dir,
                                             f'{utils.get_file_name(ranked.path)}_{num}-{i}_deltas.dat')
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
                        logging.error('RMSD minimum requirements not met in order to merge the results in mosaic mode')
                        break

                    inf_cut = int(best_list[1][3])
                    inm_cut = int(best_list[2][1])

                    delete_residues = template_modifications.TemplateModifications()
                    delete_residues.append_modification(chains=['A'],
                                                        delete_residues=[*range(inf_cut + 1, len_sequence + 1, 1)])
                    delete_residues.modify_template(pdb_in_path=inf_path, pdb_out_path=inf_path, type_modify=['delete'])
                    delete_residues = template_modifications.TemplateModifications()
                    delete_residues.append_modification(chains=['A'], delete_residues=[*range(1, inm_cut, 1)])
                    delete_residues.modify_template(pdb_in_path=pdb_out, pdb_out_path=pdb_out, type_modify=['delete'])

                    merge_pdbs_list.append(pdb_out)
                    inf_path = pdb_out

                if len(merge_pdbs_list) == len(combination):
                    bioutils.merge_pdbs_in_one_chain(list_of_paths_of_pdbs_to_merge=merge_pdbs_list,
                                                     pdb_out_path=os.path.join(self.results_dir, f'ranked_{num}.pdb'))

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
        templates_cluster, _ = bioutils.cc_and_hinges_analysis(pdbs=self.output.templates_list,
                                                               binaries_path=self.binaries_paths,
                                                               output_dir=self.results_dir)
        if templates_cluster:
            logging.error(f'Templates can be grouped in {len(templates_cluster)} clusters')
            for cluster in templates_cluster:
                cluster_paths = [pdb.path for pdb in cluster]
                name_job = f'cluster_{counter}'
                label_job = f'Cluster {counter}'
                new_path = os.path.join(self.cluster_path, name_job)
                logging.error(f'Launching an VAIRO job in {new_path} with the following templates:')
                if len(cluster_paths) > 1:
                    logging.error(', '.join([utils.get_file_name(template_in) for template_in in cluster_paths]))
                else:
                    logging.error('Without templates')
                counter += 1
                yml_path = self.create_cluster(job_path=new_path, templates=cluster_paths)
                bioutils.run_vairo(yml_path=yml_path, input_path=new_path)
                rankeds = utils.read_rankeds(input_path=new_path)

                results_path = os.path.join(new_path, os.path.basename(self.run_dir),
                                            os.path.basename(self.results_dir))
                rankeds_path_list = []

                for ranked in rankeds:
                    rankeds_path_list.append(ranked.path)
                    nonsplit_filename = f'{ranked.name}.pdb'
                    if len(templates_cluster) > 1:
                        new_name = f'{name_job}_{ranked.name}.pdb'
                        shutil.copy2(os.path.join(results_path, nonsplit_filename),
                                     os.path.join(self.results_dir, new_name))
                    else:
                        shutil.copy2(os.path.join(results_path, nonsplit_filename), self.results_dir)

                if len(templates_cluster) <= 1:
                    logging.error(
                        'Only one cluster has been created, so all information will appear in the same output file.')
                    features_file_path = os.path.join(results_path, 'features.pkl')
                    self.set_feature(features.create_features_from_file(features_file_path))

                self.cluster_list.append(structures.Cluster(
                    name=name_job,
                    label=label_job,
                    path=new_path,
                    relative_path=os.path.join(os.path.basename(self.output_dir),
                                               os.path.relpath(new_path, self.output_dir),
                                               os.path.basename(self.output.html_path)),
                    encoded_path=utils.encode_data(os.path.join(new_path, os.path.basename(self.output.html_path))),
                    rankeds={utils.get_file_name(ranked_path): ranked_path for ranked_path in rankeds_path_list},
                    templates={pdb.name: pdb.path for pdb in cluster}
                ))

    def create_plot_gantt(self, reduced: bool):
        gantt_plots_both, legend_both = plots.plot_gantt(plot_type='both', plot_path=self.output.plots_path,
                                                         a_air=self, reduced=reduced)
        gantt_plots_template, legend_template = plots.plot_gantt(plot_type='templates',
                                                                 plot_path=self.output.plots_path,
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

    def extract_results(self, region_predicted: bool):
        self.output.extract_results(vairo_struct=self, region_predicted=region_predicted)

    def analyse_output(self):
        self.output.analyse_output(
            sequence_assembled=self.sequence_assembled,
            binaries_paths=self.binaries_paths,
            experimental_pdbs=self.experimental_pdbs
        )

    def align_experimental_pdbs(self):
        aligned_experimental_pdbs_list = []
        sequence_list = [sequence.fasta_path for sequence in self.sequence_assembled.sequence_list_expanded]
        for experimental in self.experimental_pdbs:
            try:
                pdb_out_path = os.path.join(self.experimental_dir,
                                            f'{utils.get_file_name(experimental)}_experimental.pdb')
                experimental_aligned_path = bioutils.align_pdb(pdb_in_path=experimental, pdb_out_path=pdb_out_path,
                                                               sequences_list=sequence_list,
                                                               databases=self.alphafold_paths)
                if experimental_aligned_path is None:
                    raise Exception()
                aligned_experimental_pdbs_list.append(experimental_aligned_path)
            except Exception as e:
                logging.error(f'Not possible to align experimental pdb {experimental}')
                aligned_experimental_pdbs_list.append(experimental)
                pass
        self.experimental_pdbs = aligned_experimental_pdbs_list

    def delete_mutations(self) -> str:
        logging.error('Proceding to launch ARICMBOLDO_AIR in order to delete the mutations')
        if not self.output.ranked_list:
            return
        mutations_dir = os.path.join(self.run_dir, 'delete_mutations')
        utils.create_dir(dir_path=mutations_dir, delete_if_exists=False)
        mutations_run_dir = os.path.join(mutations_dir, os.path.basename(self.run_dir))
        mutations_results_dir = os.path.join(mutations_run_dir, 'results')
        old_results_dir = os.path.join(self.run_dir, 'old_results_dir')
        yml_path = os.path.join(mutations_dir, 'config.yml')
        best_ranked = self.output.ranked_list[0]
        pdb_path = shutil.copy2(best_ranked.path, os.path.join(mutations_dir, 'selected.pdb'))

        utils.delete_old_rankeds(self.output_dir)
        with open(yml_path, 'w') as f_out:
            f_out.write(f'mode: guided\n')
            f_out.write(f'output_dir: {mutations_dir}\n')
            f_out.write(f'run_dir: {mutations_run_dir}\n')
            f_out.write(f'af2_dbs_path: {self.af2_dbs_path}\n')
            f_out.write(f'glycines: {self.glycines}\n')
            f_out.write(f'run_af2: {self.run_af2}\n')
            f_out.write(f'\nsequences:\n')
            for sequence_in in self.sequence_assembled.sequence_list:
                f_out.write('-')
                f_out.write(f' fasta_path: {sequence_in.fasta_path}\n')
                f_out.write(f'  num_of_copies: {sequence_in.num_of_copies}\n')
                new_positions = [position + 1 if position != -1 else position for position in sequence_in.positions]
                f_out.write(f'  positions: {",".join(map(str, new_positions))}\n')
            f_out.write(f'\ntemplates:\n')
            f_out.write('-')
            f_out.write(f' pdb: {pdb_path}\n')
            f_out.write(f'  legacy: True\n')
            f_out.write(f'  modifications:\n')
            f_out.write(f'  -  chain: all\n')
            f_out.write(f'     mutations:\n')
            f_out.write(f'     - numbering_residues: 1-100000\n')
            f_out.write(f'       mutate_with: ALA\n')

        bioutils.run_vairo(yml_path=yml_path, input_path=mutations_dir)
        if os.path.exists(old_results_dir):
            shutil.rmtree(old_results_dir)
        shutil.move(self.results_dir, old_results_dir)
        if os.path.exists(self.results_dir):
            shutil.rmtree(self.results_dir)
        shutil.copytree(mutations_results_dir, self.results_dir)

    def create_cluster(self, job_path: str, templates: List[str]) -> str:
        yml_path = os.path.join(job_path, 'config.yml')
        features_path = os.path.join(job_path, 'features.pkl')
        utils.create_dir(dir_path=job_path, delete_if_exists=False)
        new_features = features.Features(self.sequence_assembled.sequence_assembled)
        for template_in in templates:
            index = self.feature.get_index_by_name(utils.get_file_name(template_in))
            template_dict = self.feature.get_template_by_index(index)
            new_features.set_template_features(new_templates=template_dict)
        total_msa = self.feature.get_msa_length() if self.cluster_templates_msa == -1 else self.cluster_templates_msa + 1
        if self.cluster_templates_msa != 0:
            new_features.set_msa_features(new_msa=self.feature.msa_features, start=1, finish=total_msa,
                                          delete_positions=self.cluster_templates_msa_mask)

        new_features.write_pkl(features_path)

        with open(yml_path, 'w') as f_out:
            f_out.write(f'mode: guided\n')
            f_out.write(f'output_dir: {job_path}\n')
            f_out.write(f'run_dir: {os.path.join(job_path, os.path.basename(self.run_dir))}\n')
            f_out.write(f'af2_dbs_path: {self.af2_dbs_path}\n')
            f_out.write(f'glycines: {self.glycines}\n')
            f_out.write(f'run_af2: {self.run_af2}\n')
            f_out.write(f'mosaic: {self.mosaic}\n')
            if self.mosaic_partition:
                txt_aux = []
                for partition in self.mosaic_partition:
                    txt_aux.append("-".join(map(str, partition)))
                f_out.write(f'mosaic_partition: {",".join(map(str, txt_aux))}\n')
            f_out.write(f'\nsequences:\n')
            for sequence_in in self.sequence_predicted_assembled.sequence_list:
                f_out.write(f'- fasta_path: {sequence_in.fasta_path}\n')
                f_out.write(f'  num_of_copies: {sequence_in.num_of_copies}\n')
                new_positions = [position + 1 if position != -1 else position for position in sequence_in.positions]
                f_out.write(f'  positions: {",".join(map(str, new_positions))}\n')
                if sequence_in.mutations_dict.items():
                    f_out.write(f'  mutations:\n')
                    for residue, values in sequence_in.mutations_dict.items():
                        f_out.write(f'  - {residue}: {",".join(map(str, values))}\n')
            f_out.write(f'\nfeatures:\n')
            f_out.write(f'- path: {features_path}\n')
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
            if self.pymol_show_list:
                f_out.write(f'show_pymol: {",".join(map(str, self.pymol_show_list))}\n')
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
                if self.cluster_templates_msa_mask:
                    f_out.write(
                        f'cluster_templates_msa_mask: {",".join(map(str, self.cluster_templates_msa_mask))}\n')
                if self.cluster_templates_sequence is not None:
                    f_out.write(f'cluster_templates_sequence: {self.cluster_templates_sequence}\n')
            if self.library_list:
                f_out.write(f'\nappend_library:\n')
                for library in self.library_list:
                    f_out.write(f'- path: {library.path}\n')
                    f_out.write(f'  aligned: {library.aligned}\n')
                    if library.add_to_msa:
                        f_out.write(f'  add_to_msa: {library.add_to_msa}\n')
                    if library.add_to_templates:
                        f_out.write(f'  add_to_templates: {library.add_to_templates}\n')
                    if library.numbering_query and library.numbering_library:
                        f_out.write(f'  numbering_query: {library.numbering_query}\n')
                        f_out.write(f'  numbering_library: {library.numbering_library}\n')
            if self.features_input:
                f_out.write(f'\nfeatures:\n')
                for feat in self.features_input:
                    f_out.write(f'- path: {feat.path}\n')
                    f_out.write(f'  keep_msa: {feat.keep_msa}\n')
                    f_out.write(f'  keep_templates: {feat.keep_templates}\n')
                    if feat.msa_mask:
                        f_out.write(f'  msa_mask: {",".join(map(str, feat.msa_mask))}\n')
                    f_out.write(f'  numbering_query: {feat.numbering_query}\n')
                    f_out.write(
                        f'  numbering_features: {feat.numbering_features}\n')
                    if feat.replace_sequence is not None:
                        f_out.write(f'  sequence: {feat.replace_sequence}\n')

                    if feat.mutate_residues:
                        f_out.write(f'  mutations:\n')
                        for residue, values in feat.mutate_residues.items():
                            f_out.write(f'  - {residue}: {",".join(map(str, values))}\n')
            f_out.write(f'\nsequences:\n')
            for sequence_in in self.sequence_predicted_assembled.sequence_list:
                f_out.write(f'- fasta_path: {sequence_in.fasta_path}\n')
                f_out.write(f'  num_of_copies: {sequence_in.num_of_copies}\n')
                new_positions = [position + 1 if position != -1 else position for position in sequence_in.positions]
                f_out.write(f'  positions: {",".join(map(str, new_positions))}\n')
                if sequence_in.predict_region:
                    f_out.write(f'  predict_region: {"-".join(map(str, sequence_in.predict_region))}\n')
                if sequence_in.mutations_dict.items():
                    f_out.write(f'  mutations:\n')
                    for residue, values in sequence_in.mutations_dict.items():
                        f_out.write(f'  - {residue}: {",".join(map(str, values))}\n')
            if self.templates_list:
                f_out.write(f'\ntemplates:\n')
                for template_in in self.templates_list:
                    f_out.write(f'- pdb: {template_in.pdb_path}\n')
                    f_out.write(f'  add_to_msa: {template_in.add_to_msa}\n')
                    f_out.write(f'  add_to_templates: {template_in.add_to_templates}\n')
                    f_out.write(f'  generate_multimer: {template_in.generate_multimer}\n')
                    f_out.write(f'  aligned: {template_in.aligned}\n')
                    f_out.write(f'  legacy: {template_in.legacy}\n')
                    f_out.write(f'  strict: {template_in.strict}\n')

                    if template_in.modifications_struct.modifications_list:
                        f_out.write(f'  modifications:\n')
                        for modification in template_in.modifications_struct.modifications_list:
                            f_out.write(f'  - chain: {modification.chain}\n')
                            if modification.position + 1:
                                f_out.write(f'    position: {modification.position + 1}\n')
                            if modification.maintain_residues:
                                f_out.write(
                                    f'    maintain_residues: {", ".join(map(str, modification.maintain_residues))}\n')
                            if modification.delete_residues:
                                f_out.write(
                                    f'    delete_residues: {", ".join(map(str, modification.delete_residues))}\n')
                            f_out.write(f'      when: {modification.when}\n')

                            if modification.mutations:
                                f_out.write(f'      mutations:\n')
                                for mutation in modification.mutations:
                                    f_out.write(f'      - numbering_residues: {", ".join(map(str, mutation.mutate_residues_number))}\n')
                                    f_out.write(f'        mutate_with: {mutation.mutate_with}\n')

    def __repr__(self) -> str:
        return f' \
        output_dir: {self.output_dir} \n \
        run_af2: {self.run_af2}'
