import logging
import os
import shutil
from typing import List, Dict, Union
from libs import alphafold_classes, bioutils, change_res, output_air, template, utils, features, sequence
from jinja2 import Environment, FileSystemLoader


MIN_RMSD_SPLIT = 5

class StructureAir:

    def __init__(self, parameters_dict: Dict):

        self.output_dir: str
        self.run_dir: str
        self.input_dir: str
        self.input_path: str
        self.log_path: str
        self.af2_dbs_path: str
        self.sequence_assembled = sequence.SequenceAssembled
        self.afrun_list: List[alphafold_classes.AlphaFoldRun] = []
        self.alphafold_paths: alphafold_classes.AlphaFoldPaths
        self.templates_list: List[template.Template] = []
        self.run_af2: bool = True
        self.verbose: bool = True
        self.small_bfd: bool = False
        self.glycines: int = 50
        self.template_positions_list: List[List] = []
        self.reference: Union[template.Template, None] = None
        self.custom_features: bool = True
        self.experimental_pdb: Union[str, None] = None
        self.mosaic: Union[int, None] = None
        self.mosaic_overlap: int = 150
        self.feature: Union[features.Features, None] = None
        self.output: output_air.OutputAir
        self.state: int = 0

        self.output_dir = utils.get_mandatory_value(input_load=parameters_dict, value='output_dir')
        self.run_dir = parameters_dict.get('run_dir', os.path.join(self.output_dir, 'run'))
        self.input_dir = os.path.join(self.run_dir, 'input')
        self.log_path = os.path.join(self.output_dir, 'output.log')
        self.input_path = os.path.join(self.input_dir, 'config.yml')
        self.output = output_air.OutputAir(output_dir=self.output_dir)

        utils.create_dir(self.output_dir)
        utils.create_dir(self.run_dir)
        utils.create_dir(self.input_dir)

        self.af2_dbs_path = utils.get_mandatory_value(input_load=parameters_dict, value='af2_dbs_path')
        self.run_af2 = parameters_dict.get('run_alphafold', self.run_af2)
        self.verbose = parameters_dict.get('verbose', self.verbose)
        self.custom_features = parameters_dict.get('custom_features', self.custom_features)
        self.mosaic = parameters_dict.get('mosaic', self.mosaic)
        self.small_bfd = parameters_dict.get('small_bfd', self.small_bfd)

        experimental_pdb = parameters_dict.get('experimental_pdb', self.experimental_pdb)
        if experimental_pdb is not None:
            experimental_pdb = bioutils.check_pdb(experimental_pdb, self.input_dir)
            self.experimental_pdb = os.path.join(self.run_dir, os.path.basename(experimental_pdb))
            try:
                bioutils.generate_multimer_from_pdb(experimental_pdb, self.experimental_pdb)
            except:
                logging.info('Not possible to generate the multimer for the experimental pdb')

        sequence_list = []
        if 'sequences' not in parameters_dict:
            raise Exception('No sequences detected. Mandatory input')
        else:
            for parameters_sequence in parameters_dict.get('sequences'):
                new_sequence = sequence.Sequence(parameters_sequence, self.input_dir)
                sequence_list.append(new_sequence)
        self.sequence_assembled = sequence.SequenceAssembled(sequence_list, self.glycines)

        if self.mosaic is None:
            self.mosaic = 1

        if not os.path.exists(self.af2_dbs_path):
            raise Exception('af2_dbs_path does not exist')
        if 'templates' not in parameters_dict:
            logging.info('No templates detected')
        else:
            reference = parameters_dict.get('reference')
            for parameters_template in parameters_dict.get('templates'):
                new_template = template.Template(parameters_dict=parameters_template, output_dir=self.run_dir,
                                                 input_dir=self.input_dir, num_of_copies=self.sequence_assembled.total_copies)
                self.templates_list.append(new_template)
                if new_template.pdb_id == reference:
                    self.reference = new_template

            for element in self.templates_list:
                element.set_reference_templates(self)

            self.order_templates_with_restrictions()

            if self.reference is None:
                self.reference = self.templates_list[0]

        self.alphafold_paths = alphafold_classes.AlphaFoldPaths(af2_dbs_path=self.af2_dbs_path)

    def generate_output(self):
        render_dict = {}

        template_str = open(f'{utils.get_main_path()}/templates/output.html', 'r').read()
        jinja_template = Environment(loader=FileSystemLoader(f'{utils.get_main_path()}/templates/')).from_string(
            template_str)

        render_dict['frobenius_equation'] = utils.encode_data(input_data=f'{utils.get_main_path()}/templates/frobenius_equation.png')
        render_dict['frobenius_equation2'] = utils.encode_data(input_data=f'{utils.get_main_path()}/templates/frobenius_equation2.png')
        render_dict['custom_features'] = self.custom_features
        render_dict['mosaic'] = self.mosaic
        render_dict['total_copies'] = self.sequence_assembled.total_copies
        render_dict['number_templates'] = len(self.templates_list)
        render_dict['number_alignments'] = len([template_path for template_list in self.template_positions_list for template_path in template_list if template_path is not None])

        with open(self.input_path, 'r') as f_in:
            render_dict['bor_text'] = f_in.read()

        with open(self.log_path, 'r') as f_in:
            render_dict['log_text'] = f_in.read()

        if self.feature is not None:
            self.output.create_plot_gantt(self)
            render_dict['gantt'] = [utils.encode_data(plot) for plot in self.output.gantt_plots_path]

        if os.path.exists(self.output.plddt_plot_path):
            render_dict['plddt'] = utils.encode_data(self.output.plddt_plot_path)

        if self.output.ranked_list:
            render_dict['table'] = {}
            plddt_dict = {}
            secondary_dict = {}
            rmsd_dict = {}
            energies_dict = {}
            bests_dict = {}
            interfaces_dict = {}
            frobenius_dict = {}
            template_dict = {}
            filtered_dict = {}
            green_color = 40
            blue_color = 40

            for ranked in self.output.ranked_list:
                if ranked.best:
                    bests_dict[ranked.name] = {
                        'name': ranked.name,
                        'path': ranked.path,
                        'template': ranked.superposition_templates[0],
                        'rmsd': ranked.rmsd,
                        'encoded': utils.encode_data(ranked.split_path),
                        'color': f'hsl(120, 100% ,{green_color}%)'
                    }
                    if ranked.superposition_templates:
                        template_dict.setdefault(ranked.superposition_templates[0].template, []).append(ranked.name)
                    green_color += 10
                elif ranked.filtered:
                    filtered_dict[ranked.name] = {
                        'name': ranked.name,
                        'path': ranked.path,
                        'encoded': utils.encode_data(ranked.split_path),
                        'color': f'hsl(229, 100% ,{blue_color}%)'
                    }
                    blue_color += 10

                plddt_dict[ranked.name] = round(ranked.plddt)
                secondary_dict[ranked.name] = {'ah': ranked.ah, 'bs': ranked.bs,
                                                'number_total_residues': ranked.total_residues
                                               }
                if ranked.energies is not None:
                    energies_dict[ranked.name] = {'kinetic': ranked.energies.kinetic,
                                                  'potential': ranked.energies.potential
                                                  }
                rmsd_dict[ranked.name] = {}
                for ranked_template in ranked.superposition_templates:
                    rmsd_dict[ranked.name][ranked_template.template] = {'rmsd': ranked_template.rmsd,
                                                                        'aligned_residues': ranked_template.aligned_residues,
                                                                        'total_residues': ranked_template.total_residues
                                                                        }
                if ranked.filtered and ranked.interfaces:
                    interfaces_dict[ranked.name] = {
                        'bests': [interface for interface in ranked.interfaces if interface.core > 0],
                        'others': [interface for interface in ranked.interfaces if interface.core == 0]
                    }

                if ranked.frobenius_plots:
                    ordered_list = sorted(ranked.frobenius_plots, key=lambda x: x.core, reverse=True)
                    frobenius_dict[ranked.name] = {}
                    frobenius_dict[ranked.name]['best'] = ordered_list.pop(0)
                    if ordered_list:
                        frobenius_dict[ranked.name]['worst'] = ordered_list.pop()
                    if ordered_list:
                        frobenius_dict[ranked.name]['others'] = ordered_list

            if self.templates_list:
                render_dict['templates_list'] = self.templates_list
            if self.output.ranked_list:
                render_dict['bests_list'] = self.output.ranked_list
            if bests_dict:
                render_dict['bests_dict'] = bests_dict
            if filtered_dict:
                render_dict['filtered_dict'] = filtered_dict
            if template_dict:
                render_dict['template_dict'] = template_dict
            if plddt_dict:
                render_dict['table']['plddt_dict'] = plddt_dict
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
                render_dict['table']['experimental_dict'] = self.output.experimental_dict

            self.output.write_tables(rmsd_dict=rmsd_dict, secondary_dict=secondary_dict, plddt_dict=plddt_dict,
                                     energies_dict=energies_dict)        
        
        render_dict['state'] = self.get_state_text()

        with open(self.output.html_path, 'w') as f_out:
            f_out.write(jinja_template.render(data=render_dict))

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
        partitions = utils.chunk_string(len(self.sequence_assembled.sequence_assembled), self.mosaic,
                                        overlap=self.mosaic_overlap)
        for i, feature in enumerate(features_list):
            name = f'results_{i}'
            path = os.path.join(self.run_dir, name)
            sequence_chunk = self.sequence_assembled.sequence_assembled[partitions[i][0]:partitions[i][1]]
            afrun = alphafold_classes.AlphaFoldRun(output_dir=path,
                                                   sequence=sequence_chunk,
                                                   custom_features=self.custom_features,
                                                   small_bfd=self.small_bfd,
                                                   start_chunk=partitions[i][0],
                                                   finish_chunk=partitions[i][1],
                                                   feature=feature)
            self.afrun_list.append(afrun)
            #afrun.run_af2(alphafold_paths=self.alphafold_paths)

    def merge_results(self):
        if len(self.afrun_list) == 1:
            self.run_dir = self.afrun_list[0].results_dir
        else:
            results_dir = os.path.join(self.run_dir, 'results')
            best_rankeds_dir = os.path.join(results_dir, 'best_rankeds')

            utils.create_dir(results_dir, delete_if_exists=True)
            utils.create_dir(best_rankeds_dir, delete_if_exists=True)

            best_ranked_list = []
            for afrun in self.afrun_list:
                ranked_list = output_air.read_rankeds(input_path=afrun.results_dir)
                if not ranked_list:
                    logging.info('No ranked PDBs found')
                    return
                plot_path = os.path.join(afrun.results_dir, 'plddt.png')
                _, best_ranked = output_air.plot_plddt(plot_path=plot_path, ranked_list=ranked_list)
                new_ranked_path = os.path.join(best_rankeds_dir,
                                               f'ranked_{afrun.start_chunk + 1}-{afrun.finish_chunk}.pdb')
                shutil.copy2(best_ranked.path, new_ranked_path)
                best_ranked.set_path(path=new_ranked_path)
                best_ranked_list.append(best_ranked)

            inf_path = best_ranked_list[0].path
            merge_pdbs_list = [inf_path]
            for i, ranked in enumerate(best_ranked_list[1:]):
                len_sequence = len(bioutils.extract_sequence(self.afrun_list[i].fasta_path))
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
                    lines = [line.replace('CA','').split() for line in lines]
                    for deltas in zip(lines, lines[1:], lines[2:], lines[3:]):
                        deltas_sum = sum([float(delta[0]) for delta in deltas])
                        if deltas_sum <= best_min:
                            best_list = deltas
                            best_min = deltas_sum

                if not best_list:
                    raise Exception('RMSD minimum requirements not met in order to merge the results in mosaic mode.')

                inf_cut = int(best_list[1][3])
                inm_cut = int(best_list[2][1])
                delete_residues = change_res.ChangeResidues(chain_res_dict={'A': [*range(inf_cut + 1, len_sequence + 1, 1)]})
                delete_residues.delete_residues(pdb_in_path=inf_path, pdb_out_path=inf_path)
                delete_residues = change_res.ChangeResidues(chain_res_dict={'A': [*range(1, inm_cut, 1)]})
                delete_residues.delete_residues(pdb_in_path=pdb_out, pdb_out_path=pdb_out)
                merge_pdbs_list.append(pdb_out)
                inf_path = pdb_out
            
            bioutils.merge_pdbs_in_one_chain(list_of_paths_of_pdbs_to_merge=merge_pdbs_list,
                                             pdb_out_path=os.path.join(results_dir, 'ranked_0.pdb'))
            self.run_dir = results_dir

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

    def write_input_file(self):
        with open(self.input_path, 'w') as f_out:
            f_out.write(f'output_dir: {self.output_dir}\n')
            f_out.write(f'run_dir: {self.run_dir}\n')
            f_out.write(f'af2_dbs_path: {self.af2_dbs_path}\n')
            f_out.write(f'verbose: {self.verbose}\n')
            f_out.write(f'glycines: {self.glycines}\n')
            f_out.write(f'run_af2: {self.run_af2}\n')
            if self.reference is not None:
                f_out.write(f'reference: {self.reference.pdb_path}\n')
            if self.experimental_pdb is not None:
                f_out.write(f'experimental_pdb: {self.experimental_pdb}\n')
            f_out.write(f'custom_features: {self.custom_features}\n')
            f_out.write(f'small_bfd: {self.small_bfd}\n')
            f_out.write(f'mosaic: {self.mosaic}\n')
            f_out.write(f'\nsequences:\n')
            for sequence_in in self.sequence_assembled.sequence_list:
                f_out.write('-')
                f_out.write(f' fasta_path: {sequence_in.fasta_path}\n')
                f_out.write(f'  num_of_copies: {sequence_in.num_of_copies}\n')
                new_positions = [position + 1 if position != -1 else position for position in sequence_in.positions]
                f_out.write(f'  positions: {",".join(map(str, new_positions))}\n')

            if self.templates_list:
                f_out.write(f'\ntemplates:\n')
                for template_in in self.templates_list:
                    f_out.write('-')
                    f_out.write(f' pdb: {template_in.pdb_path}\n')
                    f_out.write(f'  add_to_msa: {template_in.add_to_msa}\n')
                    f_out.write(f'  add_to_templates: {template_in.add_to_templates}\n')
                    f_out.write(f'  generate_multimer: {template_in.generate_multimer}\n')

                    if template_in.change_res_list:
                        f_out.write(f'  change_res:\n')
                        for change in template_in.change_res_list:
                            f_out.write('  -')
                            if change.resname is not None:
                                f_out.write(f' resname: {change.resname}\n')
                            elif change.sequence is not None:
                                f_out.write(f' fasta_path: {change.fasta_path}\n')
                            f_out.write(f'    when: {change.when}\n')
                            for key, value in change.chain_res_dict.items():
                                f_out.write(f'    {key}: {",".join(map(str, value))}\n')

                    if template_in.match_restrict_list:
                        f_out.write(f'  match:\n')
                        for match in template_in.match_restrict_list:
                            f_out.write('  -')
                            f_out.write(f' chain: {match.chain}\n')
                            if match.position != '':
                                f_out.write(f'    position: {match.position + 1}\n')
                            if match.residues is not None:
                                f_out.write(f'    residues: {",".join(map(str, match.residues))}\n')
                            if match.reference is not None:
                                f_out.write(f'    reference: {match.reference}\n')
                            if match.reference_chain is not None:
                                f_out.write(f'    reference_chain: {match.reference_chain}\n')

    def __repr__(self) -> str:
        return f' \
        output_dir: {self.output_dir} \n \
        run_af2: {self.run_af2}'
