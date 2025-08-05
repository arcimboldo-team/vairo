import logging
import os
import shutil
import sys
from itertools import combinations
from typing import Dict, List
import pandas as pd
from ALEPH.aleph.core import ALEPH
from libs import bioutils, utils, sequence, structures, plots

PERCENTAGE_FILTER = 0.8
QSCORE_MINIMUM = 0.3


def get_best_ranked_by_template(cluster_list: List, ranked_list: List) -> Dict:
    return_dict = {}
    for i, cluster in enumerate(cluster_list):
        ranked = next((ranked.split_path for ranked in ranked_list if f'cluster_{i}' in ranked.name), None)
        if ranked is None:
            ranked = ranked_list[0].split_path
        return_dict.update(dict.fromkeys(cluster, ranked))
    return return_dict


class OutputStructure:

    def __init__(self, output_dir: str):
        self.plots_path: str = f'{output_dir}/plots'
        self.frobenius_path: str = f'{output_dir}/frobenius'
        self.sequence_path: str = os.path.join(self.frobenius_path, 'sequence.fasta')
        self.templates_path: str = f'{output_dir}/templates'
        self.interfaces_path: str = f'{output_dir}/interfaces'
        self.analysis_path: str = f'{self.plots_path}/analysis.txt'
        self.plddt_plot_path: str = f'{self.plots_path}/plddt.png'
        self.sequence_plot_path: str = f'{self.plots_path}/sequence_plot.png'
        self.analysis_plot_path: str = f'{self.plots_path}/cc_analysis_plot.png'
        self.dendogram_plot_path: str = f'{self.plots_path}/frobenius_dendogram.png'
        self.analysis_ranked_plot_path: str = f'{self.plots_path}/cc_analysis_ranked_plot.png'
        self.html_path: str = f'{output_dir}/output.html'
        self.html_complete_path: str = f'{output_dir}/output_complete.html'
        self.pymol_session_path: str = f'{output_dir}/pymol_session.pse'
        self.gantt_plots: structures.GanttPlot = None
        self.gantt_complete_plots: structures.GanttPlot = None
        self.ranked_list: List[structures.Ranked] = []
        self.ranked_filtered_list: List[structures.Ranked] = []
        self.templates_list: List[structures.TemplateExtracted] = []
        self.output_dir: str = output_dir
        self.experimental_dict = {}
        self.experimental_list: List[structures.ExperimentalPdb] = []
        self.best_experimental: str = None
        self.results_dir: str = ''
        self.templates_nonsplit_dir: str = ''
        self.rankeds_split_dir: str = ''
        self.rankeds_without_mutations_dir: str = ''
        self.tmp_dir: str = ''
        self.conservation_ranked_path: str = ''
        self.group_ranked_by_qscore_dict: dict = {}
        self.templates_selected: List = []
        self.dendogram_struct: structures.Dendogram = None
        self.num_interfaces: int = 0
        self.filtered_ranked_reason_dict: dict = {}

        utils.create_dir(dir_path=self.plots_path, delete_if_exists=True)
        utils.create_dir(dir_path=self.templates_path, delete_if_exists=True)
        utils.create_dir(dir_path=self.interfaces_path, delete_if_exists=True)
        utils.create_dir(dir_path=self.frobenius_path, delete_if_exists=True)

    def extract_results(self, vairo_struct, region_predicted):
        if region_predicted:
            assembled_seq = vairo_struct.sequence_predicted_assembled
        else:
            assembled_seq = vairo_struct.sequence_assembled

        # Read all templates and rankeds, if there are no ranked, raise an error
        self.results_dir = vairo_struct.results_dir
        self.templates_nonsplit_dir = f'{self.results_dir}/templates_nonsplit'
        self.templates_split_originalseq_dir = f'{self.results_dir}/templates_split_originalseq'
        self.rankeds_split_dir = f'{self.results_dir}/rankeds_split'
        self.rankeds_nonsplit_dir = f'{self.results_dir}/rankeds_nonsplit'

        utils.create_dir(dir_path=self.templates_nonsplit_dir, delete_if_exists=True)
        utils.create_dir(dir_path=self.templates_split_originalseq_dir, delete_if_exists=True)
        utils.create_dir(dir_path=self.rankeds_split_dir, delete_if_exists=True)
        utils.create_dir(dir_path=self.rankeds_nonsplit_dir, delete_if_exists=True)

        logging.error('Extracting templates from the features file')
        if vairo_struct.feature is not None:
            templates_nonsplit_dict = vairo_struct.feature.write_all_templates_in_features(
                output_dir=self.templates_nonsplit_dir,
                print_number=False)
            self.templates_list = [structures.TemplateExtracted(path=path) for path in templates_nonsplit_dict.values()]

        # Split the templates with chains
        for template in self.templates_list:
            template.add_percentage(assembled_seq.get_percentages(template.path))
            if sum(template.percentage_list) == 0:
                logging.info('Template {template.name} does not have any sequence coverage. Skipping')
                continue
            template.set_split_path(os.path.join(self.templates_path, f'{template.name}.pdb'))
            template.set_sequence_msa(list(bioutils.extract_sequence_msa_from_pdb(template.path).values())[0])
            template.set_identity(
                bioutils.sequence_identity(template.sequence_msa, assembled_seq.sequence_assembled))
            bioutils.split_chains_assembly(pdb_in_path=template.path,
                                           pdb_out_path=template.split_path,
                                           sequence_assembled=assembled_seq)
            _, compactness = bioutils.run_spong(pdb_in_path=template.path,
                                                spong_path=vairo_struct.binaries_paths.spong_path)
            template.set_compactness(compactness)
            _, perc = bioutils.generate_ramachandran(pdb_path=template.split_path)
            template.set_ramachandran(perc)
            template_struct = vairo_struct.get_template_by_id(template.name)
            template.set_template(template=template_struct,
                                  originalseq_path=os.path.join(self.templates_split_originalseq_dir,
                                                                f'{template.name}.pdb'))

        # Delete templates with percentage to 0
        self.templates_list = [template for template in self.templates_list if sum(template.percentage_list) != 0]

        logging.error('Reading predictions from the results folder')
        self.ranked_list = utils.read_rankeds(input_path=self.results_dir)
        self.select_templates()
        if not self.ranked_list:
            logging.error('No predictions found')
            return

        # Create a plot with the ranked pLDDTs, also, calculate the maximum pLDDT
        bioutils.write_sequence(sequence_name=utils.get_file_name(self.sequence_path),
                                sequence_amino=assembled_seq.sequence_assembled,
                                sequence_path=self.sequence_path)


        # Copy the rankeds to the without mutations directory and remove the query sequences mutations from them
        for ranked in self.ranked_list:
            ranked.set_path(shutil.copy2(ranked.path, self.rankeds_nonsplit_dir))

            if not region_predicted:
                bioutils.shift_pdb(pdb_in_path=ranked.path, sequence_predicted_assembled=vairo_struct.sequence_predicted_assembled, 
                                   sequence_assembled=vairo_struct.sequence_assembled)


            ranked.set_split_path(os.path.join(self.rankeds_split_dir, os.path.basename(ranked.path)))
            bioutils.split_chains_assembly(pdb_in_path=ranked.path,
                                           pdb_out_path=ranked.split_path,
                                           sequence_assembled=assembled_seq)
            ranked.set_plddt()

        plots.plot_plddt(plot_path=self.plddt_plot_path, ranked_list=self.ranked_list)
        max_plddt = max([ranked.plddt for ranked in self.ranked_list])

        for ranked in self.ranked_list:
            accepted_ramachandran, perc = bioutils.generate_ramachandran(pdb_path=ranked.split_path,
                                                                         output_dir=self.plots_path)
            if perc is not None:
                perc = round(perc, 2)
            ranked.set_ramachandran(perc)
            accepted_compactness, compactness = bioutils.run_spong(pdb_in_path=ranked.path,
                                                                   spong_path=vairo_struct.binaries_paths.spong_path)
            ranked.set_compactness(compactness)
            bioutils.remove_hydrogens(ranked.split_path, ranked.split_path)
            ranked.set_encoded(ranked.split_path)
            if accepted_ramachandran and accepted_compactness and ranked.plddt >= (PERCENTAGE_FILTER * max_plddt):
                ranked.set_filtered(True)
                logging.error(f'Prediction {ranked.name} has been accepted')
                ranked.set_split_path(
                    shutil.copy2(ranked.split_path, os.path.join(self.output_dir, os.path.basename(ranked.path))))
            else:
                ranked.set_filtered(False)
                logging.error(f'Prediction {ranked.name} has been filtered:')
                conditions = [
                    (not accepted_ramachandran, 'Ramachandran above limit'),
                    (not accepted_compactness, 'Compactness below limit'),
                    (ranked.plddt < (PERCENTAGE_FILTER * max_plddt), 'PLDDT too low')
                ]
                error_str = [msg for condition, msg in conditions if condition]
                error_str = ', '.join(error_str)
                self.filtered_ranked_reason_dict[ranked.name] = error_str
                logging.error(f'    {error_str}')

        # Superpose the experimental pdb with all the rankeds and templates
        logging.error('Superposing experimental pdbs with predictions and templates')
        for experimental in vairo_struct.experimental_pdbs:
            self.experimental_list.append(structures.ExperimentalPdb(path=experimental))
            aux_dict = {}
            for pdb in self.ranked_list + self.templates_list:
                rmsd, aligned_residues, quality_q = bioutils.gesamt_pdbs(pdb_reference=pdb.split_path,
                                                                         pdb_superposed=experimental)
                if rmsd is not None:
                    rmsd = round(rmsd, 2)
                    total_residues = bioutils.get_number_residues(pdb.split_path)
                    aux_dict[pdb.name] = structures.PdbRanked(pdb.path, rmsd, aligned_residues,
                                                              total_residues, quality_q)
                    if pdb in self.ranked_list:
                        strct = structures.PdbRanked(experimental, rmsd, aligned_residues, total_residues, quality_q)
                        pdb.add_experimental(strct)
                else:
                    aux_dict[pdb.name] = structures.PdbRanked(pdb.path, None, None, None, None)
            self.experimental_dict[utils.get_file_name(experimental)] = aux_dict

        # Select the best ranked
        if vairo_struct.experimental_pdbs:
            [pdb.sort_experimental_rankeds() for pdb in self.ranked_list]
            logging.error(
                'Experimental pdbs found. Selecting the best prediction taking into account the qscore with the experimental pdbs')
            sorted_ranked_list = sorted(self.ranked_list, key=lambda ranked: (
                ranked.filtered, ranked.superposition_experimental[0].qscore), reverse=True)
        else:
            logging.error('No experimental pdbs found. Selecting best prediction by PLDDT')
            sorted_ranked_list = sorted(self.ranked_list, key=lambda ranked: (ranked.filtered, ranked.plddt),
                                        reverse=True)
        if not sorted_ranked_list:
            self.ranked_list.sort(key=lambda x: x.plddt, reverse=True)
            logging.error(
                'There are no predictions that meet the minimum quality requirements. All predictions were filtered. Check the tables')
        else:
            self.ranked_list = sorted_ranked_list

        if vairo_struct.feature is not None:
            self.conservation_ranked_path = os.path.join(self.results_dir, f'{ranked.name}_conservation.pdb')
            bioutils.conservation_pdb(self.ranked_list[0].path, self.conservation_ranked_path,
                                      vairo_struct.feature.get_msa_sequences())
            bioutils.split_chains_assembly(pdb_in_path=self.conservation_ranked_path,
                                           pdb_out_path=self.conservation_ranked_path,
                                           sequence_assembled=assembled_seq)

    def analyse_output(self, sequence_assembled: sequence.SequenceAssembled, experimental_pdbs: List[str],
                       binaries_paths):

        if not self.ranked_list:
            return

        self.tmp_dir = os.path.join(self.results_dir, 'temp')
        utils.create_dir(dir_path=self.tmp_dir, delete_if_exists=True)

        store_old_dir = os.getcwd()
        os.chdir(self.tmp_dir)

        reference_superpose = self.ranked_list[0].path
        if self.ranked_list[0].superposition_experimental:
            self.best_experimental = utils.get_file_name(self.ranked_list[0].superposition_experimental[0].pdb)

        # Store the superposition of the experimental with the best ranked
        for experimental in experimental_pdbs:
            bioutils.gesamt_pdbs(pdb_reference=reference_superpose, pdb_superposed=experimental,
                                 output_path=experimental)

        # Superpose rankeds and store the superposition with the best one
        logging.error(f'Best prediction is {self.ranked_list[0].name}')
        logging.error('Superposing predictions and templates with the best prediction')
        results = [items for items in combinations(self.ranked_list, r=2)]
        for result in results:
            if result[0].name == self.ranked_list[0].name:
                rmsd, _, qscore = bioutils.gesamt_pdbs(pdb_reference=result[0].split_path,
                                                       pdb_superposed=result[1].split_path,
                                                       output_path=result[1].split_path)
            else:
                rmsd, _, qscore = bioutils.gesamt_pdbs(pdb_reference=result[0].split_path,
                                                       pdb_superposed=result[1].split_path)

            result[0].set_ranked_to_qscore_dict(qscore=qscore, ranked_name=result[1].name)
            result[1].set_ranked_to_qscore_dict(qscore=qscore, ranked_name=result[0].name)

        # Group rankeds by how close they are between them
        for ranked in self.ranked_list:
            if ranked.filtered:
                ranked.set_minimized_path(os.path.join(self.results_dir, f'{ranked.name}_minimized.pdb'))
                try:
                    ranked.set_potential_energy(
                        bioutils.run_openmm(pdb_in_path=ranked.split_path, pdb_out_path=ranked.minimized_path))
                except:
                    logging.info(f'Not possible to calculate the energies for pdb {ranked.path}')
                found = False
                for ranked2 in self.ranked_list:
                    if ranked2.filtered and ranked2.name != ranked.name \
                            and ranked2.name in self.group_ranked_by_qscore_dict \
                            and ranked.qscore_dict.get(ranked2.name, float('inf')) >= QSCORE_MINIMUM:
                        self.group_ranked_by_qscore_dict[ranked2.name].append(ranked)
                        found = True
                        ranked.set_qscore(ranked2.qscore_dict[ranked.name])
                        if self.ranked_list[0].name == ranked2.name:
                            ranked.set_best(True)
                        break

                if not found:
                    self.group_ranked_by_qscore_dict[ranked.name] = [ranked]
                    ranked.set_qscore(0)
                    if self.ranked_list[0].name == ranked.name:
                        ranked.set_best(True)

        self.ranked_filtered_list = [ranked for ranked in self.ranked_list if ranked.filtered]

        # Use frobenius
        #templates_nonsplit_paths_list = [template.path for template in self.templates_list]
        #dendogram_file = os.path.join(self.tmp_dir, 'dendogram.txt')
        #dendogram_plot = os.path.join(self.tmp_dir, 'clustering_dendogram_angles.png')
        #if 1 < len(self.templates_list) <= 15:
        #    logging.error('Creating dendogram and clusters with ALEPH')
        #    with open(dendogram_file, 'w') as sys.stdout:
        #        _, _, _, _, _, _, _, _, dendogram_list = ALEPH.frobenius(references=templates_nonsplit_paths_list,
        #                                                                 targets=templates_nonsplit_paths_list,
        #                                                                 write_plot=True,
        #                                                                 write_matrix=True)
        #    sys.stdout = sys.__stdout__
        #    if dendogram_list:
        #        shutil.copy2(dendogram_plot, self.dendogram_plot_path)
        #        logging.error('The groups generated by frobenius are the following ones:')
        #        for i, templates in enumerate(dendogram_list):
        #            logging.error(f'    Group {i}: {" ".join(templates)}')

        #        self.dendogram_struct = structures.Dendogram(dendogram_list=dendogram_list,
        #                                                     dendogram_plot=self.dendogram_plot_path,
        #                                                    encoded_dendogram_plot=utils.encode_data(
        #                                                         self.dendogram_plot_path))
        #else:
        #    logging.error('Not possible to calculate the dendrogram with just one sample')

        # Generate CCANALYSIS plots, one without rankeds and another one with rankeds.

        logging.error('Analysing results with hinges and ccanalysis')
        templates_cluster_list, analysis_dict = bioutils.cc_and_hinges_analysis(pdbs=self.templates_list,
                                                                                binaries_path=binaries_paths,
                                                                                output_dir=self.results_dir)
        if analysis_dict:
            plots.plot_cc_analysis(plot_path=self.analysis_plot_path,
                                   analysis_dict=analysis_dict,
                                   clusters=templates_cluster_list)

        templates_cluster_ranked_list, analysis_dict_ranked = bioutils.cc_and_hinges_analysis(
            pdbs=self.templates_list + self.ranked_list,
            binaries_path=binaries_paths,
            output_dir=self.results_dir)

        if analysis_dict_ranked:
            plots.plot_cc_analysis(plot_path=self.analysis_ranked_plot_path, analysis_dict=analysis_dict_ranked,
                                   clusters=templates_cluster_ranked_list, predictions=True)

        # Superpose each template with all the rankeds.
        if self.templates_list:
            for i, ranked in enumerate(self.ranked_list):
                for template in self.templates_list:
                    total_residues = bioutils.get_number_residues(template.split_path)
                    rmsd, aligned_residues, quality_q = bioutils.gesamt_pdbs(pdb_reference=ranked.split_path,
                                                                             pdb_superposed=template.split_path)
                    if rmsd is not None:
                        rmsd = round(rmsd, 2)
                    ranked.add_template(
                        structures.PdbRanked(template.name, rmsd, aligned_residues, total_residues, quality_q))

                ranked.sort_template_rankeds()

        best_ranked_dict = get_best_ranked_by_template(templates_cluster_list, self.ranked_list)

        for template in self.templates_list:
            if best_ranked_dict and template.split_path in best_ranked_dict:
                bioutils.gesamt_pdbs(pdb_reference=best_ranked_dict[template.split_path],
                                     pdb_superposed=template.split_path, output_path=template.split_path)
            else:
                bioutils.gesamt_pdbs(pdb_reference=self.ranked_list[0].split_path, pdb_superposed=template.split_path,
                                     output_path=template.split_path)

        logging.error(
            'Analysing energies with openMM, interfaces with PISA and secondary structure information with ALEPH')
        if sequence_assembled.total_copies == 1:
            logging.error('Skipping interfaces generation. There is only one chain in the predictions')

        # Use aleph to generate domains and calculate secondary structure percentage
        for pdb_in in self.ranked_list + self.experimental_list + self.templates_list:
            results_dict, domains_dict = bioutils.aleph_annotate(output_path=self.tmp_dir, pdb_path=pdb_in.split_path)
            if results_dict is not None:
                pdb_in.set_secondary_structure(ah=results_dict['ah'], bs=results_dict['bs'],
                                               total_residues=results_dict['number_total_residues'])
            else:
                pdb_in.set_secondary_structure(ah=None, bs=None, total_residues=None)

            if sequence_assembled.total_copies > 1:
                if pdb_in in self.experimental_list + self.templates_list + self.ranked_filtered_list:
                    if isinstance(pdb_in, structures.TemplateExtracted):
                        interfaces_data_list = bioutils.find_interface_from_pisa(pdb_in.originalseq_path,
                                                                                 self.interfaces_path)
                    else:
                        interfaces_data_list = bioutils.find_interface_from_pisa(pdb_in.split_path,
                                                                                 self.interfaces_path)
                    if interfaces_data_list and domains_dict is not None:
                        for i, interface in enumerate(interfaces_data_list):
                            if not interface.chain1 in domains_dict or not interface.chain2 in domains_dict:
                                continue
                            code = f'{interface.chain1}-{interface.chain2}'
                            dimers_path = os.path.join(self.interfaces_path, f'{pdb_in.name}_{code}.pdb')
                            bioutils.create_interface_domain(pdb_in_path=pdb_in.split_path,
                                                             pdb_out_path=dimers_path,
                                                             interface=interface,
                                                             domains_dict=domains_dict)
                            interface.set_structure(dimers_path)
                        pdb_in.set_interfaces(interfaces_data_list)

        if sequence_assembled.total_copies > 1:
            if self.best_experimental is not None:
                exp = [experimental for experimental in self.experimental_list if
                       experimental.name == self.best_experimental][0]
                self.num_interfaces = len(exp.interfaces)
                exp.set_accepted_interfaces(True)
            else:
                self.num_interfaces = len(self.ranked_list[0].interfaces)

            self.ranked_list[0].set_accepted_interfaces(True)
            for ranked in self.ranked_list[1:]:
                if len(ranked.interfaces) >= self.num_interfaces:
                    ranked.set_accepted_interfaces(True)

        self.select_templates()
        os.chdir(store_old_dir)

    def select_templates(self):
        if len(self.templates_list) > 20:
            sorted_percentages = sorted(self.templates_list, key=lambda x: sum(x.percentage_list), reverse=True)[:20]
            self.templates_selected = [template.name for template in sorted_percentages]
            change_pos = -1
            for ranked in self.ranked_list:
                if ranked.superposition_templates:
                    if not ranked.superposition_templates[0].pdb in self.templates_selected:
                        self.templates_selected[change_pos] = ranked.superposition_templates[0].pdb
                        change_pos -= 1
        else:
            self.templates_selected = [template.name for template in self.templates_list]

    def write_tables(self, rmsd_dict: Dict, ranked_qscore_dict: Dict, secondary_dict: Dict, plddt_dict: Dict,
                     energies_dict: Dict):
        with open(self.analysis_path, 'w') as f_in:
            if bool(rmsd_dict):
                f_in.write('\n\n')
                f_in.write('Superpositions of rankeds and templates\n')
                data = {'ranked': rmsd_dict.keys()}
                for templates in rmsd_dict.values():
                    for template, value in templates.items():
                        data.setdefault(template, []).append(
                            f'{value["rmsd"]} {value["aligned_residues"]} ({value["total_residues"]})')
                df = pd.DataFrame(data)
                f_in.write(df.to_markdown())

            if bool(ranked_qscore_dict):
                f_in.write('\n\n')
                f_in.write('Superposition between predictions (QSCORE)\n')
                data = {'ranked': ranked_qscore_dict.keys()}
                for ranked in ranked_qscore_dict.values():
                    for key, value in ranked.items():
                        data.setdefault(key, []).append(value)
                df = pd.DataFrame(data)
                f_in.write(df.to_markdown())

            if bool(secondary_dict):
                f_in.write('\n\n')
                f_in.write('Secondary structure percentages calculated with ALEPH\n')
                data = {'ranked': secondary_dict.keys(),
                        'ah': [value['ah'] for value in secondary_dict.values()],
                        'bs': [value['bs'] for value in secondary_dict.values()],
                        'number_total_residues': [value['number_total_residues'] for value in secondary_dict.values()]
                        }
                df = pd.DataFrame(data)
                f_in.write(df.to_markdown())

            if bool(plddt_dict):
                f_in.write('\n\n')
                f_in.write('ranked information \n')
                data = {'ranked': plddt_dict.keys(),
                        'plddt': [value['plddt'] for value in plddt_dict.values()],
                        'compactness': [value['compactness'] for value in plddt_dict.values()],
                        'ramachandran': [value['ramachandran'] for value in plddt_dict.values()]
                        }
                df = pd.DataFrame(data)
                f_in.write(df.to_markdown())

            if bool(energies_dict):
                f_in.write('\n\n')
                f_in.write('OPENMM Energies\n')
                data = {'ranked': energies_dict.keys(),
                        'potential': energies_dict.values()
                        }
                df = pd.DataFrame(data)
                f_in.write(df.to_markdown())

            if bool(self.experimental_dict):
                f_in.write('\n\n')
                f_in.write(f'Superposition with experimental structures\n')
                data = {'experimental': self.experimental_dict.keys()}
                for keys_pdbs in self.experimental_dict.values():
                    for key, value in keys_pdbs.items():
                        if value.rmsd is not None:
                            data.setdefault(key, []).append(
                                f'{value.rmsd} ({value.aligned_residues} of {value.total_residues}), {value.qscore}')
                        else:
                            data.setdefault(key, []).append('None')
                df = pd.DataFrame(data)
                f_in.write(df.to_markdown())

            f_in.write('\n\n')
