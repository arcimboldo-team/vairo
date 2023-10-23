import logging
import os
import shutil
import sys
from itertools import combinations
from typing import Dict, List
import pandas as pd
from ALEPH.aleph.core import ALEPH
from libs import bioutils, features, utils, sequence, structures, plots

PERCENTAGE_FILTER = 0.8
PERCENTAGE_MAX_RMSD = 1


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
        self.templates_list: List[structures.TemplateExtracted] = []
        self.output_dir: str = output_dir
        self.experimental_dict = {}
        self.results_dir: str = ''
        self.templates_nonsplit_dir: str = ''
        self.rankeds_split_dir: str = ''
        self.rankeds_without_mutations_dir: str = ''
        self.tmp_dir: str = ''
        self.group_ranked_by_rmsd_dict: dict = {}
        self.template_interfaces: dict = {}
        self.templates_selected: List = []
        self.dendogram_struct: structures.Dendogram = None

        utils.create_dir(dir_path=self.plots_path, delete_if_exists=True)
        utils.create_dir(dir_path=self.templates_path, delete_if_exists=True)
        utils.create_dir(dir_path=self.interfaces_path, delete_if_exists=True)
        utils.create_dir(dir_path=self.frobenius_path, delete_if_exists=True)

    def extract_results(self, vairo_struct):

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
            template.add_percentage(vairo_struct.sequence_assembled.get_percentages(template.path))
            template.set_split_path(os.path.join(self.templates_path, f'{template.name}.pdb'))
            template.set_sequence_msa(list(bioutils.extract_sequence_msa_from_pdb(template.path).values())[0])
            template.set_identity(
                bioutils.sequence_identity(template.sequence_msa, vairo_struct.sequence_assembled.sequence_assembled))
            bioutils.split_chains_assembly(pdb_in_path=template.path,
                                           pdb_out_path=template.split_path,
                                           sequence_assembled=vairo_struct.sequence_assembled)
            _, compactness = bioutils.run_spong(pdb_in_path=template.path, spong_path=vairo_struct.binaries_paths.spong_path)
            template.set_compactness(compactness)
            _, perc = bioutils.generate_ramachandran(pdb_path=template.split_path)
            template.set_ramachandran(perc)
            template_struct = vairo_struct.get_template_by_id(template.name)
            template.set_template(template=template_struct, originalseq_path=os.path.join(self.templates_split_originalseq_dir, f'{template.name}.pdb'))

        logging.error('Reading predictions from the results folder')
        self.ranked_list = utils.read_rankeds(input_path=self.results_dir)
        self.select_templates()
        if not self.ranked_list:
            logging.error('No predictions found')
            return

        # Create a plot with the ranked pLDDTs, also, calculate the maximum pLDDT
        max_plddt = plots.plot_plddt(plot_path=self.plddt_plot_path, ranked_list=self.ranked_list)
        bioutils.write_sequence(sequence_name=utils.get_file_name(self.sequence_path),
                                sequence_amino=vairo_struct.sequence_assembled.sequence_assembled, sequence_path=self.sequence_path)

        # Copy the rankeds to the without mutations directory and remove the query sequences mutations from them
        for ranked in self.ranked_list:
            ranked.set_path(shutil.copy2(ranked.path, self.rankeds_nonsplit_dir))
            accepted_compactness, compactness = bioutils.run_spong(pdb_in_path=ranked.path,
                                                                   spong_path=vairo_struct.binaries_paths.spong_path)
            ranked.set_compactness(compactness)
            ranked.set_split_path(os.path.join(self.rankeds_split_dir, os.path.basename(ranked.path)))
            mapping = bioutils.split_chains_assembly(pdb_in_path=ranked.path,
                                                     pdb_out_path=ranked.split_path,
                                                     sequence_assembled=vairo_struct.sequence_assembled)
            accepted_ramachandran, perc = bioutils.generate_ramachandran(pdb_path=ranked.split_path,
                                                                         output_dir=self.plots_path)
            if perc is not None:
                perc = round(perc, 2)
            ranked.set_ramachandran(perc)
            ranked.set_mapping(mapping)
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
                if not accepted_ramachandran:
                    logging.error(f'    Ramachandran above limit')
                if not accepted_compactness:
                    logging.error(f'    Compactness below limit')
                if ranked.plddt < (PERCENTAGE_FILTER * max_plddt):
                    logging.error(f'    PLDDT too low')
        # Superpose the experimental pdb with all the rankeds and templates
        logging.error('Superposing experimental pdbs with predictions and templates')
        for experimental in vairo_struct.experimental_pdbs:
            aux_dict = {}
            for pdb in self.ranked_list + self.templates_list:
                rmsd, aligned_residues, quality_q = bioutils.gesamt_pdbs([pdb.split_path, experimental])
                if rmsd is not None:
                    rmsd = round(rmsd, 2)
                    total_residues = bioutils.get_number_residues(pdb.path)
                    aux_dict[pdb.name] = structures.PdbRanked(pdb.path, rmsd, aligned_residues,
                                                              total_residues, quality_q)
                    if pdb in self.ranked_list:
                        strct = structures.PdbRanked(experimental, rmsd, aligned_residues, total_residues, quality_q)
                        pdb.add_experimental(strct)
            self.experimental_dict[utils.get_file_name(experimental)] = aux_dict

            # Select the best ranked
        sorted_ranked_list = []
        if vairo_struct.experimental_pdbs:
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

    def analyse_output(self, sequence_assembled: sequence.SequenceAssembled, experimental_pdbs: List[str],
                       binaries_paths):

        if not self.ranked_list:
            return

        self.tmp_dir = f'{self.results_dir}/tmp'
        utils.create_dir(dir_path=self.tmp_dir, delete_if_exists=True)

        store_old_dir = os.getcwd()
        os.chdir(self.tmp_dir)

        reference_superpose = self.ranked_list[0].path

        # Store the superposition of the experimental with the best ranked
        for experimental in experimental_pdbs:
            bioutils.gesamt_pdbs([reference_superpose, experimental], experimental)

        # Superpose rankeds and store the superposition with the best one
        logging.error(f'Best prediction is {self.ranked_list[0].name}')
        logging.error('Superposing predictions and templates with the best prediction')
        results = [items for items in combinations(self.ranked_list, r=2)]
        for result in results:
            if result[0].name == self.ranked_list[0].name:
                rmsd, _, _ = bioutils.gesamt_pdbs([result[0].split_path, result[1].split_path], result[1].split_path)
            else:
                rmsd, _, _ = bioutils.gesamt_pdbs([result[0].split_path, result[1].split_path])

            result[0].set_ranked_to_rmsd_dict(rmsd=rmsd, ranked_name=result[1].name)
            result[1].set_ranked_to_rmsd_dict(rmsd=rmsd, ranked_name=result[0].name)

        # Group rankeds by how close they are between them
        for ranked in self.ranked_list:
            if ranked.filtered:
                found = False
                for ranked2 in self.ranked_list:
                    if ranked2.filtered and ranked2.name != ranked.name \
                            and ranked2.name in self.group_ranked_by_rmsd_dict \
                            and ranked.rmsd_dict.get(ranked2.name, float('inf')) <= PERCENTAGE_MAX_RMSD:
                        self.group_ranked_by_rmsd_dict[ranked2.name].append(ranked)
                        found = True
                        ranked.set_rmsd(ranked2.rmsd_dict[ranked.name])
                        if self.ranked_list[0].name == ranked2.name:
                            ranked.set_best(True)
                        break

                if not found:
                    self.group_ranked_by_rmsd_dict[ranked.name] = [ranked]
                    ranked.set_rmsd(0)
                    if self.ranked_list[0].name == ranked.name:
                        ranked.set_best(True)

        # Use frobenius
        templates_nonsplit_paths_list = [template.path for template in self.templates_list]
        dendogram_file = os.path.join(self.tmp_dir, 'dendogram.txt')
        dendogram_plot = os.path.join(self.tmp_dir, 'clustering_dendogram_angles.png')
        if len(self.templates_list) > 1:
            logging.error('Creating dendogram and clusters with ALEPH')
            with open(dendogram_file, 'w') as sys.stdout:
                _, _, _, _, _, _, _, _, dendogram_list = ALEPH.frobenius(references=templates_nonsplit_paths_list,
                                                                         targets=templates_nonsplit_paths_list,
                                                                         write_plot=True,
                                                                         write_matrix=True)
            sys.stdout = sys.__stdout__
            if dendogram_list:
                shutil.copy2(dendogram_plot, self.dendogram_plot_path)
                logging.error('The groups generated by frobenius are the following ones:')
            for i, templates in enumerate(dendogram_list):
                logging.error(f'Group {i}: {" ".join(templates)}')

            self.dendogram_struct = structures.Dendogram(dendogram_list=dendogram_list,
                                                         dendogram_plot=self.dendogram_plot_path,
                                                         encoded_dendogram_plot=utils.encode_data(
                                                             self.dendogram_plot_path))
        else:
            logging.error('Not possible to calculate the dendrogram with just one sample')

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
                    total_residues = bioutils.get_number_residues(template.path)
                    rmsd, aligned_residues, quality_q = bioutils.gesamt_pdbs([ranked.split_path, template.split_path])
                    if rmsd is not None:
                        rmsd = round(rmsd, 2)
                        ranked.add_template(
                            structures.PdbRanked(template.name, rmsd, aligned_residues, total_residues, quality_q))
                ranked.sort_template_rankeds()

        best_ranked_dict = get_best_ranked_by_template(templates_cluster_list, self.ranked_list)

        for template in self.templates_list:
            if best_ranked_dict and template.split_path in best_ranked_dict:
                bioutils.gesamt_pdbs([best_ranked_dict[template.split_path], template.split_path], template.split_path)
            else:
                bioutils.gesamt_pdbs([self.ranked_list[0].split_path, template.split_path], template.split_path)

        logging.error(
            'Analysing energies with openMM, interfaces with PISA and secondary structure information with ALEPH')
        if sequence_assembled.total_copies == 1:
            logging.error('Skipping interfaces generation. There is only one chain in the predictions')

        # Use aleph to generate domains and calculate secondary structure percentage
        for ranked in self.ranked_list:
            results_dict, domains_dict = bioutils.aleph_annotate(output_path=self.tmp_dir, pdb_path=ranked.split_path)
            if domains_dict is None or results_dict is None:
                break
            ranked.set_secondary_structure(ah=results_dict['ah'], bs=results_dict['bs'],
                                           total_residues=results_dict['number_total_residues'])
            if ranked.filtered:
                ranked.set_minimized_path(os.path.join(self.results_dir, f'{ranked.name}_minimized.pdb'))
                try:
                    ranked.set_potential_energy(
                        bioutils.run_openmm(pdb_in_path=ranked.path, pdb_out_path=ranked.minimized_path))
                except:
                    logging.debug(f'Not possible to calculate the energies for pdb {ranked.path}')

                ranked_chains_list = bioutils.get_chains(ranked.split_path)
                interfaces_data_list = bioutils.find_interface_from_pisa(ranked.split_path, self.interfaces_path)
                if interfaces_data_list:
                    deltas_list = [interface['deltaG'] for interface in interfaces_data_list]
                    deltas_list = utils.normalize_list([deltas_list])
                    for i, interface in enumerate(interfaces_data_list):
                        if not interface["chain1"] in domains_dict or not interface["chain2"] in domains_dict:
                            continue
                        seq1_name = sequence_assembled.get_sequence_name(ranked_chains_list.index(interface["chain1"]))
                        seq2_name = sequence_assembled.get_sequence_name(ranked_chains_list.index(interface["chain2"]))
                        code = f'{seq1_name}{interface["chain1"]}-{seq2_name}{interface["chain2"]}'
                        dimers_path = os.path.join(self.interfaces_path, f'{ranked.name}_{code}.pdb')
                        interface['bfactor'] = deltas_list[i]
                        if not ((float(interface['se_gain1']) < 0) and (float(interface['se_gain2']) < 0)):
                            interface['bfactor'] = abs(max(deltas_list)) * 2
                        extended_res_dict = bioutils.create_interface_domain(pdb_in_path=ranked.split_path,
                                                                             pdb_out_path=dimers_path,
                                                                             interface=interface,
                                                                             domains_dict=domains_dict)
                        renum_residues_list = []
                        renum_residues_list.extend(utils.renum_residues(extended_res_dict[interface['chain1']],
                                                                        mapping=ranked.mapping[interface['chain1']]))
                        renum_residues_list.extend(utils.renum_residues(extended_res_dict[interface['chain2']],
                                                                        mapping=ranked.mapping[interface['chain2']]))
                        ranked.add_interface(structures.Interface(name=code,
                                                                  res_list=renum_residues_list,
                                                                  chain1=interface["chain1"],
                                                                  chain2=interface["chain2"],
                                                                  se_gain1=float(interface['se_gain1']),
                                                                  se_gain2=float(interface['se_gain2']),
                                                                  solvation1=float(interface['solvation1']),
                                                                  solvation2=float(interface['solvation2'])
                                                                  ))

        for template in self.templates_list:
            frobenius_file = os.path.join(self.frobenius_path, f'frobenius_{template.name}.txt')
            matrices = os.path.join(self.tmp_dir, 'matrices')
            template_matrix = os.path.join(matrices, f'{utils.get_file_name(template.path)}_ang.npy')
            path_list = [ranked.path for ranked in self.ranked_list if ranked.filtered]
            with open(frobenius_file, 'w') as sys.stdout:
                _, list_targets, list_core, _, list_frobenius_angles, list_frobenius_distances, list_plot_ang, list_plot_dist, _ = ALEPH.frobenius(
                    references=[template.path],
                    targets=list(path_list), write_plot=True,
                    write_matrix=True)
            sys.stdout = sys.__stdout__
            ranked_filtered = [ranked for ranked in self.ranked_list if ranked.filtered]
            template_chains_list = bioutils.get_chains(template.split_path)
            template_interface_list = []
            if len(template_chains_list) > 1:
                interfaces_data_list = bioutils.find_interface_from_pisa(template.originalseq_path,
                                                                         self.interfaces_path)
                for interface in interfaces_data_list:
                    seq1_name = sequence_assembled.get_sequence_name(template_chains_list.index(interface["chain1"]))
                    seq2_name = sequence_assembled.get_sequence_name(template_chains_list.index(interface["chain2"]))
                    template_interface_list.append(f'{seq1_name}{interface["chain1"]}-{seq2_name}{interface["chain2"]}')
                self.template_interfaces[template.name] = template_interface_list
            else:
                logging.debug(f'Skipping interface search as there is only one chain in pdb {template.name}')

            for ranked in ranked_filtered:
                index = list_targets.index(ranked.name)
                ranked.add_frobenius_plot(
                    template=template.name,
                    dist_plot=[shutil.copy2(plot, self.frobenius_path) for plot in list_plot_dist if
                               ranked.name in os.path.basename(plot)][0],
                    ang_plot=[shutil.copy2(plot, self.frobenius_path) for plot in list_plot_ang if
                              ranked.name in os.path.basename(plot)][0],
                    dist_coverage=list_frobenius_distances[index],
                    ang_coverage=list_frobenius_angles[index],
                    core=list_core[index]
                )
                ranked_matrix = os.path.join(matrices, f'{ranked.name}_ang.npy')
                if not template_interface_list:
                    continue
                for interface in ranked.interfaces:
                    frobenius_file = os.path.join(self.frobenius_path, f'frobenius_{ranked.name}_{interface.name}.txt')
                    with open(frobenius_file, 'w') as sys.stdout:
                        fro_distance, fro_core, plot = ALEPH.frobenius_submatrices(path_ref=template_matrix,
                                                                                   path_tar=ranked_matrix,
                                                                                   residues_tar=interface.res_list,
                                                                                   write_plot=True,
                                                                                   title=f'Interface: {interface.name}')
                    sys.stdout = sys.__stdout__
                    new_name = os.path.join(self.frobenius_path, f'{template.name}_{ranked.name}_{interface.name}.png')
                    plot_path = os.path.join(self.tmp_dir, os.path.basename(plot))
                    interface.add_frobenius_information(
                        template=template.name,
                        dist_coverage=fro_distance,
                        core=fro_core,
                        dist_plot=shutil.copy2(plot_path, new_name)
                    )

        self.select_templates()
        os.chdir(store_old_dir)

    def select_templates(self):
        if len(self.templates_list) > 20:
            sorted_percentages = sorted(self.templates_list, key=lambda x: sum(x.percentage_list), reverse=True)[:20]
            self.templates_selected = [template.name for template in sorted_percentages]
            change_pos = -1
            for ranked in self.ranked_list:
                if ranked.superposition_templates:
                    if not ranked.superposition_templates[0].template in self.templates_selected:
                        self.templates_selected[change_pos] = ranked.superposition_templates[0].template
                        change_pos -= 1
        else:
            self.templates_selected = [template.name for template in self.templates_list]


    def write_tables(self, rmsd_dict: Dict, ranked_rmsd_dict: Dict, secondary_dict: Dict, plddt_dict: Dict,
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

            if bool(ranked_rmsd_dict):
                f_in.write('\n\n')
                f_in.write('Superposition between predictions\n')
                data = {'ranked': ranked_rmsd_dict.keys()}
                for ranked in ranked_rmsd_dict.values():
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
                        data.setdefault(key, []).append(
                            f'{value.rmsd} ({value.aligned_residues} of {value.total_residues}), {value.qscore}')
                df = pd.DataFrame(data)
                f_in.write(df.to_markdown())

            f_in.write('\n\n')
