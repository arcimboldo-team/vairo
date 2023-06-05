import logging
import os
import shutil
import sys
from itertools import combinations
from typing import Dict, List
import pandas as pd
from Bio.PDB import PDBParser, Selection

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


class OutputAir:

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
        self.analysis_ranked_plot_path: str = f'{self.plots_path}/cc_analysis_ranked_plot.png'
        self.html_path: str = f'{output_dir}/output.html'
        self.gantt_plots: structures.GanttPlot = None
        self.ranked_list: List[structures.Ranked] = []
        self.output_dir: str = output_dir
        self.experimental_dict = {}
        self.results_dir: str = ''
        self.templates_nonsplit_dir: str = ''
        self.rankeds_nonsplit_dir: str = ''
        self.rankeds_split_dir: str = ''
        self.tmp_dir: str = ''
        self.group_ranked_by_rmsd_dict: dict = {}
        self.template_interfaces: dict = {}
        self.templates_dict = {}
        self.templates_nonsplit_dict = {}
        self.percentage_sequences = {}

        utils.create_dir(dir_path=self.plots_path, delete_if_exists=True)
        utils.create_dir(dir_path=self.templates_path, delete_if_exists=True)
        utils.create_dir(dir_path=self.interfaces_path, delete_if_exists=True)
        utils.create_dir(dir_path=self.frobenius_path, delete_if_exists=True)

    def create_plot_gantt(self, a_air):
        gantt_plots_both, legend_both = plots.plot_gantt(plot_type='both', plot_path=self.plots_path,
                                                         a_air=a_air)
        gantt_plots_template, legend_template = plots.plot_gantt(plot_type='templates', plot_path=self.plots_path,
                                                                 a_air=a_air)
        gantt_plots_msa, legend_msa = plots.plot_gantt(plot_type='msa', plot_path=self.plots_path, a_air=a_air)
        self.gantt_plots = structures.GanttPlot(plot_both=utils.encode_data(gantt_plots_both),
                                                legend_both=legend_both,
                                                plot_template=utils.encode_data(gantt_plots_template),
                                                legend_template=legend_template,
                                                plot_msa=utils.encode_data(gantt_plots_msa),
                                                legend_msa=legend_msa)
        plots.plot_sequence(plot_path=self.sequence_plot_path, a_air=a_air)

    def analyse_output(self, results_dir: str, sequence_assembled: sequence.SequenceAssembled,
                       feature: features.Features, experimental_pdbs: List[str], cc_analysis_paths,
                       cluster_templates: bool = False):
        # Read all templates and rankeds, if there are no ranked, raise an error

        self.results_dir = results_dir
        self.templates_nonsplit_dir = f'{self.results_dir}/templates_nonsplit'
        self.rankeds_split_dir = f'{self.results_dir}/rankeds_split'
        self.tmp_dir = f'{self.results_dir}/tmp'

        utils.create_dir(dir_path=self.templates_nonsplit_dir, delete_if_exists=True)
        utils.create_dir(dir_path=self.rankeds_split_dir, delete_if_exists=True)
        utils.create_dir(dir_path=self.tmp_dir, delete_if_exists=True)

        store_old_dir = os.getcwd()
        os.chdir(self.tmp_dir)

        if feature is not None:
            self.templates_nonsplit_dict = feature.write_all_templates_in_features(
                output_dir=self.templates_nonsplit_dir,
                print_number=False)
        # Split the templates with chains
        for template, template_path in self.templates_nonsplit_dict.items():
            self.percentage_sequences[template] = sequence_assembled.get_percentages(template_path)
            new_pdb_path = os.path.join(self.templates_path, f'{template}.pdb')
            self.templates_dict[template] = new_pdb_path
            shutil.copy2(template_path, new_pdb_path)
            bioutils.split_chains_assembly(pdb_in_path=new_pdb_path,
                                           pdb_out_path=new_pdb_path,
                                           sequence_assembled=sequence_assembled)

        self.ranked_list = utils.read_rankeds(input_path=self.results_dir)

        # dendogram_file = os.path.join(self.run_dir, 'dendogram.txt')
        # dendogram_plot = os.path.join(self.run_dir, 'clustering_dendogram_angles.png')
        # with open(dendogram_file, 'w') as sys.stdout:
        #    _, _, _, _, _, _, _, _, dendogram_list = ALEPH.frobenius(references=list(  .values()),
        #                                            targets=list(template_nonsplit.values()), write_plot=True,
        #                                            write_matrix=True)
        # sys.stdout = sys.__stdout__
        # if dendogram_list:
        #    shutil.copy2(dendogram_plot, self.template_dendogram)
        #    if not custom_features:
        #        for templates in dendogram_list:
        #            self.templates_cluster.append([template_nonsplit[template] for template in templates])

        if not self.ranked_list or cluster_templates:
            logging.info('No ranked PDBs found')
            os.chdir(store_old_dir)
            return

        # Create a plot with the ranked pLDDTs, also, calculate the maximum pLDDT
        max_plddt = plots.plot_plddt(plot_path=self.plddt_plot_path, ranked_list=self.ranked_list)
        bioutils.write_sequence(sequence_name=utils.get_file_name(self.sequence_path),
                                sequence_amino=sequence_assembled.sequence_assembled, sequence_path=self.sequence_path)

        # Save split path of all rankeds, taking into account the split dir
        [ranked.set_split_path(os.path.join(self.rankeds_split_dir, os.path.basename(ranked.path))) for ranked in
         self.ranked_list]

        # Sort list of ranked by pLDDT
        self.ranked_list.sort(key=lambda x: x.plddt, reverse=True)
        shutil.copy2(self.ranked_list[0].path, self.ranked_list[0].split_path)
        results = [items for items in combinations(self.ranked_list, r=2)]
        for result in results:
            if result[0].name == self.ranked_list[0].name:
                rmsd, _, _ = bioutils.gesamt_pdbs([result[0].path, result[1].path], result[1].split_path)

            else:
                rmsd, _, _ = bioutils.gesamt_pdbs([result[0].path, result[1].path])
            if not os.path.exists(result[1].split_path):
                shutil.copy2(result[1].path, result[1].split_path)
            result[0].set_ranked_to_rmsd_dict(rmsd=rmsd, ranked_name=result[1].name)
            result[1].set_ranked_to_rmsd_dict(rmsd=rmsd, ranked_name=result[0].name)

        reference_superpose = self.ranked_list[0].path
        green_color = 55

        # Filter rankeds, split them in chains.
        for ranked in self.ranked_list:
            if ranked.plddt >= (PERCENTAGE_FILTER * max_plddt):
                ranked.set_filtered(True)
                ranked.set_split_path(
                    shutil.copy2(ranked.split_path, os.path.join(self.output_dir, os.path.basename(ranked.path))))

            mapping = bioutils.split_chains_assembly(pdb_in_path=ranked.split_path,
                                                     pdb_out_path=ranked.split_path,
                                                     sequence_assembled=sequence_assembled)
            ranked.set_mapping(mapping)
            ranked.set_encoded(ranked.split_path)
            bioutils.remove_hydrogens(ranked.split_path, ranked.split_path)

        for ranked in self.ranked_list:
            if ranked.filtered:
                found = False
                for ranked2 in self.ranked_list:
                    if ranked2.filtered and ranked2.name != ranked.name \
                            and ranked2.name in self.group_ranked_by_rmsd_dict \
                            and ranked.rmsd_dict[ranked2.name] is not None \
                            and ranked.rmsd_dict[ranked2.name] <= PERCENTAGE_MAX_RMSD:
                        self.group_ranked_by_rmsd_dict[ranked2.name].append(ranked)
                        found = True
                        ranked.set_rmsd(ranked2.rmsd_dict[ranked.name])
                        if self.ranked_list[0].name == ranked2.name:
                            ranked.set_green_color(green_color)
                            ranked.set_best(True)
                            green_color -= 5
                        break

                if not found:
                    self.group_ranked_by_rmsd_dict[ranked.name] = [ranked]
                    ranked.set_rmsd(0)
                    if self.ranked_list[0].name == ranked.name:
                        ranked.set_green_color(green_color)
                        ranked.set_best(True)
                        green_color -= 5

        # Generate CCANALYSIS plots, one without rankeds and another one with rankeds.

        templates_cluster_list, analysis_dict = bioutils.cc_and_hinges_analysis(paths_in=self.templates_dict,
                                                                                binaries_path=cc_analysis_paths,
                                                                                output_path=self.results_dir,
                                                                                length_sequences=self.percentage_sequences)
        if analysis_dict:
            plots.plot_cc_analysis(plot_path=self.analysis_plot_path,
                                   analysis_dict=analysis_dict,
                                   clusters=templates_cluster_list)
        aux_dict = dict({ranked.name: ranked.split_path for ranked in self.ranked_list}, **self.templates_dict)

        templates_cluster_ranked_list, analysis_dict_ranked = bioutils.cc_and_hinges_analysis(paths_in=aux_dict,
                                                                                              binaries_path=cc_analysis_paths,
                                                                                              output_path=self.results_dir,
                                                                                              length_sequences=self.percentage_sequences)

        if analysis_dict_ranked:
            plots.plot_cc_analysis(plot_path=self.analysis_ranked_plot_path, analysis_dict=analysis_dict_ranked,
                                   clusters=templates_cluster_ranked_list, predictions=True)

        # Superpose each template with all the rankeds.
        if self.templates_dict:
            for i, ranked in enumerate(self.ranked_list):
                for template, template_path in self.templates_dict.items():
                    total_residues = bioutils.get_number_residues(template_path)
                    rmsd, aligned_residues, quality_q = bioutils.gesamt_pdbs([ranked.split_path, template_path])
                    rmsd = round(rmsd, 2) if rmsd is not None else rmsd
                    ranked.add_template(structures.TemplateRanked(template, rmsd, aligned_residues, total_residues))
                ranked.sort_template_rankeds()

        best_ranked_dict = get_best_ranked_by_template(templates_cluster_list, self.ranked_list)

        for template_path in self.templates_dict.values():
            if best_ranked_dict and template_path in best_ranked_dict:
                bioutils.gesamt_pdbs([best_ranked_dict[template_path], template_path], template_path)
            else:
                bioutils.gesamt_pdbs([self.ranked_list[0].split_path, template_path], template_path)

        # Use aleph to generate domains and calculate secondary structure percentage

        for ranked in self.ranked_list:
            results_dict, domains_dict = bioutils.aleph_annotate(output_path=self.tmp_dir, pdb_path=ranked.split_path)
            if domains_dict is None or results_dict is None:
                break
            ranked.set_secondary_structure(ah=results_dict['ah'], bs=results_dict['bs'],
                                           total_residues=results_dict['number_total_residues'])
            if ranked.filtered:
                ranked.set_minimized_path(os.path.join(self.tmp_dir, f'{ranked.name}_minimized.pdb'))
                ranked.set_energies(bioutils.run_openmm(pdb_in_path=ranked.path, pdb_out_path=ranked.minimized_path))
                interfaces_data_list = bioutils.find_interface_from_pisa(ranked.split_path, self.interfaces_path)
                if interfaces_data_list:
                    deltas_list = [interface['deltaG'] for interface in interfaces_data_list]
                    deltas_list = utils.normalize_list([deltas_list])
                    for i, interface in enumerate(interfaces_data_list):
                        if not interface["chain1"] in domains_dict or not interface["chain2"] in domains_dict:
                            continue
                        code = f'{interface["chain1"]}{interface["chain2"]}'
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

        # Superpose the experimental pdb with all the rankeds and templates
        for experimental in experimental_pdbs:
            aux_dict = {}
            for pdb in [ranked.split_path for ranked in self.ranked_list] + list(self.templates_dict.values()):
                rmsd, aligned_residues, quality_q = bioutils.gesamt_pdbs([pdb, experimental])
                total_residues = bioutils.get_number_residues(pdb)
                rmsd = round(rmsd, 2) if rmsd is not None else str(rmsd)
                aux_dict[utils.get_file_name(pdb)] = structures.TemplateRanked(pdb, rmsd, aligned_residues,
                                                                               total_residues)
            output_pdb = os.path.join(self.output_dir, os.path.basename(experimental))
            bioutils.gesamt_pdbs([reference_superpose, experimental], output_pdb)
            self.experimental_dict[utils.get_file_name(experimental)] = aux_dict

        for template, template_path in self.templates_nonsplit_dict.items():
            frobenius_file = os.path.join(self.frobenius_path, f'frobenius_{template}.txt')
            matrices = os.path.join(self.tmp_dir, 'matrices')
            template_matrix = os.path.join(matrices, f'{utils.get_file_name(template_path)}_ang.npy')
            path_list = [ranked.path for ranked in self.ranked_list if ranked.filtered]
            with open(frobenius_file, 'w') as sys.stdout:
                _, list_targets, list_core, _, list_frobenius_angles, list_frobenius_distances, list_plot_ang, list_plot_dist, _ = ALEPH.frobenius(
                    references=[template_path],
                    targets=list(path_list), write_plot=True,
                    write_matrix=True)
            sys.stdout = sys.__stdout__
            ranked_filtered = [ranked for ranked in self.ranked_list if ranked.filtered]

            interfaces_data_list = bioutils.find_interface_from_pisa(self.templates_dict[template],
                                                                     self.interfaces_path)
            template_interface_list = []
            for interface in interfaces_data_list:
                template_interface_list.append(f'{interface["chain1"]}{interface["chain2"]}')
            self.template_interfaces[template] = template_interface_list

            for ranked in ranked_filtered:
                index = list_targets.index(ranked.name)
                ranked.add_frobenius_plot(
                    template=template,
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
                    new_name = os.path.join(self.frobenius_path, f'{template}_{ranked.name}_{interface.name}.png')
                    plot_path = os.path.join(self.tmp_dir, os.path.basename(plot))
                    interface.add_frobenius_information(
                        template=template,
                        dist_coverage=fro_distance,
                        core=fro_core,
                        dist_plot=shutil.copy2(plot_path, new_name)
                    )

        os.chdir(store_old_dir)

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
                data = {'ranked': rmsd_dict.keys()}
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
                f_in.write('rankeds PLDDT\n')
                data = {'ranked': plddt_dict.keys(),
                        'plddt': [value for value in plddt_dict.values()],
                        }
                df = pd.DataFrame(data)
                f_in.write(df.to_markdown())

            if bool(energies_dict):
                f_in.write('\n\n')
                f_in.write('OPENMM Energies\n')
                data = {'ranked': energies_dict.keys(),
                        'kinetic': [value['kinetic'] for value in energies_dict.values()],
                        'potential': [value['potential'] for value in energies_dict.values()]
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
                            f'{value.rmsd} {value.aligned_residues} ({value.total_residues})')
                df = pd.DataFrame(data)
                f_in.write(df.to_markdown())

            f_in.write('\n\n')
