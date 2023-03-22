import logging
import os
import re
import shutil
import sys
import statistics
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from typing import Dict, List
from Bio.PDB import PDBParser, Selection
from ALEPH.aleph.core import ALEPH
from itertools import combinations
from libs import bioutils, change_res, features, utils, sequence, structures

PERCENTAGE_FILTER = 0.8
PERCENTAGE_MAX_RMSD = 1
GROUPS = ['GAVLI', 'FYW', 'CM', 'ST', 'KRH', 'DENQ', 'P']
MATPLOTLIB_FONT = 14
plt.set_loglevel('WARNING')


def check_ranked(path: str) -> bool:
    return re.match('ranked_[0-9]+.pdb', path) or re.match('cluster_[0-9]+_ranked_[0-9]+.pdb', path)


def read_rankeds(input_path: str) -> List[str]:
    ranked_paths = [path for path in os.listdir(input_path) if check_ranked(path)]
    return [structures.Ranked(os.path.join(input_path, path)) for path in utils.sort_by_digit(ranked_paths)]


def get_group(res: str) -> str:
    group = [s for s in GROUPS if res in s]
    if group:
        return group[0]
    return res


def compare_sequences(sequence1: str, sequence2: str) -> List[str]:
    # Given two sequences with same length, return a list showing
    # if there is a match, a group match, they are different, or
    # they are not aligned
    return_list = []

    for i in range(0, len(sequence1)):
        if i < len(sequence2):
            res1 = sequence1[i]
            res2 = sequence2[i]
            if res1 == res2:
                return_list.append('0')
            elif res2 == '-':
                return_list.append('-')
            elif get_group(res1) == get_group(res2):
                return_list.append('0.4')
            else:
                return_list.append('0.7')
        else:
            return_list.append('-')

    return return_list


def convert_residues(residues_list: List[List], sequence_assembled):
    for i in range(0, len(residues_list)):
        if residues_list[i] is not None:
            for residue in residues_list[i]:
                result = sequence_assembled.get_real_residue_number(i, residue)
                if result is not None:
                    residues_list.append(result)
    return residues_list


def plot_plddt(plot_path: str, ranked_list: List) -> float:
    plt.figure(figsize=(18, 6))
    plt.rcParams.update({'font.size': MATPLOTLIB_FONT})
    for ranked in ranked_list:
        return_dict = bioutils.read_bfactors_from_residues(pdb_path=ranked.path)
        plddt_list = [value for value in list(return_dict.values())[0] if value is not None]
        res_list = [int(item) for item in range(1, len(plddt_list) + 1)]
        ranked.set_plddt(round(statistics.median(map(float, plddt_list)), 2))
        plt.plot(res_list, plddt_list, label=ranked.name)
    plt.legend()
    plt.xlabel('residue number')
    plt.ylabel('pLDDT')
    plt.savefig(plot_path, dpi=100)
    plt.cla()
    max_plddt = max([ranked.plddt for ranked in ranked_list])
    return max_plddt


def plot_cc_analysis(plot_path: str, analysis_dict: Dict, clusters: List, predictions: bool = False):
    plt.figure(figsize=(8, 8))
    plt.rcParams.update({'font.size': MATPLOTLIB_FONT})
    for key, values in analysis_dict.items():
        if key.startswith('cluster_'):
            color = 'red'
        else:
            color = 'blue'
        if key in [utils.get_file_name(clus) for clus in clusters[0]]:
            plt.scatter(values.x, values.y, marker='*', color=color, label='Cluster 0')
        else:
            plt.scatter(values.x, values.y, color=color, label='Cluster 1')
        plt.annotate(key, (values.x, values.y), horizontalalignment='right', verticalalignment='top',)
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())    
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    if predictions:
        plt.title('TEMPLATES AND PREDICTIONS CLUSTERING')
    else:
        plt.title('TEMPLATES CLUSTERING')
    plt.savefig(plot_path, dpi=100)
    plt.cla()


def plot_sequence(plot_path: str, a_air):
    plt.rcParams.update({'font.size': MATPLOTLIB_FONT})
    fig, ax = plt.subplots(1, figsize=(16, 0.5))
    legend_seq = [Patch(label='Sequence', color='tab:cyan'),Patch(label='Linker', color='tab:blue')]
    ax.legend(handles=legend_seq[::-1], loc='upper center', bbox_to_anchor=(0.5, -0.4), fancybox=True, framealpha=0.5, ncol=2)
    for i in range(a_air.sequence_assembled.total_copies):
        ax.barh('sequence', a_air.sequence_assembled.get_sequence_length(i),
                left=a_air.sequence_assembled.get_starting_length(i) + 1, color='tab:cyan')
        ax.barh('sequence', a_air.sequence_assembled.glycines,
                left=a_air.sequence_assembled.get_starting_length(i) + a_air.sequence_assembled.get_sequence_length(
                    i) + 1, color='tab:blue')

        xcenters = (a_air.sequence_assembled.get_starting_length(i) + 1) + a_air.sequence_assembled.get_sequence_length(i) / 2
        ax.text(xcenters, 0, a_air.sequence_assembled.get_sequence_name(i), ha='center', va='center')

    ax_secondary = ax.secondary_xaxis('top')
    ax_secondary.set_xticks([a_air.sequence_assembled.get_starting_length(i) + 1 for i in range(a_air.sequence_assembled.total_copies)], rotation=45)
    ax_secondary.set_xticks(list(ax_secondary.get_xticks())+[a_air.sequence_assembled.get_starting_length(i) + a_air.sequence_assembled.get_sequence_length(i) + 1 for i in range(a_air.sequence_assembled.total_copies)], rotation=45)
    ax_secondary.set_xticklabels([1]*a_air.sequence_assembled.total_copies+[a_air.sequence_assembled.get_sequence_length(i)+1 for i in range(a_air.sequence_assembled.total_copies)], rotation = 45)
    ax.set_xticks([a_air.sequence_assembled.get_starting_length(i) + 1 for i in range(a_air.sequence_assembled.total_copies)], rotation=45)
    ax.set_xticks(list(ax.get_xticks())+[a_air.sequence_assembled.get_starting_length(i) + a_air.sequence_assembled.get_sequence_length(i) + 1 for i in range(a_air.sequence_assembled.total_copies)], rotation=45)
    ax.set_xticklabels(ax.get_xticks(), rotation = 45)
    ax.set_xlim(0, len(a_air.sequence_assembled.sequence_assembled) + a_air.sequence_assembled.glycines)
    ax.set_yticks([])
    fig.tight_layout()
    fig.subplots_adjust(top=.95)
    plt.savefig(plot_path, bbox_inches='tight', dpi=100)
    plt.cla()


def plot_gantt(plot_type: str, plot_path: str, a_air) -> str:
    plt.rcParams.update({'font.size': MATPLOTLIB_FONT})
    fig, ax = plt.subplots(1, figsize=(16, 2))
    legend_elements = []
    legend_seq = [Patch(label='Sequence', color='tab:cyan'),Patch(label='Linker', color='tab:blue')]
    number_of_templates = 1

    total_length = len(a_air.sequence_assembled.sequence_assembled) + a_air.sequence_assembled.glycines
    for i in range(a_air.sequence_assembled.total_copies):
        ax.barh('sequence', a_air.sequence_assembled.get_sequence_length(i),
                left=a_air.sequence_assembled.get_starting_length(i) + 1, color='tab:cyan')
        ax.barh('sequence', a_air.sequence_assembled.glycines,
                left=a_air.sequence_assembled.get_starting_length(i) + a_air.sequence_assembled.get_sequence_length(
                    i) + 1, color='tab:blue')

    if plot_type == 'msa':
        title = 'MSA'
        file = os.path.join(plot_path, 'msa_gantt.png')
    else:  # should be template:
        title = 'TEMPLATES'
        file = os.path.join(plot_path, 'template_gantt.png')

    if plot_type == 'msa':
        names = a_air.feature.get_names_msa()
    else:
        names = a_air.feature.get_names_templates()

    names = [name for name in names if name != '']
    if len(names) > 30:
        number_of_templates += 1
        aligned_sequences = [0] * len(a_air.sequence_assembled.sequence_assembled)
        for name in names:
            if plot_type == 'msa':
                features_search = a_air.feature.get_msa_by_name(name)
            else:
                features_search = a_air.feature.get_sequence_by_name(name)
            aligned_sequence = compare_sequences(a_air.sequence_assembled.sequence_assembled, features_search)
            new_array = [0 if align == '0' else 1 for align in aligned_sequence]
            aligned_sequences = np.add(new_array, aligned_sequences)
        aligned_sequences = [aligned / len(names) for aligned in aligned_sequences]
        for i in range(1, len(aligned_sequences)):
            ax.barh('Percentage', 1, left=i, height=0.5, color=str(aligned_sequences[i - 1]))
    else:
        pdb_hits_path = os.path.join(a_air.results_dir, 'msas/pdb_hits.hhr')
        hhr_text = ''
        if os.path.exists(pdb_hits_path):
            hhr_text = open(pdb_hits_path, 'r').read()
        for j, name in reversed(list(enumerate(names))):
            number_of_templates += 1
            template = a_air.get_template_by_id(name)
            changed_residues = []
            changed_fasta = []
            deleted_residues = []
            if template is not None:
                changed_residues, changed_fasta, _ = template.get_changes()
                changed_residues = convert_residues(changed_residues, a_air.sequence_assembled)
                changed_fasta = convert_residues(changed_fasta, a_air.sequence_assembled)
                text = ''
                if len(name) > 4:
                    template_name = f'T{j + 1}'
                    text = f'\n{template_name} ({name}):'
                else:
                    template_name = name
                    text =  f'\n{template_name}:'
                for alignment in template.get_results_alignment():
                    if alignment is not None:
                        text += f'\n   Chain {alignment.database.chain}: Aligned={alignment.aligned_columns}({alignment.total_columns}) Evalue={alignment.evalue} Identities={alignment.identities}'
                    else:
                        text += f'\n   No alignment'
                legend_elements.append(text)
            else:
                template_name = name

            if plot_type == 'msa':
                features_search = a_air.feature.get_msa_by_name(name)
            else:
                features_search = a_air.feature.get_sequence_by_name(name)
                if hhr_text != '':
                    match = re.findall(rf'(\d+ {name.upper()}.*$)', hhr_text, re.M)
                    if match:
                        match_split = match[0].split()[-9:]
                        legend_elements.append(f'{template_name}: Aligned={match_split[5]}({match_split[8].replace("(","").replace(")","")}) Evalue={match_split[2]}')
                    else:
                        legend_elements.append(f'{template_name}')

            if features_search is not None:
                aligned_sequence = compare_sequences(a_air.sequence_assembled.sequence_assembled, features_search)
                for i in range(1, len(features_search)):
                    if aligned_sequence[i - 1] != '-':
                        ax.barh(template_name, 1, left=i, height=0.5, color=str(aligned_sequence[i - 1]))
                        if i in changed_residues:
                            ax.barh(template_name, 1, left=i, height=0.25, align='edge', color='yellow')
                        elif i in changed_fasta:
                            ax.barh(template_name, 1, left=i, height=0.25, align='edge', color='red')
                        else:  
                            ax.barh(template_name, 1, left=i, height=0.25, align='edge', color='white')
                        ax.barh(template_name, 1, left=i, height=0.1, align='edge', color=str(aligned_sequence[i - 1]))


    ax.xaxis.grid(color='k', linestyle='dashed', alpha=0.4, which='both')
    plt.setp([ax.get_xticklines()], color='k')
    ax.set_xlim(0, total_length)
    ax.set_ylim(0, number_of_templates)
    if number_of_templates > 10:
        fig.set_size_inches(16, number_of_templates*0.7)
    else:
        fig.set_size_inches(16, number_of_templates)
    legend_elements.append('The templates gray scale shows the similarity between the aligned template sequence and the input sequence.\n'
                           'The darker parts indicate that the residues are the same or belong to the same group.\n')
    legend_elements.append('Yellow shows which residues have been changed to another specific residue\n'
                            'whereas the red shows which residues have been changed from another query sequence.\n'
                            'No information (white) implies that no modifications have been done.')
    legend_elements.reverse()
    plt.figtext(0.05, -0.05, '\n'.join(legend_elements), va='top')

    ax.set_xticks([a_air.sequence_assembled.get_starting_length(i) + 1 for i in range(a_air.sequence_assembled.total_copies)])
    ax.set_xticks(list(ax.get_xticks())+[a_air.sequence_assembled.get_starting_length(i) + a_air.sequence_assembled.get_sequence_length(i) + 1 for i in range(a_air.sequence_assembled.total_copies)])
    
    cut_chunk = [list(tup) for tup in a_air.chunk_list]
    cut_chunk = utils.remove_list_layer(cut_chunk)
    ax.set_xticks(list(ax.get_xticks())+[cut+1 for cut in cut_chunk])
    ax.set_xticklabels(ax.get_xticks(), rotation = 45)

    ax.set_xlabel('Residue number')
    ax.set_ylabel('Sequences')
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_color('k')

    ax.legend(handles=legend_seq[::-1], loc='upper right', framealpha=0.5, fancybox=True, ncol=2)
    fig.tight_layout()
    fig.subplots_adjust(top=.95)
    plt.title(title)
    plt.savefig(file, bbox_inches='tight', dpi=100)
    plt.cla()
    return file


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
        self.gantt_plots_path: List[str] = []
        self.ranked_list: List[structures.Ranked] = []
        self.output_dir: str = output_dir
        self.experimental_dict = {}
        self.results_dir: str
        self.templates_nonsplit_dir: str
        self.rankeds_nonsplit_dir: str
        self.rankeds_split_dir: str
        self.aleph_results_path: str
        self.group_ranked_by_rmsd_dict: dict = {}
        self.template_interfaces: dict = {}
        self.templates_cluster: List = [[],[]]
        self.templates_predictions_cluster: List = [[],[]]

        utils.create_dir(dir_path=self.plots_path, delete_if_exists=True)
        utils.create_dir(dir_path=self.templates_path, delete_if_exists=True)
        utils.create_dir(dir_path=self.interfaces_path, delete_if_exists=True)
        utils.create_dir(dir_path=self.frobenius_path, delete_if_exists=True)


    def create_plot_gantt(self, a_air):
        self.gantt_plots_path = [plot_gantt(plot_type='template', plot_path=self.plots_path, a_air=a_air)]
        self.gantt_plots_path.append(plot_gantt(plot_type='msa', plot_path=self.plots_path, a_air=a_air))
        plot_sequence(plot_path=self.sequence_plot_path, a_air=a_air)


    def cc_analysis(self, paths_in: Dict, cc_analysis_paths: structures.CCAnalysis):
        cc_path = os.path.join(self.results_dir, 'cc_analysis')
        ranked = False
        utils.create_dir(cc_path, delete_if_exists=True)
        paths = [shutil.copy2(path, cc_path) for path in paths_in.values()]
        trans_dict = {}
        for index, path in enumerate(paths):
            if check_ranked(os.path.basename(path)):
                ranked = True
                bfactors_dict = bioutils.read_bfactors_from_residues(path)
                for chain, residues in bfactors_dict.items():
                    for i in range(len(residues)):
                        if bfactors_dict[chain][i] is not None:
                            bfactors_dict[chain][i] = round(bfactors_dict[chain][i] - 70.0, 2)
                residues_dict = bioutils.read_residues_from_pdb(path)
                change = change_res.ChangeResidues(chain_res_dict=residues_dict, chain_bfactors_dict=bfactors_dict)
                change.change_bfactors(path, path)
            new_path = os.path.join(cc_path, f'orig.{str(index)}.pdb')
            os.rename(os.path.join(cc_path, path), new_path)
            trans_dict[index] = utils.get_file_name(path)
        if trans_dict:
            logging.info('Running pdb2cc and ccanalysis with the following templates:')
            logging.info(' ,'.join([f'{key} = {value}' for key, value in trans_dict.items()]))
            output_pdb2cc = bioutils.run_pdb2cc(templates_path=cc_path, pdb2cc_path=cc_analysis_paths.pd2cc_path)
            if os.path.exists(output_pdb2cc):
                output_cc = bioutils.run_cc_analysis(input_path=output_pdb2cc,
                                                     cc_analysis_path=cc_analysis_paths.cc_analysis_path)
                if os.path.exists(output_cc):
                    cc_analysis_dict = utils.parse_cc_analysis(file_path=output_cc)
                    clean_dict = {}
                    for key, values in cc_analysis_dict.items():
                        if values.module > 0.1 or values.module < -0.1:
                            clean_dict[trans_dict[int(key) - 1]] = values
                    points = np.array([[values.x, values.y] for values in clean_dict.values()])
                    kmeans = KMeans(n_clusters=2)
                    kmeans.fit(points)
                    aux_templates_cluster = [[],[]]
                    for i, label in enumerate(kmeans.labels_):
                        aux_templates_cluster[int(label)].append(paths_in[list(clean_dict.keys())[i]])
                    if self.templates_cluster[0]:
                        first_element = self.templates_cluster[0][0]
                        if first_element in aux_templates_cluster[1]:
                            aux_templates_cluster.append(aux_templates_cluster.pop(0))
                    if not ranked:
                        self.templates_cluster = aux_templates_cluster
                        plot_path = self.analysis_plot_path
                    else:
                        self.templates_predictions_cluster = aux_templates_cluster
                        plot_path = self.analysis_ranked_plot_path
                    plot_cc_analysis(plot_path=plot_path, analysis_dict=clean_dict,
                                clusters=aux_templates_cluster, predictions=ranked)

    def analyse_output(self, results_dir: str, sequence_assembled: sequence.SequenceAssembled, feature: features.Features,
                       experimental_pdb: str, cc_analysis_paths, cluster_templates: bool = False):
        # Read all templates and rankeds, if there are no ranked, raise an error

        self.results_dir = results_dir
        self.templates_nonsplit_dir = f'{self.results_dir}/templates_nonsplit'
        self.rankeds_split_dir = f'{self.results_dir}/rankeds_split'
        self.tmp_dir = f'{self.results_dir}/tmp'
        self.aleph_results_path = f'{self.tmp_dir}/output.json'

        utils.create_dir(dir_path=self.templates_nonsplit_dir, delete_if_exists=True)
        utils.create_dir(dir_path=self.rankeds_split_dir, delete_if_exists=True)
        utils.create_dir(dir_path=self.tmp_dir, delete_if_exists=True)

        store_old_dir = os.getcwd()
        os.chdir(self.tmp_dir)

        templates_dict = {}
        templates_nonsplit_dict = {}

        if feature is not None:
            templates_nonsplit_dict = feature.write_all_templates_in_features(output_dir=self.templates_nonsplit_dir,
                                                                        print_number=False)
        # Split the templates with chains
        for template, template_path in templates_nonsplit_dict.items():
            new_pdb_path = os.path.join(self.templates_path, f'{template}.pdb')
            templates_dict[template] = new_pdb_path

            if not cluster_templates:
                shutil.copy2(template_path, new_pdb_path)
                bioutils.split_chains_assembly(pdb_in_path=new_pdb_path,
                                           pdb_out_path=new_pdb_path,
                                           sequence_assembled=sequence_assembled)

        self.ranked_list = read_rankeds(input_path=self.results_dir)

        #dendogram_file = os.path.join(self.run_dir, 'dendogram.txt')
        #dendogram_plot = os.path.join(self.run_dir, 'clustering_dendogram_angles.png') 
        #with open(dendogram_file, 'w') as sys.stdout:
        #    _, _, _, _, _, _, _, _, dendogram_list = ALEPH.frobenius(references=list(  .values()),
        #                                            targets=list(template_nonsplit.values()), write_plot=True,
        #                                            write_matrix=True)
        #sys.stdout = sys.__stdout__
        #if dendogram_list:
        #    shutil.copy2(dendogram_plot, self.template_dendogram)
        #    if not custom_features:
        #        for templates in dendogram_list:
        #            self.templates_cluster.append([template_nonsplit[template] for template in templates])
        
        self.cc_analysis(paths_in=templates_dict, cc_analysis_paths=cc_analysis_paths)

        if not self.ranked_list:
            logging.info('No ranked PDBs found')
            os.chdir(store_old_dir)
            return
        
        # Create a plot with the ranked pLDDTs, also, calculate the maximum pLDDT
        max_plddt = plot_plddt(plot_path=self.plddt_plot_path, ranked_list=self.ranked_list)
        bioutils.write_sequence(sequence_name=utils.get_file_name(self.sequence_path),
                                sequence_amino=sequence_assembled.sequence_assembled, sequence_path=self.sequence_path)

        # Save split path of all rankeds, taking into account the split dir
        [ranked.set_split_path(os.path.join(self.rankeds_split_dir, os.path.basename(ranked.path))) for ranked in self.ranked_list]

        # Sort list of ranked by pLDDT
        self.ranked_list.sort(key=lambda x: x.plddt, reverse=True)
        shutil.copy2(self.ranked_list[0].path, self.ranked_list[0].split_path)
        results = [items for items in combinations(self.ranked_list, r=2)]
        for result in results:
            if result[0].name == self.ranked_list[0].name:
                if not cluster_templates:
                    rmsd, _, _ = bioutils.superpose_pdbs([result[1].path, result[0].path], result[1].split_path)
                else:
                    rmsd, _, _ = bioutils.superpose_pdbs([result[1].path, result[0].path])
                    shutil.copy2(result[1].path, result[1].split_path)
            else:
                rmsd, _, _ = bioutils.superpose_pdbs([result[1].path, result[0].path])
            result[0].set_ranked_to_rmsd_dict(rmsd=rmsd, ranked_name=result[1].name)
            result[1].set_ranked_to_rmsd_dict(rmsd=rmsd, ranked_name=result[0].name)

        reference_superpose = self.ranked_list[0].path
        green_color = 40

        # Filter rankeds, split them in chains.
        for ranked in self.ranked_list:
            if ranked.plddt >= (PERCENTAGE_FILTER * max_plddt):
                ranked.set_filtered(True)
                ranked.set_split_path(shutil.copy2(ranked.split_path, os.path.join(self.output_dir, os.path.basename(ranked.path))))

            mapping = bioutils.split_chains_assembly(pdb_in_path=ranked.split_path,
                                                     pdb_out_path=ranked.split_path,
                                                     sequence_assembled=sequence_assembled)
            ranked.set_mapping(mapping)
            ranked.set_encoded(ranked.split_path)

        for ranked in self.ranked_list:
            if ranked.filtered:
                found = False
                for ranked2 in self.ranked_list:
                    if ranked2.filtered and ranked2.name != ranked.name and ranked2.name in self.group_ranked_by_rmsd_dict and ranked.rmsd_dict[ranked2.name] <= PERCENTAGE_MAX_RMSD:
                        self.group_ranked_by_rmsd_dict[ranked2.name].append(ranked)
                        found = True
                        ranked.set_rmsd(ranked2.rmsd_dict[ranked.name])
                        if self.ranked_list[0].name == ranked2.name:
                            ranked.set_green_color(green_color)
                            ranked.set_best(True)
                            green_color += 10
                        break

                if not found:
                    self.group_ranked_by_rmsd_dict[ranked.name] = [ranked]
                    ranked.set_rmsd(0)
                    if self.ranked_list[0].name == ranked.name:
                        ranked.set_green_color(green_color)
                        ranked.set_best(True)
                        green_color += 10

        aux_dict = dict({ranked.name:ranked.split_path for ranked in self.ranked_list}, **templates_dict)
        self.cc_analysis(paths_in=aux_dict, cc_analysis_paths=cc_analysis_paths)

        # Superpose each template with all the rankeds.
        if templates_dict:
            for i, ranked in enumerate(self.ranked_list):
                for template, template_path in templates_dict.items():
                    total_residues = len(
                        [res for res in Selection.unfold_entities(PDBParser().get_structure(template, template_path), 'R')])
                    if i == 0 and not cluster_templates:
                        rmsd, aligned_residues, quality_q = bioutils.superpose_pdbs([template_path, ranked.split_path], template_path)
                    else:
                        rmsd, aligned_residues, quality_q = bioutils.superpose_pdbs([template_path, ranked.split_path])
                    if rmsd is not None:
                        rmsd = round(rmsd, 2)
                    ranked.add_template(structures.TemplateRanked(template, rmsd, aligned_residues, total_residues))
                ranked.sort_template_rankeds()

        # Use aleph to generate domains and calculate secondary structure percentage
        for ranked in self.ranked_list:
            aleph_file = os.path.join(self.tmp_dir, f'aleph_{ranked.name}.txt')
            with open(aleph_file, 'w') as sys.stdout:
                try:
                    ALEPH.annotate_pdb_model(reference=ranked.split_path, strictness_ah=0.45, strictness_bs=0.2,
                                            peptide_length=3, width_pic=1, height_pic=1, write_graphml=False, write_pdb=True)
                except:
                    pass
            sys.stdout = sys.__stdout__
            if os.path.exists(self.aleph_results_path):
                result_dict = utils.parse_aleph_annotate(file_path=self.aleph_results_path)
                ranked.set_secondary_structure(ah=result_dict['ah'], bs=result_dict['bs'],
                                               total_residues=result_dict['number_total_residues'])
                aleph_txt_path = f'{self.tmp_dir}/aleph_{ranked.name}.txt'
                domains_dict = utils.parse_aleph_ss(aleph_txt_path)
            else:
                break
            if ranked.filtered and os.path.exists(aleph_txt_path):
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
        if experimental_pdb is not None:
            for ranked in self.ranked_list:
                rmsd, nalign, quality_q = bioutils.superpose_pdbs([experimental_pdb, ranked.split_path])
                self.experimental_dict[ranked.name] = round(rmsd, 2)
            for template, template_path in templates_dict.items():
                rmsd, nalign, quality_q = bioutils.superpose_pdbs([experimental_pdb, template_path])
                self.experimental_dict[template] = round(rmsd, 2)
            output_pdb = os.path.join(self.output_dir, os.path.basename(experimental_pdb))
            bioutils.superpose_pdbs([experimental_pdb, reference_superpose], output_pdb)

        for template, template_path in templates_nonsplit_dict.items():
            frobenius_file = os.path.join(self.frobenius_path, f'frobenius_{template}.txt')
            matrices = os.path.join(self.tmp_dir, 'matrices')
            template_matrix = os.path.join(matrices, f'{utils.get_file_name(template_path)}_ang.npy')
            path_list = [ranked.path for ranked in self.ranked_list if ranked.filtered]
            with open(frobenius_file, 'w') as sys.stdout:
                _, list_targets, list_core, _, list_frobenius_angles, list_frobenius_distances, list_plot_ang, list_plot_dist, _ = ALEPH.frobenius(references=[template_path],
                                                                                                                                                targets=list(path_list), write_plot=True,
                                                                                                                                                write_matrix=True)
            sys.stdout = sys.__stdout__
            ranked_filtered = [ranked for ranked in self.ranked_list if ranked.filtered]

            interfaces_data_list = bioutils.find_interface_from_pisa(templates_dict[template], self.interfaces_path)
            template_interface_list = []
            for interface in interfaces_data_list:
                template_interface_list.append(f'{interface["chain1"]}{interface["chain2"]}')
            self.template_interfaces[template] = template_interface_list

            for ranked in ranked_filtered:
                index = list_targets.index(ranked.name)
                ranked.add_frobenius_plot(
                    template=template, 
                    dist_plot=[shutil.copy2(plot, self.frobenius_path) for plot in list_plot_dist if ranked.name in os.path.basename(plot)][0],
                    ang_plot=[shutil.copy2(plot, self.frobenius_path) for plot in list_plot_ang if ranked.name in os.path.basename(plot)][0], 
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
                        fro_distance, fro_core, plot = ALEPH.frobenius_submatrices(path_ref=template_matrix, path_tar=ranked_matrix,
                                                                                    residues_tar=interface.res_list, write_plot=True,
                                                                                    title=f'Interface: {interface.name}')
                    sys.stdout = sys.__stdout__
                    new_name = os.path.join(self.frobenius_path, f'{template}_{ranked.name}_{interface.name}.png')
                    plot_path = os.path.join(self.tmp_dir, os.path.basename(plot))
                    interface.add_frobenius_information(
                        template = template,
                        dist_coverage = fro_distance,
                        core = fro_core,
                        dist_plot = shutil.copy2(plot_path, new_name)
                    )

        os.chdir(store_old_dir)


    def write_tables(self, rmsd_dict: Dict, ranked_rmsd_dict: Dict, secondary_dict: Dict, plddt_dict: Dict, energies_dict: Dict):
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
                f_in.write(f'Superposition with experimental structure\n')
                data = {'pdb': self.experimental_dict.keys(),
                        'rmsd': self.experimental_dict.values()
                        }
                df = pd.DataFrame(data)
                f_in.write(df.to_markdown())

            f_in.write('\n\n')
