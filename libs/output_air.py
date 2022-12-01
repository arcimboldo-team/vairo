import logging
import os
import re
import shutil
import sys
import statistics
import matplotlib.pyplot as plt

from matplotlib.patches import Patch
import pandas as pd
from typing import Dict, List
from Bio.PDB import PDBParser, Selection
from ALEPH.aleph.core import ALEPH
from libs import bioutils, features, utils, sequence, structures

PERCENTAGE_FILTER = 0.8
GROUPS = ['GAVLI', 'FYW', 'CM', 'ST', 'KRH', 'DENQ', 'P']
plt.set_loglevel('WARNING')


def plot_plddt(plot_path: str, ranked_models_dict: Dict) -> Dict:
    return_plddt_dict = {}

    plt.clf()
    for ranked, ranked_path in utils.sort_by_digit(ranked_models_dict):
        plddt_list = []
        with open(ranked_path) as f:
            for line in f.readlines():
                if line[:4] == 'ATOM' and line[13:16] == 'CA ':
                    plddt_list.append(float(line[60:66].replace(" ", "")))
        res_list = [int(item) for item in range(1, len(plddt_list) + 1)]
        return_plddt_dict[ranked] = round(statistics.median(map(float, plddt_list)), 2)
        plt.plot(res_list, plddt_list, label=ranked)
    plt.legend()
    plt.xlabel('residue number')
    plt.ylabel('pLDDT')
    plt.savefig(plot_path)

    return return_plddt_dict

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

def plot_gantt(plot_type: str, plot_path: str, a_air):
    fig , (ax, ax1) = plt.subplots(2, figsize=(16, 6), gridspec_kw={'height_ratios':[6, 1]})
    total_length = len(a_air.sequence_assembled.sequence_assembled) + a_air.sequence_assembled.glycines
    height = 0.3
    for i in range(a_air.sequence_assembled.total_copies):
        ax.barh('sequence', a_air.sequence_assembled.get_sequence_length(i), height=height,
                left=a_air.sequence_assembled.get_starting_length(i) + 1, color='tab:cyan')
        ax.barh('sequence', a_air.sequence_assembled.glycines, height=height,
                left=a_air.sequence_assembled.get_starting_length(i) + a_air.sequence_assembled.get_sequence_length(i) + 1,
                color='tab:blue')

    if plot_type == 'msa':
        title = 'MSA'
        file = os.path.join(plot_path, 'msa_gantt.png')
    else:  # should be template:
        title = 'TEMPLATES'
        file = os.path.join(plot_path, 'template_gantt.png')

    names = a_air.feature.get_names()
    legend_elements = []
    for j, name in reversed(list(enumerate(names))):
        template = a_air.get_template_by_id(name)
        results_alignment_text = template.get_results_alignment_text()
        template_name = f'T{j+1}'
        legend_elements.append(Patch(label=f'{template_name} ({name}): {results_alignment_text}'))
        if plot_type == 'msa':
            features_search = a_air.feature.get_msa_by_name(name)
        else:  # should be template:
            features_search = a_air.feature.get_sequence_by_name(name)
        if features_search is not None:
            aligned_sequence = compare_sequences(a_air.sequence_assembled.sequence_assembled, features_search)
            for i in range(1, len(features_search)):
                if aligned_sequence[i - 1] != '-':
                    ax.barh(template_name, 1, height=height, left=i, color=str(aligned_sequence[i - 1]))
    ax.xaxis.grid(color='k', linestyle='dashed', alpha=0.4, which='both')
    plt.setp([ax.get_xticklines()], color='k')
    ax.set_xlim(0, total_length)
    ax.set_xticks([a_air.sequence_assembled.get_starting_length(i) + 1 for i in range(a_air.sequence_assembled.total_copies)])
    ax.set_xlabel('Residue number')
    ax.set_ylabel('Sequences')
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_color('k')

    # Put a legend below current axis
    legend = ax1.legend(handles=legend_elements[::-1], loc='upper left', handlelength=0, handletextpad=0, fancybox=True)
    plt.setp(legend.get_texts(), color='k')
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.set_xticks([])
    ax1.set_yticks([])
    plt.suptitle(title)
    plt.savefig(file, bbox_inches='tight', dpi=100)
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
        self.html_path: str = f'{output_dir}/output.html'
        self.gantt_plots_path: List[str] = []
        self.output_dir: str = output_dir
        self.plddt_dict: Dict = {}
        self.experimental_dict: Dict = {}
        self.secondary_dict: Dict = {}
        self.rmsd_dict: Dict = {}
        self.frobenius_plots: Dict = {}
        self.ranked_filtered: Dict = {}
        self.best_ranked: str = ''


        self.run_dir: str = None
        self.nonsplit_path: str = None
        self.aleph_results_path: str = None

        utils.create_dir(dir_path=self.plots_path, delete_if_exists=True)
        utils.create_dir(dir_path=self.templates_path, delete_if_exists=True)
        utils.create_dir(dir_path=self.interfaces_path, delete_if_exists=True)
        utils.create_dir(dir_path=self.frobenius_path, delete_if_exists=True)


    def create_plot_gantt(self, a_air):
        self.gantt_plots_path = [plot_gantt(plot_type='template', plot_path=self.plots_path, a_air=a_air)]
        self.gantt_plots_path.append(plot_gantt(plot_type='msa', plot_path=self.plots_path, a_air=a_air))


    def set_run_dir(self, run_dir: str):
        self.run_dir = run_dir
        self.nonsplit_path = f'{run_dir}/nonsplit'
        self.aleph_results_path = f'{run_dir}/output.json'
        utils.create_dir(dir_path=self.nonsplit_path, delete_if_exists=True)


    def analyse_output(self, sequence_assembled: sequence.SequenceAssembled, feature: features.Features, experimental_pdb: str):
        # Read all templates and rankeds, if there are no ranked, raise an error
        template_dict = {}
        template_nonsplit = {}
        if feature is not None:
            template_nonsplit = feature.write_all_templates_in_features(output_dir=self.nonsplit_path, print_number=False)

        # Split the templates with chains
        for template, template_path in template_nonsplit.items():
            new_pdb_path = os.path.join(self.templates_path, f'{template}.pdb')
            shutil.copy2(template_path, new_pdb_path)
            template_dict[template] = new_pdb_path
            bioutils.split_chains_assembly(pdb_in_path=new_pdb_path,
                                        pdb_out_path=new_pdb_path,
                                        sequence_assembled=sequence_assembled)

        ranked_models_dict = {utils.get_file_name(ranked): os.path.join(self.run_dir, ranked) for ranked in
                            os.listdir(self.run_dir) if re.match('ranked_[0-9]+.pdb', ranked)}
        if not bool(ranked_models_dict):
            logging.info('No ranked PDBs found')
            return

        # Create a plot with the ranked plddts, also, calculate the maximum plddt
        self.plddt_dict = plot_plddt(plot_path=self.plddt_plot_path, ranked_models_dict=ranked_models_dict)
        max_plddt = max(self.plddt_dict.values())
        self.best_ranked = utils.get_key_for_value(max_plddt, self.plddt_dict)
        bioutils.write_sequence(sequence_name=utils.get_file_name(self.sequence_path), sequence_amino=sequence_assembled.sequence_assembled, sequence_path=self.sequence_path)

        # Filter rankeds, split them in chains.
        mappings = {}
        ranked_nonsplit_filtered = {}
        for ranked, ranked_path in ranked_models_dict.items():
            if self.plddt_dict[ranked] >= (PERCENTAGE_FILTER * max_plddt):
                ranked_nonsplit_filtered[ranked] = ranked_path
                self.ranked_filtered[ranked] = ranked_path
                new_pdb_path = os.path.join(self.output_dir, os.path.basename(ranked_path))
            else:
                new_pdb_path = os.path.join(self.run_dir, f'split_{os.path.basename(ranked_path)}')

            shutil.copy2(ranked_path, new_pdb_path)
            mapping = bioutils.split_chains_assembly(pdb_in_path=new_pdb_path,
                                                    pdb_out_path=new_pdb_path,
                                                    sequence_assembled=sequence_assembled)
            mappings[ranked] = mapping
            ranked_models_dict[ranked] = new_pdb_path

        # Save superpositions of rankeds and templates
        reference_superpose = utils.sort_by_digit(list(ranked_models_dict.values()))[0]
        for path in utils.sort_by_digit(list(ranked_models_dict.values()))[1:] + list(template_dict.values()):
            bioutils.superpose_pdbs([path, reference_superpose], path)

        # Superpose each template with all the rankeds.
        if template_dict:
            for ranked, ranked_path in utils.sort_by_digit(ranked_models_dict):
                for template, template_path in template_dict.items():
                    if not template in self.rmsd_dict:
                        self.rmsd_dict[template] = {}

                    res_list_length = len(
                        [res for res in Selection.unfold_entities(PDBParser().get_structure(template, template_path), 'R')])
                    rmsd, nalign, quality_q = bioutils.superpose_pdbs([template_path, ranked_path])
                    self.rmsd_dict[template][ranked] = {'rmsd': rmsd, 'nalign': nalign, 'res_list_length': res_list_length}


        # Use aleph to generate domains and calculate secondary structure percentage
        interfaces_dict = {}
        for ranked, ranked_path in utils.sort_by_digit(ranked_models_dict):
            aleph_file = os.path.join(self.run_dir, f'aleph_{ranked}.txt')
            with open(aleph_file, 'w') as sys.stdout:
                ALEPH.annotate_pdb_model(reference=ranked_path, strictness_ah=0.45, strictness_bs=0.2, peptide_length=3,
                                        width_pic=1, height_pic=1, write_graphml=False, write_pdb=True)
            sys.stdout = sys.__stdout__
            if os.path.exists(self.aleph_results_path):
                self.secondary_dict[ranked] = utils.parse_aleph_annotate(file_path=self.aleph_results_path)
                aleph_txt_path = f'{self.run_dir}/aleph_{ranked}.txt'
                domains_dict = utils.parse_aleph_ss(aleph_txt_path)
            else:
                break
            if ranked in self.ranked_filtered and os.path.exists(aleph_txt_path):
                # print(bioutils.run_openmm(ranked_path))
                interfaces_data_list = bioutils.find_interface_from_pisa(ranked_path, self.interfaces_path)
                if interfaces_data_list:
                    interfaces_dict[ranked] = []
                    deltas_list = [interface['deltaG'] for interface in interfaces_data_list]
                    deltas_list = utils.normalize_list([deltas_list])
                    for i, interface in enumerate(interfaces_data_list):
                        code = f'{utils.get_file_name(ranked_path)}_{interface["chain1"]}{interface["chain2"]}'
                        dimers_path = os.path.join(self.interfaces_path, f'{code}.pdb')
                        interface['bfactor'] = deltas_list[i]
                        if not ((float(interface['se_gain1']) < 0) and (float(interface['se_gain2']) < 0)):
                            interface['bfactor'] = abs(max(deltas_list)) * 2
                        extended_res_dict = bioutils.create_interface_domain(pdb_in_path=ranked_path,
                                                                            pdb_out_path=dimers_path,
                                                                            interface=interface,
                                                                            domains_dict=domains_dict)
                        renum_residues_list = []
                        renum_residues_list.extend(utils.renum_residues(extended_res_dict[interface['chain1']],
                                                                        mapping=mappings[ranked][interface['chain1']]))
                        renum_residues_list.extend(utils.renum_residues(extended_res_dict[interface['chain2']],
                                                                        mapping=mappings[ranked][interface['chain2']]))
                        inter_name = structures.InterfaceName(name=code, res_list=renum_residues_list)
                        interfaces_dict[ranked].append(inter_name)

        # Superpose the experimental pdb with all the rankeds and templates
        if experimental_pdb is not None:
            for ranked, ranked_path in ranked_models_dict.items():
                rmsd, nalign, quality_q = bioutils.superpose_pdbs([experimental_pdb, ranked_path])
                self.experimental_dict[ranked] = rmsd
            for template, template_path in template_dict.items():
                rmsd, nalign, quality_q = bioutils.superpose_pdbs([experimental_pdb, template_path])
                self.experimental_dict[template] = rmsd
            output_pdb = os.path.join(self.output_dir, os.path.basename(experimental_pdb))
            bioutils.superpose_pdbs([experimental_pdb, reference_superpose], output_pdb)

        for template, template_path in template_nonsplit.items():
            frobenius_file = os.path.join(self.frobenius_path, f'frobenius_{template}.txt')
            matrices = os.path.join(self.run_dir, 'matrices')
            template_matrix = os.path.join(matrices, f'{utils.get_file_name(template_path)}_ang.npy')
            with open(frobenius_file, 'w') as sys.stdout:
                _, _, _, _, _, list_plot_ang, list_plot_dist = ALEPH.frobenius(references=[template_path], targets=list(ranked_nonsplit_filtered.values()), write_plot=True, write_matrix=True)
            sys.stdout = sys.__stdout__
            self.frobenius_plots['general'] = [shutil.copy2(plot, self.frobenius_path) for plot in (list_plot_ang+list_plot_dist)]
            for ranked, ranked_path in utils.sort_by_digit(ranked_nonsplit_filtered):
                ranked_matrix = os.path.join(matrices, f'{ranked}_ang.npy')
                interfaces_plots_list = []
                for interface_list in interfaces_dict[ranked]:
                    frobenius_file = os.path.join(self.frobenius_path, f'frobenius_{interface_list.name}.txt')
                    with open(frobenius_file, 'w') as sys.stdout:
                        _,_, plot = ALEPH.frobenius_submatrices(path_ref=template_matrix, path_tar=ranked_matrix, residues_tar=interface_list.res_list, write_plot=True, title=f'Interface: {interface_list.name}')
                    sys.stdout = sys.__stdout__
                    new_name = os.path.join(self.frobenius_path, f'{interface_list.name}.png')
                    plot_path = os.path.join(self.run_dir, os.path.basename(plot))
                    interfaces_plots_list.append(shutil.copy2(plot_path, new_name))
                self.frobenius_plots[ranked] = interfaces_plots_list
