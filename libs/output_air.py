import logging
import os
import re
import shutil
import sys
import statistics
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from typing import Dict, List
from Bio.PDB import PDBParser, Selection
from ALEPH.aleph.core import ALEPH
from libs import bioutils, features, utils, sequence, structures

PERCENTAGE_FILTER = 0.8
PERCENTAGE_BEST = 0.9
GROUPS = ['GAVLI', 'FYW', 'CM', 'ST', 'KRH', 'DENQ', 'P']
plt.set_loglevel('WARNING')


def read_rankeds(input_path: str):
    ranked_paths = [path for path in os.listdir(input_path) if re.match('ranked_[0-9]+.pdb', path)]
    return [structures.Ranked(os.path.join(input_path, path)) for path in utils.sort_by_digit(ranked_paths)]


def plot_plddt(plot_path: str, ranked_list: List) -> float:
    plt.clf()
    for ranked in ranked_list:
        plddt_list = []
        with open(ranked.path) as f:
            for line in f.readlines():
                if line[:4] == 'ATOM' and line[13:16] == 'CA ':
                    plddt_list.append(float(line[60:66].replace(" ", "")))
        res_list = [int(item) for item in range(1, len(plddt_list) + 1)]
        ranked.set_plddt(round(statistics.median(map(float, plddt_list)), 2))
        plt.plot(res_list, plddt_list, label=ranked.name)
    plt.legend()
    plt.xlabel('residue number')
    plt.ylabel('pLDDT')
    plt.savefig(plot_path)

    max_plddt = max([ranked.plddt for ranked in ranked_list])
    return max_plddt


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
    fig, (ax, ax1) = plt.subplots(2, figsize=(16, 6), gridspec_kw={'height_ratios': [6, 1]})
    total_length = len(a_air.sequence_assembled.sequence_assembled) + a_air.sequence_assembled.glycines
    height = 0.3
    for i in range(a_air.sequence_assembled.total_copies):
        ax.barh('sequence', a_air.sequence_assembled.get_sequence_length(i), height=height,
                left=a_air.sequence_assembled.get_starting_length(i) + 1, color='tab:cyan')
        ax.barh('sequence', a_air.sequence_assembled.glycines, height=height,
                left=a_air.sequence_assembled.get_starting_length(i) + a_air.sequence_assembled.get_sequence_length(
                    i) + 1,
                color='tab:blue')

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
    legend_elements = []

    names = [name for name in names if name != '']
    if len(names) > 30:
        aligned_sequences = [0] * len(a_air.sequence_assembled.sequence_assembled)
        for name in names:
            features_search = a_air.feature.get_msa_by_name(name)
            aligned_sequence = compare_sequences(a_air.sequence_assembled.sequence_assembled, features_search)
            new_array = [0 if align == '0' else 1 for align in aligned_sequence]
            aligned_sequences = np.add(new_array, aligned_sequences)
        aligned_sequences = [aligned / len(names) for aligned in aligned_sequences]
        for i in range(1, len(aligned_sequences)):
            ax.barh('Percentage', 1, height=height, left=i, color=str(aligned_sequences[i - 1]))

    else:
        pdb_hits_path = os.path.join(a_air.run_dir, 'msas/pdb_hits.hhr')
        hhr_text = ''
        if os.path.exists(pdb_hits_path):
            hhr_text = open(pdb_hits_path, 'r').read()
        for j, name in reversed(list(enumerate(names))):
            template = a_air.get_template_by_id(name)
            if template is not None:
                results_alignment_text = template.get_results_alignment_text()
                if len(name) > 4:
                    template_name = f'T{j + 1}'
                    legend_elements.append(Patch(label=f'{template_name} ({name}): {results_alignment_text}'))
                else:
                    template_name = name
                    legend_elements.append(Patch(label=f'{template_name}: {results_alignment_text}'))
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
                        legend_elements.append(Patch(label=f'{template_name}: Aligned={match_split[5]}({match_split[8].replace("(","").replace(")","")}) Evalue={match_split[2]}'))
                    else:
                        legend_elements.append(Patch(label=f'{template_name}'))


            if features_search is not None:
                aligned_sequence = compare_sequences(a_air.sequence_assembled.sequence_assembled, features_search)
                for i in range(1, len(features_search)):
                    if aligned_sequence[i - 1] != '-':
                        ax.barh(template_name, 1, height=height, left=i, color=str(aligned_sequence[i - 1]))

    ax.xaxis.grid(color='k', linestyle='dashed', alpha=0.4, which='both')
    plt.setp([ax.get_xticklines()], color='k')
    ax.set_xlim(0, total_length)
    ax.set_xticks(
        [a_air.sequence_assembled.get_starting_length(i) + 1 for i in range(a_air.sequence_assembled.total_copies)])

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
        self.ranked_list: List[structures.Ranked] = []
        self.output_dir: str = output_dir
        self.experimental_dict = {}
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

    def analyse_output(self, sequence_assembled: sequence.SequenceAssembled, feature: features.Features,
                       experimental_pdb: str):
        # Read all templates and rankeds, if there are no ranked, raise an error
        template_dict = {}
        template_nonsplit = {}

        if feature is not None:
            template_nonsplit = feature.write_all_templates_in_features(output_dir=self.nonsplit_path,
                                                                        print_number=False)

        # Split the templates with chains
        for template, template_path in template_nonsplit.items():
            new_pdb_path = os.path.join(self.templates_path, f'{template}.pdb')
            shutil.copy2(template_path, new_pdb_path)
            template_dict[template] = new_pdb_path
            bioutils.split_chains_assembly(pdb_in_path=new_pdb_path,
                                           pdb_out_path=new_pdb_path,
                                           sequence_assembled=sequence_assembled)

        self.ranked_list = read_rankeds(input_path=self.run_dir)

        if not self.ranked_list:
            logging.info('No ranked PDBs found')
            return

        # Create a plot with the ranked pLDDTs, also, calculate the maximum pLDDT
        max_plddt = plot_plddt(plot_path=self.plddt_plot_path, ranked_list=self.ranked_list)
        bioutils.write_sequence(sequence_name=utils.get_file_name(self.sequence_path),
                                sequence_amino=sequence_assembled.sequence_assembled, sequence_path=self.sequence_path)

        # Sort list of ranked by pLDDT
        self.ranked_list.sort(key=lambda x: x.plddt, reverse=True)

        # Filter rankeds, split them in chains.
        for ranked in self.ranked_list:
            if ranked.plddt >= (PERCENTAGE_FILTER * max_plddt):
                ranked.set_filtered(True)
                ranked.set_split_path(os.path.join(self.output_dir, os.path.basename(ranked.path)))
                if ranked.plddt >= (PERCENTAGE_BEST * max_plddt):
                    ranked.set_best(True)
            else:
                ranked.set_split_path(os.path.join(self.run_dir, f'split_{os.path.basename(ranked.path)}'))

            shutil.copy2(ranked.path, ranked.split_path)
            mapping = bioutils.split_chains_assembly(pdb_in_path=ranked.split_path,
                                                     pdb_out_path=ranked.split_path,
                                                     sequence_assembled=sequence_assembled)
            ranked.set_mapping(mapping)

        # Save superpositions of rankeds and templates
        reference_superpose = self.ranked_list[0].split_path

        for path in [ranked.split_path for ranked in self.ranked_list[1:]] + list(template_dict.values()):
            bioutils.superpose_pdbs([path, reference_superpose], path)

        # Superpose each template with all the rankeds.
        if template_dict:
            for ranked in self.ranked_list:
                for template, template_path in template_dict.items():
                    total_residues = len(
                        [res for res in
                         Selection.unfold_entities(PDBParser().get_structure(template, template_path), 'R')])
                    rmsd, aligned_residues, quality_q = bioutils.superpose_pdbs([template_path, ranked.split_path])
                    if rmsd is not None:
                        rmsd = round(rmsd, 2)
                    ranked.add_template(structures.TemplateRanked(template, rmsd, aligned_residues, total_residues))
                ranked.sort_template_rankeds()

        # Use aleph to generate domains and calculate secondary structure percentage
        for ranked in self.ranked_list:
            aleph_file = os.path.join(self.run_dir, f'aleph_{ranked.name}.txt')
            with open(aleph_file, 'w') as sys.stdout:
                try:
                    ALEPH.annotate_pdb_model(reference=ranked.split_path, strictness_ah=0.45, strictness_bs=0.2,
                                            peptide_length=3,
                                            width_pic=1, height_pic=1, write_graphml=False, write_pdb=True)
                except:
                    pass
            sys.stdout = sys.__stdout__
            if os.path.exists(self.aleph_results_path):
                result_dict = utils.parse_aleph_annotate(file_path=self.aleph_results_path)
                ranked.set_secondary_structure(ah=result_dict['ah'], bs=result_dict['bs'],
                                               total_residues=result_dict['number_total_residues'])
                aleph_txt_path = f'{self.run_dir}/aleph_{ranked.name}.txt'
                domains_dict = utils.parse_aleph_ss(aleph_txt_path)
            else:
                break
            if ranked.filtered and os.path.exists(aleph_txt_path):
                ranked.set_minimized_path(os.path.join(self.run_dir, f'{ranked.name}_minimized.pdb'))
                #ranked.set_energies(bioutils.run_openmm(pdb_in_path=ranked.path, pdb_out_path=ranked.minimized_path))
                interfaces_data_list = bioutils.find_interface_from_pisa(ranked.split_path, self.interfaces_path)
                if interfaces_data_list:
                    deltas_list = [interface['deltaG'] for interface in interfaces_data_list]
                    deltas_list = utils.normalize_list([deltas_list])
                    for i, interface in enumerate(interfaces_data_list):
                        code = f'{utils.get_file_name(ranked.split_path)}_{interface["chain1"]}{interface["chain2"]}'
                        dimers_path = os.path.join(self.interfaces_path, f'{code}.pdb')
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
                        ranked.add_interface(structures.Interface(name=code, res_list=renum_residues_list))

        # Superpose the experimental pdb with all the rankeds and templates
        if experimental_pdb is not None:
            for ranked in self.ranked_list:
                rmsd, nalign, quality_q = bioutils.superpose_pdbs([experimental_pdb, ranked.split_path])
                self.experimental_dict[ranked.name] = round(rmsd, 2)
            for template, template_path in template_dict.items():
                rmsd, nalign, quality_q = bioutils.superpose_pdbs([experimental_pdb, template_path])
                self.experimental_dict[template] = round(rmsd, 2)
            output_pdb = os.path.join(self.output_dir, os.path.basename(experimental_pdb))
            bioutils.superpose_pdbs([experimental_pdb, reference_superpose], output_pdb)

        for template, template_path in template_nonsplit.items():
            frobenius_file = os.path.join(self.frobenius_path, f'frobenius_{template}.txt')
            matrices = os.path.join(self.run_dir, 'matrices')
            template_matrix = os.path.join(matrices, f'{utils.get_file_name(template_path)}_ang.npy')
            path_list = [ranked.path for ranked in self.ranked_list if ranked.filtered]
            with open(frobenius_file, 'w') as sys.stdout:
                _, _, _, _, _, list_plot_ang, list_plot_dist = ALEPH.frobenius(references=[template_path],
                                                                               targets=list(path_list), write_plot=True,
                                                                               write_matrix=True)
            sys.stdout = sys.__stdout__
            frobenius_file_dict = utils.parse_frobenius(frobenius_file)
            ranked_filtered = [ranked for ranked in self.ranked_list if ranked.filtered]
            for ranked in ranked_filtered:
                ranked.add_frobenius_plot(
                    template=template, 
                    dist_plot=[shutil.copy2(plot, self.frobenius_path) for plot in list_plot_dist if ranked.name in os.path.basename(plot)][0],
                    ang_plot=[shutil.copy2(plot, self.frobenius_path) for plot in list_plot_ang if ranked.name in os.path.basename(plot)][0], 
                    dist_coverage=frobenius_file_dict['dist'][ranked.name], 
                    ang_coverage=frobenius_file_dict['ang'][ranked.name], 
                )
                ranked_matrix = os.path.join(matrices, f'{ranked.name}_ang.npy')
                for interface in ranked.interfaces:
                    frobenius_file = os.path.join(self.frobenius_path, f'frobenius_{interface.name}.txt')
                    with open(frobenius_file, 'w') as sys.stdout:
                        fro_distance, fro_core, plot = ALEPH.frobenius_submatrices(path_ref=template_matrix, path_tar=ranked_matrix,
                                                                 residues_tar=interface.res_list, write_plot=True,
                                                                 title=f'Interface: {interface.name}')
                    sys.stdout = sys.__stdout__
                    new_name = os.path.join(self.frobenius_path, f'{interface.name}.png')
                    plot_path = os.path.join(self.run_dir, os.path.basename(plot))

                    interface.add_frobenius_information(
                        dist_coverage = fro_distance,
                        core = fro_core,
                        dist_plot = shutil.copy2(plot_path, new_name)
                    )


    def write_tables(self, rmsd_dict: Dict, secondary_dict: Dict, plddt_dict: Dict, energies_dict: Dict):
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
