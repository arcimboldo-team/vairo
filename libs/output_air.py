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
from libs import bioutils, features, structure_air, utils, sequence, structures

PERCENTAGE_FILTER = 0.8
GROUPS = ['GAVLI', 'FYW', 'CM', 'ST', 'KRH', 'DENQ', 'P']
plt.set_loglevel('WARNING')



def plot_plddt(plot_path: str, ranked_models_dict: Dict) -> Dict:
    return_plddt_dict = {}

    plt.clf()
    for ranked, ranked_path in ranked_models_dict.items():
        plddt_list = []
        with open(ranked_path) as f:
            for line in f.readlines():
                if line[:4] == 'ATOM' and line[13:16] == 'CA ':
                    plddt_list.append(float(line[60:66].replace(" ", "")))
        res_list = [int(item) for item in range(1, len(plddt_list) + 1)]
        return_plddt_dict[ranked] = statistics.median(map(float, plddt_list))
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

def plot_gantt(plot_type: str, plot_path: str, sequence_assembled: sequence.SequenceAssembled,
               feature: features.Features, a_air: structure_air.StructureAir):
    fig , (ax, ax1) = plt.subplots(2, figsize=(16, 6), gridspec_kw={'height_ratios':[6, 1]})
    total_length = len(sequence_assembled.sequence_assembled) + sequence_assembled.glycines
    height = 0.3
    for i in range(sequence_assembled.total_copies):
        ax.barh('sequence', sequence_assembled.get_sequence_length(i), height=height,
                left=sequence_assembled.get_starting_length(i) + 1, color='tab:cyan')
        ax.barh('sequence', sequence_assembled.glycines, height=height,
                left=sequence_assembled.get_starting_length(i) + sequence_assembled.get_sequence_length(i) + 1,
                color='tab:blue')

    if plot_type == 'msa':
        title = 'MSA'
        file = os.path.join(plot_path, 'msa_gantt.png')
    else:  # should be template:
        title = 'TEMPLATES'
        file = os.path.join(plot_path, 'template_gantt.png')

    names = feature.get_names()
    legend_elements = []
    for j, name in reversed(list(enumerate(names))):
        template = a_air.get_template_by_id(name)
        results_alignment_text = template.get_results_alignment_text()
        template_name = f'T{j+1}'
        legend_elements.append(Patch(label=f'{template_name} ({name}): {results_alignment_text}'))
        if plot_type == 'msa':
            features_search = feature.get_msa_by_name(name)
        else:  # should be template:
            features_search = feature.get_sequence_by_name(name)
        if features_search is not None:
            aligned_sequence = compare_sequences(sequence_assembled.sequence_assembled, features_search)
            for i in range(1, len(features_search)):
                if aligned_sequence[i - 1] != '-':
                    ax.barh(template_name, 1, height=height, left=i, color=str(aligned_sequence[i - 1]))
    ax.xaxis.grid(color='k', linestyle='dashed', alpha=0.4, which='both')
    plt.setp([ax.get_xticklines()], color='k')
    ax.set_xlim(0, total_length)
    ax.set_xticks([sequence_assembled.get_starting_length(i) + 1 for i in range(sequence_assembled.total_copies)])
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

class OutputAir:

    def __init__(self, output_dir, run_dir):
        self.plots_path = f'{output_dir}/plots'
        self.frobenius_path = f'{output_dir}/frobenius'
        self.sequence_path = os.path.join(self.frobenius_path, 'sequence.fasta')
        self.templates_path = f'{output_dir}/templates'
        self.interfaces_path = f'{output_dir}/interfaces'
        self.analysis_path = f'{self.plots_path}/analysis.txt'
        self.aleph_results_path = f'{run_dir}/output.json'
        self.plddt_plot_path = f'{self.plots_path}/plddt.png'
        self.nonsplit_path = f'{run_dir}/nonsplit'
        self.html_path = f'{output_dir}/output.html'
        self.run_dir = run_dir
        self.output_dir=output_dir

        utils.create_dir(dir_path=self.plots_path, delete_if_exists=True)
        utils.create_dir(dir_path=self.templates_path, delete_if_exists=True)
        utils.create_dir(dir_path=self.interfaces_path, delete_if_exists=True)
        utils.create_dir(dir_path=self.frobenius_path, delete_if_exists=True)
        utils.create_dir(dir_path=self.nonsplit_path, delete_if_exists=True)


    def create_plot_gantt(self, sequence_assembled: sequence.SequenceAssembled, feature: features.Features, a_air: structure_air.StructureAir):
        
        plot_gantt(plot_type='template', plot_path=self.plots_path, sequence_assembled=sequence_assembled, feature=feature, a_air=a_air)
        plot_gantt(plot_type='msa', plot_path=self.plots_path, sequence_assembled=sequence_assembled, feature=feature, a_air=a_air)


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
        plddt_dict = plot_plddt(plot_path=self.plddt_plot_path, ranked_models_dict=ranked_models_dict)
        max_plddt = max(plddt_dict.values())
        bioutils.write_sequence(sequence_name=utils.get_file_name(self.sequence_path), sequence_amino=sequence_assembled.sequence_assembled, sequence_path=self.sequence_path)

        # Filter rankeds, split them in chains.
        ranked_filtered = []
        mappings = {}
        ranked_nonsplit_filtered = {}
        for ranked, ranked_path in ranked_models_dict.items():
            if plddt_dict[ranked] >= (PERCENTAGE_FILTER * max_plddt):
                ranked_nonsplit_filtered[ranked] = ranked_path
                ranked_filtered.append(ranked)
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
        rmsd_dict = {}
        if template_dict:
            for ranked, ranked_path in utils.sort_by_digit(ranked_models_dict):
                for template, template_path in template_dict.items():
                    res_list_length = len(
                        [res for res in Selection.unfold_entities(PDBParser().get_structure(template, template_path), 'R')])
                    rmsd, nalign, quality_q = bioutils.superpose_pdbs([template_path, ranked_path])
                    try:
                        rmsd_dict[f'{template}'].append(f'{rmsd}, {nalign} ({res_list_length})')
                    except KeyError:
                        rmsd_dict[f'{template}'] = [f'{rmsd}, {nalign} ({res_list_length})']

        # Use aleph to generate domains and calculate secondary structure percentage
        interfaces_dict = {}
        secondary_dict = {}
        openmm_dict = {}
        for ranked, ranked_path in utils.sort_by_digit(ranked_models_dict):
            aleph_file = os.path.join(self.run_dir, f'aleph_{ranked}.txt')
            with open(aleph_file, 'w') as sys.stdout:
                ALEPH.annotate_pdb_model(reference=ranked_path, strictness_ah=0.45, strictness_bs=0.2, peptide_length=3,
                                        width_pic=1, height_pic=1, write_graphml=False, write_pdb=True)
            sys.stdout = sys.__stdout__
            if os.path.exists(self.aleph_results_path):
                secondary_dict[ranked] = utils.parse_aleph_annotate(file_path=self.aleph_results_path)
                aleph_txt_path = f'{self.run_dir}/aleph_{ranked}.txt'
                domains_dict = utils.parse_aleph_ss(aleph_txt_path)
            else:
                break
            if ranked in ranked_filtered and os.path.exists(aleph_txt_path):
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
        experimental_dict = {}
        if experimental_pdb is not None:
            for ranked, ranked_path in ranked_models_dict.items():
                rmsd, nalign, quality_q = bioutils.superpose_pdbs([experimental_pdb, ranked_path])
                experimental_dict[ranked] = rmsd
            for template, template_path in template_dict.items():
                rmsd, nalign, quality_q = bioutils.superpose_pdbs([experimental_pdb, template_path])
                experimental_dict[template] = rmsd
            output_pdb = os.path.join(self.output_dir, os.path.basename(experimental_pdb))
            bioutils.superpose_pdbs([experimental_pdb, reference_superpose], output_pdb)

        for template, template_path in template_nonsplit.items():
            frobenius_file = os.path.join(self.frobenius_path, f'frobenius_{template}.txt')
            matrices = os.path.join(self.run_dir, 'matrices')
            template_matrix = os.path.join(matrices, f'{utils.get_file_name(template_path)}_ang.npy')
            with open(frobenius_file, 'w') as sys.stdout:
                _, _, _, _, _, list_plot_ang, list_plot_dist = ALEPH.frobenius(references=[template_path], targets=list(ranked_nonsplit_filtered.values()), write_plot=True, write_matrix=True)
            sys.stdout = sys.__stdout__
            [shutil.copy2(plot, self.frobenius_path) for plot in (list_plot_ang+list_plot_dist)]
            for ranked, ranked_path in ranked_nonsplit_filtered.items():
                ranked_matrix = os.path.join(matrices, f'{ranked}_ang.npy')
                for interface_list in interfaces_dict[ranked]:
                    frobenius_file = os.path.join(self.frobenius_path, f'frobenius_{interface_list.name}.txt')
                    with open(frobenius_file, 'w') as sys.stdout:
                        _,_, plot = ALEPH.frobenius_submatrices(path_ref=template_matrix, path_tar=ranked_matrix, residues_tar=interface_list.res_list, write_plot=True)
                    sys.stdout = sys.__stdout__
                    new_name = os.path.join(self.frobenius_path, f'{interface_list.name}.png')
                    plot_path = os.path.join(self.run_dir, os.path.basename(plot))
                    shutil.copy2(plot_path, new_name)

        with open(self.analysis_path, 'w') as f_in:

            if bool(rmsd_dict):
                f_in.write('\n')
                f_in.write('Superpositions of rankeds and templates\n')
                rows = []
                for key in rmsd_dict.keys():
                    rows.append([key] + rmsd_dict[key])
                df = pd.DataFrame(rows, columns=['template'] + utils.sort_by_digit(list(ranked_models_dict)))
                f_in.write(df.to_markdown())

            if bool(secondary_dict):
                f_in.write('\n\n')
                f_in.write('Secondary structure percentages calculated with ALEPH\n')
                rows = []
                for key in secondary_dict.keys():
                    rows.append([key] + list(secondary_dict[key].values()))
                df = pd.DataFrame(rows, columns=['ranked'] + list(secondary_dict[key]))
                f_in.write(df.to_markdown())

            if bool(openmm_dict):
                f_in.write('\n\n')
                f_in.write('OPENMM\n')
                rows = []
                for key in openmm_dict.keys():
                    rows.append([key] + list(openmm_dict[key].values()))
                df = pd.DataFrame(rows, columns=['ranked'] + list(openmm_dict[key]))
                f_in.write(df.to_markdown())

            if bool(plddt_dict):
                f_in.write('\n\n')
                f_in.write('rankeds PLDDT\n')
                rows = []
                for key in plddt_dict.keys():
                    rows.append([key, plddt_dict[key]])
                df = pd.DataFrame(rows, columns=['ranked', 'plddt'])
                f_in.write(df.to_markdown())

            if bool(experimental_dict):
                f_in.write('\n\n')
                f_in.write(f'Superposition with experimental structure {experimental_pdb}\n')
                rows = []
                for key in experimental_dict.keys():
                    rows.append([key, experimental_dict[key]])
                df = pd.DataFrame(rows, columns=['pdb', 'rmsd'])
                f_in.write(df.to_markdown())

            f_in.write('\n\n')
