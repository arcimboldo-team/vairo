import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
import re
import statistics
from libs import bioutils, utils
from typing import Dict, List

MATPLOTLIB_FONT = 14
plt.set_loglevel('WARNING')


def plot_ramachandran(plot_path: str, phi_psi_angles: List[List[float]]):
    fig, ax = plt.subplots(figsize=(8, 8))
    phi = [x[0] for x in phi_psi_angles]
    psi = [x[1] for x in phi_psi_angles]
    ax.plot(phi, psi, 'o', markersize=3, alpha=0.5)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xlabel(r'$\phi$')
    ax.set_ylabel(r'$\psi$')
    ax.set_title(f'Ramachandran plot for {utils.get_file_name(plot_path)}')
    plt.savefig(plot_path, dpi=100, bbox_inches='tight')
    plt.cla()


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
    text = []
    markers = ['.', '*', 's', 'P']
    for i, cluster in enumerate(clusters):
        text_cluster = f'Cluster {i}:'
        for path in cluster:
            name = utils.get_file_name(path)
            params = analysis_dict[name]
            if name.startswith('cluster_'):
                color = 'red'
            else:
                color = 'blue'
            plt.scatter(params.coord[0], params.coord[1], marker=markers[i], color=color, label=f'Cluster {i}')
            plt.annotate(name, (params.coord[0], params.coord[1]), horizontalalignment='right',
                         verticalalignment='top')

            if len(text_cluster) < 60:
                text_cluster += f' {name},'
            else:
                text.append(text_cluster)
                text_cluster = f' {name},'
        text.append(text_cluster[:-1] + '\n')

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    plt.xlim([-1, 1])
    plt.ylim([-1, 1])
    if predictions:
        plt.title('TEMPLATES AND PREDICTIONS CLUSTERING')
    else:
        plt.title('TEMPLATES CLUSTERING')
    plt.figtext(0.05, -0.15, '\n'.join(text))
    plt.savefig(plot_path, dpi=100, bbox_inches='tight')
    plt.cla()


def plot_sequence(plot_path: str, a_air):
    plt.rcParams.update({'font.size': MATPLOTLIB_FONT})
    fig, ax = plt.subplots(1, figsize=(16, 0.5))
    legend_seq = [Patch(label='Sequence', color='tab:cyan'), Patch(label='Linker', color='tab:blue')]
    ax.legend(handles=legend_seq[::-1], loc='upper center', bbox_to_anchor=(0.5, -0.4), fancybox=True, framealpha=0.5,
              ncol=2)
    for i in range(a_air.sequence_assembled.total_copies):
        ax.barh('sequence', a_air.sequence_assembled.get_sequence_length(i),
                left=a_air.sequence_assembled.get_starting_length(i) + 1, color='tab:cyan')
        ax.barh('sequence', a_air.sequence_assembled.glycines,
                left=a_air.sequence_assembled.get_starting_length(i) + a_air.sequence_assembled.get_sequence_length(
                    i) + 1, color='tab:blue')

        xcenters = (a_air.sequence_assembled.get_starting_length(i) + 1) + a_air.sequence_assembled.get_sequence_length(
            i) / 2
        ax.text(xcenters, 0, a_air.sequence_assembled.get_sequence_name(i), ha='center', va='center')

    ax_secondary = ax.secondary_xaxis('top')
    ax_secondary.set_xticks(
        ticks=[a_air.sequence_assembled.get_starting_length(i) + 1 for i in
               range(a_air.sequence_assembled.total_copies)],
        rotation=45)
    ax_secondary.set_xticks(
        ticks=list(ax_secondary.get_xticks()) + [
            a_air.sequence_assembled.get_starting_length(i) + a_air.sequence_assembled.get_sequence_length(i) + 1 for i
            in range(a_air.sequence_assembled.total_copies)],
        rotation=45)
    ax_secondary.set_xticklabels(
        labels=[1] * a_air.sequence_assembled.total_copies + [a_air.sequence_assembled.get_sequence_length(i) + 1 for i
                                                              in range(a_air.sequence_assembled.total_copies)],
        rotation=45)
    ax.set_xticks(
        ticks=[a_air.sequence_assembled.get_starting_length(i) + 1 for i in
               range(a_air.sequence_assembled.total_copies)],
        rotation=45)
    ax.set_xticks(
        ticks=list(ax.get_xticks()) + [
            a_air.sequence_assembled.get_starting_length(i) + a_air.sequence_assembled.get_sequence_length(i) + 1 for i
            in range(a_air.sequence_assembled.total_copies)],
        rotation=45)
    ax.set_xticklabels(labels=ax.get_xticks(), rotation=45)
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
    legend_seq = [Patch(label='Sequence', color='tab:cyan'), Patch(label='Linker', color='tab:blue')]
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
        add_sequences = [0] * len(a_air.sequence_assembled.sequence_assembled)
        for name in names:
            if plot_type == 'msa':
                features_search = a_air.feature.get_msa_by_name(name)
            else:
                features_search = a_air.feature.get_sequence_by_name(name)
            aligned_sequence = bioutils.compare_sequences(a_air.sequence_assembled.sequence_assembled, features_search)
            aligned_sequence = [1 if align == '-' else float(align) for align in aligned_sequence]
            add_sequences = np.add(aligned_sequence, add_sequences)
        add_sequences = [aligned / len(names) for aligned in add_sequences]
        for i in range(1, len(add_sequences)):
            ax.barh('Percentage', 1, left=i, height=0.5, color=str(add_sequences[i - 1]))
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
            if template is not None:
                changed_residues, changed_fasta, _ = template.get_changes()
                changed_residues = bioutils.convert_residues(changed_residues, a_air.sequence_assembled)
                changed_fasta = bioutils.convert_residues(changed_fasta, a_air.sequence_assembled)
                if len(name) > 4:
                    template_name = f'T{j + 1}'
                    text = f'\n{template_name} ({name}):'
                else:
                    template_name = name
                    text = f'\n{template_name}:'
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
                        legend_elements.append(
                            f'{template_name}: Aligned={match_split[5]}({match_split[8].replace("(", "").replace(")", "")}) Evalue={match_split[2]}')
                    else:
                        legend_elements.append(f'{template_name}')

            if features_search is not None:
                aligned_sequence = bioutils.compare_sequences(a_air.sequence_assembled.sequence_assembled,
                                                              features_search)
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
        fig.set_size_inches(16, number_of_templates * 0.7)
    else:
        fig.set_size_inches(16, number_of_templates)
    legend_elements.append(
        'The templates gray scale shows the similarity between the aligned template sequence and the input sequence.\n'
        'The darker parts indicate that the residues are the same or belong to the same group.\n')
    legend_elements.append('Yellow shows which residues have been changed to another specific residue\n'
                           'whereas the red shows which residues have been changed from another query sequence.\n'
                           'No information (white) implies that no modifications have been done.')
    legend_elements.reverse()
    plt.figtext(0.05, -0.05, '\n'.join(legend_elements), va='top')

    ax.set_xticks(
        [a_air.sequence_assembled.get_starting_length(i) + 1 for i in range(a_air.sequence_assembled.total_copies)])
    ax.set_xticks(list(ax.get_xticks()) + [
        a_air.sequence_assembled.get_starting_length(i) + a_air.sequence_assembled.get_sequence_length(i) + 1 for i in
        range(a_air.sequence_assembled.total_copies)])

    cut_chunk = [list(tup) for tup in a_air.chunk_list]
    cut_chunk = utils.remove_list_layer(cut_chunk)
    ax.set_xticks(list(ax.get_xticks()) + [cut + 1 for cut in cut_chunk])
    ax.set_xticklabels(ax.get_xticks(), rotation=45)

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
