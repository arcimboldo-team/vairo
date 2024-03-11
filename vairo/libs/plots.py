import os
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from libs import bioutils, utils
from typing import Dict, List

MATPLOTLIB_FONT = 14
plt.set_loglevel('WARNING')


def generate_minor_ticks(ax_list: List[int], step: int) -> List[int]:
    ax_list = sorted(ax_list)
    minor_ticks = []
    for i in range(0, len(ax_list) - 1, step):
        current_num = ax_list[i]
        next_num = ax_list[i + 1]
        step = (next_num - current_num) / 6  # Divide by 6 to get 5 numbers between each pair
        minor_ticks.extend([current_num + j * step for j in range(1, 6)])
    return minor_ticks


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


def plot_plddt(plot_path: str, ranked_list: List):
    plt.figure(figsize=(18, 6))
    plt.rcParams.update({'font.size': MATPLOTLIB_FONT})
    for ranked in ranked_list:
        return_dict = bioutils.read_bfactors_from_residues(pdb_path=ranked.path)
        plddt_list = [value for value in list(return_dict.values())[0] if value is not None]
        res_list = [int(item) for item in range(1, len(plddt_list) + 1)]
        plt.plot(res_list, plddt_list, label=ranked.name)
    plt.legend(loc='upper right')
    plt.xlabel('residue number')
    plt.ylabel('pLDDT')
    plt.ylim(0, 100)
    plt.savefig(plot_path, dpi=100)
    plt.cla()


def plot_cc_analysis(plot_path: str, analysis_dict: Dict, clusters: List, predictions: bool = False):
    plt.figure(figsize=(8, 8))
    plt.rcParams.update({'font.size': MATPLOTLIB_FONT})
    text = []
    markers = ['.', '*', 's', 'P']
    for i, cluster in enumerate(clusters):
        text_cluster = f'Cluster {i}:'
        for pdb in cluster:
            name = pdb.name
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

    if len(text) > 6:
        plt.figtext(0.05, -0.22, '\n'.join(text))
    else:
        plt.figtext(0.05, -0.15, '\n'.join(text))
    plt.savefig(plot_path, dpi=100, bbox_inches='tight')
    plt.cla()


def plot_sequence(plot_path: str, a_air):
    plt.rcParams.update({'font.size': MATPLOTLIB_FONT})
    color_seq = '#2e75b6'
    color_link = '#7f7f7f'
    fig, ax = plt.subplots(1, figsize=(16, 0.5))
    lines_leg = [Line2D([0], [0], color=color_seq, linewidth=3),
                 Line2D([0], [0], color=color_link, linewidth=3, linestyle='dashed')]
    lines_leg_lab = ['Query sequence', 'Linker']
    fig.legend(lines_leg, lines_leg_lab, loc='upper center', bbox_to_anchor=(0.5, -0.4), ncol=2, frameon=False)

    for i in range(a_air.sequence_assembled.total_copies):
        ax.barh('sequence', a_air.sequence_assembled.get_sequence_length(i),
                left=a_air.sequence_assembled.get_starting_length(i) + 1, color=color_seq)
        if i < a_air.sequence_assembled.total_copies - 1:
            for num in range(0, a_air.sequence_assembled.glycines, 8):
                ax.barh('sequence', 4,
                        left=a_air.sequence_assembled.get_finishing_length(i) + 2 + num, color=color_link, height=0.2,
                        zorder=2)

        xcenters = (a_air.sequence_assembled.get_starting_length(i) + 1) + a_air.sequence_assembled.get_sequence_length(
            i) / 2
        ax.text(xcenters, 0, a_air.sequence_assembled.get_sequence_name(i), ha='center', va='center', color='white')

    ax_secondary = ax.secondary_xaxis('top')
    ax_secondary.set_xticks(
        ticks=[a_air.sequence_assembled.get_starting_length(i) + 1 for i in
               range(a_air.sequence_assembled.total_copies - 1)], rotation=45)
    ax_secondary.set_xticks(
        ticks=list(ax_secondary.get_xticks()) + [a_air.sequence_assembled.get_finishing_length(i) + 2 for i
                                                 in range(a_air.sequence_assembled.total_copies)], rotation=45)
    ax_secondary.set_xticklabels(
        labels=[1] * a_air.sequence_assembled.total_copies + [a_air.sequence_assembled.get_sequence_length(i) + 2 for i
                                                              in range(a_air.sequence_assembled.total_copies - 1)],
        rotation=45)
    ax.set_xticks(
        ticks=[a_air.sequence_assembled.get_starting_length(i) + 1 for i in
               range(a_air.sequence_assembled.total_copies)],
        rotation=45)
    ax.set_xticks(
        ticks=list(ax.get_xticks()) + [a_air.sequence_assembled.get_finishing_length(i) + 2 for i
                                       in range(a_air.sequence_assembled.total_copies)],
        rotation=45)

    ax.tick_params('both', length=10, which='major')
    ax.tick_params('both', length=5, which='minor')

    ax.set_xticks(generate_minor_ticks(list(ax.get_xticks()), step=2), minor=True)
    ax.set_xticklabels(labels=ax.get_xticks(), rotation=45)
    ax.set_xlim(0, len(a_air.sequence_assembled.sequence_assembled))
    ax.set_yticks([])
    fig.tight_layout()
    fig.subplots_adjust(top=.95)
    plt.savefig(plot_path, bbox_inches='tight', dpi=100)
    plt.cla()


def plot_gantt(plot_type: str, plot_path: str, a_air, reduced: bool = False) -> str:
    plt.rcParams.update({'font.size': MATPLOTLIB_FONT})
    color_seq = '#2e75b6'
    color_link = '#7f7f7f'
    fig, ax = plt.subplots(1, figsize=(16, 2))

    ax1 = ax.twiny()
    ax1.xaxis.set_ticks_position('bottom')
    ax1.spines[['right', 'top', 'left']].set_visible(False)

    legend_elements = []

    number_of_templates = 1
    total_length = len(a_air.sequence_assembled.sequence_assembled)
    msa_found = False
    templates_found = False

    mutated_residues = a_air.sequence_assembled.get_mutated_residues_list()
    for i in mutated_residues:
        ax.barh('sequence', 1, left=i + 1, align='edge', color='yellow', height=0.35, zorder=3)

    for i in range(a_air.sequence_assembled.total_copies):
        ax.barh('sequence', a_air.sequence_assembled.get_sequence_length(i),
                left=a_air.sequence_assembled.get_starting_length(i) + 1, color=color_seq, height=0.7, zorder=2)
        ax.barh('sequence', a_air.sequence_assembled.get_sequence_length(i),
                left=a_air.sequence_assembled.get_starting_length(i) + 1, align='edge', color=color_seq, height=0.23,
                zorder=3)
        if i < a_air.sequence_assembled.total_copies - 1:
            for num in range(0, a_air.sequence_assembled.glycines, 8):
                ax.barh('sequence', 4,
                        left=a_air.sequence_assembled.get_finishing_length(i) + 2 + num, color=color_link, height=0.2,
                        zorder=2)

        xcenters = (a_air.sequence_assembled.get_starting_length(i) + 1) + a_air.sequence_assembled.get_sequence_length(
            i) / 2
        ax.text(xcenters, 0, a_air.sequence_assembled.get_sequence_name(i), fontsize='small', ha='center', va='center',
                color='white')

    if plot_type == 'msa':
        title = 'MSA'
        file = os.path.join(plot_path, 'msa_gantt.png')
    elif plot_type == 'templates':  # should be template:
        title = 'TEMPLATES'
        file = os.path.join(plot_path, 'template_gantt.png')
    else:
        title = 'TEMPLATES and ALIGNED SEQUENCES (MSA)'
        file = os.path.join(plot_path, 'template_msa_gantt.png')

    names = a_air.feature.get_names_msa()
    names = [name for name in names if name != '']
    if ((len(names) > 20 and plot_type == 'msa') or plot_type == 'both') and len(names) > 0:
        number_of_templates += 1
        sequences_msa = [a_air.feature.get_msa_by_name(name) for name in names]
        new_sequences = bioutils.calculate_coverage_scaled(query_seq=a_air.sequence_assembled.sequence_mutated_assembled, sequences=sequences_msa)

        if plot_type == 'both':
            name = 'MSA'
        else:
            name = 'Percentage'
        for i in range(len(new_sequences)):
            msa_found = True
            ax.barh(name, 1, left=i + 1, height=0.5, color=str(1-new_sequences[i]), zorder=2)

    if plot_type != 'msa' or len(names) <= 20:
        if plot_type != 'msa':
            names = a_air.feature.get_names_templates()
        if reduced:
            names_selected = [name for name in names if name in a_air.output.templates_selected]
        if not reduced or not names_selected:
            names_selected = names
        pdb_hits_path = os.path.join(a_air.results_dir, 'msas/pdb_hits.hhr')
        hhr_text = ''
        if os.path.exists(pdb_hits_path):
            hhr_text = open(pdb_hits_path, 'r').read()

        if reduced and len(names) > 20:
            number_of_templates += 1
            sequences_templates = [a_air.feature.get_sequence_by_name(name[1]) for name in reversed(list(enumerate(names))) if a_air.feature.get_sequence_by_name(name[1]) is not None]
            new_sequences = bioutils.calculate_coverage_scaled(query_seq=a_air.sequence_assembled.sequence_mutated_assembled, sequences=sequences_templates)
            add_sequences = [0] * len(a_air.sequence_assembled.sequence_assembled)
            for i in range(len(add_sequences)):
                ax.barh('Templates', 1, left=i + 1, height=0.5, color=str(1-new_sequences[i]), zorder=2)

        long_names = any([name for name in names_selected if len(name) > 7])
        for j, name in reversed(list(enumerate(names_selected))):
            templates_found = True
            number_of_templates += 1
            template = a_air.get_template_by_id(name)
            changed_residues = []
            changed_fasta = []

            if long_names:
                template_name = f"M{j + 1}" if plot_type == "msa" else f"T{j + 1}"
                text = f'\n{template_name} ({name})'
            else:
                template_name = name
                text = f'\n{template_name}'

            if template is not None:
                changed_residues, changed_fasta, _, _ = template.get_changes()
                changed_residues = bioutils.convert_residues(changed_residues, a_air.sequence_assembled)
                changed_fasta = bioutils.convert_residues(changed_fasta, a_air.sequence_assembled)
                chains = template.get_chain_by_position()
                for i, alignment in enumerate(template.get_results_alignment()):
                    if alignment is not None:
                        text += f'\n\tChain {alignment.chain}: Aligned={alignment.aligned_columns}({alignment.total_columns}) Evalue={alignment.evalue} Identities={alignment.identities}'
                    elif alignment is None and template.aligned and chains[i] is not None:
                        text += f'\n\tChain {chains[i]}: Prealigned'
                    else:
                        text += f'\n\tNot used'

            if plot_type == 'msa':
                features_search = a_air.feature.get_msa_by_name(name)
            else:
                features_search = a_air.feature.get_sequence_by_name(name)
                if hhr_text != '':
                    evalue, aligned, identity, total_residues = bioutils.parse_pdb_hits_hhr(hhr_text, name.upper())
                    if evalue is not None:
                        text += f' Aligned={aligned}({total_residues}) Evalue={evalue} Identity={identity}'

            legend_elements.append(text)

            if features_search is not None:
                aligned_sequence, _ = bioutils.compare_sequences(a_air.sequence_assembled.sequence_mutated_assembled,
                                                                 features_search)
                for i in range(len(features_search)):
                    if aligned_sequence[i] != 0:
                        aux_aligned = 1 - aligned_sequence[i] / 6
                        if i + 1 in changed_residues:
                            ax.barh(template_name, 1, left=i + 1, height=0.25, align='edge', color='yellow', zorder=3)
                        elif i + 1 in changed_fasta:
                            ax.barh(template_name, 1, left=i + 1, height=0.25, align='edge', color='red', zorder=3)
                        else:
                            ax.barh(template_name, 1, left=i + 1, height=0.25, align='edge', zorder=3,
                                    color=str(aux_aligned))
                        ax.barh(template_name, 1, left=i + 1, height=0.1, align='edge', zorder=3,
                                color=str(aux_aligned))
                        ax.barh(template_name, 1, left=i + 1, height=0.5, zorder=2, color=str(aux_aligned))

    if number_of_templates == 1:
        index = 2.1
    elif number_of_templates == 2:
        index = number_of_templates * 1.2
    elif number_of_templates == 3:
        index = number_of_templates * 1
    elif number_of_templates == 4:
        index = number_of_templates * 0.8
    elif number_of_templates == 5:
        index = number_of_templates * 0.75
    elif number_of_templates == 6:
        index = number_of_templates * 0.7
    else:
        index = number_of_templates * 0.5

    lines_leg = [Line2D([0], [0], color=color_seq, linewidth=3),
                 Line2D([0], [0], color=color_link, linewidth=3, linestyle='dashed')]
    lines_leg_lab = ['Query Sequence', 'Linker']
    fig.legend(lines_leg, lines_leg_lab, loc="lower left", bbox_to_anchor=(0.70, 0), ncol=2, frameon=False)

    fig.set_size_inches(16, index)
    plt.setp([ax.get_xticklines()], color='k')
    ax.set_xlim(-round(ax.get_xlim()[1] - total_length) / 6, total_length + round(ax.get_xlim()[1] - total_length) / 6)
    ax.set_ylim(-0.5, number_of_templates + 0.5)

    color_template = '#dcede9'
    color_msa = '#d7e9cb'
    color_sequence = '#fbf7e6'
    if number_of_templates == 1:
        ax.axhspan(-0.5, ax.get_ylim()[1], facecolor='0.5', zorder=1, color=color_sequence)
    else:
        if plot_type == 'both':
            ax.axhspan(-0.5, 0.5, facecolor='0.5', zorder=1, color=color_sequence)
            if msa_found:
                if templates_found:
                    ax.axhspan(0.5, 1.5, facecolor='0.5', zorder=1, color=color_msa)
                    ax.axhspan(1.5, ax.get_ylim()[1], facecolor='0.5', zorder=1, color=color_template)
                else:
                    ax.axhspan(0.5, ax.get_ylim()[1], facecolor='0.5', zorder=1, color=color_msa)
            elif templates_found:
                ax.axhspan(0.5, ax.get_ylim()[1], facecolor='0.5', zorder=1, color=color_template)
        elif plot_type == 'msa':
            ax.axhspan(-0.5, 0.5, facecolor='0.5', zorder=1, color=color_sequence)
            ax.axhspan(0.5, ax.get_ylim()[1], facecolor='0.5', zorder=1, color=color_msa)
        else:
            ax.axhspan(-0.5, 0.5, facecolor='0.5', zorder=1, color=color_sequence)
            ax.axhspan(0.5, ax.get_ylim()[1], facecolor='0.5', zorder=1, color=color_template)

    legend_elements.append('Yellow flags residues substituted by another type.\n'
                           'Red flags residues substituted following a sequence ('
                           'typically to match query sequence).\n'
                           'The gray scale on the template bar expresses similarity to the query sequence (black '
                           'identical, the lighter the more dissimilar).\n')
    legend_elements.reverse()

    ax.xaxis.grid(color='k', linestyle='dashed', alpha=0.4, which='major')
    ax.set_xticks(
        [a_air.sequence_assembled.get_starting_length(i) + 1 for i in range(a_air.sequence_assembled.total_copies)])
    ax.set_xticks(list(ax.get_xticks()) + [a_air.sequence_assembled.get_finishing_length(i) + 2 for i in
                                           range(a_air.sequence_assembled.total_copies)])

    length = a_air.sequence_assembled.length
    if length < 200:
        partition = 50
    elif length < 3000:
        partition = 100
    else:
        partition = 200

    length = a_air.sequence_assembled.length
    rounded_partitions_labels = []
    rounded_partitions_aux = []
    rounded_partitions_num = []
    threshold = length * 0.03
    x_ticks = list(ax.get_xticks())
    for i in range(partition, length, partition):
        close_to_tick = any(i > num - threshold and i < num + threshold for num in x_ticks)
        if not close_to_tick:
            rounded_partitions_labels.append(i)
            rounded_partitions_aux.append(i)
        rounded_partitions_num.append(i)

    # cut_chunk = [list(tup) for tup in a_air.chunk_list]
    # cut_chunk = utils.remove_list_layer(cut_chunk)
    # ax.set_xticks(list(ax.get_xticks()) + [cut + 1 for cut in cut_chunk])
    ax.set_xticklabels(ax.get_xticks(), rotation=45)

    ax1.set_xticks(rounded_partitions_aux, major=True)
    ax1.set_xticklabels(rounded_partitions_labels, rotation=45)
    ax1.set_xticks(list(ax1.get_xticks()) + rounded_partitions_num, major=True)
    ax1.set_xticks(generate_minor_ticks([1] + list(ax1.get_xticks()), step=1), minor=True)
    ax1.tick_params('x', length=7, which='major', labelsize='medium')
    ax1.tick_params('x', length=3, which='minor')
    ax1.set_xlim(-round(ax.get_xlim()[1] - total_length) / 6, total_length + round(ax.get_xlim()[1] - total_length) / 6)

    ax.tick_params('x', length=10, which='major', width=2, color='red', labelsize='large')
    ax.tick_params('x', length=5, which='minor')

    ax.set_xlabel('Residue number')
    ax.set_ylabel('Information')
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_color('k')

    fig.tight_layout()
    fig.subplots_adjust(top=.95)
    plt.title(title)
    plt.savefig(file, bbox_inches='tight', dpi=100)
    plt.cla()
    return file, ''.join(legend_elements)
