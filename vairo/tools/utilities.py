#! /usr/bin/env python3
import sys
sys.path.append('../')

import pickle
import os
import logging
from libs import features, bioutils, output, utils, structures


def write_features(features_path: str, output_dir: str = None):
    with open(os.path.abspath(features_path), 'rb') as f:
        data = pickle.load(f)
    if output_dir is None:
        output_dir = os.getcwd()
    features.write_templates_in_features(data, output_dir)


def print_features(features_path: str):
    logging.error = print
    features.print_features_from_file(features_path)


def generate_features(query_path: str, fasta_path: str):
    path = os.path.join(os.getcwd(), 'features.pkl')
    query = bioutils.extract_sequence(query_path)
    sequences = bioutils.extract_sequences(fasta_path)
    feature = features.Features(query)
    [feature.append_row_in_msa(sequence=seq, sequence_id=seq_id) for seq_id, seq in sequences.items()]
    write_features(path)


def hinges(template_path: str):
    output_path = os.path.join(template_path, 'hinges')
    os.listdir(template_path)
    templates_dict = {utils.get_file_name(path): os.path.join(template_path, path) for path in os.listdir(template_path)
                      if path.endswith('.pdb')}

    binaries_path = structures.CCAnalysis(os.path.join(utils.get_main_path(), 'binaries'))
    templates_cluster = bioutils.hinges(paths_in=templates_dict,
                                        binaries_path=binaries_path,
                                        output_path=output_path)

    for i, values in enumerate(templates_cluster):
        print(f'Group {i}: {",".join(values)}')


def ccanalysis(template_path: str):
    output_path = os.path.join(template_path, 'ccanalysis')
    os.listdir(template_path)
    templates_dict = {utils.get_file_name(path): os.path.join(template_path, path) for path in os.listdir(template_path)
                      if path.endswith('.pdb')}
    cc_analysis = structures.CCAnalysis(os.path.join(utils.get_main_path(), 'binaries'))
    templates_cluster_list, analysis_dict = bioutils.cc_analysis(paths_in=templates_dict, cc_analysis_paths=cc_analysis,
                                                                 cc_path=output_path, n_clusters=2)
    if analysis_dict:
        output.plot_cc_analysis(plot_path=os.path.join(output_path, 'plot.png'), analysis_dict=analysis_dict,
                                    clusters=templates_cluster_list)

def superposition_chains(pdb1_path: str, pdb2_path: str):
    ret_dict = bioutils.superposition_by_chains(pdb1_in_path=pdb1_path, pdb2_in_path=pdb2_path)
    for key3, i3 in ret_dict.items():
        print(key3)
        for key2, i2 in i3.items():
            print(key2)
            for key1, i1 in i2.items():
                print(key1, i1)

def superpose_chains(pdb1_path: str, pdb2_path: str, tmp_dir: str):
    utils.create_dir(tmp_dir)
    pdb1_chains_dict = bioutils.split_pdb_in_chains(pdb_path=pdb1_path, output_dir=tmp_dir)
    pdb2_chains_dict = bioutils.split_pdb_in_chains(pdb_path=pdb2_path, output_dir=tmp_dir)

def renumber():
    def check_consecutive(numbers):
        # Check if the difference between each pair of consecutive numbers is equal to 1
        for i in range(len(numbers) - 1):
            if numbers[i + 1] - numbers[i] != 1:
                return False
        return True


    # Specify the folder path containing the PDB files
    folder_path = "/Users/pep/work/transfers/clusters_lib"
    # Get a list of all PDB files in the folder
    pdb_files = [os.path.join(folder_path, file) for file in os.listdir(folder_path) if file.endswith(".pdb")]

    list_pdbs = []

    # Loop through each PDB file
    for pdb_file in pdb_files:

        structure = bioutils.get_structure(pdb_file)

        # Initialize counters for CYS residues and consecutive positions
        cys_count = 0
        save_residues = []
        save_pdb = False

        # Iterate over all residues in the structure
        for model in structure:
            for chain in model:
                residues = list(chain.get_residues())
                residues = sorted(residues, key=lambda x: bioutils.get_resseq(x))
                for j, residue in enumerate(residues):
                    # Check if the residue is CYS
                    if residue.get_resname() == 'CYS':
                        cys_count += 1
                        if cys_count == 1:
                            try:
                                list_cys = [residues[j+i] for i in range(-5, 2)]
                                list_cys = [bioutils.get_resseq(res)-1 for res in list_cys]
                                if check_consecutive(list_cys):
                                    save_residues.extend(list_cys)
                                else:
                                    raise Exception
                            except:
                                cys_count = 3
                                pass
                        if cys_count == 2:
                            try:
                                list_cys = [residues[j+i] for i in range(-5, 3)]
                                list_cys = [bioutils.get_resseq(res)-1 for res in list_cys]
                                if check_consecutive(list_cys):
                                    save_residues.extend(list_cys)
                                    if utils.get_file_name(pdb_file)[:4] not in list_pdbs:
                                        list_pdbs.append(utils.get_file_name(pdb_file)[:4])
                                        save_pdb = True
                                else:
                                    raise Exception
                            except:
                                pass
                

        if save_pdb:
            print(save_residues)
            if len(save_residues) != 15:
                raise Exception
            bioutils.copy_positions_of_pdb(pdb_file, os.path.join("/Users/pep/work/transfers/library", utils.get_file_name(pdb_file))+'.pdb', save_residues)
            print(f"Renumbering complete for {pdb_file}. Renumbered file saved as {utils.get_file_name(pdb_file)}.")



if __name__ == "__main__":
    print('Usage: utilities.py function input')
    print('Functions: write_features, print_features')
    logging.error = print
    args = sys.argv
    globals()[args[1]](*args[2:])
