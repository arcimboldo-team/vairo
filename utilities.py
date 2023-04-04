#! /usr/bin/env python3

import pickle
import os
import sys
from libs import features, bioutils, utils, structures
import logging


def write_features(features_path: str, output_dir: str = None):
    with open(os.path.abspath(features_path), 'rb') as f:
        data = pickle.load(f)
    if output_dir is None:
        output_dir = os.getcwd()
    features.write_templates_in_features(data, output_dir)


def print_features(features_path: str):
    logging.info = print
    features.print_features_from_file(features_path)


def generate_features(query_path: str, fasta_path: str):
    path = os.path.join(os.getcwd(), 'features.pkl')
    query = bioutils.extract_sequence(query_path)
    sequences = bioutils.extract_sequences(fasta_path)
    feature = features.Features(query)
    [feature.append_row_in_msa(sequence=seq, sequence_id=seq_id) for seq_id, seq in sequences.items()]
    write_features(path)


def ccanalysis(template_path: str):
    output_path = os.path.join(template_path, 'ccanalysis')
    os.listdir(template_path)
    templates_dict = {utils.get_file_name(path):os.path.join(template_path, path) for path in os.listdir(template_path) if path.endswith('.pdb')}
    cc_analysis = structures.CCAnalysis(os.path.join(utils.get_main_path(), 'binaries'))
    bioutils.cc_analysis(paths_in=templates_dict, cc_analysis_paths=cc_analysis, cc_path=output_path)


if __name__ == "__main__":
    print('Usage: utilities.py function input')
    print('Functions: write_features, print_features')
    args = sys.argv
    globals()[args[1]](*args[2:])