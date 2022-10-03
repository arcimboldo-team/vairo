#! /usr/bin/env python3

import pickle
import os
import sys
from libs import features
import logging

def write_features(features_path: str, output_dir: str = None):
    with open(os.path.abspath(features_path), 'rb') as f:
        data = pickle.load(f)
    if output_dir is None:
        output_dir = os.getcwd()
    features.write_templates_in_features(data, output_dir)

def print_features(features_path: str):
    logging.info = print
    features.print_features(features_path)

if __name__ == "__main__":
    args = sys.argv
    globals()[args[1]](*args[2:])