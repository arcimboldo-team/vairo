import pickle as pickle
import numpy as np
import re
import string
import subprocess

def read_features_from_pkl(pickle_path):

    features_file = open(f'{pickle_path}', 'rb')
    features = pickle.load(features_file)
    features_file.close()

    return features


def show_keys_in_features(features):

    for key in features.keys():
        print(key, features[key].shape)


def read__sequence(fasta_path):

    with open(fasta_path, 'r') as f:
        fasta_lines = f.readlines()
        sequence = fasta_lines[1].split('\n')[0]

    return sequence


def find_mon_seq_in_assembly_seq(monomer_seq, assembly_seq):

    results_dict = {}

    seq_range_in_chains = [(m.start() + 1, m.end()) for m in re.finditer(monomer_seq, assembly_seq)]

    chain_list = string.ascii_uppercase
    for num, seq_range in enumerate(seq_range_in_chains):
        results_dict[chain_list[num]] = seq_range

    return results_dict

def split_chains_in_pdb(pdb_in_path, pdb_out_path, chain_res_range_dict):

    output_file = open(pdb_out_path, 'w')

    with open(pdb_in_path, 'r') as file:
        lines = file.readlines()

        for num, key in enumerate(chain_res_range_dict.keys()):
            for line in lines:
                if line.startswith('ATOM'):
                    # print(line)
                    if int(line[22:26]) in range(seq_range_in_chains[key][0], seq_range_in_chains[key][1] + 1):
                        new_res_num = str(int(line[22:26]) - seq_range_in_chains[key][0] + 1)
                        output_file.write(line[:21] + key + new_res_num.rjust(4) + line[26:])
                else:
                    pass


def find_interface_from_pisa(input_pdb_path):

    subprocess.Popen(['pisa', 'temp', '-analyse', input_pdb_path],
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE).communicate()

    pisa_output = subprocess.Popen(['pisa', 'temp', '-list', 'interfaces'], stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE).communicate()[0].decode('utf-8')

    match1 = [m.start() for m in re.finditer(" LIST OF INTERFACES", pisa_output)][0]
    match2 = [m.start() for m in re.finditer(" ##:  serial number", pisa_output)][0]

    interfaces_dict = {}
    for line in pisa_output[match1:match2].split('\n')[4:-2]:
        area = line.split('|')[3][:8].replace(' ', '')
        deltaG = line.split('|')[3][8:15].replace(' ', '')
        chain1 = line.split('|')[1].replace(' ', '')
        chain2 = line.split('|')[2].split()[0].replace(' ', '')
        interfaces_dict[f'{chain1}{chain2}'] =  (area, deltaG)

    return interfaces_dict


def split_dimers_in_pdb(input_pdb_path, output_pdb_path, chain1, chain2):

    output = open(output_pdb_path, 'w')

    input_lines = open(input_pdb_path, 'r').readlines()
    for line in input_lines:
        if line[21:22] in [chain1, chain2]:
            output.write(line)
    output.close()




FEATURES_PKL_PATH = '/cri4/albert/Desktop/ATZR_PAPER/features.pkl'
INPUT_PDB = '/cri4/albert/Desktop/ATZR_PAPER/ranked_0.pdb'
MONOMER_FASTA = '/cri4/albert/Desktop/atzr.fasta'

features = read_features_from_pkl(pickle_path=FEATURES_PKL_PATH)
show_keys_in_features(features=features)

monomer_seq = read__sequence(fasta_path=MONOMER_FASTA)
assembly_seq = features['sequence'][0].decode('utf-8')

seq_range_in_chains = find_mon_seq_in_assembly_seq(monomer_seq=monomer_seq, assembly_seq=assembly_seq)

split_chains_in_pdb(pdb_in_path=INPUT_PDB,
                    pdb_out_path=f'{INPUT_PDB[:-4]}_split.pdb',
                    chain_res_range_dict=seq_range_in_chains)

interfaces_dict = find_interface_from_pisa(input_pdb_path=f'{INPUT_PDB[:-4]}_split.pdb')


output_path = '/'.join(f'{INPUT_PDB[:-4]}_split.pdb'.split('/')[:-1])
for interfaces in interfaces_dict:
    chain1, chain2 = list(interfaces)

    split_dimers_in_pdb(input_pdb_path=f'{INPUT_PDB[:-4]}_split.pdb',
                        output_pdb_path=f'{output_path}/{chain1}{chain2}.pdb',
                        chain1=chain1,
                        chain2=chain2)




# for interface in interfaces_dict:
#     chain1, chain2 = list(interface)
#     for line in input_lines:
#         print(line)
