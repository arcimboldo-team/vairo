import pickle5 as pickle
import sys
sys.path.append('/cri4/albert/repos/arcimboldo_air_v4/alphafold')
from alphafold.common import residue_constants
import numpy as np

# features_pkl_path = '/home/albert/Desktop/JOBS/4/features.pkl'
features_pkl_path = '/localdata/albert/SEQUENCE_DEPENDENCE_PAPER/JOBS/ATZR/2_fixed_Nter_polyalaCter/features.pkl'

with open(f"{features_pkl_path}", "rb") as input_file:
    features = pickle.load(input_file)

for key in features.keys():
    print(key, features[key].shape)

print('\n')
print('MSA:')
for num, name in enumerate(features['msa_uniprot_accession_identifiers']):
    print(name.decode('utf-8'),'\n')
    print(''.join([residue_constants.ID_TO_HHBLITS_AA[res] for res in features['msa'][num].tolist()]))
    print('\n')

print('TEMPLATES:')
print('\n')
for num, seq in enumerate(features['template_sequence']):
    print(f'{features["template_domain_names"][num].decode("utf-8")}:\n')
    for i in range(4):
        print('\t'+''.join(np.array_split(list(seq.decode('utf-8')),4)[i].tolist()))
    print('\n')

