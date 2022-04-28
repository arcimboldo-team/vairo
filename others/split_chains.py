#!/cri4/albert/anaconda3/envs/af2/bin/python

import os
import string


list_res_ranges = '1-300,351-650,701-1000,1051-1350'

pdb_path = '/localdata/albert/ATZR_PAPER/POLYALA_PREDICTIONS/atzr_polyala_predictions_from_assambled_monomers/SUMMARIES/atzr-1ixc/experimental_ranked_0.pdb'


dict = {}

initial_res_num, final_res_num = [], []


num_first_chain = list_res_ranges.split(',')[0].split('-')
res_range_first_chain = list(range(int(num_first_chain[0]), int(num_first_chain[1])+1))

chain_list = list(string.ascii_uppercase[:len(list_res_ranges.split(','))])
for num, chain in enumerate(chain_list):
    res_list = list_res_ranges.split(',')[num].split('-')
    final_res_num += res_range_first_chain
    initial_res_num += list(range(int(res_list[0]),int(res_list[1])+1))
    for res in list(range(int(res_list[0]),int(res_list[1])+1)):
        dict[res] = chain


with open(f'{"/".join(pdb_path.split("/")[:-1])}/temp.pdb', 'w') as f_out:

    with open(pdb_path) as f:
        for line in f.readlines():
            if 'ATOM'in line:
                chain, res_num = line[21:26][0], int(line[21:26][1:])
                if res_num in dict:
                    new_res_num = final_res_num[initial_res_num.index(res_num)]
                    f_out.write(line[:21] + dict[res_num] + str(new_res_num).rjust(4) + line[26:])   
            elif 'TER' in line:
                chain, res_num = line[21:26][0], int(line[21:26][1:])
                if res_num in dict:
                    f_out.write(line[:21] + dict[res_num] + str(new_res_num).rjust(4) + line[26:]) 
            else:
                f_out.write(line)

os.system(f'mv {"/".join(pdb_path.split("/")[:-1])}/temp.pdb {pdb_path}')
