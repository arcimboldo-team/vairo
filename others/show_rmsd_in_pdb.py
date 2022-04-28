#!/cri4/albert/anaconda3/envs/af2/bin/python

import os
from Bio.PDB import PDBParser, Selection
import numpy as np


pdb1_path = '/localdata/albert/ATZR_PAPER/POLYALA_PREDICTIONS/atzr_gio/atzr_gio_in_polyala/summary/templates/BEST.pdb'
pdb2_path = '/localdata/albert/ATZR_PAPER/POLYALA_PREDICTIONS/atzr_gio/atzr_gio_in_polyala/summary/templates/newA_template.pdb'

structure1 = PDBParser(QUIET=1).get_structure('pdb1', pdb1_path)
structure1_res_id_list = [res.id[1] for res in Selection.unfold_entities(structure1, 'R')]

structure2 = PDBParser(QUIET=1).get_structure('pdb2', pdb2_path)
structure2_res_id_list = [res.id[1] for res in Selection.unfold_entities(structure2, 'R')]

common_res_between_structure1_and_structure2 = list(set(structure1_res_id_list).intersection(structure2_res_id_list))
common_res_rmsd = []
for res_id in common_res_between_structure1_and_structure2:
    
    CA1 = [res['CA'] for res in Selection.unfold_entities(structure1, 'R') if res.id[1] == res_id][0]
    CA2 = [res['CA'] for res in Selection.unfold_entities(structure2, 'R') if res.id[1] == res_id][0]
    common_res_rmsd.append(format(round(np.linalg.norm(CA1.coord - CA2.coord), 2),'.2f'))
    
    
pdb1_folder_path = '/'.join(pdb1_path.split('/')[:-1])

with open(f'{pdb1_folder_path}/temp.pdb','w') as out_f:
    
    with open(pdb1_path, 'r') as read_f:
        for line in read_f.readlines():
            if line[:4] == 'ATOM':
                if int(line[22:26]) in common_res_between_structure1_and_structure2:
                    match_index = common_res_between_structure1_and_structure2.index(int(line[22:26]))
                    out_f.write(line[:56] + common_res_rmsd[match_index].rjust(4) + line[60:])
                else:
                    out_f.write(line[:56] + format(0,'.2f').rjust(4) + line[60:])
            else:
                out_f.write(line)

os.system(f'mv {pdb1_folder_path}/temp.pdb {pdb1_path}')

     
