import os

def remove_hydrogens(pdb_path):

    output_file = open(f'{"/".join(pdb_path.split("/")[:-1])}/temp.pdb', 'w')

    counter = 0
    with open(f'{pdb_path}', 'r') as f:
        for line in f.readlines():
            if line.split()[-1] in ['N', 'C', 'O', 'S']:
                counter = counter + 1
                output_file.write(line[:6] + str(counter).rjust(5) + line[11:])


pdb_path = '/localdata/albert/ATZR_PAPER/POLYALA_PREDICTIONS/atzr_polyala_predictions_from_assambled_monomers/SUMMARIES/test/ranked_0.pdb'                
remove_hydrogens(pdb_path)
os.system(f'mv {"/".join(pdb_path.split("/")[:-1])}/temp.pdb {pdb_path}')

