import os

def convert_template_to_polyala(pdb_in_path, pdb_out_path, list_of_res_ranges):

    ala_res_list = []

    if len(list_of_res_ranges) > 1:
        for item in list_of_res_ranges:
            ala_res_list.extend(list(range(int(item[0]),int(item[1]))))
    else:
        ala_res_list = list(range(int(list_of_res_ranges[0][0]), int(list_of_res_ranges[0][1])))

    ala_atoms_list = ['N', 'CA', 'C', 'CB', 'O']

    with open(f'{pdb_out_path}', 'w') as f_out:
        with open(f'{pdb_in_path}') as f_in:
            lines = f_in.readlines()
            num = 0
            for line in lines:
                if line[:4] == 'ATOM':
                    if int(line[22:26].replace(' ', '')) in ala_res_list:
                        if line[13:16].replace(' ', '') in ala_atoms_list:
                            num = num + 1
                            f_out.write(line[:7] + str(num).rjust(4) + line[11:17] + 'ALA' + line[20:])
                    else:
                        num = num + 1
                        f_out.write(line[:7] + str(num).rjust(4) + line[11:])

def remove_hydrogens(pdb_path):

    path = "/".join(pdb_path.split('/')[:-1])

    output_file = open(f'{path}/output.pdb', 'w')

    counter = 0
    with open(f'{pdb_path}', 'r') as f:
        for line in f.readlines():
            if line.split()[-1] in ['N', 'C', 'O', 'S']:
                counter = counter + 1
                output_file.write(line[:6] + str(counter).rjust(5) + line[11:])
    os.system(f'mv {path}/output.pdb {path}/{pdb_path.split("/")[-1]}')

remove_hydrogens(pdb_path='/cri4/albert/repos/arcimboldo_air/polyala_vs_sequence/ranked_0.pdb')
#convert_template_to_polyala(pdb_in_path='/home/albert/repos/arcimboldo_air/output/3fzvA_template.pdb',
#                            pdb_out_path='/home/albert/repos/arcimboldo_air/output/templates/3fzvA_Polyala.pdb',
#                            list_of_res_ranges=[[-1000, 1000]])
