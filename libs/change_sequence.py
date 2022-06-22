from Bio.PDB import PDBParser, Selection


def change_seq(pdb_target, pdb_to_change, output_pdb):


    target_structure = PDBParser(QUIET=True).get_structure('test', pdb_target)
    target_res_list = [res for res in Selection.unfold_entities(target_structure, "R")]

    to_change_structure = PDBParser(QUIET=True).get_structure('test', pdb_to_change)
    to_change_resid_list = [res.id[1] for res in Selection.unfold_entities(to_change_structure, "R")]

    modifications_dict = {}

    for res in target_res_list:
        if res.id[1] in to_change_resid_list:
            match_res_num = res.id[1]
            modifications_dict[res.id[1]] = (res.resname, [res for res in Selection.unfold_entities(to_change_structure, "R") if res.id[1] == match_res_num][0].resname)
            print(res.id[1], res.resname, [res for res in Selection.unfold_entities(to_change_structure, "R") if res.id[1] == match_res_num][0].resname)
        else:
            modifications_dict[res.id[1]] = (res.resname, 'ALA')
            print(res.id[1], res.resname, 'ALA', '****')


    with open(output_pdb, 'w') as f:

        pdb_target_lines = open(pdb_target, 'r').readlines()
        num = 0
        for line in pdb_target_lines:
            if line[13:17].replace(' ', '') in ['C', 'O', 'N', 'CA', 'CB']:
                num = num + 1
                new_res = modifications_dict[int(line[22:26].replace(' ', ''))][1]
                f.write(line[:7] + f'{num}'.rjust(4) + line[11:17] + new_res + line[20:])


change_seq(pdb_target='/cri4/albert/repos/arcimboldo_air/polyala_vs_sequence/ranked_0.pdb',
           pdb_to_change='/cri4/albert/repos/arcimboldo_air/polyala_vs_sequence/ranked_0.pdb',
           output_pdb='/cri4/albert/repos/arcimboldo_air/polyala_vs_sequence/templates/ranked_rankedWithoutSidechains.pdb')
