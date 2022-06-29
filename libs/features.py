import sys
import alphafold
from Bio.PDB import PDBParser, MMCIFIO, PDBList, PDBIO
from Bio.PDB import MMCIFParser, Selection
import os
import re
from pathlib import Path
from ALPHAFOLD.alphafold.data import parsers, templates, mmcif_parsing, pipeline, msa_identifiers
from ALPHAFOLD.alphafold.common import residue_constants
import numpy as np
import subprocess
import string
import pickle
import logging


three_to_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P',
                'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R',
                'TRP': 'W', 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
AA_TO_ID_TO_HHBLITS = {v: k for k, v in residue_constants.ID_TO_HHBLITS_AA.items()}
ID_TO_HHBLITS_AA_3LETTER_CODE = {0: 'ALA', 1: 'CYS', 2: 'ASP', 3: 'GLU', 4: 'PHE', 5: 'GLY', 6: 'HIS',
                                 7: 'ILE', 8: 'LYS', 9: 'LEU', 10: 'MET', 11: 'ASN', 12: 'PRO', 13: 'GLN',
                                 14: 'ARG', 15: 'SER', 16: 'THR', 17: 'VAL', 18: 'TRP', 19: 'TYR', 20: 'X',
                                 21: '-'}
atom_types = ['N', 'CA', 'C', 'CB', 'O', 'CG', 'CG1', 'CG2', 'OG', 'OG1', 'SG', 'CD', 'CD1', 'CD2',
              'ND1', 'ND2', 'OD1', 'OD2', 'SD', 'CE', 'CE1', 'CE2', 'CE3', 'NE', 'NE1', 'NE2', 'OE1',
              'OE2', 'CH2', 'NH1', 'NH2', 'OH', 'CZ', 'CZ2', 'CZ3', 'NZ', 'OXT']
atom_order = {atom_type: i for i, atom_type in enumerate(atom_types)}
order_atom = {v: k for k, v in atom_order.items()}


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


def download_pdb(pdb_id, out_dir):

    logging.info(f'Downloading PDB {pdb_id}')

    pdbl = PDBList()
    pdbl.retrieve_pdb_file(pdb_id, pdir=f'{out_dir}', file_format='pdb', obsolete=False)
    os.system('rm -r obsolete')


def pdb2mmcif(pdb_in_path, cif_out_path):

    with open('/tmp/maxit.sh','w') as f:
        f.write(f'export RCSBROOT="/opt/maxit"\n')
        f.write(f'/opt/maxit/bin/maxit -input {pdb_in_path} -output {cif_out_path} -o 1\n')
    os.system('chmod +x /tmp/maxit.sh')
    os.system('bash /tmp/maxit.sh')


def cif2pdb(cif_in_path, pdb_out_path):

    p = MMCIFParser(QUIET=True)
    struc = p.get_structure('', f'{cif_in_path}')
    io = PDBIO()
    io.set_structure(struc)
    io.save(f'{pdb_out_path}')


def extract_query_sequence(query_fasta_path):

    logging.info(f'Extracting query sequence from {query_fasta_path}')

    with open(query_fasta_path, 'r') as f:
        fasta_lines = f.readlines()
        fasta_name, query_sequence = fasta_lines[0][1:-1], fasta_lines[1].split('\n')[0]

    return fasta_name, query_sequence


def empty_msa_features(query_sequence):

    msa = {'a3m': f'>query\n{query_sequence}'}
    custom_msa = parsers.parse_a3m(msa['a3m'])

    msas = [custom_msa]  # ACT: it is needed in order to introduce MSA inside a list in the code

    int_msa = []
    deletion_matrix = []
    uniprot_accession_ids = []
    species_ids = []
    seen_sequences = set()
    for msa_index, msa in enumerate(msas):
        if not msa:
            raise ValueError(f'MSA {msa_index} must contain at least one sequence.')
        for sequence_index, sequence in enumerate(msa.sequences):
            if sequence in seen_sequences:
                continue
            seen_sequences.add(sequence)
            int_msa.append(
                [residue_constants.HHBLITS_AA_TO_ID[res] for res in sequence])
            deletion_matrix.append(msa.deletion_matrix[sequence_index])
            identifiers = msa_identifiers.get_identifiers(
                msa.descriptions[sequence_index])
            uniprot_accession_ids.append(
                identifiers.uniprot_accession_id.encode('utf-8'))
            species_ids.append(identifiers.species_id.encode('utf-8'))

    num_res = len(msas[0].sequences[0])
    num_alignments = len(int_msa)
    features = {}
    features['deletion_matrix_int'] = np.array(deletion_matrix, dtype=np.int32)
    features['msa'] = np.array(int_msa, dtype=np.int32)
    features['num_alignments'] = np.array(
        [num_alignments] * num_res, dtype=np.int32)
    features['msa_uniprot_accession_identifiers'] = np.array(
        uniprot_accession_ids, dtype=np.object_)
    features['msa_species_identifiers'] = np.array(species_ids, dtype=np.object_)
    return features


def empty_template_features(query_sequence):

    ln = (len(query_sequence) if isinstance(query_sequence, str) else sum(len(s) for s in query_sequence))
    output_templates_sequence = "A" * ln

    templates_all_atom_positions = np.zeros((ln, residue_constants.atom_type_num, 3))
    templates_all_atom_masks = np.zeros((ln, residue_constants.atom_type_num))
    templates_aatype = residue_constants.sequence_to_onehot(output_templates_sequence, residue_constants.HHBLITS_AA_TO_ID)
    template_sum_probs = f'None'
    template_features = {
        "template_all_atom_positions": np.tile(templates_all_atom_positions[None], [0, 1, 1, 1]),
        "template_all_atom_masks": np.tile(templates_all_atom_masks[None], [0, 1, 1]),
        "template_sequence": [f"None".encode()] * 0,
        "template_aatype": np.tile(np.array(templates_aatype)[None], [0, 1, 1]),
        "template_domain_names": [f"None".encode()] * 0,
        "template_sum_probs": np.tile(template_sum_probs, [0, 1])
    }
    return template_features


def generate_hhsearch_db(template_cif_path, output_dir):

    with open(f"{output_dir}/pdb70_a3m.ffdata", "w") as a3m, \
         open(f"{output_dir}/pdb70_cs219.ffindex", "w") as cs219_index, \
         open(f"{output_dir}/pdb70_a3m.ffindex", "w") as a3m_index, \
         open(f"{output_dir}/pdb70_cs219.ffdata", "w") as cs219:
        id = 1000000
        index_offset = 0

        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("test", template_cif_path)
        models = list(structure.get_models())
        if len(models) != 1:
            raise ValueError(f"Only single model PDBs are supported. Found {len(models)} models.")
        model = models[0]
        for chain in model:
            amino_acid_res = []
            for res in chain:
                # if res.id[2] != " ":
                #     raise ValueError(f"PDB contains an insertion code at chain {chain.id} and residue index"
                #                      f" {res.id[1]}. These are not supported.")
                amino_acid_res.append(residue_constants.restype_3to1.get(res.resname, "X"))

            protein_str = "".join(amino_acid_res)
            a3m_str = f">{template_cif_path.split('/')[-1][:-4]}_{chain.id}\n{protein_str}\n\0"
            a3m_str_len = len(a3m_str)
            a3m_index.write(f"{id}\t{index_offset}\t{a3m_str_len}\n")
            cs219_index.write(f"{id}\t{index_offset}\t{len(protein_str)}\n")
            index_offset += a3m_str_len
            a3m.write(a3m_str)
            cs219.write("\n\0")
            id += 1


def run_hhsearch(query_fasta_path, pdb70_db, output_path):

    print(f'Running hhsearch using {pdb70_db} as database.')

    out = subprocess.Popen(['hhsearch', '-i', query_fasta_path, '-o', output_path, '-oa3m', '/home/albert/Desktop/test.a3m', '-maxseq',
                            '1000000', '-d', pdb70_db, '-glob'],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = out.communicate()
    hhr = stdout.decode('utf-8')

    return hhr


def read_remark_350(pdb_path, use_pisa=False):

    pdb_id = pdb_path

    if not use_pisa:
        logging.info(f'Reading REMARK 350 from {pdb_path}.')
        pdb_text = open(pdb_id, 'r').read()
    else:
        print(f'Generating REMARK 350 for {pdb_path} with PISA.')
        subprocess.Popen(['pisa', 'temp', '-analyse', pdb_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        pisa_output = subprocess.Popen(['pisa', 'temp', '-350'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0]
        pdb_text = pisa_output.decode('utf-8')

    match_biomolecules = [m.start() for m in re.finditer(r'REMARK 350 BIOMOLECULE:', pdb_text)] # to know how many biomolecules there are.
    if len(match_biomolecules) == 1:
        match_last_350 = [m.start() for m in re.finditer(r'REMARK 350', pdb_text)][-1]
        match_end_in_last_350 = [m.end() for m in re.finditer(r'\n', pdb_text[match_last_350:])][-1]
        remark_350_text = pdb_text[match_biomolecules[0]:(match_last_350+match_end_in_last_350)]
    else:
        print('(It seem there is more than one biological assembly from REMARK 350. Only'
              ' "BIOMOLECULE 1" will be considered for the assembly generation)')
        remark_350_text = pdb_text[match_biomolecules[0]:match_biomolecules[1]-1]

    match_biomt1 = [m.start() for m in re.finditer(r'REMARK 350   BIOMT1', remark_350_text)]
    match_biomt3 = [m.end() for m in re.finditer(r'REMARK 350   BIOMT3', remark_350_text)]

    end_remark_350_block = [m.start() for m in re.finditer('\n', remark_350_text[match_biomt3[-1]:])]

    transformation_blocks_indices = match_biomt1 + [match_biomt3[-1] + end_remark_350_block[0] + 1]

    transformations_list = []
    for index in range(len(transformation_blocks_indices) - 1):
        block = remark_350_text[transformation_blocks_indices[index]:transformation_blocks_indices[index+1]]
        matrix = [item.split()[4:8] for item in block.split('\n')[:-1]]
        r11, r12, r13, t1 = matrix[0]
        r21, r22, r23, t2 = matrix[1]
        r31, r32, r33, t3 = matrix[2]
        transformation = [[r11, r12, r13], [r21, r22, r23], [r31, r32, r33], [t1, t2, t3]]
        transformations_list.append(transformation)

    chain_list = [line.split(':')[-1].replace(" ", "").split(',') for line in remark_350_text.split('\n')
                  if 'REMARK 350 APPLY THE FOLLOWING TO CHAINS:' in line][0]

    return chain_list, transformations_list


def rot_and_trans(pdb_path, out_pdb_path, rot_tra_matrix):

    r11, r12, r13 = rot_tra_matrix[0]
    r21, r22, r23 = rot_tra_matrix[1]
    r31, r32, r33 = rot_tra_matrix[2]
    t1, t2, t3 = rot_tra_matrix[3]

    with open('/tmp/pdbset.sh', 'w') as f:
        f.write(f'pdbset xyzin {pdb_path} xyzout {out_pdb_path} << eof\n')
        f.write(
            f'rotate {float(r11)} {float(r12)} {float(r13)} {float(r21)} {float(r22)} {float(r23)} {float(r31)} {float(r32)} {float(r33)}\n')
        f.write(f'shift {float(t1)} {float(t2)} {float(t3)}\n')
        f.write('end\n')
        f.write('eof')
    subprocess.Popen(['bash', '/tmp/pdbset.sh'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
    os.system(f'rm /tmp/pdbset.sh')


def change_chain_and_apply_offset_in_single_chain(pdb_in_path, pdb_out_path, offset=None, chain=None):

    with open('/tmp/pdbset.sh', 'w') as f:
        f.write(f'pdbset xyzin {pdb_in_path} xyzout {pdb_out_path} << eof\n')
        if offset:
            f.write(f'renumber increment {offset}\n')
        if chain:
            f.write(f'chain {chain}\n')
        f.write('end\n')
        f.write('eof')
    subprocess.Popen(['bash', '/tmp/pdbset.sh'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
    os.system(f'rm /tmp/pdbset.sh')


def merge_pdbs(list_of_paths_of_pdbs_to_merge, merged_pdb_path):

    with open(merged_pdb_path, 'w') as f:
        counter = 0
        for pdb_path in list_of_paths_of_pdbs_to_merge:
            for line in open(pdb_path,'r').readlines():
                if line[:4] == 'ATOM':
                    counter += 1
                    f.write(line[:4] + str(counter).rjust(7) + line[11:])


def merge_features(sequence_features, msa_features, template_features):

    logging.info(f'Merging sequence, msa and template features!')

    return {**sequence_features, **msa_features, **template_features}


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


def write_pkl_from_features(features, out_path):

    logging.info(f'Writting all input features in {out_path}')

    with open(out_path, 'wb') as handle:
        pickle.dump(features, handle, protocol=pickle.HIGHEST_PROTOCOL)





class Features:

    def __init__(self, query_sequence):

        self.query_sequence = query_sequence
        self.sequence_features = pipeline.make_sequence_features(sequence=self.query_sequence,
                                                                 description='Query',
                                                                 num_res=len(self.query_sequence))
        self.msa_features = empty_msa_features(query_sequence=self.query_sequence)
        self.template_features = empty_template_features(query_sequence=self.query_sequence)

    def extract_template_features_from_pdb(self, hhr_path, pdb_id, chain_id, mmcif_db): # TODO: template names must be in lowercase

        hhr_text = open(hhr_path, 'r').read()

        matches = re.finditer(r'No\s+\d+', hhr_text)
        matches_positions = [match.start() for match in matches] + [len(hhr_text)]

        detailed_lines_list = []
        for i in range(len(matches_positions) - 1):
            detailed_lines_list.append(hhr_text[matches_positions[i]:matches_positions[i + 1]].split('\n')[:-3])

        if 'aligned_' in pdb_id:

            hits_list = [detailed_lines for detailed_lines in detailed_lines_list if
                     detailed_lines[1].split('_')[1] == f'{pdb_id.split("_")[-1]}']
        else:
            hits_list = [detailed_lines for detailed_lines in detailed_lines_list if
                         detailed_lines[1].split('_')[0][1:] == f'{pdb_id}']

        detailed_lines = hits_list[0]

        file_id = f'{pdb_id.lower()}'
        hit = parsers._parse_hhr_hit(detailed_lines)

        template_sequence = hit.hit_sequence.replace('-', '')
        mapping = templates._build_query_to_hit_index_mapping(
            hit.query, hit.hit_sequence, hit.indices_hit, hit.indices_query,
            self.query_sequence)

        mmcif_string = open(f'{mmcif_db}/{file_id}.cif').read()
        parsing_result = mmcif_parsing.parse(file_id=file_id, mmcif_string=mmcif_string)

        template_features, realign_warning = templates._extract_template_features(
            mmcif_object=parsing_result.mmcif_object,
            pdb_id=file_id,
            mapping=mapping,
            template_sequence=template_sequence,
            query_sequence=self.query_sequence,
            template_chain_id=chain_id,
            kalign_binary_path='kalign')
        template_features['template_sum_probs'] = np.array([[hit.sum_probs]])
        template_features['template_aatype'] = np.array([template_features['template_aatype']])
        template_features['template_all_atom_masks'] = np.array([template_features['template_all_atom_masks']])
        template_features['template_all_atom_positions'] = np.array([template_features['template_all_atom_positions']])
        template_features['template_domain_names'] = np.array([template_features['template_domain_names']])
        template_features['template_sequence'] = np.array([template_features['template_sequence']])

        return template_features

    def extract_template_features_from_aligned_pdb_and_sequence(self, pdb_path, chain_ID):

        # WARNING: input PDB must be aligned to the MSA part in features #

        name = pdb_path.split('/')[-1][:-4]

        seq_length = len(self.query_sequence)

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('test', pdb_path)

        template_sequence = '-' * (seq_length)
        template_res_list = [res for res in Selection.unfold_entities(structure, "R")
                             if res.get_parent().id == chain_ID and res.id[0] != 'W']
        for res in template_res_list:
            template_sequence = template_sequence[:res.id[1]] + three_to_one[res.resname] + template_sequence[
                                                                                            res.id[1]:]
        template_sequence = np.array([template_sequence[:seq_length + 1]])[0]

        atom_masks = []
        for i, res in enumerate(template_sequence):
            if res == '-':
                a37_in_res = [0] * 37
                atom_masks.append(a37_in_res)
            else:
                a37_in_res = [0] * 37
                list_of_atoms_in_res = [atom.id for atom in [resi for resi in template_res_list if resi.id[1] == i][0]]
                for atom in list_of_atoms_in_res:
                    index = atom_types.index(atom)
                    a37_in_res[index] = 1.
                atom_masks.append(a37_in_res)
        template_all_atom_masks = np.array([atom_masks[1:]])

        template_container = []
        for i, res in enumerate(template_all_atom_masks[0]):
            res_container = []
            for j, atom in enumerate(res):
                if atom == 1.:
                    resi = [res for res in Selection.unfold_entities(structure, "R")
                            if res.get_parent().id == chain_ID and res.id[0] != 'W' and res.id[1] == (i + 1)][0]
                    res_container.append(resi[atom_types[j]].coord)
                else:
                    res_container.append(np.array([0.] * 3))
            template_container.append(res_container)
        template_all_atom_positions = np.array([template_container])

        template_domain_names = np.array([(f'{name}_' + chain_ID).encode('ascii')])

        template_aatype_container = []
        for res in template_sequence[1:]:
            aa_container = [0] * 22
            aa_container[AA_TO_ID_TO_HHBLITS[res]] = 1
            template_aatype_container.append(aa_container)
        template_aatype = np.array([template_aatype_container])

        template_sum_probs = np.array([100.])

        template_sequence_to_add = np.array([template_sequence[1:].encode('ascii')])
        template_all_atom_masks_to_add = template_all_atom_masks
        template_all_atom_positions_to_add = template_all_atom_positions
        template_domain_names_to_add = template_domain_names
        template_aatype_to_add = template_aatype
        template_sum_probs_to_add = template_sum_probs

        template_features = {}

        template_features['template_sequence'] = template_sequence_to_add
        template_features['template_all_atom_masks'] = template_all_atom_masks_to_add
        template_features['template_all_atom_positions'] = template_all_atom_positions_to_add
        template_features['template_domain_names'] = template_domain_names_to_add
        template_features['template_aatype'] = template_aatype_to_add
        template_features['template_sum_probs'] = np.array([template_sum_probs_to_add])

        return template_features

    def append_new_template_features(self, new_template_features, custom_sum_prob = None):

        self.template_features['template_all_atom_positions'] = np.vstack([self.template_features['template_all_atom_positions'], new_template_features['template_all_atom_positions']])
        self.template_features['template_all_atom_masks'] = np.vstack([self.template_features['template_all_atom_masks'], new_template_features['template_all_atom_masks']])
        self.template_features['template_aatype'] = np.vstack([self.template_features['template_aatype'], new_template_features['template_aatype']])
        self.template_features['template_sequence'] = np.hstack([self.template_features['template_sequence'], new_template_features['template_sequence']])
        self.template_features['template_domain_names'] = np.hstack([self.template_features['template_domain_names'], new_template_features['template_domain_names']])
        if not custom_sum_prob:
            self.template_features['template_sum_probs'] = np.vstack([self.template_features['template_sum_probs'], new_template_features['template_sum_probs']])
        else:
            self.template_features['template_sum_probs'] = np.vstack([self.template_features['template_sum_probs'], custom_sum_prob])

        return self.template_features

    def append_row_in_msa(self, sequence, msa_uniprot_accession_identifiers):

        sequence_array = np.array([AA_TO_ID_TO_HHBLITS[res] for res in sequence])
        self.msa_features['msa'] = np.vstack([self.msa_features['msa'], sequence_array])
        self.msa_features['msa_uniprot_accession_identifiers'] = np.hstack([self.msa_features['msa_uniprot_accession_identifiers'], msa_uniprot_accession_identifiers.encode()])
        self.msa_features['deletion_matrix_int'] = np.vstack([self.msa_features['deletion_matrix_int'], np.zeros(self.msa_features['msa'].shape[1])])
        self.msa_features['msa_species_identifiers'] = np.hstack([self.msa_features['msa_species_identifiers'], ''])
        self.msa_features['num_alignments'] = np.full(self.msa_features['num_alignments'].shape, len(self.msa_features['msa']))

    def write_all_templates_in_features(self, output_path):

        for i, pdb_name in enumerate(self.template_features['template_domain_names']):
            pdb, chain = pdb_name.decode('utf-8').split('_')[:-1][0], pdb_name.decode('utf-8').split('_')[-1]

            # pdb, chain = pdb_name.decode('utf-8').split('_')
            output_pdb = open(output_path + '/' + pdb + chain + '_template.pdb', 'w')
            template_domain_index = np.where(self.template_features['template_domain_names'] == pdb_name)[0][0]
            true_seq = ''
            atom_num_int = 0  # AQUI
            for index, atoms_mask in enumerate(self.template_features['template_all_atom_masks'][template_domain_index][:]):
                template_residue_masks = self.template_features['template_aatype'][template_domain_index][index]
                template_residue_masks_index = np.where(template_residue_masks == 1)[0][0]
                res_type = ID_TO_HHBLITS_AA_3LETTER_CODE[template_residue_masks_index]
                list_of_atoms_in_residue = [order_atom[i] for i, atom in enumerate(atoms_mask) if atom == 1]
                for atom in list_of_atoms_in_residue:
                    atom_num_int = atom_num_int + 1
                    atom_remark = 'ATOM'.ljust(6)
                    atom_num = str(atom_num_int).rjust(5)
                    atom_name = atom.ljust(4)
                    res_name = res_type.ljust(3)
                    chain_id = pdb_name.decode('utf-8').split('_')[-1].rjust(1)
                    # chain_id = pdb_name.decode('utf-8').split('_')[1].rjust(1)
                    res_num = str(index + 1).rjust(4)
                    x_coord = str('%8.3f' % (float(str(
                        self.template_features['template_all_atom_positions'][template_domain_index][index][
                            atom_types.index(atom)][
                            0])))).rjust(8)
                    y_coord = str('%8.3f' % (float(str(
                        self.template_features['template_all_atom_positions'][template_domain_index][index][
                            atom_types.index(atom)][
                            1])))).rjust(8)
                    z_coord = str('%8.3f' % (float(str(
                        self.template_features['template_all_atom_positions'][template_domain_index][index][
                            atom_types.index(atom)][
                            2])))).rjust(8)
                    occ = str('%6.2f' % (float('1.0'))).rjust(6)
                    bfact = str('%6.2f' % (float('25.0'))).ljust(6)
                    atom_type = atom[0].rjust(12)
                    output_pdb.write(
                        f'{atom_remark}{atom_num}  {atom_name}{res_name} {chain_id}{res_num}    {x_coord}{y_coord}{z_coord}{occ}{bfact}{atom_type}\n')

            output_pdb.close()

    def complete_msa_from_template_features(self, template_features):

        msa_from_templates_list = [(''.join(residue_constants.ID_TO_HHBLITS_AA[res] for res in self.msa_features['msa'][0]), 'Query')]
        for num, seq in enumerate(template_features['template_sequence']):
            msa = ''.join([f'>{seq[1]}\n{seq[0]}\n' for seq in msa_from_templates_list])
            if seq.decode('utf-8') not in msa:
                msa_from_templates_list.append((seq.decode('utf-8'), template_features['template_domain_names'][num].decode('utf-8')))

        for seq in msa_from_templates_list[1:]:
            self.append_row_in_msa(sequence=seq[0], msa_uniprot_accession_identifiers=seq[1])

