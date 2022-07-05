from Bio.PDB import PDBParser, Selection
from libs.features import *
from Bio.PDB.PDBIO import PDBIO
from alphafold.data import parsers
import glob


def extract_query_sequence(query_fasta_path):

    logging.info(f'Extracting query sequence from {query_fasta_path}')

    with open(query_fasta_path, 'r') as f:
        fasta_lines = f.readlines()
        fasta_name, query_sequence = fasta_lines[0][1:-1], fasta_lines[1].split('\n')[0]

    return fasta_name, query_sequence


def split_pdb_by_chains(pdb_in, output_dir):

    parser = PDBParser(QUIET=1)
    io = PDBIO()

    name = pdb_in.split('/')[-1][:-4]

    structure = parser.get_structure(name, pdb_in)
    pdb_chains = structure.get_chains()

    splitted_pdbs = []
    for chain in pdb_chains:
        io.set_structure(chain)
        io.save(f'{output_dir}/{name}_{chain.get_id()}.pdb')
        splitted_pdbs.append(f'{name}_{chain.get_id()}')

    return splitted_pdbs


def generate_hhsearch_db(pdb_path, output_dir):

    with open(f"{output_dir}/pdb70_a3m.ffdata", "w") as a3m, \
         open(f"{output_dir}/pdb70_cs219.ffindex", "w") as cs219_index, \
         open(f"{output_dir}/pdb70_a3m.ffindex", "w") as a3m_index, \
         open(f"{output_dir}/pdb70_cs219.ffdata", "w") as cs219:
        id = 1000000
        index_offset = 0

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("test", pdb_path)
        models = list(structure.get_models())
        if len(models) != 1:
            raise ValueError(f"Only single model PDBs are supported. Found {len(models)} models.")
        model = models[0]
        for chain in model:
            amino_acid_res = []
            for res in chain:
                amino_acid_res.append(residue_constants.restype_3to1.get(res.resname, "X"))
            protein_str = "".join(amino_acid_res)
            a3m_str = f">{pdb_path.split('/')[-1][:-4]}_{chain.id}\n{protein_str}\n\0"
            a3m_str_len = len(a3m_str)
            a3m_index.write(f"{id}\t{index_offset}\t{a3m_str_len}\n")
            cs219_index.write(f"{id}\t{index_offset}\t{len(protein_str)}\n")
            index_offset += a3m_str_len
            a3m.write(a3m_str)
            cs219.write("\n\0")
            id += 1


def run_hhsearch(query_fasta_path, pdb70_db, output_path):

    print(f'Running hhsearch using {pdb70_db} as database.')

    out = subprocess.Popen(['hhsearch', '-i', query_fasta_path, '-o', output_path, '-maxseq',
                            '1000000', '-d', pdb70_db, '-glob'],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = out.communicate()
    hhr = stdout.decode('utf-8')

    return hhr


def mapping_between_template_and_query_sequence(hhr_path, pdb_id, chain):

    hhr_text = open(hhr_path, 'r').read()

    matches = re.finditer(r'No\s+\d+', hhr_text)
    matches_positions = [match.start() for match in matches] + [len(hhr_text)]

    detailed_lines_list = []
    for i in range(len(matches_positions) - 1):
        detailed_lines_list.append(hhr_text[matches_positions[i]:matches_positions[i + 1]].split('\n')[:-3])

    hits_list = [detailed_lines for detailed_lines in detailed_lines_list if
                  detailed_lines[1][1:] == f'{pdb_id}_{chain}']

    best_hit = hits_list[0]

    hit = parsers._parse_hhr_hit(best_hit)
    mapping = templates._build_query_to_hit_index_mapping(hit.query, hit.hit_sequence, hit.indices_hit,
                                                          hit.indices_query, hit.query.replace('-', ''))

    inv_mapping = {v: k for k, v in mapping.items()}

    return inv_mapping


def aligned_template_to_query_sequence(mapping, pdb_in_path, pdb_out_path):

    with open(pdb_out_path, 'w') as f_out:

        with open(pdb_in_path, 'r') as f_in:
            lines = f_in.readlines()
            for line in lines:
                if line.startswith('ATOM'):
                    old_res_num = int(line[22:26])
                    if old_res_num in mapping:
                        f_out.write(line[:22] + str(mapping[old_res_num]).rjust(4) + line[26:])


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

    transformations_dict = {}
    for chain in chain_list:
         transformations_dict[chain] = transformations_list

    return transformations_dict


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
        f.write(f'chain A\n')
        f.write('end\n')
        f.write('eof')
    subprocess.Popen(['bash', '/tmp/pdbset.sh'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
    os.system(f'rm /tmp/pdbset.sh')


def superpose_pdbs(query_pdb, target_pdb, output_superposition=True):

    if output_superposition:
        superpose_output = subprocess.Popen(['superpose', f'{query_pdb}', '-s', '-all', f'{target_pdb}', '-s', '-all',
                                             '-o', f'{query_pdb[:-4]}_superposed.pdb'],
                                            stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
    else:
        superpose_output = subprocess.Popen(['superpose', f'{query_pdb}', '-s', '-all', f'{target_pdb}', '-s', '-all'],
                                            stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
    for line in superpose_output.split('\n'):
        if 'r.m.s.d:' in line:
            rmsd = float(line.split()[1])
        if 'quality Q:' in line:
            quality_q = line.split()[2]
        if 'Nalign:' in line:
            nalign = line.split()[1]

    try:
        match1 = [m.start() for m in re.finditer("TEXT:Residue alignment:", superpose_output)][0]
        match2 = [m.start() for m in re.finditer("`-------------'----------'-------------'", superpose_output)][0]
        alignment_output = superpose_output[match1:match2].split('\n')[5:]
        aligned_res_list = []
        for line in alignment_output:
            if line[23:25] == '**':
                aligned_res_list.append(int(line[36:39].replace(' ', '')))
    except:

        rmsd, nalign, quality_q, aligned_res_list = None, None, None, None


    return rmsd, nalign, quality_q, aligned_res_list







OUTPUT_DIR = '/cri4/albert/Desktop/OUTPUT'
TEMPLATES_DIR = '/cri4/albert/Desktop/TEMPLATES'
QUERY_FASTA_PATH = '/cri4/albert/Desktop/atzr.fasta'
COPIES = 4

isExist = os.path.exists(OUTPUT_DIR)
if not isExist:
    os.makedirs(OUTPUT_DIR)



list_of_templates = os.listdir(TEMPLATES_DIR)
# for template_name in list_of_templates:
#     isExist = os.path.exists(f'{OUTPUT_DIR}/{template_name[:-4]}')
#     if not isExist:
#         os.makedirs(f'{OUTPUT_DIR}/{template_name[:-4]}')
#
#     input_pdb_path = f'{TEMPLATES_DIR}/{template_name}'
#
#
#     transformations_dict = read_remark_350(pdb_path=input_pdb_path, use_pisa=False)
#
#     splitted_pdbs = split_pdb_by_chains(pdb_in=input_pdb_path, output_dir=f'{OUTPUT_DIR}/{template_name[:-4]}')
#     generate_hhsearch_db(pdb_path=input_pdb_path, output_dir=f'{OUTPUT_DIR}/{template_name[:-4]}')
#     run_hhsearch(query_fasta_path=QUERY_FASTA_PATH, pdb70_db=f'{OUTPUT_DIR}/{template_name[:-4]}/pdb70', output_path=f'{OUTPUT_DIR}/{template_name[:-4]}/output.hhr')
#
#
#     for pdb in splitted_pdbs:
#         pdb_id, chain = pdb.split('_')
#         mapping = mapping_between_template_and_query_sequence(hhr_path=f'{OUTPUT_DIR}/{template_name[:-4]}/output.hhr', pdb_id=pdb_id, chain=chain)
#         aligned_template_to_query_sequence(mapping=mapping,
#                                            pdb_in_path=f'{OUTPUT_DIR}/{template_name[:-4]}/{pdb}.pdb',
#                                            pdb_out_path=f'{OUTPUT_DIR}/{template_name[:-4]}/{pdb}_aligned.pdb')
#
#     counter = 0
#     for chain in transformations_dict:
#         for trans in transformations_dict[chain]:
#             counter += 1
#             rot_and_trans(pdb_path=f'{OUTPUT_DIR}/{template_name[:-4]}/{template_name[:-4]}_{chain}_aligned.pdb', out_pdb_path=f'{OUTPUT_DIR}/{template_name[:-4]}/{counter}.pdb', rot_tra_matrix=trans)

reference_template = list_of_templates[0][:-4]

for template_name in list_of_templates[1:]:
    print('\n')
    for query_pdb in glob.glob(f'{OUTPUT_DIR}/{template_name[:-4]}/[0-90].pdb'):
        for target_pdb in glob.glob(f'{OUTPUT_DIR}/{reference_template}/[0-90].pdb'):
            rmsd, nalign, quality_q, aligned_res_list = superpose_pdbs(query_pdb=query_pdb, target_pdb=target_pdb, output_superposition=False)
            print(query_pdb, target_pdb, rmsd, nalign)





