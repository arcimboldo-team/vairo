import logging
import subprocess
from Bio import PDB
from ALPHAFOLD.alphafold.common import residue_constants

def generate_hhsearch_db(template_cif_path: str, output_dir: str):

    with open(f"{output_dir}/pdb70_a3m.ffdata", "w") as a3m, \
         open(f"{output_dir}/pdb70_cs219.ffindex", "w") as cs219_index, \
         open(f"{output_dir}/pdb70_a3m.ffindex", "w") as a3m_index, \
         open(f"{output_dir}/pdb70_cs219.ffdata", "w") as cs219:
        id = 1000000
        index_offset = 0

        parser = PDB.MMCIFParser(QUIET=True)
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

def run_hhsearch(fasta_path: str, pdb70_db: str, output_path: str) -> str:

    logging.info(f'Running hhsearch using {pdb70_db} as database.')

    out = subprocess.Popen(['hhsearch', '-i', fasta_path, '-o', output_path, '-maxseq',
                            '1000000', '-d', pdb70_db, '-glob'],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = out.communicate()
    hhr = stdout.decode('utf-8')

    return hhr