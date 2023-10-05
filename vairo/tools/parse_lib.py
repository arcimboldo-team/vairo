from Bio.PDB import PDBParser, PDBIO, Selection
from io import StringIO
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue


def get_coordinates(struct):
    for chain in struct[0]:
        for residue in chain:
            return residue['CA']

def read_file(path):
    parser = PDBParser()
    content = open(path, 'r').readlines()
    filtered_lines = []
    splitted_file = []
    atom_started = False
    for line in content:
        if line.startswith("ATOM"):
            atom_started = True
            filtered_lines.append(line)
        elif atom_started and line.startswith("TER"):
            splitted_file.append(parser.get_structure("pdb", StringIO(''.join(filtered_lines))))           
            filtered_lines = []
            atom_started = False

    return splitted_file

# Your PDB data
pdb_path_reference = "/Users/pep/work/transfers/BETA_lib_uud/1/0/1vdv_0_1.pdb"
pdb_data_reference = read_file(pdb_path_reference)

pdb_path_obj = "/Users/pep/work/transfers/BETA_lib_uud/1/0/3tmr_0_66.pdb"
pdb_data_obj = read_file(pdb_path_obj)


# Initialize a PDB parser

correct_order = []
for segment in pdb_data_obj:
    min_distance = float('inf')  # Initialize with a large value
    closest_coordinate = None
    coord_obj = get_coordinates(segment)
    for num, reference_struct in enumerate(pdb_data_reference):
        distance = coord_obj - get_coordinates(reference_struct)
        if distance < min_distance:
            min_distance = distance
            closest_coordinate = num
    correct_order.append(closest_coordinate)

merged_structure = Structure("merged")
merged_model = Model(0)
merged_chain = Chain('A')


residue_id = 1
for order in correct_order:
    add_struct = pdb_data_obj[order]
    for residue in Selection.unfold_entities(add_struct, 'R'):
        new_residue = Residue((' ', residue_id, ' '), residue.get_resname(), '')
        
        # Copy atoms from the original residue to the new residue
        for atom in residue:
            new_atom = atom.copy()
            new_atom.set_parent(new_residue)
            new_residue.add(new_atom)
        # Add the new residue to the merged chain
        merged_chain.add(new_residue)
        # Increment residue ID
        residue_id += 1

# Add the merged chain to the model and model to the structure
merged_model.add(merged_chain)
merged_structure.add(merged_model)

# Save the merged and renumbered structure to a new PDB file
io = PDBIO()
io.set_structure(merged_structure)
io.save("merged_structure.pdb")
    