Requires maxit, pisa, superpose and AlphaFold2 installation

Usage: arcimboldo_air.py configuration.bor

The configuration.bor has to be a file with YAML format and can contain the following
paramaters:

output_dir (mandatory, string): Path to the results directory.
fasta_path (mandatory, string): Path to the sequence file.
num_of_copies (mandatory, integer): Number of sequence repetitions.
af2_dbs_path (mandatory, string): Path to the AlphaFold2 databases, they have to be downloaded.
run_dir (optional, string, 'run'): Path to the directory where AlphaFold2 will be run.
verbose (optional, bool, true): Enable debugging.
glycines (optional, integer, 50): Number of glycines between sequence copies.
run_alphafold (optional, bool, true): Run AlphaFold2, it can be used to generate a features.pkl without going further.
reference (optional, string, ''): Existing pdbid or path to a pdb
experimental_pdb (optional, string, ''): Existing pdbid or path to a pdb
use_features (optional, bool, false): Run AlphaFold2 without any modification and with AlphaFold2 generated features.pkl


Templates can be added to the templates section inside the features.pkl, as well as their sequence to the msa, 
it is possible to add as many templates as the user requiers:

templates:
- pdb (mandatory, string): Existing pdbid or path to a pdb
  add_to_msa (optional, bool, false): Add tempalte to the msa.
  add_to_templates (optional, bool, true): Add template to the features.pkl
  generate_multimer (optional, bool, true if num_of_copies > 1 else false):
  sum_prob (optional, integer, ''):
  aligned (optional, bool, false): If the template has already been aligned with the sequence.
  reference (optional, string, ''): Existing pdbid or path to a pdb
  
  change_res: -> It can be a chain or 'ALL' so this can be applied to all chains
    - All: int -> separated with ','. It can also contain ranges: 10-20
      B:
      resname: str -> string with residue name

  match:
    - chain A
      position: None, Any, X
      residues: 100-120
      reference:
      reference_chain:




Example of a configuration.bor:


templates:
- pdb: 3fxq
  add_to_msa: false
  add_to_templates: true
  generate_multimer: false
  change_res:
    - resname: 'ALA'
      B: 10-50
  match:
      - chain: B
        position: 1
        residues: 30
      - chain: A
        position: 2
      - chain: A
        position: 3

- pdb: 3fzv
  add_to_msa: false
  add_to_templates: true
  generate_multimer: true
  match:
    - chain: B
      position: 1
    - chain: B
      position: 1

- pdb: 4x6g
  add_to_msa: false
  add_to_templates: true
  generate_multimer: true

