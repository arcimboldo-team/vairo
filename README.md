Requires HHsearch, maxit, pisa and superpose installation.
Using AlphaFold2 and ALEPH.

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
  custom_features (optional, bool, true): Run AlphaFold2 without any modification and with AlphaFold2 generated features.pkl
  mosaic (optional, integer, None): Split the sequence in several pieces. Each piece will have mosaic size.

Templates can be added to the templates section inside the features.pkl, as well as their sequence to the msa, 
it is possible to add as many templates as the user requiers:

templates:
- pdb (mandatory, string): Existing pdbid or path to a pdb
  add_to_msa (optional, bool, false): Add tempalte to the msa.
  add_to_templates (optional, bool, true): Add template to the features.pkl
  generate_multimer (optional, bool, true if num_of_copies > 1 else false):
  sum_prob (optional, integer, None):
  aligned (optional, bool, false): If the template has already been aligned with the sequence.
  reference (optional, string, ''): Existing pdbid or path to a pdb
  
  change_res: -> Change residues of the template. It can be a chain or 'ALL' so this can be applied to all chains
    - {chain_name or ALL}: (mandatory, range): Selection of residues to apply the modification
      resname (mandatory, string): Residue name

  match: -> Set restrictions in order to insert the template into the sequence copies
    - chain (mandatory, string): Set the position of the chain
      position: (optional, string, None, X): Set an specific position
      residues: (optional, int range, ''): Selection of residues to set in a position. Can be a range or an integer (Ex: 100-120, 100), otherwise, the whole chain is going to be selected.
      reference:  (optional, string, ''): Existing pdbid or path to a pdb
      reference_chain: (optional, string, ''): Existing reference chain. The match chain will be fixed in the same position as the reference chain.

Example of a configuration.bor:

templates:
- pdb: 3fxq
  add_to_msa: false
  add_to_templates: true
  generate_multimer: false
  change_res:
    - resname: 'ALA'
      B: 10-50
      A: 1-10
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

- pdb: 4x6g
  add_to_msa: false
  add_to_templates: true
  generate_multimer: true

