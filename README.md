Requires HHsearch, maxit, pisa and superpose installation.
Using AlphaFold2 and ALEPH.

Usage: arcimboldo_air.py configuration.bor

The configuration.bor has to be a file with YAML format and can contain the following
parameters:

  output_dir (mandatory, string): Path to the results' directory.
  af2_dbs_path (mandatory, string): Path to the AlphaFold2 databases, they have to be downloaded.
  run_dir (optional, string, 'run'): Path to the directory where AlphaFold2 will be run.
  verbose (optional, bool, true): Enable debugging.
  glycines (optional, integer, 50): Number of glycines between sequence copies.
  small_bfd (optional, bool, false): Use reduced bfd library.
  run_af2 (optional, bool, true): Run AlphaFold2, it can be used to generate a features.pkl without going further.
  stop_after_msa (optional, bool, false): Run AlphaFold2 and stop it after generating the msa.
  reference (optional, string, ''): Existing pdbid or path to a pdb
  experimental_pdb (optional, string, ''): Existing pdbid or path to a pdb
  custom_features (optional, bool, true): Run AlphaFold2 without any modification and with AlphaFold2 generated features.pkl
  mosaic (optional, integer, None): Split the sequence in X partitions.
  cluster_templates (optional, bool, false): Group templates by distance and relaunch arcimboldo_air with them.

features:
- path:
  keep_msa:
  keep_templates:

sequences:
- fasta_path:
  num_of_copies:
  positions:
  name:

- fasta_path:
  num_of_copies:

Templates can be added to the templates section inside the features.pkl, as well as their sequence to the msa, it is possible to add as many templates as the user requires:

templates:
- pdb (mandatory, string): Existing pdbid or path to a pdb
  add_to_msa (optional, bool, false): Add template to the msa.
  add_to_templates (optional, bool, true): Add template to the features.pkl
  generate_multimer (optional, bool, true if num_of_copies > 1 else false):
  sum_prob (optional, integer, None):
  aligned (optional, bool, false): If the template has already been aligned with the sequence.
  legacy (optional, bool, false): If the template has been prepared (aligned, one chain)
  reference (optional, string, ''): Existing pdbid or path to a pdb
  
  change_res: -> Change residues of the template. It can be a chain or 'ALL' so this can be applied to all chains
    - {chain_name or ALL}: (mandatory, range): Selection of residues to apply the modification
      resname (mandatory, string): Residue name

  match: -> Set restrictions in order to insert the template into the sequence copies
    - chain (mandatory, string): Set the position of the chain
      position: (optional, string, None, X): Set a specific position
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
      fasta_path:
      when: after_alignment/before_alignment
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
