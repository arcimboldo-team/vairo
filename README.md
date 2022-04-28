# ARCIMBOLDO_AIR

ARCIMBOLDO_AIR code

## Current modes:

1. **Load features from an existing features.pkl file** (`--feature_mode=from_pkl`):

    Mandatory flags: `--feature_mode`

2. **Generate features for a specific chain in PDB** (`--feature_mode=from_chain_in_pdb`):

    Mandatory flags: `--feature_mode`, `--output_dir`, `--fasta_path`, `--template_pdb_id`, `--template_pdb_chain_id`

    Optional flags: `--hhr_path`

3. **Generate features for assembly in specific PDB defined through its REMARK 350** (`--feature_mode=from_assembly_in_pdb`): 

    Number of glycines used as linkers: 50.

    Mandatory flags: `--feature_mode`, `--output_dir`, `--fasta_path`, `--template_pdb_id`

    Optional flags: `--hhr_path`

4. **Generate features for custom PDB** (`--feature_mode=from_custom_pdb`):

    Sequence from custom PDB must be aligned with query sequence.

    Mandatory flags: `--feature_mode`, `--output_dir`, `--fasta_path`, `--custom_pdb_path`

    Optional flags: `--poly_ala_list_of_res_ranges`

    



 
