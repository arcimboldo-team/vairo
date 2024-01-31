Requires HH-suite and CCP4 suite installation.
Using AlphaFold2 and ALEPH.

Usage: vairo.py configuration.bor

The configuration.bor has to be a file with YAML format and can contain the following
parameters:

  mode (mandatory, string): {naive, guided}
  output_dir (mandatory, string): Path to the results' directory.
  af2_dbs_path (mandatory, string): Path to the AlphaFold2 databases, they have to be downloaded.
  run_dir (optional, string, 'run'): Path to the directory where AlphaFold2 will be run.
  glycines (optional, integer, 50): Number of glycines between sequence copies.
  small_bfd (optional, bool, false): Use reduced bfd library.
  run_af2 (optional, bool, true): Run AlphaFold2, it can be used to generate a features.pkl without going further.
  stop_after_msa (optional, bool, false): Run AlphaFold2 and stop it after generating the msa.
  reference (optional, string, ''): Existing pdbid or path to a pdb
  experimental_pdbs (optional, List, []): List of existing pdbid or path to a pdb
  mosaic (optional, integer, None): Split the sequence in X partitions.
  mosaic_partition (optional, range, None): Split the sequence by the number of residues.
  mosaic_seq_partition (optional, range, None): Split the sequence by the number of sequences.
  cluster_templates (optional, bool, false, true if naive): Cluster templates from preprocessed features.pkl
  cluster_templates_msa (optional, int, -1): Select a specific number of sequences to add to the new features.pkl MSA (-1 to add all the sequences) from the preprocessed features.pkl MSA
  cluster_templates_msa_mask: (optional, range, None): Delete specific residues from the MSA sequences.
  cluster_templates_sequence: (optional, path, None): Replace templates sequence with sequence in path.
  show_pymol: (optional, range, None): Select regions to zoom in the pymol session.

Several sequences can be added to the query sequence. Each Fasta can have several copies and can be inserted in a specific position
inside the query sequence. All the sequences will be concatenated using Glycines.
sequences:
- fasta_path:
  num_of_copies:
  positions:
  name:
  mutations:
    - 'G': 10, 20

Several features.pkl can be merged to the final features.pkl. It is possible to select the amount of templates and sequences to be inserted.
Delete specific regions of the MSA, replace the sequence of the templates and insert the features in specific regions in the final features.

append_library:
  - path: (mandatory, string): Path to a folder, pdb or fasta. The folder can contain fastas or pdbs too.
    add_to_msa: (optional, bool, True): Append the sequences of the pdbs or fastas to the MSA
    add_to_templates: (optional, bool, False): Only if they are pdbs. Append the pdbs to the templates.
    positions: Dictionary with the positions of the input pdbs to the query sequence
               {Select positions of the pdb or fasta (1-4)}: {Put those residues in the selected positions (5-8)}

features:
- path: (mandatory, string): Path to features.pkl
  keep_msa: (optional, int, -1): Keep the msa of the features.pkl 
  keep_templates: (optional, int, -1): Keep the templates of the features.pkl
  msa_mask: (optional, range, None): Delete specific residues from the MSA sequences.
  sequence: (optional, path, None): Replace template sequence with sequence in path.
  positions: (optional, range, None): Insert the features.pkl in a specific region inside the query sequence.

Templates can be added to the templates section inside the features.pkl, as well as their sequence to the msa, it is possible to add as many templates as the user requires:

templates:
- pdb (mandatory, string): Existing pdbid or path to a pdb
  add_to_msa (optional, bool, False): Add template to the msa.
  add_to_templates (optional, bool, True): Add template to the features.pkl
  generate_multimer (optional, bool, True):
  sum_prob (optional, integer, False):
  strict (optional, bool True): Check the evalues of the alignment
  aligned (optional, bool, false): If the template has already been aligned with the sequence.
  legacy (optional, bool, false): If the template has been prepared (aligned, one chain)
  reference (optional, string, ''): Existing pdbid or path to a pdb
  
  change_res: -> Change residues of the template. It can be a chain or 'ALL' so this can be applied to all chains
    - {chain_name or ALL}: (mandatory, range): Selection of residues to apply the modification
      resname (optional, string): Residue name
      fasta_path (optional, string): Fasta path to replace the sequence
      when (optional, string, after_alignment): When change the chain, before_alignment or after_alignment

  match: -> Set restrictions in order to insert the template into the sequence copies
    - chain (mandatory, string): Set the position of the chain
      position: (optional, string, None, X): Set a specific position
      residues: (optional, int range, ''): Selection of residues to set in a position. Can be a range or an integer (Ex: 100-120, 100), otherwise, the whole chain is going to be selected.
      reference:  (optional, string, ''): Existing pdbid or path to a pdb
      reference_chain: (optional, string, ''): Existing reference chain. The match chain will be fixed in the same position as the reference chain.


PATHS:
All information can be found in the output_dir directory, which is an input parameter in the configuration file. Inside the output_dir
we can find the following folders and files:
- html: The results html, with all the plots, statistics about the run and the analysis of the predictions.
- log: The log file, containing all the info from the execution.
- plots: All the plots generated by the output analysis.
- frobenius: All the plots generated by aleph.
- interfaces: The results of interfaces analysis done by PISA. The output txt file and the interfaces PDBs.
- clustering: If clustering has been activated, the clusters jobs can be found inside the directory.
- input: All the input files used in the run.
- run: All the run information.
- templates: Extracted templates from the features.pkl. split in chains
- rankeds: Generated rankeds by AlphaFold2, split in chains (if necessary)
Inside run directory we can find the following folders:
- results: Results of the AlphaFold2 run.
- Templates folder: Named by each template name. It has the databases generated to align the template.
- Sequences folder: Named by each sequence. It has the alignments of the templates with that sequence.

Inside the results directory we can find the following folders:
- tmp: Information generated by running external programs, like aleph.
- ccanalysis/ccnalaysis_ranked: PDBs used for the ccanalysis and the results.
- msas: Information generated by AlphaFold2, it has the extracted sequences and the template alignments
- templates_nonsplit: Extracted templates from the features.pkl. not split in chains
- rankeds_split: Generated rankeds by AlphaFold2, split in chains (if necessary)
- rankeds: Generated rankeds by AlphaFold2, not split in chains.