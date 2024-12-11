Requires HH-suite and CCP4 suite installation.
Using AlphaFold2 and ALEPH.

Usage: run_vairo.py [-h] [-check] input_path

Flags:
 [-check]: Check the .bor information and if is well parsed.


Configuration File
-------------------
The configuration file must be in YAML format and can contain the following parameters:

Mandatory parameters:
  mode (string): Options are {naive, guided}.
  output_dir (string): Path to the results' directory.
  af2_dbs_path (string): Path to the AlphaFold2 databases, which must be downloaded.

Optional parameters:
    run_dir (string, default 'run'): Path to the directory where AlphaFold2 will be run.
    glycines (integer, default 50): Number of glycines between sequence copies.
    small_bfd (bool, default False): Use reduced bfd library.
    run_af2 (bool, default True): Run AlphaFold2; can be used to generate a features.pkl without going further.
    stop_after_msa (bool, default False): Run AlphaFold2 and stop it after generating the MSA.
    reference (string, default ''): Existing pdbid or path to a pdb.
    experimental_pdbs (List, default []): List of existing pdbid or paths to pdbs.
    mosaic (integer, default None): Split the sequence into X partitions.
    mosaic_partition (range, default None): Split the sequence by the number of residues.
    mosaic_seq_partition (range, default None): Split the sequence by the number of sequences.
    cluster_templates (bool, default False, True if naive): Cluster templates from preprocessed features.pkl.
    cluster_templates_msa (int, default -1): Select a specific number of sequences to add to the new features.pkl MSA (-1 to add all sequences).
    cluster_templates_msa_mask (range, default None): Delete specific residues from the MSA sequences.
    cluster_templates_sequence (path, default None): Replace template sequence with sequence in path.
    show_pymol (str, default None): Select regions to zoom in the pymol session using selection language in pymol, separated by commas. 

Several sequences can be added to the query sequence. Each FASTA file can contain multiple copies of a sequence, which can be inserted at specific positions within the query sequence. All sequences will be concatenated using glycines.
sequences: List of sequences.
    fasta_path (str, default None): Path to the FASTA file.
    num_of_copies (int, default None): Number of copies of the same sequence.
    positions (List[int], default None): Positions in the query sequence.
    name (str, default file name from fasta_path): Name of the sequence.
    mutations: List of mutations.
        Format: 'aminoacid': residues. One-letter code of the amino acid, and the residues that will be mutated to this amino acid.
        Example: 'G': 10, 20 
    predict_region (range, default None): Predict a certain region of the sequence instead of the entire sequence.

An already extracted library can be used as input. It can contain FASTA files or PDBs. It is possible to insert it into specific regions of the query sequence.
append_library:
  - path: (string): Path to a folder, PDB, or FASTA. The folder can contain FASTA or PDB files as well.
    add_to_msa: (bool, default True): Append the sequences of the PDB or FASTA files to the MSA.
    add_to_templates: (bool, default False): Append the PDB files to the templates (only if they are PDBs).
    numbering_query (List[int], default None): Positions of the query sequence where the library will be inserted.
    numbering_library (range, default None): Positions of the library; those residues will be selected and inserted in the positions specified by the parameter numbering_query

Multiple features.pkl files can be combined into a final features.pkl file. You can select the number of templates and sequences to be inserted, delete specific regions of the MSA, replace the sequences of the templates, and insert the feature.pkl file into specific regions of the query sequence.
features:
- path: (mandatory, string): Path to a features.pkl
  keep_msa: (int, default -1): Indicates the number of sequences to select from the MSA. If -1, it will use all the sequences. For any other value, it will use the X sequences that have more coverage with the query sequence.
  keep_templates: (int, default -1): Indicates the number of templates to select from the features.pkl. If -1, it will use all the templates. For any other value, it will use the X templates that have more coverage with the query sequence.
  msa_mask: (range, default None): Delete specific residues from the MSA.
  sequence: (str, default None): Replace all the template sequences with the sequence in the specified path.
  numbering_query: (List[int], default 1): Positions where the features will be placed. By default, it starts at position 1, but multiple positions can be specified (e.g., 1, 50, 100) and split the features.pkl using the numbering_features keyword.
  numbering_features: (List[ranges], default None): The specified ranges will be placed in the positions specified by numbering_query.
  positions: (range, default None): Insert the features.pkl in the query sequence. The position refers to the sequence position, whereas in numbering_query and numbering_features, it refers to the residue positions in the entire query sequence.


Templates can be added to the templates section inside the features.pkl, as well as their sequences to the MSA. It is possible to add as many templates as the user requires.
templates:
- pdb (str): Path to a PDB file or an existing PDB ID.
  add_to_msa (bool, default False): Add the template's sequence to the MSA.
  add_to_templates (bool, default True): Add the template to the features.pkl.
  generate_multimer (bool, default True): Generate the multimer from the PDB.
  strict (bool, default True): Check the E-values of the alignment. If activated, it will discard the template if the alignment E-values are below threshold.
  aligned (bool, default False): If the template has already been aligned with the sequence, activate this to skip the alignment
  legacy (bool, default False): If the template has been prepared (aligned and in one chain) to match the whole query sequence.
  reference (str, default None): Reference to be used in order to insert it into the query sequence.
  
  modifications: List of modifications for the template. Modify a chain/chains in a PDB file with various options.
     - chain (str, default None): Select different chains, a single chain or use the keyword All in order to modify all the chains in the PDB.
       position (int, default None): Insert the chain (mandatory a single chain) into a specific sequence in the query sequence.
       maintain_residues (List[int], default None): The selected residues will be kept, and the rest will be discarded.
       delete_residues (List[int], default None): The selected residues will be discarded, the rest will be kept.
       when (str, default 'after_alignment'): When the modification will be applied, 'before_alignment' or 'after_alignment'
       mutations: List of mutations.
          - numbering_residues (List[int], default None): Residues where the mutations will be applied
            mutate_with (str, default None): The amino acid to mutate with, specified as a three-letter code or a FASTA file path.


PATHS:
All information can be found in the output_dir directory, which is an input parameter in the configuration file. Inside the output_dir, we can find the following folders and files:
- html: The results html, with all the plots, statistics about the run and the analysis of the predictions.
- log: The log file, containing all the info from the execution.
- plots: All the plots generated by the output analysis.
- frobenius: All the plots generated by aleph.
- interfaces: The results of interfaces analysis done by PISA. The output txt file and the interfaces PDBs.
- clustering: If clustering has been activated, the clusters jobs can be found inside the directory.
- input: All the input files used in the run.
- run: All the run information.
- templates: Extracted templates from the features.pkl, split into chains.
- rankeds: Generated rankeds by AlphaFold2, split into chains (if necessary).
Inside run directory, We can find the following folders:
- results: Results of the AlphaFold2 run.
- Templates folder: Named by each template name. It contains the databases generated to align the template.
- Sequences folder: Named by each sequence. It contains the alignments of the templates with that sequence.

Inside the results directory, We can find the following folders:
- tmp: Information generated by running external programs, like Aleph.
- ccanalysis/ccnalaysis_ranked: PDBs used for the ccanalysis and the results.
- msas: Information generated by AlphaFold2. It contains the extracted sequences and the template alignments.
- templates_nonsplit: Extracted templates from the features.pkl, not split into chains.
- rankeds_split: Generated rankeds by AlphaFold2, split in chains (if necessary)
- rankeds: Generated rankeds by AlphaFold2, not split in chains.
