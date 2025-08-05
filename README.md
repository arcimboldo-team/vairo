# VAIRO
Guiding structural model predictions with experimental information

-------------------

## Prerequisites
* AlphaFold2
* HH-suite
* CCP4 suite
* ALEPH
* MAXIT

-------------------

## Installation

In order to install VAIRO and its interface VAIROGUI, you need to run the installer script located in tools/install_vairo.sh. This script handles conda setup and installs all VAIRO dependencies within a dedicated environment.

Execute the installer script:
```
bash tools/install_vairo.sh
```

The script will:
1. Check for an existing conda installation and install it if missing.
2. Create and activate a conda environment for VAIRO.
3. Install all Python and system dependencies required by VAIRO.
4. Verify that system libraries like CUDA drivers or MAXIT are already present.

-------------------

## Usage

To run the command-line program:

```
vairo [-h] [-check] <config.yaml>
```

| Flag      | Description                                    |
|-----------|------------------------------------------------|
| `-h`      | Show help and exit                             |
| `-check`  | Validate configuration (*.yml*) file parsing   |    


To launch the graphical interface:
```
vairogui
```

-------------------

## Configuration File (YAML)
-------------------
The configuration file must be in valid YAML. Below are all supported sections and parameters.

### 1. Mandatory keys

    mode (string) Choose one of: naive, guided.
    output_dir (string) Directory where results will be saved.
    af2_dbs_path (string) Path to the AlphaFold2 databases (must be pre-downloaded).

### 2. Common optional keys

    run_dir (string, default: "run") Directory where AlphaFold2 jobs will run.
    glycines (integer, default: 50) Number of glycine residues to insert between concatenated sequences.
    small_bfd (boolean, default: false) Use reduced BFD library.
    run_af2 (boolean, default: true) Run AlphaFold2 (otherwise stop after generating features.pkl file).
    stop_after_msa (boolean, default: false) Run AlphaFold2 up to MSA generation, then exit.
    reference (string, default: "") PDB ID or path to PDB file to be used as global reference.
    experimental_pdbs (list of strings, default: []) List of PDB IDs or paths to PDB files for result comparison.
    mosaic (integer, default: null) Split the sequence into X partitions.
    mosaic_partition (range, default: null) Residue based partitioning.
    mosaic_seq_partition (range, default: null) Sequence numbering partitioning.
    cluster_templates (boolean, default: false - becomes true if mode: naive) Cluster templates from preprocessed features.pkl.
    cluster_templates_msa (integer, default: -1) Number of sequences to add to the MSA (-1 = all).
    cluster_templates_msa_mask (sequence range, default: null) Remove specific residues from MSA sequences.
    cluster_templates_sequence (string path, default: null) Replace templates sequences using FASTA at given path.
    show_pymol (string, default: null) Pymol selection string (comma-separated regions) to zoom into.


### 3. Query sequence
Define one or more sequences to generate the query sequence. All sequences will be concatenated using glycine linkers.
```
sequences:
    - fasta_path (string, mandatory) Path to the FASTA file.
      num_of_copies (integer, default: 1) Number of copies of the sequence.
      positions (list of integers, default: [], any position) Insertion position in the query.
      name (string, default: file name from fasta_path) Sequence name.
      predict_region (range, default: null) Predict only this subsequence instead of the full length.
      mutations (map) Map three-letter amino acid codes to residue indices. Example:
        - 'ALA': 10, 20
```

### 4. Add templates
Customize PDB templates for insertion into features.pkl.
```
templates:
    - pdb (string, mandatory) Path to a PDB file or existing PDB ID.
      add_to_msa (boolean, default: false) Add the template’s sequence to the MSA.
      add_to_templates (boolean, default: true) Include the template in features.pkl.
      generate_multimer (boolean, default: true) Generate a multimeric assembly from the PDB.
      strict (boolean, default: true) Discard templates with E-values below threshold.
      aligned (boolean, default: false) Skip alignment if already aligned.
      legacy (boolean, default: false) Use pre-aligned, single-chain template for the full query.
      reference (string, default: null) Reference to be used in order to insert it into the query sequence.
      modifications (List) Chain-level edits before/after alignment. Each modification can include:
         - chain (string, mandatory) chain ID or All.
           position (integer, default: null) Insertion position in query (if single chain).
           maintain_residues (list of integers, default: null) Selected residues will be kept, and the rest will be deleted.
           delete_residues (list of integers, default: null) Selected residues will be deleted, the rest will be kept.
           when (string, default: after_alignment) before_alignment or after_alignment.
           mutations (List) Modifications in the residues:
              - numbering_residues (list of integers, mandatory) Residue positions where the mutations will be applied.
                mutate_with (string, mandatory) The amino acid to mutate to, specified as a three‑letter code or as a FASTA file path.
```

### 5. Add features
Merge or slice existing features.pkl files from other AlphaFold2 runs into your run.
```
features:
    - path (string, mandatory) Path to an existing features.pkl file.
      keep_msa (integer, default: -1) -1 = all sequences; otherwise top X by coverage.
      keep_templates (integer, default: -1) -1 = all templates; otherwise top X by coverage.
      msa_mask (range, default: null) Remove this residue range from the MSA.
      sequence (string, default: null) FASTA file to replace all template sequences.
      numbering_query (list of integers, default: null) Insertion positions in the query sequence.
      numbering_features (list of ranges, default: null) Map feature blocks into the positions given by numbering_query.
      positions (range, default: null) Inserts the features.pkl into the query sequence. The position refers to the sequence index, whereas in numbering_query and numbering_features, it refers to the residue positions in the entire query sequence.
      mutations (map) Map three-letter amino acid codes to residue indices. Example:
        - 'ALA': 10, 20
```

### 6. Append library
Append existing FASTA/PDB files from a library into your run.
```
append_library:
    - path: (string, mandatory) Path to a directory, PDB, or FASTA file.
      add_to_msa (boolean, default: true) Append sequences to the MSA.
      add_to_templates (boolean, default: false) Append PDBs to the templates.
      numbering_query (list of integers, default: null) Insertion positions in the query.
      numbering_library (list of ranges, default: null) Residue range from the library entry to insert.
```

### 7. Configuration file example
```
mode: guided
output_dir: /path/to/output
af2_dbs_path: /path/to/af2_dbs
run_af2: True
experimental_pdbs: /path/to/references/experimental.pdb

sequences:
- fasta_path: /path/to/data/seq1.fasta
  num_of_copies: 1
- fasta_path: /path/to/data/seq2.fasta
  num_of_copies: 1
- fasta_path: /path/to/data/seq3.fasta
  num_of_copies: 1
- fasta_path: /path/to/data/seq4.fasta
  num_of_copies: 1

features:
- path: /path/to/features1.pkl
  keep_msa: 30
  keep_templates: 0
  numbering_query: 1

- path: /path/to/features2.pkl
  keep_msa: 30
  keep_templates: 0
  msa_mask: 276-477, 652-857
  numbering_query: 1

- path: /path/to/features3.pkl
  keep_msa: 30
  keep_templates: 0
  msa_mask: 8-250
  numbering_query: 4

templates:
- pdb: /path/to/templates/template.pdb
  add_to_msa: true
  add_to_templates: True
  generate_multimer: False
  aligned: true
  modifications:
  - chain: A
    position: 1
    mutations:
    - numbering_residues: 276-477
      mutate_with: /path/to/data/seq1.fasta

```

-------------------

## Output information

All information is located in the output_dir directory, which is specified as an input parameter in the configuration file. Inside output_dir, you will find the following folders and files:
- output.html: Contains the results in HTML format, including all plots, run statistics, and prediction analyses.
- output.log: The log file with detailed information from the execution.
- plots/: All plots generated by the output analysis.
- frobenius/: Plots generated by ALEPH.
- interfaces/: Results of the interface analysis performed by PISA.
- clustering/: (If clustering is enabled) Contains the results related to clustering jobs.
- input/: All input files used in the run.
- run/: Stores runtime information and outputs (see below for details).
- templates/: Templates extracted from the features.pkl, split by chains.
- rankeds: Ranked models generated by AlphaFold2, split by chains.

Inside the run/ directory, you will find:
- results: Results of the AlphaFold2 run (see below for details).
- Templates folder: Subfolders named after each template, containing the databases generated to align each template.
- Sequences folder: Subfolders named after each sequence, containing alignments of the templates with the corresponding sequence.

Inside the run/results/ directory, you will find:
- tmp/: Contains intermediate files generated by external programs (e.g., Aleph).
- ccanalysis/ and ccanalysis_ranked/: PDB files used for the cc_analysis run.
- msas/: Information generated by AlphaFold2. It contains the extracted sequences and the template alignments.
- templates_nonsplit/: Templates extracted from features.pkl, not split by chains.
- rankeds_split/: Ranked models generated by AlphaFold2, split by chains.
- rankeds/: Ranked models generated by AlphaFold2, not split by chains.
