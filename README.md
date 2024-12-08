# FDTranSearcher

FDTranSearcher is a tool designed for detecting functional DNA transposons in genomes. Currently, it supports Linux and macOS systems. For Windows users, BLAST+ must be manually downloaded because conda does not provide BLAST for Windows.

## Table of Contents

- [Background](#background)
- [Installation](#installation)
- [Modes](#modes)
- [Parameters](#parameters)
- [Usage](#usage)
- [Citation](#citation)

## Background

Transposon elements (TEs) are elements that can move around the genome and can be broadly classified into two main categories:  copy-and-paste retrotransposons and cut-and-paste DNA transposons. Studies have shown that DNA transposons play important roles in plant genomes for gene regulation and evolution.

However, the search and verification for functional DNA transposons with transposase in genomes is still a major challenge. To facilitate further research, we plan to develop tools for finding functional DNA transposons in the maize genome, and try to expand to other plant genomes. In order to address different application scenarios, we intend to design three modules: 
- pipeline construction based on existing databases to find known functional DNA transposons 
- de novo algorithms based on structural features to find unknown functional DNA transposons
- build a ensemble tool for comprehensive functional DNA transposon annotation

## Installation

### Using Conda (Recommended)

After downloading the source code, we recommend creating a virtual environment that meets the requirements of this tool. Run the following commands in your terminal:

```bash
cd FDTranSearcher
conda env create -f environment.yml
conda activate FDTranSearcher
```

### Manual Installation

Requirements:

- Python >=3.9.18
- BLAST >= 2.5.0
- Ubuntu  20.04.2 LTS(Recommended)
- Required Python packages:
  - biopython >= 1.78
  - PyYAML >=6.0.2
  - pandas
 
### Additional Tool Installation

Additionally, you need to manually install the required open-source tool [miniprot](https://github.com/lh3/miniprot). If you have problems installing miniprot, you can also [click](https://github.com/lh3/miniprot) to jump to miniprot.:

```bash
git clone https://github.com/lh3/miniprot.git
cd miniprot && make
```
## Modes

FDTranSearcher has three modes:

- Reference-based (r)
- De-novo (d)
- Ensemble (a)

The usage of these modes is integrated into a bash script controller (FDTranSearcher.sh). Start by checking the parameter usage guidance for each mode:

```bash
bash FDTranSearcher.sh r -h
bash FDTranSearcher.sh d -h
```
## Parameters

The tool accepts the following parameters:

- `-h, --help`
  - Show this help message and exit

- `-g, --genome_fasta_path GENOME_FASTA_PATH`
  - Path to the genome FASTA file that you want to search

- `-ref_dna, --reference_transposon_file REFERENCE_TRANSPOSON_FILE`
  - Path to the reference transposon sequence file

- `-miniprot, --miniprot_path MINIPROT_PATH`
  - Path to the miniprot executable

- `-ref_tp, --Tpase TPASE`
  - Path to the transposase amino acid sequences file

- `-or, --output_report_file OUTPUT_REPORT_FILE`
  - Path to the output report file

- `-ot, --output_csv_file OUTPUT_CSV_FILE`
  - Path to the output CSV file

- `-gff, --output_gff3_file OUTPUT_GFF3_FILE`
  - Path to the output GFF3 file

- `--evalue_threshold EVALUE_THRESHOLD` (default: 1e-5)
  - E-value threshold for BLASTN

- `--identity_threshold IDENTITY_THRESHOLD` (default: 90)
  - Identity threshold for BLASTN

- `--alignment_length_threshold ALIGNMENT_LENGTH_THRESHOLD` (default: 100)
  - Alignment length threshold for BLASTN

- `--max_target_seqs MAX_TARGET_SEQS` (default: 500)
  - Maximum number of target sequences for BLASTN

- `--max_hsps MAX_HSPS` (default: 10)
  - Maximum number of HSPs for BLASTN

- `--threads THREADS` (default: all available CPU cores)
  - Number of threads to use for parallel operations

- `--pattern_size PATTERN_SIZE` (default: 20)
  - Size of the TSD (Target Site Duplication) pattern

- `--gap_size GAP_SIZE` (default: 200)
  - Maximum allowed gap between TSDs

- `--tir_size TIR_SIZE` (default: 100)
  - Size of the TIR (Terminal Inverted Repeat) pattern

- `--mismatch_allowed MISMATCH_ALLOWED` (default: 2)
  - Number of allowed mismatches in TIRs

- `--mini_size MINI_SIZE` (default: 500)
  - Minimum size of the structure

- `--max_size MAX_SIZE` (default: 5000)
  - Maximum size of the structure

- `--report_mode REPORT_MODE` (default: "all")
  - Report mode (options: "all", "human_readable", "gff3", "table")


## Usage

### Reference-based Module
To use the reference-based module, make sure you have two reference files:

- Consensus transposon DNA sequences
- Transposase protein sequences

The consensus datasets for [DNA transposons](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-018-0144-1) and [transposases](https://www.nature.com/articles/nature22971) for the five DNA TE families in maize, used for testing, will be automatically downloaded as built-in references. 

Here is an example usage:

```bash
bash FDTranSearcher.sh r [-h] -g GENOME_FASTA_PATH -ref_dna REFERENCE_TRANSPOSON_FILE -miniprot MINIPROT_PATH [-ref_tp TPASE] [-or OUTPUT_REPORT_FILE] [-ot OUTPUT_CSV_FILE]
                        [-gff OUTPUT_GFF3_FILE] [--evalue_threshold EVALUE_THRESHOLD] [--identity_threshold IDENTITY_THRESHOLD]
                        [--alignment_length_threshold ALIGNMENT_LENGTH_THRESHOLD] [--max_target_seqs MAX_TARGET_SEQS] [--max_hsps MAX_HSPS] [--threads THREADS]
                        [--pattern_size PATTERN_SIZE] [--gap_size GAP_SIZE] [--tir_size TIR_SIZE] [--mismatch_allowed MISMATCH_ALLOWED] [--mini_size MINI_SIZE]
                        [--max_size MAX_SIZE] [--report_mode REPORT_MODE]
```

### De-novo Module

You don't need any reference files for the de-novo module. Here is an example usage:

```bash
bash FDTranSearcher.sh d --input INPUT --output-prefix OUTPUT_PREFIX
```

### Ensemble Module

You must enter the parameters required for the reference-based module and the de-novo module:

```bash
bash FDTranSearcher.sh a
```

Enter the appropriate parameters when you see the following:

```bash
Running both reference-based and de novo modules...
Please enter parameters for the reference-based module, e.g., -g file -o output:

Please enter parameters for the de novo module, e.g., -i file -o output:

Please enter the path to save the final merged results, e.g., /path/to/output.gff:
```

## Citation

[miniprot](https://github.com/lh3/miniprot) for finding transposase sequences

The consensus datasets for [DNA transposons](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-018-0144-1) and [transposases](https://www.nature.com/articles/nature22971) for the five DNA TE families in maize.

