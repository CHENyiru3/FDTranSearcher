# FDTranSearcher

FDTranSearcher is a bioinformatics tool for detecting functional DNA transposons in plants genomic sequences. It is a part of ZJE BMI3 ICA.  

### Supported Operating Systems
- Linux ✅
- macOS ✅
- Windows (with limitations) ⚠️

> **Important**: Windows users must manually download and install BLAST+ from the [NCBI website](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.5.0/), as conda does not provide BLAST packages for Windows.

## Table of Contents

- [Background](#background)
- [Installation](#installation)
- [Modes](#modes)
- [Parameters](#parameters)
- [Usage](#usage)
- [Citation](#citation)

## Background

Transposon elements (TEs) are elements that can move around the genome and can be broadly classified into two main categories:  copy-and-paste retrotransposons and cut-and-paste DNA transposons. Studies have shown that DNA transposons play important roles in plant genomes for gene regulation and evolution.

However, the search and verification for  **functional DNA transposons with transposase activity**  in genomes is still a major challenge. To facilitate further research, we plan to develop tools for finding functional DNA transposons in the maize genome, and try to expand to other plant genomes. In order to address different application scenarios, we designed three modules: 
1. **Reference-based Module**  
   A pipeline leveraging existing databases to identify known functional DNA transposons.
2. **De Novo Module**  
   An algorithm that uses structural features to predict unknown functional DNA transposons.
3. **Ensemble Tool mode**  
   A comprehensive module that integrates the outputs of both pipelines for robust functional DNA transposon annotation.


<img width="1261" alt="b56a72c1fb789708e50bff028954b4c" src="https://github.com/user-attachments/assets/eb9a2094-3bd2-4590-a390-5c90f50d5f1e">



## Installation

### Using Conda (Recommended)

After downloading the source code, we recommend creating a virtual environment that meets the requirements of this tool. Follow the commands below:

```bash
git clone https://github.com/CHENyiru3/FDTranSearcher
cd FDTranSearcher
# Create and activate the Conda environment:
conda env create -f environment.yml
conda activate FDTranSearcher
```


### Manual Installation
If you prefer a manual setup, ensure you meet the following system and package requirements:

#### Environment Requirements
- **Operating System**: Ubuntu 20.04.2 LTS (Recommended)
  
- **Python**: Version >= 3.9.18

----conda packages----
- **BLAST**: Version >= 2.5.0  
- biopython >= 1.78
- psutil=6.1.0
- tqdm=4.66.2
- pandas
- numpy

### Additional Tool Installation

You also need to manually install the open-source tool miniprot. Use the following commands to clone and build it:
```bash
git clone https://github.com/lh3/miniprot.git
cd miniprot && make
```
If you encounter any issues, refer to the [official miniprot GitHub page](https://github.com/lh3/miniprot) for additional instructions or troubleshooting.


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

## Usage

### Reference-based Module
To use the reference-based module, make sure you have two reference files:

- Consensus transposon DNA sequences
- Transposase protein sequences

The consensus datasets for [DNA transposons](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-018-0144-1) and [transposases](https://www.nature.com/articles/nature22971) for the five DNA TE families in maize, used for testing, will be automatically downloaded as built-in references. 

The built-in transposase reference can be used in most plants; for easily downloading vaild consensus transposon DNA sequences of other species, we recommend you search in [PlantRep](http://www.plantrep.cn/).

Here is an example usage:

```bash
usage: reference_based_main.py [-h] -g GENOME_FASTA_PATH
                               -ref_dna TARGET_SEQUENCE_FILE
                               -miniprot MINIPROT_PATH
                               [-ref_tp TPASE]
                               [-or OUTPUT_REPORT_FILE]
                               [-ot OUTPUT_CSV_FILE]
                               [-gff OUTPUT_GFF3_FILE]
                               [--blast_evalue_threshold BLAST_EVALUE_THRESHOLD]
                               [--blast_identity_threshold BLAST_IDENTITY_THRESHOLD]
                               [--blast_alignment_length_threshold BLAST_ALIGNMENT_LENGTH_THRESHOLD]
                               [--blast_max_target_seqs BLAST_MAX_TARGET_SEQS]
                               [--blast_max_hsps BLAST_MAX_HSPS]
                               [--threads THREADS]
                               [--pattern_size PATTERN_SIZE]
                               [--gap_size GAP_SIZE]
                               [--tir_size TIR_SIZE]
                               [--mismatch_allowed MISMATCH_ALLOWED]
                               [--mini_size MINI_SIZE]
                               [--max_size MAX_SIZE]
                               [--report_mode {all,human_readable,gff3,table}]
                               [--extension EXTENSION]
                               [--miniprot_threads MINIPROT_THREADS]
                               [--miniprot_c MINIPROT_C]
                               [--miniprot_m MINIPROT_M]
                               [--miniprot_p MINIPROT_P]
                               [--miniprot_N MINIPROT_N]
                               [--miniprot_O MINIPROT_O]
                               [--miniprot_J MINIPROT_J]
                               [--miniprot_F MINIPROT_F]
                               [--miniprot_K MINIPROT_K]
                               [--miniprot_outn MINIPROT_OUTN]
                               [--miniprot_outs MINIPROT_OUTS]
                               [--miniprot_outc MINIPROT_OUTC]

Reference-based Transposable Element Analysis

optional arguments:
  -h, --help            show this help message and exit
  -g GENOME_FASTA_PATH, --genome_fasta_path GENOME_FASTA_PATH
                        Path to the genome FASTA file (REQUIRED)
  -ref_dna TARGET_SEQUENCE_FILE, --target_sequence_file TARGET_SEQUENCE_FILE
                        Path to the target sequence file (reference TEs) (REQUIRED)
  -miniprot MINIPROT_PATH, --miniprot_path MINIPROT_PATH
                        Path to the miniprot executable (REQUIRED)
  -ref_tp TPASE, --Tpase TPASE
                        Path to the Tpase sequences file
  -or OUTPUT_REPORT_FILE, --output_report_file OUTPUT_REPORT_FILE
                        Path to the output report file
  -ot OUTPUT_CSV_FILE, --output_csv_file OUTPUT_CSV_FILE
                        Path to the output CSV file
  -gff OUTPUT_GFF3_FILE, --output_gff3_file OUTPUT_GFF3_FILE
                        Path to the output GFF3 file
  --blast_evalue_threshold BLAST_EVALUE_THRESHOLD
                        E-value threshold for BLASTN (default: 1e-5)
  --blast_identity_threshold BLAST_IDENTITY_THRESHOLD
                        Identity threshold for BLASTN (default: 60.0)
  --blast_alignment_length_threshold BLAST_ALIGNMENT_LENGTH_THRESHOLD
                        Alignment length threshold for BLASTN (default: 30)
  --blast_max_target_seqs BLAST_MAX_TARGET_SEQS
                        Maximum number of target sequences for BLASTN (default: 1000)
  --blast_max_hsps BLAST_MAX_HSPS
                        Maximum number of HSPs for BLASTN (default: 10)
  --threads THREADS     Number of threads to use (default: 16)
  --pattern_size PATTERN_SIZE
                        Size of the TSD pattern (default: 8)
  --gap_size GAP_SIZE   Maximum allowed gap between TSDs (default: 1500)
  --tir_size TIR_SIZE   Size of the TIR pattern (default: 5)
  --mismatch_allowed MISMATCH_ALLOWED
                        Number of allowed mismatches in TIRs (default: 2)
  --mini_size MINI_SIZE
                        Minimum size of the structure (default: 2000)
  --max_size MAX_SIZE   Maximum size of the structure (default: 15000)
  --report_mode {all,human_readable,gff3,table}
                        Report mode:
                                'all': Generate detailed report, CSV, and GFF3.
                                'human_readable': Generate only a summary CSV report.
                                'gff3': Generate only a GFF3 report.
                                'table': Generate only a detailed report.
                                (default: all)
  --extension EXTENSION
                        Extension size for TE sequence extraction (default: 3000)
  --miniprot_threads MINIPROT_THREADS
                        Number of threads for miniprot (default: 16)
  --miniprot_c MINIPROT_C
                        Minimum spanning chain score for miniprot (default: 50000)
  --miniprot_m MINIPROT_M
                        Minimal alignment score for miniprot (default: 10)
  --miniprot_p MINIPROT_P
                        Minimal alignment score fraction for miniprot (default: 0.2)
  --miniprot_N MINIPROT_N
                        Maximum number of intron hints for miniprot (default: 200)
  --miniprot_O MINIPROT_O
                        Maximum number of introns for miniprot (default: 3)
  --miniprot_J MINIPROT_J
                        Maximum intron length for miniprot (default: 8)
  --miniprot_F MINIPROT_F
                        Maximum number of exons for miniprot (default: 8)
  --miniprot_K MINIPROT_K
                        Maximum memory to use for miniprot (default: 5M)
  --miniprot_outn MINIPROT_OUTN
                        Maximum number of alignments to output for miniprot (default: 5000)
  --miniprot_outs MINIPROT_OUTS
                        Output score threshold for miniprot (default: 0.5)
  --miniprot_outc MINIPROT_OUTC
                        Output coverage threshold for miniprot (default: 0.03)
```

### De-novo Module

You don't need any reference files for the de-novo module. Here is an example usage:

```bash
usage: transposon_analyzer.py [-h] --input INPUT
                              --output-prefix OUTPUT_PREFIX
                              [--min-tsd-pattern-size MIN_TSD_PATTERN_SIZE]
                              [--max-tsd-pattern-size MAX_TSD_PATTERN_SIZE]
                              [--gap-size GAP_SIZE]
                              [--min-tir-size MIN_TIR_SIZE]
                              [--max-tir-size MAX_TIR_SIZE]
                              [--max-tir-mismatch MAX_TIR_MISMATCH]
                              [--conserve-site-range CONSERVE_SITE_RANGE]
                              [--subterminal-threshold SUBTERMINAL_THRESHOLD]
                              [--subterminal-length SUBTERMINAL_LENGTH]
                              [--motif MOTIF]
                              [--chunk-size CHUNK_SIZE]
                              [--num-processes NUM_PROCESSES]
                              [--conserve-sites CONSERVE_SITES]

Analyze genome for transposons

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        Input genome file in FASTA format
  --output-prefix OUTPUT_PREFIX, -o OUTPUT_PREFIX
                        Prefix for output files
  --min-tsd-pattern-size MIN_TSD_PATTERN_SIZE
                        Minimum TSD pattern size for search
                        (default: 5)
  --max-tsd-pattern-size MAX_TSD_PATTERN_SIZE
                        Maximum TSD pattern size for search
                        (default: 10)
  --gap-size GAP_SIZE   Expected transposable element search
                        gap (default: 2500 means length
                        between 2500-5000bp)
  --min-tir-size MIN_TIR_SIZE
                        Minimum TIR size (default: 5)
  --max-tir-size MAX_TIR_SIZE
                        Maximum TIR size (default: 15)
  --max-tir-mismatch MAX_TIR_MISMATCH
                        Maximum allowed mismatches in TIR
                        (default: 1)
  --conserve-site-range CONSERVE_SITE_RANGE
                        Search range around conserved reserve
                        sites positions (default: 50)
  --subterminal-threshold SUBTERMINAL_THRESHOLD
                        Threshold for subterminal motif
                        enrichment percentage (default: 5.0)
  --subterminal-length SUBTERMINAL_LENGTH
                        Length of subterminal region to
                        search for motifs (default: 250)
  --motif MOTIF         User-specified subterminal motif
                        (default: AAAGGG)
  --chunk-size CHUNK_SIZE
                        Size of sequence chunks for parallel
                        processing (default: 25000)
  --num-processes NUM_PROCESSES
                        Number of parallel processes
                        (default: 10)
  --conserve-sites CONSERVE_SITES
                        Comma-separated conserve sites in
                        format [AminoAcid][Position] (e.g.,
                        D301,E719,Q500). Default:
                        D301,D376,E719
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

