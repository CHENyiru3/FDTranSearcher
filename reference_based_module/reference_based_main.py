import os
import subprocess
from dataclasses import dataclass
from typing import List, Dict, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
from io import StringIO
import tempfile
import datetime
import csv
import math
import argparse
import yaml
import os
from tqdm import tqdm
import multiprocessing
from multiprocessing import Pool

from miniprot_filter import run_miniprot_analysis
from GFFEntry_class import GFFEntry, parse_gff_attributes
from architecture_identification import analyze_structures
from result_reporter import generate_csv_report, generate_gff3_report, generate_report
from blast_search import run_blastn_search, extract_transposable_elements
def run_analysis(
    genome_fasta_path: str,
    target_sequence_file: str,
    miniprot_path: str,
    Tpase: str,
    output_report_file: str = None,
    output_csv_file: str = None,
    output_gff3_file: str = None,
    evalue_threshold: float = 1e-5,
    identity_threshold: float = 60.0,
    alignment_length_threshold: int = 30,
    max_target_seqs: int = 1000,
    max_hsps: int = 10,
    threads: int = 16,
    # additional optional parameters
    extension = 3000, mini_size = 2000, max_size = 15000, pattern_size = 8, gap_size = 1500,
    tir_size = 5, mismatch_allowed = 2, report_mode = 'all'
):
    try:
        # First step: use blastn to search for transposable elements
        print("Step 1: Running BLASTN search...")
        gff_entries = run_blastn_search(
            target_sequence_file,
            genome_fasta_path,
            evalue_threshold,
            identity_threshold,
            alignment_length_threshold,
            max_target_seqs,
            max_hsps,
            threads,
        )

        print(f"Found {len(gff_entries)} transposable elements.")

        # Step 2: extract the sequences of transposable elements
        print("Step 2: Extracting transposable element sequences...")
        transposable_elements = extract_transposable_elements(
            gff_entries,
            genome_fasta_path,
            extension
        )
        print(f"Extracted {len(transposable_elements)} transposable element sequences.")

        # Step 3: run miniprot analysis to identify protein regions
        print("Step 3: Running miniprot analysis...")
        miniprot_results = run_miniprot_analysis(
            miniprot_path,
            Tpase,
            transposable_elements,
            threads
        )
        print(f"miniprot analysis found {len(miniprot_results)} protein regions.")

        analysis_results = []
        for te_record in transposable_elements:
            # Directly get position information from annotations
            blast_match = te_record.annotations['blast_match']

            # Use a more flexible matching method
            protein_regions = []
            for miniprot_entry in miniprot_results:
                # Modify matching logic
                if (miniprot_entry.seqid == te_record.id or
                    te_record.id in miniprot_entry.seqid or
                    miniprot_entry.seqid.split('_')[-1] in te_record.id):
                    protein_regions.append((miniprot_entry.start, miniprot_entry.end))

            # Only keep transposons that have found transposase proteins
            if protein_regions:
                try:
                    structures = analyze_structures(
                        str(te_record.seq),
                        protein_regions
                    )

                    # Check if structures are empty
                    if not structures:
                        print(f"Filter: No structures found for TE {te_record.id}")
                        continue

                    # Calculate the precise genomic position of the transposon structure
                    te_structure_start = min(
                        structures[0]['TSD1'][0],
                        structures[0]['TIR1'][0],
                        structures[0]['TIR2'][0],
                        structures[0]['TSD2'][0]
                    )
                    te_structure_end = max(
                        structures[0]['TSD1'][1],
                        structures[0]['TIR1'][1],
                        structures[0]['TIR2'][1],
                        structures[0]['TSD2'][1]
                    )
                    architecture_score = structures[0]['scores']['total_score']
                    # Calculate the precise position in the reference genome
                    blast_position = f"{blast_match['chromosome']}:{blast_match['start']}-{blast_match['end']}"
                    element_position = f"{blast_match['chromosome']}:{blast_match['start'] + te_structure_start}-{blast_match['start'] + te_structure_end}"

                    analysis_results.append({
                        'te_id': te_record.id,
                        'blast_position': blast_position,
                        'element_position': element_position,
                        'blast_match': blast_match,
                        'te_structure_position': {
                            'start': te_structure_start,
                            'end': te_structure_end
                        },
                        'protein_regions': protein_regions,
                        'structures': structures,
                        'source': 'FDT_reference_pipeline',
                        'architecture_score': architecture_score,
                        'type': 'functional_transposable_element',
                        'strand': blast_match['strand'],  # Use BLAST's strand information
                        'score': 0.0,
                        'miniprot_strand': miniprot_entry.strand,  # Optionally record miniprot's strand in attributes
                        'phase': '.',
                    })

                except Exception as e:
                    print(f"Error processing TE {te_record.id}: {str(e)}")
                    continue

        print(f"Successfully analyzed {len(analysis_results)} transposable elements.")

        # Step 5: generate reports
        print("Step 5: Generating reports...")
        if report_mode == 'all':
            generate_report(analysis_results, output_report_file)
            print(f"Detailed report saved to {output_report_file}.")

            generate_csv_report(analysis_results, output_csv_file)
            print(f"Detailed CSV report saved to {output_csv_file}.")

            generate_gff3_report(analysis_results, output_gff3_file)
            print(f"GFF3 report saved to {output_gff3_file}.")
        elif report_mode == 'human_readable':
            generate_csv_report(analysis_results, output_csv_file)
            print(f"Summary CSV report saved to {output_csv_file}.")
        elif report_mode == 'gff3':
            generate_gff3_report(analysis_results, output_gff3_file)
            print(f"GFF3 report saved to {output_gff3_file}.")
        elif report_mode == "table":
            generate_report(analysis_results, output_report_file)
            print(f"Detailed report saved to {output_report_file}.")
    finally:
        for file in [f for f in os.listdir() if f.endswith('.tmp')]:
            os.remove(file)

def load_config(config_file):
    """Load configuration from a YAML file."""
    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)
    return config

def merge_args_with_config(args, config):
    """Override config dictionary with command-line arguments if provided."""
    for key, value in vars(args).items():
        if value is not None:  # Command-line argument takes precedence if provided
            config[key] = value
    return config
import argparse

def main():
    parser = argparse.ArgumentParser(description="Reference-based Transposable Element Analysis")

    # Add all arguments directly in the main function
    parser.add_argument("-g", "--genome_fasta_path", required=True, help="Path to the genome FASTA file")
    parser.add_argument("-ref", "--target_sequence_file", required=True, help="Path to the target sequence file")
    parser.add_argument("-miniprot", "--miniprot_path", required=True, help="Path to the miniprot executable")
    parser.add_argument("-tp", "--Tpase", help="Path to the hAT sequences file")
    parser.add_argument("-or", "--output_report_file", help="Path to the output report file")
    parser.add_argument("-ot", "--output_csv_file",  help="Path to the output CSV file")
    parser.add_argument("-gff", "--output_gff3_file",  help="Path to the output GFF3 file")
    parser.add_argument("--evalue_threshold", type=float, default=1e-5, help="E-value threshold for BLASTN")
    parser.add_argument("--identity_threshold", type=float, default=60.0, help="Identity threshold for BLASTN")
    parser.add_argument("--alignment_length_threshold", type=int, default=30, help="Alignment length threshold for BLASTN")
    parser.add_argument("--max_target_seqs", type=int, default=1000, help="Maximum number of target sequences for BLASTN")
    parser.add_argument("--max_hsps", type=int, default=10, help="Maximum number of HSPs for BLASTN")
    parser.add_argument("--threads", type=int, default=16, help="Number of threads to use")
    parser.add_argument("--pattern_size", type=int, default=8, help="Size of the TSD pattern")
    parser.add_argument("--gap_size", type=int, default=1500, help="Maximum allowed gap between TSDs")
    parser.add_argument("--tir_size", type=int, default=5, help="Size of the TIR pattern")
    parser.add_argument("--mismatch_allowed", type=int, default=2, help="Number of allowed mismatches in TIRs")
    parser.add_argument("--mini_size", type=int, default=2000, help="Minimum size of the structure")
    parser.add_argument("--max_size", type=int, default=15000, help="Maximum size of the structure")
    parser.add_argument("--report_mode", type=str, default="all", help="Report mode (all, human_readable, gff3, table)")

    args = parser.parse_args()

    # Now you can use args to access all the parsed arguments
    run_analysis(
        genome_fasta_path=args.genome_fasta_path,
        target_sequence_file=args.target_sequence_file,
        miniprot_path=args.miniprot_path,
        Tpase=args.Tpase,
        output_report_file=args.output_report_file,
        output_csv_file=args.output_csv_file,
        output_gff3_file=args.output_gff3_file,
        evalue_threshold=args.evalue_threshold,
        identity_threshold=args.identity_threshold,
        alignment_length_threshold=args.alignment_length_threshold,
        max_target_seqs=args.max_target_seqs,
        max_hsps=args.max_hsps,
        threads=args.threads,
        pattern_size=args.pattern_size,
        gap_size=args.gap_size,
        tir_size=args.tir_size,
        mismatch_allowed=args.mismatch_allowed,
        mini_size=args.mini_size,
        max_size=args.max_size,
        report_mode=args.report_mode
        
    )

if __name__ == "__main__":
    main()