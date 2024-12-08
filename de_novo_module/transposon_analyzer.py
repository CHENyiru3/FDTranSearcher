from Bio import SeqIO
from Bio.Seq import Seq
import csv
import pandas as pd
from typing import List, Tuple, Dict, Optional, Generator, Set
from datetime import datetime
from multiprocessing import Pool
import re
from tqdm import tqdm
import psutil
import os
from collections import defaultdict
import numpy as np
import logging

# Import custom modules
from element_search import process_chunk
from sequence_tools import extract_sequence, parse_reserve_sites
from structure_verification import StructureVerification
from result_processor import process_results

# Configure logging
logging.basicConfig(
    filename='transposon_analyzer.log',
    filemode='a',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

class TransposonAnalyzer:
    """
    Main class for transposon analysis pipeline.
    Handles the complete workflow of identifying and analyzing transposons in genomic sequences.
    """
    def __init__(self, **params):
        """
        Initialize TransposonAnalyzer with analysis parameters.
        
        Args:
            **params: Dictionary containing analysis parameters
        """
        # Set default parameters
        self.params = {
            'min_tsd_pattern_size': 5,
            'max_tsd_pattern_size': 10,
            'gap_size': 2500,
            'min_tir_size': 5,
            'max_tir_size': 15,
            'max_tir_mismatch': 1,
            'conserve_site_range': 50,
            'chunk_size': 25000,
            'num_processes': 10,
            'subterminal_threshold': 5.0,
            'motif': 'AAAGGG',
            'subterminal_length': 250,
            'min_cds_distance': params['min_cds_distance']  # Default value
        }
        # Update with user-provided parameters
        self.params.update(params)
        
        # Initialize reserve_sites
        self.reserve_sites = []
        
        # Set min_cds_distance from params
        self.min_cds_distance = self.params['min_cds_distance']
        
        # Validate parameters
        self._validate_params()
    
    def _validate_params(self):
        """
        Validate parameter values.
        
        Raises:
            ValueError: If any parameter values are invalid
        """
        if self.params['conserve_site_range'] < 0:
            raise ValueError("Conserve site search range must be non-negative")
        if self.params['gap_size'] < 0: 
            raise ValueError("Gap size must be non-negative")
        if self.params['min_tsd_pattern_size'] > self.params['max_tsd_pattern_size']:
            raise ValueError("min_tsd_pattern_size cannot be greater than max_tsd_pattern_size")
        if self.params['min_tir_size'] > self.params['max_tir_size']:
            raise ValueError("min_tir_size cannot be greater than max_tir_size")
        if self.params['subterminal_threshold'] < 0: 
            raise ValueError("Subterminal threshold must be non-negative")
        
        # Validate motif
        motif = self.params['motif']
        if not re.fullmatch(r'[ATCG]+', motif):
            raise ValueError(f"Invalid motif '{motif}'. Motif must consist of A, T, C, G only.")
        
        # Validate subterminal_length
        max_subterminal = (2500 + self.params['gap_size']) / 2
        if not isinstance(self.params['subterminal_length'], int):
            raise ValueError("subterminal_length must be an integer")
        if self.params['subterminal_length'] >= max_subterminal:
            raise ValueError(f"subterminal_length must be less than {max_subterminal}")

    def run_analysis(self, input_file: str, output_prefix: str, reserve_sites: List[str]):
        """
        Run complete transposon analysis pipeline.
    
        Args:
            input_file: Path to input FASTA file
            output_prefix: Prefix for output files
            reserve_sites: List of user-specified reserve sites (e.g., ['D301', 'E719'])
        """
        # Process reserve sites
        self.reserve_sites = parse_reserve_sites(reserve_sites)
    
        # Phase 1: Search transposons
        initial_output = f"{output_prefix}_raw.csv"
        self._search_transposons(input_file, initial_output)
    
        # Phase 2: Process results
        final_csv = f"{output_prefix}_processed.csv"
        final_gff = f"{output_prefix}.gff3"
        process_results(initial_output, final_csv, final_gff, self.reserve_sites, self.params)
    
        print(f"Analysis complete. Results written to:")
        print(f"  Raw CSV: {initial_output}")    
        print(f"  Processed CSV: {final_csv}")   
        print(f"  GFF3: {final_gff}")

    def _search_transposons(self, input_file: str, output_file: str):
        """
        Execute transposon search and save raw results.

        Args:
            input_file: Path to input FASTA file
            output_file: Path for output raw CSV file
        """
        start_time = datetime.now()
        logging.info(f"Analysis started at: {start_time}")
        print(f"Analysis started at: {start_time}")
        
        sequences = extract_sequence(input_file)
        chunk_size = self.params['chunk_size']
        num_processes = self.params['num_processes']
        
        max_memory_usage = 0.0  # Initialize max memory usage
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            # Generate reserve site column names
            reserve_site_headers = []
            for i in range(1, len(self.reserve_sites)+1):
                reserve_site_headers.extend([
                    f'conserve_site_{i}_distance'
                ])
            
            writer.writerow([
                'Sequence_ID',
                'TSD1_Start', 'TSD1_End',
                'TIR1_Start', 'TIR1_End',
                'TIR2_Start', 'TIR2_End',
                'TSD2_Start', 'TSD2_End',
                'TSD_Sequence',
                'TIR1_Sequence',
                'TIR2_Sequence',
                f'Left_Subterminal_{self.params["motif"]}_%',
                f'Right_Subterminal_{self.params["motif"]}_%',
                'Strand',
                'CDS_Start',
                'CDS_End',
                'CDS_Length',
                'Protein_Length',
                'Distance_to_Left_TIR',
                'Distance_to_Right_TIR',
                'Pattern_Match',
                *reserve_site_headers
            ])
            
            total_structures = 0
            
            for seq_id, seq in sequences:
                logging.info(f"Processing sequence: {seq_id}")
                print(f"\nProcessing sequence: {seq_id}")
                chunks = []
                for i in range(0, len(seq), chunk_size):
                    chunk = seq[i:min(i+chunk_size, len(seq))]
                    chunk_id = f"{seq_id}_chunk_{i//chunk_size}"
                    chunks.append((chunk_id, chunk, seq_id, i, self.params, self.reserve_sites))
                
                with Pool(processes=num_processes) as pool:
                    for result in tqdm(pool.imap_unordered(process_chunk, chunks), 
                                     total=len(chunks), 
                                     desc=f"Processing {seq_id}"):
                        chunk_id, matches, memory_usage = result
                        if memory_usage > max_memory_usage:
                            max_memory_usage = memory_usage
                        if matches:
                            for match in matches:
                                row = [
                                    match['seq_id'],
                                    match['tsd1_start'],
                                    match['tsd1_end'],
                                    match['tir1_start'],
                                    match['tir1_end'],
                                    match['tir2_start'],
                                    match['tir2_end'],
                                    match['tsd2_start'],
                                    match['tsd2_end'],
                                    match['tsd_seq'],
                                    match['tir1_seq'],
                                    match['tir2_seq'],
                                    f"{match['left_percentage']:.2f}",
                                    f"{match['right_percentage']:.2f}",
                                    match['strand'],
                                    match['cds_start'],
                                    match['cds_end'],
                                    match['cds_length'],
                                    match['protein_length'],
                                    match['distance_to_left_tir'],
                                    match['distance_to_right_tir'],
                                    match['pattern_match'],
                                ]
                                # Add reserve site distances
                                for i in range(1, len(self.reserve_sites)+1):
                                    distance = match.get(f'reserve_site_{i}_distance', None)
                                    row.append(distance if distance is not None else '')
                                
                                writer.writerow(row)
                                total_structures += 1
        
        runtime = datetime.now() - start_time
        with open("analysis_summary.txt", 'a') as f:
            f.write(f"Genome Analysis Summary\n")
            f.write(f"Start Time: {start_time}\n")
            f.write(f"End Time: {datetime.now()}\n")
            f.write(f"Total Runtime: {runtime}\n")
            f.write(f"Total Structures Found: {total_structures}\n")
            f.write(f"Maximum Memory Usage: {max_memory_usage:.2f} MB\n")  
            f.write("\nAnalysis Parameters:\n")
            for key, value in self.params.items():
                if key != 'motif':
                    f.write(f"{key}: {value}\n")
            f.write(f"motif: {self.params['motif']}\n")
            f.write("conserve_sites: " + ', '.join([f"{site['amino_acid']}{site['position']}" for site in self.reserve_sites]) + "\n")

def main():
    """
    Main function for the transposon analysis pipeline.
    Handles command line argument parsing and initiates the analysis process.
    """
    import argparse
    parser = argparse.ArgumentParser(description='Analyze genome for transposons')
    
    # Input/output parameters
    parser.add_argument('--input', '-i', type=str, required=True,
                      help='Input genome file in FASTA format')
    parser.add_argument('--output-prefix', '-o', type=str, required=True,
                      help='Prefix for output files')
    
    # Analysis parameters
    parser.add_argument('--min-tsd-pattern-size', type=int, default=5,
                      help='Minimum TSD pattern size for search (default: 5)')
    parser.add_argument('--max-tsd-pattern-size', type=int, default=10,
                      help='Maximum TSD pattern size for search (default: 10)')
    parser.add_argument('--gap-size', type=int, default=2500,
                      help='Expected transposable element search gap (default: 2500 means length between 2500-5000bp)')
    parser.add_argument('--min-tir-size', type=int, default=5,
                      help='Minimum TIR size (default: 5)')
    parser.add_argument('--max-tir-size', type=int, default=15,
                      help='Maximum TIR size (default: 15)')
    parser.add_argument('--max-tir-mismatch', type=int, default=1,
                      help='Maximum allowed mismatches in TIR (default: 1)')
    parser.add_argument('--conserve-site-range', type=int, default=50,
                      help='Search range around conserved reserve sites positions (default: 50)')
    parser.add_argument('--subterminal-threshold', type=float, default=5.0,
                  help='Threshold for subterminal motif enrichment percentage (default: 5.0)')
    parser.add_argument('--subterminal-length', type=int, default=250,
                      help='Length of subterminal region to search for motifs (default: 250)')
    parser.add_argument('--motif', type=str, default='AAAGGG',
                      help='User-specified subterminal motif (default: AAAGGG)')
    parser.add_argument('--chunk-size', type=int, default=25000,
                      help='Size of sequence chunks for parallel processing (default: 25000)')
    parser.add_argument('--num-processes', type=int, default=10,
                      help='Number of parallel processes (default: 10)')
    
    # Reserve sites
    parser.add_argument('--conserve-sites', type=str, default=None,
    help='Comma-separated conserve sites in format [AminoAcid][Position] (e.g., D301,E719,Q500). Default: D301,D376,E719')

    
    args = parser.parse_args()
    
    # Default Conserved Sites
    default_reserve_sites = ['D301', 'D376', 'E719']
    
    valid_amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    # Checking user input Conserved Sites
    if args.conserve_sites:
        reserve_sites = [site.strip() for site in args.conserve_sites.split(',') if site.strip()]
        print(f"Using user-provided conserve sites: {', '.join(reserve_sites)}")
        
        # Validation of each conserved site
        for site in reserve_sites:
            # Extract letters and numbers
            match = re.match(r'([A-Za-z])(\d+)', site)
            if not match:
                print(f"Error: Invalid format for reserve site: {site}. Must be in the form '[AminoAcid][Position]' (e.g., D301).")
                return  # Exit the program if format is invalid

            letter, number = match.groups()
            number = int(number)
            
            # 1. Check that the letters are valid
            if letter not in valid_amino_acids:
                print(f"Error: Invalid amino acid letter '{letter}' in conserve site: {site}. Must be one of {', '.join(valid_amino_acids)}.")
                return  # Exit the program if letter is invalid
            
            # 2. Check that the numbers are not out of range
            max_cds_distance = (2500 + args.gap_size) - (args.subterminal_length * 2)
            if number * 3 >= max_cds_distance:
                print(f"Error: Conserve site {site} is out of range. The maximum valid position is {(max_cds_distance // 3)}.")
                return  # Exit the program if number is out of range
    else:
        print("No --conserve-sites provided. Using default reserve sites: D301, D376, E719.")
        reserve_sites = default_reserve_sites

    # Calculate the minimum CDS distance
    numbers = [int(re.search(r'\d+', site).group()) for site in reserve_sites]
    min_cds_distance = max(numbers) * 3
    
    # Build parameter dictionary
    params = {
        'min_tsd_pattern_size': args.min_tsd_pattern_size,
        'max_tsd_pattern_size': args.max_tsd_pattern_size,
        'gap_size': args.gap_size,
        'min_tir_size': args.min_tir_size,
        'max_tir_size': args.max_tir_size,
        'max_tir_mismatch': args.max_tir_mismatch,
        'conserve_site_range': args.conserve_site_range,
        'subterminal_threshold': args.subterminal_threshold,
        'subterminal_length': args.subterminal_length,
        'motif': args.motif,
        'chunk_size': args.chunk_size,
        'num_processes': args.num_processes,
        'min_cds_distance': min_cds_distance  # Add min_cds_distance to params
    }
    
    # Create analyzer instance and run analysis
    analyzer = TransposonAnalyzer(**params)
    analyzer.run_analysis(args.input, args.output_prefix, reserve_sites)

if __name__ == "__main__":
    main()