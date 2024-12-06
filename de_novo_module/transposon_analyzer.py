#!/usr/bin/env python3
import csv
import psutil
import os
from datetime import datetime
from multiprocessing import Pool
from tqdm import tqdm
from Bio import SeqIO
from typing import List, Tuple
import argparse

from structure_verification import verify_transposon_structure
from result_processor import process_raw_results

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
        self.max_memory_usage = 0 
        # Set default parameters
        self.params = {
            'min_tsd_pattern_size': 5,
            'max_tsd_pattern_size': 10,
            'gap_size': 2500,
            'min_tir_size': 5,
            'max_tir_size': 15,
            'max_tir_mismatch': 1,
            'search_dde_range': 50,
            'chunk_size': 25000,
            'num_processes': 10,
            'subterminal_threshold': 5.0
        }
        # Update with user-provided parameters
        self.params.update(params)
        
        # Validate parameters
        self._validate_params()

    def _validate_params(self):
        """
        Validate parameter values.
        
        Raises:
            ValueError: If any parameter values are invalid
        """
        if self.params['search_dde_range'] < 0:
            raise ValueError("Search range must be non-negative")
        if self.params['gap_size'] < 0: 
            raise ValueError("Gap size must be non-negative")
        if self.params['min_tsd_pattern_size'] > self.params['max_tsd_pattern_size']:
            raise ValueError("min_tsd_pattern_size cannot be greater than max_tsd_pattern_size")
        if self.params['min_tir_size'] > self.params['max_tir_size']:
            raise ValueError("min_tir_size cannot be greater than max_tir_size")
        if self.params['subterminal_threshold'] < 0:
            raise ValueError("Subterminal threshold must be non-negative")

    def run_analysis(self, input_file: str, output_prefix: str):
        """
        Run complete transposon analysis pipeline.
        
        Args:
            input_file: Path to input FASTA file
            output_prefix: Prefix for output files
        """
        # Phase 1: Search transposons
        initial_output = f"{output_prefix}_raw.csv"
        self._search_transposons(input_file, initial_output)
        
        # Phase 2: Process results
        final_csv = f"{output_prefix}_processed.csv"
        final_gff = f"{output_prefix}.gff3"
        process_raw_results(initial_output, final_csv, final_gff)
        
        print(f"Analysis complete. Results written to:")
        print(f"  Raw CSV: {initial_output}")
        print(f"  Processed CSV: {final_csv}")
        print(f"  GFF3: {final_gff}")

    def _update_max_memory_usage(self):
        """
        Update maximum memory usage tracking.
        """
        current_process = psutil.Process(os.getpid())
        memory_usage = current_process.memory_info().rss / 1024 / 1024  # Convert to MB
        self.max_memory_usage = max(self.max_memory_usage, memory_usage)

    def _extract_sequence(self, input_file: str) -> List[Tuple[str, str]]:
        """
        Read sequences from FASTA file.
        
        Args:
            input_file: Path to input FASTA file
            
        Returns:
            List of tuples containing (sequence_id, sequence)
        """
        print(f"Reading sequence from {input_file}")
        sequences = []
        total_length = 0
        for record in SeqIO.parse(input_file, "fasta"):
            print(f"Found sequence: {record.id}, length: {len(record.seq)}")
            sequences.append((record.id, str(record.seq).upper()))
            total_length += len(record.seq)
        print(f"Total sequence length: {total_length} bp")
        return sequences

    def _process_chunk(self, args: Tuple) -> Tuple[str, List[dict]]:
        """
        Process a single sequence chunk.
        
        Args:
            args: Tuple containing (chunk_id, sequence, seq_id, base_position)
            
        Returns:
            Tuple of (chunk_id, list of matches)
        """
        chunk_id, sequence, seq_id, base_pos = args
        self._update_max_memory_usage()
        process = psutil.Process(os.getpid())
        mem_usage = process.memory_info().rss / 1024 / 1024
        
        print(f"Processing chunk {chunk_id} of {seq_id} (Memory usage: {mem_usage:.2f} MB)")
        
        try:
            matches = list(verify_transposon_structure(sequence, seq_id, self.params))
            
            for match in matches:
                for key in ['tsd1_start', 'tsd1_end', 'tir1_start', 'tir1_end',
                           'tir2_start', 'tir2_end', 'tsd2_start', 'tsd2_end',
                           'cds_start', 'cds_end', 'd301_position', 'd367_position',
                           'e719_position']:
                    if match[key] is not None:  # Only adjust non-None values
                        match[key] += base_pos
            
            return chunk_id, matches
        except Exception as e:
            print(f"Error in chunk {chunk_id}: {str(e)}")
            return chunk_id, []

    def _search_transposons(self, input_file: str, output_file: str):
        """
        Execute transposon search and save raw results.
        
        Args:
            input_file: Path to input FASTA file
            output_file: Path for output raw CSV file
        """
        start_time = datetime.now()
        print(f"Analysis started at: {start_time}")
        
        sequences = self._extract_sequence(input_file)
        chunk_size = self.params['chunk_size']
        num_processes = self.params['num_processes']
        
        self._update_max_memory_usage()
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([
                'Sequence_ID',
                'TSD1_Start', 'TSD1_End',
                'TIR1_Start', 'TIR1_End',
                'TIR2_Start', 'TIR2_End',
                'TSD2_Start', 'TSD2_End',
                'TSD_Sequence',
                'TIR1_Sequence',
                'TIR2_Sequence',
                'Left_Subterminal_AAAGGG_%',
                'Right_Subterminal_AAAGGG_%',
                'Strand',
                'CDS_Start',
                'CDS_End',
                'CDS_Length',
                'Protein_Length',
                'Distance_to_Left_TIR',
                'Distance_to_Right_TIR',
                'D301_Position',
                'D367_Position',
                'E719_Position'
            ])
            
            total_structures = 0
            
            for seq_id, seq in sequences:
                print(f"\nProcessing sequence: {seq_id}")
                chunks = []
                for i in range(0, len(seq), chunk_size):
                    chunk = seq[i:min(i+chunk_size, len(seq))]
                    chunk_id = f"{seq_id}_chunk_{i//chunk_size}"
                    chunks.append((chunk_id, chunk, seq_id, i))
                
                with Pool(processes=num_processes) as pool:
                    for result in tqdm(pool.imap_unordered(self._process_chunk, chunks), 
                                     total=len(chunks), 
                                     desc=f"Processing {seq_id}"):
                        chunk_id, matches = result
                        if matches:
                            for match in matches:
                                writer.writerow([
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
                                    match['d301_position'],
                                    match['d367_position'],
                                    match['e719_position']
                                ])
                                total_structures += 1
        
        runtime = datetime.now() - start_time
        with open("analysis_summary.txt", 'w') as f:
            f.write(f"Genome Analysis Summary\n")
            f.write(f"Start Time: {start_time}\n")
            f.write(f"End Time: {datetime.now()}\n")
            f.write(f"Total Runtime: {runtime}\n")
            f.write(f"Total Structures Found: {total_structures}\n")
            f.write(f"Maximum Memory Usage: {self.max_memory_usage:.2f} MB\n")
            f.write("\nAnalysis Parameters:\n")
            for key, value in self.params.items():
                f.write(f"{key}: {value}\n")
                


def main():
    parser = argparse.ArgumentParser(description='De novo Transposon Analysis Algorithm')
    
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
    parser.add_argument('--search-dde-range', type=int, default=50,
                        help='Search range around conserved DDE positions (default: 50)')
    parser.add_argument('--subterminal-threshold', type=float, default=5.0,
                        help='Threshold for subterminal AAAGGG/CCCTTT enrichment percentage (default: 5.0)')
    parser.add_argument('--chunk-size', type=int, default=25000,
                        help='Size of sequence chunks for parallel processing (default: 25000)')
    parser.add_argument('--num-processes', type=int, default=10,
                        help='Number of parallel processes (default: 10)')
    
    args = parser.parse_args()

    # Build parameter dictionary
    params = {
        'min_tsd_pattern_size': args.min_tsd_pattern_size,
        'max_tsd_pattern_size': args.max_tsd_pattern_size,
        'gap_size': args.gap_size,
        'min_tir_size': args.min_tir_size,
        'max_tir_size': args.max_tir_size,
        'max_tir_mismatch': args.max_tir_mismatch,
        'search_dde_range': args.search_dde_range,
        'subterminal_threshold': args.subterminal_threshold,
        'chunk_size': args.chunk_size,
        'num_processes': args.num_processes
    }
    
    # Create analyzer instance and run analysis
    analyzer = TransposonAnalyzer(**params)
    analyzer.run_analysis(args.input, args.output_prefix)

if __name__ == "__main__":
    main()