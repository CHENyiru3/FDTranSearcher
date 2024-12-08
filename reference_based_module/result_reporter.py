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
import os

from GFFEntry_class import GFFEntry,parse_gff_attributes

def generate_report(analysis_results: List[Dict], output_file: str):
    with open(output_file, 'w') as f:
        f.write("Transposable Element Analysis Report\n")
        f.write(f"Analysis Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Total Transposable Elements with Protein: {len(analysis_results)}\n\n")

        for result in analysis_results:
            f.write(f"TE ID: {result['te_id']}\n")

            # position information
            f.write("Genomic Positions:\n")
            f.write(f"  BLAST Match Position: {result['blast_position']}\n")
            f.write(f"  Element Structure Position: {result['element_position']}\n")

            # detials of protein regions
            f.write("Protein Regions:\n")
            for i, (start, end) in enumerate(result['protein_regions'], 1):
                f.write(f"  Region {i}: {start}-{end}\n")

            # strucure information
            f.write("\nTransposable Element Structures:\n")
            for i, structure in enumerate(result['structures'][:3], 1):
                f.write(f"  Structure {i}:\n")
                for key in ['TSD1', 'TIR1', 'TIR2', 'TSD2']:
                    start, end = structure[key]
                    f.write(f"    {key}: {start}-{end} ({end - start + 1} bp)\n")
                    f.write(f"    {key} Sequence: {structure[f'{key.lower()}_seq']}\n")

                f.write("\n")

            f.write("-" * 50 + "\n\n")


def generate_csv_report(analysis_results: List[Dict], output_file: str):
    """
    Generates a CSV file based on the analysis results.
    
    Args:
        analysis_results (list): List of dictionaries containing analysis results.
        output_file (str): Path to the output CSV file.

    Returns:
        None 
    """
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = [
            'te_id',
            'blast_chromosome',
            'blast_start',
            'blast_end',
            'element_chromosome',
            'element_start',
            'element_end',
            'protein_region_start',
            'protein_region_end'
        ]

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for result in analysis_results:
            for protein_region in result['protein_regions']:
                writer.writerow({
                    'te_id': result['te_id'],
                    'blast_chromosome': result['blast_match'].get('chromosome', 'N/A'),
                    'blast_start': result['blast_match'].get('start', 'N/A'),
                    'blast_end': result['blast_match'].get('end', 'N/A'),
                    'element_chromosome': result['blast_match'].get('chromosome', 'N/A'),
                    'element_start': result['blast_match'].get('start', 'N/A') + result['te_structure_position']['start'],
                    'element_end': result['blast_match'].get('start', 'N/A') + result['te_structure_position']['end'],
                    'protein_region_start': protein_region[0],
                    'protein_region_end': protein_region[1]
                })


def generate_gff3_report(analysis_results: List[Dict], output_file: str):
    """
    Generates a GFF3 file based on the analysis results.

    Args:
        analysis_results (list): List of dictionaries containing analysis results.
        output_file (str): Path to the output GFF3 file.
    
    """
    with open(output_file, 'w') as f:
        # Write the GFF3 header
        f.write("##gff-version 3\n")

        for result in analysis_results:
            # get the information from the analysis results
            te_id = result['te_id']
            blast_match = result['blast_match']
            te_structure_position = result['te_structure_position']
            protein_regions = result['protein_regions']
            structures = result['structures']

            # if no structures found, skip
            if not structures:
                print(f"Warning: No structures found for TE {te_id}")
                continue

            # filter the best structure
            best_structure = max(structures, key=lambda x: x['scores']['total_score'])

            # extract the attributes
            attributes = []
            attributes.append(f"ID={te_id}")
            attributes.append(f"BLAST_Match={blast_match['chromosome']}:{blast_match['start']}-{blast_match['end']}")

            # position information of the element
            if te_structure_position:
                attributes.append(f"Element_Structure={blast_match['chromosome']}:{blast_match['start'] + te_structure_position['start']}-{blast_match['start'] + te_structure_position['end']}")

            # add the protein regions
            for i, (start, end) in enumerate(protein_regions, 1):
                attributes.append(f"Protein_Region_{i}={start}-{end}")

            # add the structure information
            for key in ['TSD1', 'TIR1', 'TIR2', 'TSD2']:
                if key in best_structure:
                    start, end = best_structure[key]
                    seq = best_structure.get(f"{key.lower()}_seq", "")
                    attributes.append(f"{key}={start}-{end};{key}_Sequence={seq}")

            # write the GFF3 entry
            f.write(f"{blast_match['chromosome']}\t{result['source']}\t{result['type']}\t"
                    f"{blast_match['start']}\t{blast_match['end']}\t{result['score']}\t"
                    f"{result['strand']}\t{result['phase']}\t"
                    f"{';'.join(attributes)}\n")
