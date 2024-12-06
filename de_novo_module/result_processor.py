#!/usr/bin/env python3
import pandas as pd
import numpy as np
from collections import defaultdict
from typing import List, Dict, Tuple

def process_raw_results(input_file: str, output_csv: str, output_gff: str):
    """
    Process raw transposon analysis results into final CSV and GFF3 formats.
    
    Args:
        input_file: Path to input raw CSV file
        output_csv: Path for output processed CSV file
        output_gff: Path for output GFF3 file
    """
    print("Processing analysis results...")
    
    # Read CSV file
    df = pd.read_csv(input_file)
    
    # Group transposons by unique identifiers
    tsd_groups = create_tsd_groups(df)
    
    # Process groups and generate output files
    process_groups(tsd_groups, output_csv, output_gff)

def create_tsd_groups(df: pd.DataFrame) -> Dict:
    """
    Create groups of transposons based on their TSD/TIR positions.
    
    Args:
        df: Input DataFrame containing raw transposon data
        
    Returns:
        Dictionary of grouped transposon data
    """
    tsd_groups = defaultdict(lambda: {
        'records': [],
        'strands': [],
        'cds_starts': set(),
        'cds_ends': set(),
        'd301_positions': set(),
        'd367_positions': set(),
        'e719_positions': set()
    })
    
    for _, row in df.iterrows():
        key = (row['Sequence_ID'], row['TSD1_Start'], row['TSD1_End'],
               row['TIR1_Start'], row['TIR1_End'], row['TIR2_Start'],
               row['TIR2_End'], row['TSD2_Start'], row['TSD2_End'])
        
        group = tsd_groups[key]
        group['records'].append(row)
        group['strands'].append(row['Strand'])
        
        if row['CDS_Start'] is not None and row['CDS_End'] is not None:
            group['cds_starts'].add(row['CDS_Start'])
            group['cds_ends'].add(row['CDS_End'])
        
        if row['D301_Position'] is not None:
            group['d301_positions'].add(row['D301_Position'])
        if row['D367_Position'] is not None:
            group['d367_positions'].add(row['D367_Position'])
        if row['E719_Position'] is not None:
            group['e719_positions'].add(row['E719_Position'])
    
    return tsd_groups

def process_positions(positions: list, cds_start: float) -> tuple:
    """
    Process position information for conserved sites.
    
    Args:
        positions: List of position values
        cds_start: CDS start position for distance calculation
        
    Returns:
        Tuple of (count, mean_position, mean_distance_to_cds)
    """
    positions = [p for p in positions if p is not None]
    if not positions:
        return 0, None, None
    num = len(positions)
    mean_pos = np.mean(positions)
    mean_distance = mean_pos - cds_start if cds_start is not None else None
    return num, mean_pos, mean_distance

def process_groups(tsd_groups: Dict, output_csv: str, output_gff: str):
    """
    Process grouped transposon data and generate output files.
    
    Args:
        tsd_groups: Dictionary of grouped transposon data
        output_csv: Path for output processed CSV file
        output_gff: Path for output GFF3 file
    """
    output_records = []
    gff_records = ['##gff-version 3']
    
    for key, group in tsd_groups.items():
        base_record = group['records'][0]
        
        # Calculate dominant strand
        strand_counts = pd.Series(group['strands']).value_counts()
        dominant_strand = strand_counts.index[0]
        
        # Convert sets to lists for calculations
        unique_cds_starts = list(group['cds_starts'])
        unique_cds_ends = list(group['cds_ends'])
        unique_d301_positions = list(group['d301_positions'])
        unique_d367_positions = list(group['d367_positions'])
        unique_e719_positions = list(group['e719_positions'])
        
        # Calculate statistics
        mean_cds_start = int(np.mean(unique_cds_starts)) if unique_cds_starts else None
        mean_cds_end = int(np.mean(unique_cds_ends)) if unique_cds_ends else None
        mean_coding_length = int(mean_cds_end - mean_cds_start) if (mean_cds_start is not None and mean_cds_end is not None) else None
        mean_protein_length = int(mean_coding_length // 3) if mean_coding_length is not None else None
        
        # Calculate distances
        mean_distance_left = int(mean_cds_start - base_record['TIR1_End']) if mean_cds_start is not None else None
        mean_distance_right = int(base_record['TIR2_End'] - mean_cds_end) if mean_cds_end is not None else None
        
        # Process conserved positions
        d301_num, d301_mean, d301_distance = process_positions(unique_d301_positions, mean_cds_start)
        d367_num, d367_mean, d367_distance = process_positions(unique_d367_positions, mean_cds_start)
        e719_num, e719_mean, e719_distance = process_positions(unique_e719_positions, mean_cds_start)
        
        # Create output record
        output_record = create_output_record(
            base_record, dominant_strand, unique_cds_starts,
            mean_cds_start, mean_cds_end, mean_coding_length, mean_protein_length,
            mean_distance_left, mean_distance_right,
            d301_num, d301_mean, d301_distance,
            d367_num, d367_mean, d367_distance,
            e719_num, e719_mean, e719_distance
        )
        
        output_records.append(output_record)
        gff_records.append(create_gff_record(base_record, output_record, dominant_strand))
    
    # Save output files
    output_df = pd.DataFrame(output_records)
    output_df.to_csv(output_csv, index=False)
    
    with open(output_gff, 'w') as f:
        f.write('\n'.join(gff_records))
    
    print(f"Results processing complete.")
    print(f"Processed {len(output_records)} unique transposon structures.")

def create_output_record(base_record: Dict, dominant_strand: str, unique_cds_starts: List,
                        mean_cds_start: int, mean_cds_end: int, mean_coding_length: int,
                        mean_protein_length: int, mean_distance_left: int,
                        mean_distance_right: int, d301_num: int, d301_mean: float,
                        d301_distance: float, d367_num: int, d367_mean: float,
                        d367_distance: float, e719_num: int, e719_mean: float,
                        e719_distance: float) -> Dict:
    """
    Create a processed output record for a transposon.
    
    Args:
        Various statistical and positional information about the transposon
        
    Returns:
        Dictionary containing processed transposon information
    """
    return {
        'Sequence_ID': base_record['Sequence_ID'],
        'TSD1_Start': base_record['TSD1_Start'],
        'TSD1_End': base_record['TSD1_End'],
        'TIR1_Start': base_record['TIR1_Start'],
        'TIR1_End': base_record['TIR1_End'],
        'TIR2_Start': base_record['TIR2_Start'],
        'TIR2_End': base_record['TIR2_End'],
        'TSD2_Start': base_record['TSD2_Start'],
        'TSD2_End': base_record['TSD2_End'],
        'TSD_Sequence': base_record['TSD_Sequence'],
        'TIR1_Sequence': base_record['TIR1_Sequence'],
        'TIR2_Sequence': base_record['TIR2_Sequence'],
        'Left_Subterminal_AAAGGG_%': base_record['Left_Subterminal_AAAGGG_%'],
        'Right_Subterminal_AAAGGG_%': base_record['Right_Subterminal_AAAGGG_%'],
        'Strand': dominant_strand,
        'CDS_Number': len(unique_cds_starts),
        'Mean_CDS_Start': mean_cds_start,
        'Mean_CDS_End': mean_cds_end,
        'Mean_Coding_Sequence_Length': mean_coding_length,
        'Mean_Protein_Length': mean_protein_length,
        'Mean_Distance_to_Left_TIR': mean_distance_left,
        'Mean_Distance_to_Right_TIR': mean_distance_right,
        'D301_Number': d301_num,
        'Mean_D301_Position': int(d301_mean) if d301_mean is not None else None,
        'Mean_D301_CDS_Start_Distance': int(d301_distance) if d301_distance is not None else None,
        'D367_Number': d367_num,
        'Mean_D367_Position': int(d367_mean) if d367_mean is not None else None,
        'Mean_D367_CDS_Start_Distance': int(d367_distance) if d367_distance is not None else None,
        'E719_Number': e719_num,
        'Mean_E719_Position': int(e719_mean) if e719_mean is not None else None,
        'Mean_E719_CDS_Start_Distance': int(e719_distance) if e719_distance is not None else None
    }

def create_gff_record(base_record: Dict, output_record: Dict, dominant_strand: str) -> str:
    """
    Create a GFF3 format record for a transposon.
    
    Args:
        base_record: Original transposon record
        output_record: Processed transposon record
        dominant_strand: Dominant strand orientation
        
    Returns:
        String containing GFF3 formatted record
    """
    attributes = [
        f"ID=TE_{base_record['Sequence_ID']}_{int(base_record['TSD1_Start'])}_{int(base_record['TSD2_End'])}"
    ]
    
    for key, value in output_record.items():
        if value is not None:
            if isinstance(value, float):
                value = f"{value:.2f}"
            elif isinstance(value, int):
                value = str(value)
            attributes.append(f"{key}={value}")
    
    gff_fields = [
        base_record['Sequence_ID'],
        'FDT_de_novo_pipeline',
        'predicted_functional_transposable_element',
        str(int(base_record['TSD1_Start'])),
        str(int(base_record['TSD2_End'])),
        '0.0',
        dominant_strand,
        '.',
        ';'.join(attributes)
    ]
    
    return '\t'.join(gff_fields)