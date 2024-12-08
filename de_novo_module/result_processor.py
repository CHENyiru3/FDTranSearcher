import pandas as pd
import numpy as np
import logging
from collections import defaultdict
from typing import List, Dict, Any

def process_results(input_file: str, output_csv: str, output_gff: str, reserve_sites: List[Dict[str, Any]], params: Dict[str, Any]):
    """
    Process raw transposon analysis results into final CSV and GFF3 formats.

    Args:
        input_file: Path to input raw CSV file
        output_csv: Path for output processed CSV file
        output_gff: Path for output GFF3 file
        reserve_sites: List of reserve site dictionaries
        params: Dictionary containing analysis parameters
    """
    logging.info("Processing analysis results...")
    print("Processing analysis results...")

    try:
        # Read CSV file
        df = pd.read_csv(input_file)
    except FileNotFoundError as e:
        logging.error(f"Raw CSV file not found: {e}")
        raise
    except pd.errors.EmptyDataError as e:
        logging.error(f"Raw CSV file is empty: {e}")
        raise
    except Exception as e:
        logging.error(f"Error reading raw CSV file: {e}")
        raise

    # Group transposons by unique identifiers
    tsd_groups = defaultdict(lambda: {
        'records': [],
        'strands': [],
        'cds_starts': set(),
        'cds_ends': set(),
        'reserve_site_distances': defaultdict(set)
    })

    # First traversal: collect all data
    for _, row in df.iterrows():
        key = (row['Sequence_ID'], row['TSD1_Start'], row['TSD1_End'],
            row['TIR1_Start'], row['TIR1_End'], row['TIR2_Start'],
            row['TIR2_End'], row['TSD2_Start'], row['TSD2_End'])
    
        group = tsd_groups[key]
        group['records'].append(row)
        group['strands'].append(row['Strand'])
    
        # Storing unique CDS locations
        if not pd.isnull(row['CDS_Start']) and not pd.isnull(row['CDS_End']):
            group['cds_starts'].add(row['CDS_Start'])
            group['cds_ends'].add(row['CDS_End'])
    
        # Storing reserve site distances
        for i in range(1, len(reserve_sites)+1):
            distance = row.get(f'conserve_site_{i}_distance', None)
            if not pd.isnull(distance):
                group['reserve_site_distances'][f'reserve_site_{i}_distance'].add(distance)

    # Preparing data for output
    output_records = []
    gff_records = ['##gff-version 3']

    # Treatment of each TSD group
    for key, group in tsd_groups.items():
        # Take the basic information of the first record
        base_record = group['records'][0]
    
        # Calculate strand (take the one with the most)
        strand_counts = pd.Series(group['strands']).value_counts()
        dominant_strand = strand_counts.index[0]
    
        # Converting from a collection to a list for computation
        unique_cds_starts = list(group['cds_starts'])
        unique_cds_ends = list(group['cds_ends'])
        unique_reserve_distances = {k: list(v) for k, v in group['reserve_site_distances'].items()}

        # Calculation of CDS statistics
        mean_cds_start = np.mean(unique_cds_starts) if unique_cds_starts else None
        mean_cds_end = np.mean(unique_cds_ends) if unique_cds_ends else None
        mean_coding_length = mean_cds_end - mean_cds_start if (mean_cds_start is not None and mean_cds_end is not None) else None
        mean_protein_length = int(mean_coding_length // 3) if mean_coding_length is not None else None

        # Calculate the average distance from the TIR
        mean_distance_left = mean_cds_start - base_record['TIR1_End'] if mean_cds_start is not None else None
        mean_distance_right = base_record['TIR2_End'] - mean_cds_end if mean_cds_end is not None else None
    
        # Computationally Conserved Loci Related Information
        def process_positions(positions: list, cds_start: float) -> tuple:
            positions = [p for p in positions if p is not None]
            if not positions:
                return 0, None, None
            num = len(positions)
            mean_pos = np.mean(positions)
            mean_distance = mean_pos - cds_start if cds_start is not None else None
            return num, mean_pos, mean_distance
    
        reserve_site_info = {}
        for i, site in enumerate(reserve_sites, start=1):
            distances = unique_reserve_distances.get(f'reserve_site_{i}_distance', [])
            num, mean_pos, mean_distance = process_positions(distances, mean_cds_start)
            reserve_site_info[f'Possible_Conserve_Site_{i}_Number'] = num
            reserve_site_info[f'Mean_Possible_Conserve_Site_{i}_Distance'] = int(mean_pos) if mean_pos is not None else None
    
        # Creating Output Records
        output_record = {
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
            f'Left_Subterminal_{params["motif"]}_%': base_record[f'Left_Subterminal_{params["motif"]}_%'],
            f'Right_Subterminal_{params["motif"]}_%': base_record[f'Right_Subterminal_{params["motif"]}_%'],
            'Strand': dominant_strand,
            'CDS_Start': int(mean_cds_start) if mean_cds_start is not None else None,
            'CDS_End': int(mean_cds_end) if mean_cds_end is not None else None,
            'CDS_Length': int(mean_coding_length) if mean_coding_length is not None else None,
            'Protein_Length': mean_protein_length,
            'Mean_Distance_to_Left_TIR': int(mean_distance_left) if mean_distance_left is not None else None,
            'Mean_Distance_to_Right_TIR': int(mean_distance_right) if mean_distance_right is not None else None,
        }
    
        # Add reserve site info
        output_record.update(reserve_site_info)
    
        # Add Pattern_Matching_Score as the last column
        output_record['Pattern_Matching_Score'] = base_record['Pattern_Match']
    
        output_records.append(output_record)
    
        # Preparation of GFF3 records
        attributes = [
            f"ID=TE_{base_record['Sequence_ID']}_{int(base_record['TSD1_Start'])}_{int(base_record['TSD2_End'])}"
        ]
    
        # Add all other properties
        for key, value in output_record.items():
            if value is not None and key not in ['Sequence_ID', 'Pattern_Match']:
                if isinstance(value, float):
                    value = f"{value:.2f}"
                elif isinstance(value, int):
                    value = str(value)
                attributes.append(f"{key}={value}")
    
        # Constructing GFF3 records
        gff_record = [
            base_record['Sequence_ID'],                   # seqid
            'FDT_de_novo_pipeline',                       # source
            'predicted_functional_transposable_element',  # type
            str(int(base_record['TSD1_Start'])),          # start
            str(int(base_record['TSD2_End'])),            # end
            '0.0',                                        # score
            dominant_strand,                              # strand
            '.',                                          # phase
            ';'.join(attributes)                          # attributes
        ]
    
        gff_records.append('\t'.join(gff_record))

    # Create output DataFrame and save CSV
    output_df = pd.DataFrame(output_records)
    try:
        output_df.to_csv(output_csv, index=False)
        logging.info(f"Processed CSV saved to {output_csv}")
    except Exception as e:
        logging.error(f"Error saving processed CSV: {e}")
        raise

    # Save GFF3 file
    try:
        with open(output_gff, 'w') as f:
            f.write('\n'.join(gff_records))
        logging.info(f"GFF3 file saved to {output_gff}")
    except Exception as e:
        logging.error(f"Error saving GFF3 file: {e}")
        raise

    print(f"Results processing complete.")
    print(f"Processed {len(output_records)} unique transposon structures.")