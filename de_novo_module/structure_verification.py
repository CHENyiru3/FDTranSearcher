#!/usr/bin/env python3
from typing import List, Tuple, Dict, Optional, Generator
from sequence_tools import kmp_search, get_reverse_complement, analyze_subterminal_regions
from element_search import find_conserved_codon_positions

def find_cds_region(seq: str, window_size: int = 5000, min_cds_length: int = 2160,
                   search_range: int = 50) -> List[Tuple[int, int, int, str, Tuple[int, int, int]]]:
    """
    Search for protein-coding regions (CDS) with conserved DDE motifs.
    
    Args:
        seq: Input DNA sequence
        window_size: Maximum size of window to search for stop codons
        min_cds_length: Minimum length requirement for CDS
        search_range: Range to search for conserved codons
        
    Returns:
        List of tuples containing (start_pos, end_pos, protein_length, strand, conserved_positions)
    """
    results = []
    stop_codons = ["TAA", "TAG", "TGA"]
    
    # Search forward strand
    start_positions = kmp_search(seq, "ATG")
    for start_pos in start_positions:
        window_end = min(start_pos + window_size, len(seq)-2)
        for stop_codon in stop_codons:
            stop_positions = kmp_search(seq[start_pos:window_end], stop_codon)
            for stop_pos in stop_positions:
                actual_stop_pos = start_pos + stop_pos
                cds_length = actual_stop_pos - start_pos + 3
                if cds_length >= min_cds_length:
                    try:
                        cds_seq = seq[start_pos:actual_stop_pos+3]
                        conserved_positions = find_conserved_codon_positions(cds_seq, '+', search_range)
                        if all(pos is not None for pos in conserved_positions):
                            protein_length = cds_length // 3
                            abs_conserved_positions = tuple(pos + start_pos if pos is not None else None 
                                                         for pos in conserved_positions)
                            results.append((start_pos, actual_stop_pos + 3, protein_length, "+", abs_conserved_positions))
                    except Exception:
                        continue
    
    # Search reverse strand
    rc_seq = get_reverse_complement(seq)
    start_positions = kmp_search(rc_seq, "ATG")
    for start_pos in start_positions:
        window_end = min(start_pos + window_size, len(rc_seq)-2)
        for stop_codon in stop_codons:
            stop_positions = kmp_search(rc_seq[start_pos:window_end], stop_codon)
            for stop_pos in stop_positions:
                actual_stop_pos = start_pos + stop_pos
                cds_length = actual_stop_pos - start_pos + 3
                if cds_length >= min_cds_length:
                    try:
                        cds_seq = rc_seq[start_pos:actual_stop_pos+3]
                        conserved_positions = find_conserved_codon_positions(cds_seq, '-', search_range)
                        if all(pos is not None for pos in conserved_positions):
                            protein_length = cds_length // 3
                            orig_start = len(seq) - (actual_stop_pos + 3)
                            orig_end = len(seq) - start_pos
                            abs_conserved_positions = tuple(len(seq) - (pos + start_pos) if pos is not None else None 
                                                         for pos in conserved_positions[::-1])
                            results.append((orig_start, orig_end, protein_length, "-", abs_conserved_positions))
                    except Exception:
                        continue
    
    return results

def verify_transposon_structure(seq: str, seq_id: str, params: Dict) -> Generator[Dict, None, None]:
    """
    Verify and identify complete transposon structures.
    
    Args:
        seq: Input DNA sequence
        seq_id: Sequence identifier
        params: Dictionary of search parameters including thresholds and sizes
        
    Returns:
        Generator yielding dictionaries containing identified transposon structures
    """
    seq = seq.upper()
    
    from element_search import build_pattern_hash, find_tir_near_tsd
    
    pattern_hashes = {
        size: build_pattern_hash(seq, size)
        for size in range(params['min_tsd_pattern_size'], 
                         params['max_tsd_pattern_size'] + 1)
    }
    
    for pattern_size, pattern_hash in pattern_hashes.items():
        for pattern, positions in pattern_hash.items():
            for i, tsd1_start in enumerate(positions):
                for tsd2_start in positions[i + 1:]:
                    gap_size = params['gap_size']
                    if not (gap_size <= tsd2_start - (tsd1_start + pattern_size) <= gap_size * 2):
                        continue
                    
                    tir_result = find_tir_near_tsd(
                        seq, tsd1_start, tsd2_start, pattern_size,
                        params['min_tir_size'], params['max_tir_size'],
                        params['max_tir_mismatch']
                    )
                    
                    if tir_result:
                        tir1, tir2, tir_size = tir_result
                        
                        # Process subterminal regions
                        subterminal_left = seq[tsd1_start + pattern_size + tir_size:
                                             tsd1_start + pattern_size + tir_size + 250]
                        subterminal_right = seq[tsd2_start - tir_size - 250:tsd2_start - tir_size]
                        
                        left_percentage, left_strand = analyze_subterminal_regions(
                            subterminal_left, params['subterminal_threshold'])
                        right_percentage, right_strand = analyze_subterminal_regions(
                            subterminal_right, params['subterminal_threshold'])
                        
                        if min(left_percentage, right_percentage) >= params['subterminal_threshold']:
                            if left_strand is not None and right_strand is not None and left_strand != right_strand:
                                continue
                                
                            subterminal_strand = left_strand if left_strand is not None else right_strand
                            
                            te_start = tsd1_start + pattern_size + tir_size
                            te_end = tsd2_start - tir_size
                            te_seq = seq[te_start:te_end]
                            
                            # Find CDS regions
                            cds_results = find_cds_region(te_seq, search_range=params['search_dde_range'])
                            
                            for result in process_cds_results(cds_results, te_start, tsd1_start, tsd2_start, 
                                                            pattern_size, tir_size, pattern, tir1, tir2,
                                                            left_percentage, right_percentage, seq_id,
                                                            subterminal_strand):
                                yield result

def process_cds_results(cds_results, te_start, tsd1_start, tsd2_start, pattern_size,
                       tir_size, pattern, tir1, tir2, left_percentage, right_percentage,
                       seq_id, subterminal_strand):
    """
    Process CDS search results and generate transposon structure records.
    
    Args:
        cds_results: List of CDS search results
        te_start: Start position of transposable element
        tsd1_start, tsd2_start: Start positions of TSDs
        pattern_size: Size of TSD pattern
        tir_size: Size of TIR
        pattern: TSD sequence
        tir1, tir2: TIR sequences
        left_percentage, right_percentage: Subterminal region percentages
        seq_id: Sequence identifier
        subterminal_strand: Strand orientation from subterminal analysis
        
    Returns:
        Generator yielding processed transposon structure records
    """
    for cds_start, cds_end, protein_length, cds_strand, conserved_positions in cds_results:
        if subterminal_strand is not None and cds_strand != subterminal_strand:
            continue
        
        global_cds_start = te_start + cds_start
        global_cds_end = te_start + cds_end
        
        d301_pos, d367_pos, e719_pos = conserved_positions
        global_conserved_positions = (
            te_start + d301_pos if d301_pos is not None else None,
            te_start + d367_pos if d367_pos is not None else None,
            te_start + e719_pos if e719_pos is not None else None
        )
        
        distance_to_left_tir = abs(global_cds_start - (tsd1_start + pattern_size + tir_size))
        distance_to_right_tir = abs(global_cds_end - (tsd2_start - tir_size))
        
        yield {
            'seq_id': seq_id,
            'tsd1_start': tsd1_start,
            'tsd1_end': tsd1_start + pattern_size,
            'tir1_start': tsd1_start + pattern_size,
            'tir1_end': tsd1_start + pattern_size + tir_size,
            'tir2_start': tsd2_start - tir_size,
            'tir2_end': tsd2_start,
            'tsd2_start': tsd2_start,
            'tsd2_end': tsd2_start + pattern_size,
            'tsd_seq': pattern,
            'tir1_seq': tir1,
            'tir2_seq': tir2,
            'left_percentage': left_percentage,
            'right_percentage': right_percentage,
            'strand': cds_strand,
            'cds_start': global_cds_start,
            'cds_end': global_cds_end,
            'cds_start': global_cds_start,
            'cds_end': global_cds_end,
            'cds_length': cds_end - cds_start + 3,  # Including stop codon
            'protein_length': protein_length,
            'distance_to_left_tir': distance_to_left_tir,
            'distance_to_right_tir': distance_to_right_tir,
            'd301_position': global_conserved_positions[0],
            'd367_position': global_conserved_positions[1],
            'e719_position': global_conserved_positions[2]
        }