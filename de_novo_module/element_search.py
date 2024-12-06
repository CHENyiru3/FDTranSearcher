#!/usr/bin/env python3
from typing import Dict, List, Tuple, Optional
from collections import defaultdict
from sequence_tools import get_reverse_complement


def build_pattern_hash(seq: str, pattern_size: int) -> Dict[str, List[int]]:
    """
    Build a hash table for patterns of specified size from sequence.
    
    Args:
        seq: Input DNA sequence
        pattern_size: Size of patterns to hash
        
    Returns:
        Dictionary mapping patterns to lists of their starting positions
    """
    pattern_hash = defaultdict(list)
    for i in range(len(seq) - pattern_size + 1):
        pattern = seq[i:i + pattern_size]
        if 'N' not in pattern:
            pattern_hash[pattern].append(i)
    return pattern_hash

def find_tir_near_tsd(seq: str, tsd_start: int, tsd2_start: int, pattern_size: int,
                      min_tir_size: int, max_tir_size: int, max_mismatch: int) -> Optional[Tuple[str, str, int]]:
    """
    Search for Terminal Inverted Repeats (TIR) near Target Site Duplications (TSD).
    
    Args:
        seq: Input DNA sequence
        tsd_start: Start position of first TSD
        tsd2_start: Start position of second TSD
        pattern_size: Size of TSD pattern
        min_tir_size: Minimum size of TIR
        max_tir_size: Maximum size of TIR
        max_mismatch: Maximum allowed mismatches in TIR
        
    Returns:
        Tuple of (left_tir, right_tir, tir_size) if found, None otherwise
    """
    best_tir = None
    best_mismatch_count = max_mismatch + 1
    
    for tir_size in range(min_tir_size, max_tir_size + 1):
        left_tir = seq[tsd_start + pattern_size:tsd_start + pattern_size + tir_size]
        right_tir = seq[tsd2_start - tir_size:tsd2_start]
        
        if len(left_tir) != tir_size or len(right_tir) != tir_size:
            continue
            
        right_tir_rc = get_reverse_complement(right_tir)
        mismatches = sum(a != b for a, b in zip(left_tir, right_tir_rc))
        
        if mismatches <= max_mismatch and mismatches < best_mismatch_count:
            best_mismatch_count = mismatches
            best_tir = (left_tir, right_tir, tir_size)
    
    return best_tir

def find_conserved_codon_positions(seq: str, strand: str, search_range: int) -> Tuple[Optional[int], Optional[int], Optional[int]]:
    """
    Search for conserved DDE codons in a protein-coding sequence.
    
    Args:
        seq: Input DNA sequence
        strand: Strand orientation ('+' or '-')
        search_range: Range around expected positions to search
        
    Returns:
        Tuple of (D301_position, D367_position, E719_position),
        where positions are relative to sequence start or None if not found
    """
    d_codons = {'GAT', 'GAC'}
    e_codons = {'GAA', 'GAG'}
    
    if strand == '-':
        seq = get_reverse_complement(seq)
    
    target_positions = [
        (301 * 3 - 3, d_codons),  # D301
        (367 * 3 - 3, d_codons),  # D367
        (719 * 3 - 3, e_codons)   # E719
    ]
    
    results = []
    for target_pos, target_codons in target_positions:
        found_pos = None
        search_start = max(0, target_pos - search_range)
        search_end = min(len(seq), target_pos + search_range)
        
        for pos in range(search_start, search_end, 3):
            if pos + 3 <= len(seq):
                codon = seq[pos:pos+3]
                if codon in target_codons:
                    found_pos = pos
                    break
        
        results.append(found_pos)
    
    return tuple(results)