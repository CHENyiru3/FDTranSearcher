#!/usr/bin/env python3
from typing import List

def build_kmp_table(pattern: str) -> List[int]:
    """
    Build Knuth-Morris-Pratt partial match table for pattern matching.
    
    Args:
        pattern: The pattern string to build the table for
    
    Returns:
        A list representing the partial match table for KMP algorithm
    """
    table = [0] * len(pattern)
    j = 0
    for i in range(1, len(pattern)):
        while j > 0 and pattern[i] != pattern[j]:
            j = table[j-1]
        if pattern[i] == pattern[j]:
            j += 1
        table[i] = j
    return table

def kmp_search(text: str, pattern: str) -> List[int]:
    """
    Search for all occurrences of pattern in text using KMP algorithm.
    
    Args:
        text: The text to search in
        pattern: The pattern to search for
    
    Returns:
        List of starting positions where pattern occurs in text
    """
    if not pattern or not text:
        return []
    
    table = build_kmp_table(pattern)
    positions = []
    j = 0
    
    for i in range(len(text)):
        while j > 0 and text[i] != pattern[j]:
            j = table[j-1]
        if text[i] == pattern[j]:
            j += 1
        if j == len(pattern):
            positions.append(i - j + 1)
            j = table[j-1]
    
    return positions

def get_reverse_complement(seq: str) -> str:
    """
    Get the reverse complement of a DNA sequence.
    
    Args:
        seq: Input DNA sequence
        
    Returns:
        The reverse complement of the input sequence
    """
    complement = str.maketrans('ATCGN', 'TAGCN')
    return seq.translate(complement)[::-1]

def analyze_subterminal_regions(seq: str, threshold: float):
    """
    Analyze subterminal regions for AAAGGG/CCCTTT motif enrichment.
    
    Args:
        seq: Sequence to analyze
        threshold: Minimum percentage threshold for enrichment
        
    Returns:
        Tuple of (enrichment_percentage, strand_orientation)
        where strand_orientation is '+', '-', or None
    """
    import re
    forward_pattern = re.compile('AAAGGG|GGGAAA')
    reverse_pattern = re.compile('CCCTTT|TTTCCC')

    forward_matches = len(list(forward_pattern.finditer(seq)))
    reverse_matches = len(list(reverse_pattern.finditer(seq)))

    forward_percentage = (forward_matches * 6 / len(seq)) * 100 if seq else 0.0
    reverse_percentage = (reverse_matches * 6 / len(seq)) * 100 if seq else 0.0

    if forward_percentage >= threshold and reverse_percentage < threshold:
        return (forward_percentage, "+")
    elif reverse_percentage >= threshold and forward_percentage < threshold:
        return (reverse_percentage, "-")
    else:
        return (max(forward_percentage, reverse_percentage), None)