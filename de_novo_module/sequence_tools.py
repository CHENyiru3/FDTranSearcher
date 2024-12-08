from Bio import SeqIO
from Bio.Seq import Seq
from typing import List, Dict, Generator, Tuple
import logging
import re

def extract_sequence(input_file: str) -> Generator[Tuple[str, str], None, None]:
    """
    Read sequences from FASTA file.

    Args:
        input_file: Path to input FASTA file
            
    Yields:
        Tuple of (sequence_id, sequence)
    """
    logging.info(f"Reading sequence from {input_file}")
    total_length = 0
    try:
        for record in SeqIO.parse(input_file, "fasta"):
            seq_id = record.id
            seq = str(record.seq).upper()
            logging.info(f"Found sequence: {seq_id}, length: {len(seq)}")
            yield (seq_id, seq)
            total_length += len(seq)
        logging.info(f"Total sequence length: {total_length} bp")
    except FileNotFoundError as e:
        logging.error(f"File not found: {e}")
        raise
    except Exception as e:
        logging.error(f"Error reading FASTA file: {e}")
        raise

def parse_reserve_sites(reserve_sites: List[str]) -> List[Dict[str, str]]:
    """
    Parse user-specified reserve sites.
    
    Args:
        reserve_sites: List of reserve site strings (e.g., ['D301', 'E719'])
    
    Returns:
        List of dictionaries with 'amino_acid' and 'position'
    """
    parsed_sites = []
    pattern = re.compile(r'^([A-Z])(\d+)$')
    for site in reserve_sites:
        match = pattern.fullmatch(site)
        if not match:
            raise ValueError(f"Invalid reserve site format: '{site}'. Expected format like 'D301'.")
        amino_acid, position = match.groups()
        parsed_sites.append({
            'amino_acid': amino_acid,
            'position': int(position)
        })
    return parsed_sites

def get_reverse_complement(seq: str) -> str:
    """
    Get the reverse complement of a DNA sequence.

    Args:
        seq: Input DNA sequence
    
    Returns:
        The reverse complement of the input sequence
    """
    return str(Seq(seq).reverse_complement())

def get_codon_set(amino_acid: str) -> set:
    """
    Get the set of codons corresponding to a given amino acid.
    
    Args:
        amino_acid: Single-letter amino acid code
    
    Returns:
        Set of codons encoding the amino acid
    """
    codon_table = {
        'A': {'GCT', 'GCC', 'GCA', 'GCG'},
        'C': {'TGT', 'TGC'},
        'D': {'GAT', 'GAC'},
        'E': {'GAA', 'GAG'},
        'F': {'TTT', 'TTC'},
        'G': {'GGT', 'GGC', 'GGA', 'GGG'},
        'H': {'CAT', 'CAC'},
        'I': {'ATT', 'ATC', 'ATA'},
        'K': {'AAA', 'AAG'},
        'L': {'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'},
        'M': {'ATG'},
        'N': {'AAT', 'AAC'},
        'P': {'CCT', 'CCC', 'CCA', 'CCG'},
        'Q': {'CAA', 'CAG'},
        'R': {'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'},
        'S': {'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'},
        'T': {'ACT', 'ACC', 'ACA', 'ACG'},
        'V': {'GTT', 'GTC', 'GTA', 'GTG'},
        'W': {'TGG'},
        'Y': {'TAT', 'TAC'},
    }
    return codon_table.get(amino_acid.upper(), set())

def kmp_search(text: str, pattern: str) -> list:
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

def build_kmp_table(pattern: str) -> list:
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