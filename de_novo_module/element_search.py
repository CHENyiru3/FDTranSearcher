import logging
import os
import psutil
from typing import Tuple, List, Dict, Any
from structure_verification import StructureVerification

def process_chunk(args: Tuple[str, str, str, int, Dict[str, Any], List[Dict[str, Any]]]) -> Tuple[str, List[Dict], float]:
    """
    Process a single sequence chunk for transposon analysis.
    
    Args:
        args: Tuple containing (chunk_id, sequence, seq_id, base_position, params, reserve_sites)
        
    Returns:
        Tuple of (chunk_id, list of matches, current_memory_usage)
    """
    chunk_id, sequence, seq_id, base_pos, params, reserve_sites = args
    try:
        # Convert reserve_sites back to string format if they are dictionaries
        if isinstance(reserve_sites[0], dict):
            reserve_sites = [f"{site['amino_acid']}{site['position']}" for site in reserve_sites]
            
        verifier = StructureVerification(params, reserve_sites)  # params now includes min_cds_distance
        
        current_process = psutil.Process(os.getpid())
        memory_usage = current_process.memory_info().rss / 1024 / 1024  # Convert to MB
        logging.info(f"Processing chunk {chunk_id} of {seq_id} (Memory usage: {memory_usage:.2f} MB)")
        
        matches = list(verifier.verify_structures(sequence, seq_id))
        
        for match in matches:
            for key in ['tsd1_start', 'tsd1_end', 'tir1_start', 'tir1_end',
                       'tir2_start', 'tir2_end', 'tsd2_start', 'tsd2_end',
                       'cds_start', 'cds_end']:
                if key in match and match[key] is not None:
                    match[key] += base_pos
        
        return chunk_id, matches, memory_usage
    except Exception as e:
        logging.error(f"Error in chunk {chunk_id}: {str(e)}")
        return chunk_id, [], 0.0