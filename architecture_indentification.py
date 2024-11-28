
from typing import List, Tuple


def hamming_distance(str1, str2):
    """Calculates the Hamming distance between two strings."""
    if len(str1)!= len(str2):
        return float('inf')  # Infinite distance if lengths don't match
    distance = 0
    for i in range(len(str1)):
        if str1[i]!= str2[i]:
            distance += 1
    return distance

def structure_verification(seq, pattern_size=8, gap_size=2500, tir_size=5, mismatch_allowed=2):
    seq = seq.upper()
    seq_length = len(seq)
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}  # Simplified complement

    result = []

    tsd_hash = {}
    for i in range(seq_length - pattern_size + 1):
        pattern = seq[i:i + pattern_size]
        if 'N' in pattern:  # Still exclude patterns with N
            continue
        tsd_hash.setdefault(pattern, []).append(i)

    for pattern, positions in tsd_hash.items():
        for i in range(len(positions)):
            for j in range(i + 1, len(positions)):
                tsd_start_1 = positions[i]
                tsd_start_2 = positions[j]

                if tsd_start_2 - (tsd_start_1 + pattern_size) >= gap_size:
                    tir_1_start = tsd_start_1 + pattern_size
                    tir_2_start = tsd_start_2 - tir_size

                    if tir_2_start < 0 or tir_1_start + tir_size > seq_length:
                        continue

                    tir_1 = seq[tir_1_start:tir_1_start + tir_size]
                    tir_2 = seq[tir_2_start:tir_2_start + tir_size]

                    rev_complement_tir_2 = ''.join(complement.get(base, 'N') for base in reversed(tir_2))

                    if hamming_distance(tir_1, rev_complement_tir_2) <= mismatch_allowed:
                        result.append((tsd_start_1, tsd_start_1 + pattern_size,
                                      tir_1_start, tir_1_start + tir_size,
                                      tsd_start_2, tsd_start_2 + pattern_size,
                                      tir_2_start, tir_2_start + tir_size))

    return result

def analyze_structures(seq: str, protein_regions: List[Tuple[int, int]], pattern_size=8, gap_size=1500, tir_size=5, mismatch_allowed=2):
    matches = structure_verification(seq, pattern_size, gap_size, tir_size, mismatch_allowed)
    seq_length = len(seq)

    analyzed_structures = []

    for match in matches:
        tsd1_start, tsd1_end, tir1_start, tir1_end, tsd2_start, tsd2_end, tir2_start, tir2_end = match

        # Check for overlap with any protein region
        if any((tsd1_start < p_end and tsd1_end > p_start) or
               (tir1_start < p_end and tir1_end > p_start) or
               (tsd2_start < p_end and tsd2_end > p_start) or
               (tir2_start < p_end and tir2_end > p_start)
               for p_start, p_end in protein_regions):
            continue

        scores = score_structure_match(match, protein_regions, seq_length)

        structure_info = {
            'TSD1': (tsd1_start, tsd1_end),
            'TIR1': (tir1_start, tir1_end),
            'TIR2': (tir2_start, tir2_end),
            'TSD2': (tsd2_start, tsd2_end),
            'scores': scores,
            'tsd1_seq': seq[tsd1_start:tsd1_end],
            'tir1_seq': seq[tir1_start:tir1_end],
            'tir2_seq': seq[tir2_start:tir2_end],
            'tsd2_seq': seq[tsd2_start:tsd2_end]
        }

        analyzed_structures.append(structure_info)

    return sorted(analyzed_structures, key=lambda x: x['scores']['total_score'], reverse=True)

def score_structure_match(match_positions, protein_regions, seq_length):
    """
    Calculates a score for a genomic structure match based on several features.

    Args:
        match_positions (tuple): Tuple containing start and end positions of TSDs and TIRs.
        protein_regions (list): List of tuples, each containing start and end positions of a protein region.
        seq_length (int): Length of the full sequence.

    Returns:
        dict: A dictionary containing individual scores and the total score.
    """
    tsd1_start, tsd1_end, tir1_start, tir1_end, tsd2_start, tsd2_end, tir2_start, tir2_end = match_positions
    structure_start = min(tsd1_start, tir1_start)
    structure_end = max(tsd2_end, tir2_end)
    structure_length = structure_end - structure_start

    # 1. TSD-TIR Integrity Score (0-100) - Simplified, no similarity check
    tsd_length = (tsd1_end - tsd1_start) + (tsd2_end - tsd2_start)
    tir_length = (tir1_end - tir1_start) + (tir2_end - tir2_start)
    integrity_score = min(100, (tsd_length + tir_length) * 5)  # Increased weight back to 5 per base

    # 2. Size Appropriateness Score (0-100)
    if 2000 <= structure_length <= 20000:
        size_score = 100
    elif structure_length < 1000:
        size_score = structure_length / 10  # Linear scaling for small structures
    else:
        size_score = max(0, 100 - ((structure_length - 20000) / 1000) ** 2)  # Quadratic decrease for large structures

    # 3. Protein Proximity Score (0-100)
    protein_proximity_score = calculate_protein_proximity_score(
        tir1_start, tir1_end, tir2_start, tir2_end,
        tsd1_start, tsd1_end, tsd2_start, tsd2_end,
        protein_regions, structure_length
    )

    # 4. Symmetry Score (0-100)
    left_spacing = tir1_start - tsd1_end
    right_spacing = tsd2_start - tir2_end
    relative_spacing_diff = abs(left_spacing - right_spacing) / structure_length if structure_length > 0 else 0
    symmetry_score = max(0, 100 - (relative_spacing_diff * 500))  # Non-linear decrease based on relative spacing difference

    # Calculate total score (0-400)
    total_score = integrity_score + size_score + protein_proximity_score + symmetry_score

    return {
        'integrity_score': integrity_score,
        'size_score': size_score,
        'protein_proximity_score': protein_proximity_score,
        'symmetry_score': symmetry_score,
        'total_score': total_score,
        'structure_length': structure_length
    }

def calculate_protein_proximity_score(tir1_start, tir1_end, tir2_start, tir2_end,
                                      tsd1_start, tsd1_end, tsd2_start, tsd2_end,
                                      protein_regions, structure_length):
    """Calculates the protein proximity score based on distances between TIR/TSD and protein regions."""
    if not protein_regions:
        return 0

    tir_tsd_regions = [(tir1_start, tir1_end), (tir2_start, tir2_end),
                       (tsd1_start, tsd1_end), (tsd2_start, tsd2_end)]
    min_distances = []

    for tir_tsd_start, tir_tsd_end in tir_tsd_regions:
        tir_tsd_center = (tir_tsd_start + tir_tsd_end) / 2
        closest_distance = float('inf')
        for p_start, p_end in protein_regions:
            p_center = (p_start + p_end) / 2
            distance = abs(tir_tsd_center - p_center)
            closest_distance = min(closest_distance, distance)
        min_distances.append(closest_distance)

    avg_min_distance = sum(min_distances) / len(min_distances) if min_distances else structure_length
    proximity_score = max(0, 100 - (avg_min_distance / structure_length) * 200)
    return proximity_score