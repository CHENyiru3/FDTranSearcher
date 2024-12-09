from multiprocessing import Pool
from Bio import motifs
from Bio.Seq import Seq
from typing import List, Tuple

class CompressedTrieNode:
    def __init__(self, label=None):
        self.label = label
        self.children = {}
        self.positions = []

class CompressedTrie:
    def __init__(self):
        self.root = CompressedTrieNode()

    def insert(self, pattern, position):
        node = self.root
        while pattern:
            if not node.children:
                node.children[pattern[0]] = CompressedTrieNode(pattern)
                node.children[pattern[0]].positions.append(position)
                return

            if pattern[0] not in node.children:
                node.children[pattern[0]] = CompressedTrieNode(pattern)
                node.children[pattern[0]].positions.append(position)
                return

            child = node.children[pattern[0]]

            common_prefix = self._common_prefix(child.label, pattern)

            if len(common_prefix) == len(pattern):
                child.positions.append(position)
                return

            elif len(common_prefix) == len(child.label):
                pattern = pattern[len(common_prefix):]
                node = child
                continue

            elif len(common_prefix) < len(child.label):
                new_node = CompressedTrieNode(child.label[len(common_prefix):])
                new_node.children = child.children
                new_node.positions = child.positions

                child.label = common_prefix
                child.children = {new_node.label[0]: new_node}
                child.positions = []

                if len(common_prefix) < len(pattern):
                    pattern = pattern[len(common_prefix):]
                    if pattern[0] not in child.children:
                         child.children[pattern[0]] = CompressedTrieNode(pattern)
                    child = child.children[pattern[0]]
                    child.positions.append(position)
                    return
                else:
                    child.positions.append(position)
                    return

            elif len(common_prefix) < len(pattern):
                pattern = pattern[len(common_prefix):]
                if pattern[0] not in child.children:
                    child.children[pattern[0]] = CompressedTrieNode(pattern)
                child = child.children[pattern[0]]
                child.positions.append(position)
                return
            else:
                child.positions.append(position)
                return

    def search(self, pattern):
        node = self.root
        current_pattern = pattern
        while current_pattern:
            if not node.children:
                return []

            if current_pattern[0] not in node.children:
                return []

            child = node.children[current_pattern[0]]

            if current_pattern.startswith(child.label):
                current_pattern = current_pattern[len(child.label):]
                if not current_pattern:
                    return child.positions
                node = child
            else:
                return []

        return node.positions if not current_pattern else []

    def _common_prefix(self, str1, str2):
        i = 0
        if str1 is None or str2 is None:
            return ""
        while i < len(str1) and i < len(str2) and str1[i] == str2[i]:
            i += 1
        return str1[:i]

    def get_all_patterns(self):
        patterns = []
        self._dfs(self.root, "", patterns)
        return patterns

    def _dfs(self, node, current_pattern, patterns):
        if node.label:
            current_pattern += node.label
        if node.positions:
            patterns.append((current_pattern, node.positions))
        for child in node.children.values():
            self._dfs(child, current_pattern, patterns)

    def print_trie(self):
        self._print_trie_helper(self.root, "")

    def _print_trie_helper(self, node, prefix):
        if node.label:
            print(f"{prefix}{node.label} -> {node.positions}")
        for key, child in node.children.items():
            self._print_trie_helper(child, prefix + (node.label if node.label else "") + "-")

def hamming_distance(str1, str2):
    """Calculates the Hamming distance between two strings."""
    if len(str1) != len(str2):
        return float('inf')  # Infinite distance if lengths don't match
    distance = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            distance += 1
    return distance

def structure_verification(seq, pattern_size=8, gap_size=2500, tir_size=5, mismatch_allowed=2):
    seq = seq.upper()
    seq_length = len(seq)
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

    result = []
    tsd_trie = CompressedTrie()

    for i in range(seq_length - pattern_size + 1):
        pattern = seq[i:i + pattern_size]
        if 'N' not in pattern:
            tsd_trie.insert(pattern, i)

    for tsd_start_1 in range(seq_length - pattern_size + 1):
        pattern = seq[tsd_start_1:tsd_start_1 + pattern_size]
        if 'N' not in pattern:
            positions = tsd_trie.search(pattern)
            for tsd_start_2 in positions:
                if tsd_start_2 > tsd_start_1 and tsd_start_2 - (tsd_start_1 + pattern_size) >= gap_size:
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

# Alternative KMP implementation; but not used in the final version
# def structure_verification(seq, pattern_size=8, gap_size=2500, tir_size=5, mismatch_allowed=2):
#     def build_kmp_table(pattern):
#         m = len(pattern)
#         fail = [0] * m
#         j = 0
#         for i in range(1, m):
#             while j > 0 and pattern[j] != pattern[i]:
#                 j = fail[j-1]
#             if pattern[j] == pattern[i]:
#                 j += 1
#             fail[i] = j
#         return fail

#     def kmp_search(text, pattern):
#         if not pattern:
#             return []
#         fail = build_kmp_table(pattern)
#         positions = []
#         j = 0
#         for i in range(len(text)):
#             while j > 0 and pattern[j] != text[i]:
#                 j = fail[j-1]
#             if pattern[j] == text[i]:
#                 j += 1
#             if j == len(pattern):
#                 positions.append(i - j + 1)
#                 j = fail[j-1]
#         return positions

#     def verify_tir(tir1, tir2, mismatch_allowed):
#         complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
#         rev_complement_tir2 = ''.join(complement.get(base, 'N') 
#                                     for base in reversed(tir2))
        
#         mismatches = 0
#         for a, b in zip(tir1, rev_complement_tir2):
#             if a != b:
#                 mismatches += 1
#                 if mismatches > mismatch_allowed:
#                     return False
#         return True

#     seq = seq.upper()
#     seq_length = len(seq)
#     results = []
    
#     potential_positions = [False] * seq_length
    
#     for i in range(seq_length - pattern_size + 1): 
#         pattern = seq[i:i + pattern_size]
#         if 'N' in pattern:
#             continue
            
#         matches = kmp_search(seq[i+pattern_size:], pattern)
#         for pos in matches:
#             actual_pos = i + pattern_size + pos
#             if actual_pos - i >= gap_size:
#                 potential_positions[i] = True
#                 potential_positions[actual_pos] = True
    
#     for i in range(seq_length - pattern_size + 1):
#         if not potential_positions[i]:
#             continue
            
#         pattern = seq[i:i + pattern_size]
#         if 'N' in pattern:
#             continue
            
#         matches = kmp_search(seq[i+pattern_size:], pattern)
#         for pos in matches:
#             j = i + pattern_size + pos
#             if j - (i + pattern_size) < gap_size:
#                 continue
                
#             tir1_start = i + pattern_size
#             tir2_start = j - tir_size
            
#             if tir2_start < 0 or tir1_start + tir_size > seq_length:
#                 continue
                
#             tir1 = seq[tir1_start:tir1_start + tir_size]
#             tir2 = seq[tir2_start:tir2_start + tir_size]
            
#             if verify_tir(tir1, tir2, mismatch_allowed):
#                 results.append((
#                     i, i + pattern_size,
#                     tir1_start, tir1_start + tir_size,
#                     j, j + pattern_size,
#                     tir2_start, tir2_start + tir_size
#                 ))
    
#     return sorted(results)

def score_structure_match(match_positions, protein_regions, mini_size = 2000, max_size = 15000):
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
    if mini_size <= structure_length <= max_size:
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

def analyze_structures(seq: str, protein_regions: List[Tuple[int, int]], pattern_size=8, gap_size=1500, tir_size=5, mismatch_allowed=1, mini_size = 2000, max_size = 15000):
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

        scores = score_structure_match(match, protein_regions, seq_length,mini_size,max_size)

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