from typing import Dict, List, Optional, Generator, Tuple
import re
from collections import defaultdict
from sequence_tools import (
    get_reverse_complement, 
    get_codon_set,
    kmp_search
)

class StructureVerification:
    """
    Class for verifying transposon structures in DNA sequences.
    Handles the identification and verification of TSD, TIR, and other features.
    """
    def __init__(self, params: Dict, reserve_sites: List[str]):
        self.params = params
        self.reserve_sites = self._parse_reserve_sites(reserve_sites)

    def _parse_reserve_sites(self, reserve_sites: List[str]) -> List[Dict[str, str]]:
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

    def verify_structures(self, seq: str, seq_id: str) -> Generator[dict, None, None]:
        """Structure Validation Master Function"""
        seq = seq.upper()
    
        pattern_hashes = {
            size: self._build_pattern_hash(seq, size)
            for size in range(self.params['min_tsd_pattern_size'], 
                    self.params['max_tsd_pattern_size'] + 1)
        }
    
        motif = self.params['motif']
        subterminal_length = self.params['subterminal_length']
    
        for pattern_size, pattern_hash in pattern_hashes.items():
            for pattern, positions in pattern_hash.items():
                for i, tsd1_start in enumerate(positions):
                    for tsd2_start in positions[i + 1:]:
                        gap_size = self.params['gap_size']
                        if not (gap_size <= tsd2_start - (tsd1_start + pattern_size) <= gap_size * 2):
                            continue
                    
                        tir_result = self._find_tir_near_tsd(
                            seq, tsd1_start, tsd2_start, pattern_size
                        )
                    
                        if tir_result:
                            tir1, tir2, tir_size = tir_result
                        
                            subterminal_left = seq[tsd1_start + pattern_size + tir_size:
                                                tsd1_start + pattern_size + tir_size + subterminal_length]
                            subterminal_right = seq[tsd2_start - tir_size - subterminal_length:tsd2_start - tir_size]
                        
                            enrichment_left, strand_left, motif_scores_left = self._analyze_subterminal_regions(subterminal_left, motif)
                            enrichment_right, strand_right, motif_scores_right = self._analyze_subterminal_regions(subterminal_right, motif)
                        
                            if (enrichment_left < self.params['subterminal_threshold'] or
                                enrichment_right < self.params['subterminal_threshold']):
                                continue
                        
                            if strand_left and strand_right and strand_left == strand_right:
                                subterminal_strand = strand_left
                            elif strand_left and not strand_right:
                                subterminal_strand = strand_left
                            elif strand_right and not strand_left:
                                subterminal_strand = strand_right
                            else:
                                subterminal_strand = None
                        
                            te_start = tsd1_start + pattern_size + tir_size
                            te_end = tsd2_start - tir_size
                            te_seq = seq[te_start:te_end]
                        
                            cds_results = self._find_cds_region(te_seq)
                        
                            for cds_start, cds_end, protein_length, cds_strand, conserved_positions in cds_results:
                                if subterminal_strand is not None and cds_strand != subterminal_strand:
                                    continue
                            
                                global_cds_start = te_start + cds_start
                                global_cds_end = te_start + cds_end
                            
                                # Calculate distance to reserve sites relative to target positions
                                distance_to_reserve = {}
                                all_within_range = True
                                for idx, site in enumerate(self.reserve_sites):
                                    target_pos = (site['position'] * 3) - 3
                                    conserved_pos = self._find_conserved_codon_positions(
                                        seq=seq,
                                        strand=cds_strand,
                                        target_pos=global_cds_start + target_pos,
                                        conserve_site_range=self.params['conserve_site_range'],
                                        amino_acid=site['amino_acid']
                                    )
                                
                                    if conserved_pos is not None:
                                        distance = abs(conserved_pos - (global_cds_start + target_pos))
                                        if distance > self.params['conserve_site_range']:
                                            all_within_range = False
                                            break
                                        distance_to_reserve[f'reserve_site_{idx+1}_distance'] = distance
                                    else:
                                        all_within_range = False
                                        break
                            
                                if not all_within_range:
                                    continue
                            
                                pattern_match = self._calculate_pattern_match(
                                    pattern,
                                    self.params['motif'],
                                    list(distance_to_reserve.values()),
                                    enrichment_left,
                                    enrichment_right
                                )
                            
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
                                    'left_percentage': enrichment_left,
                                    'right_percentage': enrichment_right,
                                    'strand': cds_strand,
                                    'cds_start': global_cds_start,
                                    'cds_end': global_cds_end,
                                    'cds_length': cds_end - cds_start,
                                    'protein_length': protein_length,
                                    'distance_to_left_tir': abs(global_cds_start - (tsd1_start + pattern_size + tir_size)),
                                    'distance_to_right_tir': abs(global_cds_end - (tsd2_start - tir_size)),
                                    'pattern_match': pattern_match,
                                    **distance_to_reserve
                                }

    def _build_pattern_hash(self, seq: str, pattern_size: int) -> Dict[str, List[int]]:
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

    def _find_tir_near_tsd(self, seq: str, tsd_start: int, tsd2_start: int, pattern_size: int) -> Optional[Tuple[str, str, int]]:
        """
        Search for Terminal Inverted Repeats (TIR) near Target Site Duplications (TSD).

        Args:
            seq: Input DNA sequence
            tsd_start: Start position of first TSD
            tsd2_start: Start position of second TSD
            pattern_size: Size of TSD pattern

        Returns:
            Tuple of (left_tir, right_tir, tir_size) if found, None otherwise
        """
        best_tir = None
        best_mismatch_count = self.params['max_tir_mismatch'] + 1
        
        for tir_size in range(self.params['min_tir_size'], self.params['max_tir_size'] + 1):
            left_tir_start = tsd_start + pattern_size
            left_tir_end = left_tir_start + tir_size
            right_tir_start = tsd2_start - tir_size
            right_tir_end = tsd2_start
            
            if left_tir_end > len(seq) or right_tir_start < 0:
                continue
            
            left_tir = seq[left_tir_start:left_tir_end]
            right_tir = seq[right_tir_start:right_tir_end]
            
            right_tir_rc = get_reverse_complement(right_tir)
            mismatches = sum(a != b for a, b in zip(left_tir, right_tir_rc))
            
            if mismatches <= self.params['max_tir_mismatch'] and mismatches < best_mismatch_count:
                best_mismatch_count = mismatches
                best_tir = (left_tir, right_tir, tir_size)
        
        return best_tir

    def _analyze_subterminal_regions(self, seq: str, motif: str) -> Tuple[float, Optional[str], Dict[str, float]]:
        """
        Analyze subterminal regions for user-specified motif enrichment.

        Args:
            seq: Sequence to analyze
            motif: User-specified motif

        Returns:
            Tuple of (enrichment_percentage, strand_orientation, motif_match_scores)
            where strand_orientation is '+', '-', or None
        """
        total_seq_length = len(seq)

        if total_seq_length == 0:
            return 0.0, None, {}

        forward_pattern = re.compile(motif)
        reverse_pattern = re.compile(get_reverse_complement(motif))

        forward_matches = len(list(forward_pattern.finditer(seq)))
        reverse_matches = len(list(reverse_pattern.finditer(seq)))

        forward_percentage = (forward_matches * len(motif) / total_seq_length) * 100
        reverse_percentage = (reverse_matches * len(motif) / total_seq_length) * 100

        if forward_percentage >= self.params['subterminal_threshold'] and reverse_percentage < self.params['subterminal_threshold']:
            strand = '+'
            enrichment = forward_percentage
        elif reverse_percentage >= self.params['subterminal_threshold'] and forward_percentage < self.params['subterminal_threshold']:
            strand = '-'
            enrichment = reverse_percentage
        else:
            strand = None
            enrichment = max(forward_percentage, reverse_percentage)

        motif_match_scores = {motif: max(forward_percentage, reverse_percentage)}

        return enrichment, strand, motif_match_scores

    def _find_cds_region(self, seq: str, window_size: int = 5000) -> List[Tuple[int, int, int, str, Tuple[Optional[int], ...]]]:
        """
        Search for protein-coding regions (CDS) with conserved motif.
    
        Args:
            seq: Input DNA sequence
            window_size: Maximum size of window to search for stop codons
    
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
                    if cds_length >= self.params['min_cds_distance']:
                        conserved_positions = self._find_conserved_codon_positions_dynamic(seq[start_pos:actual_stop_pos+3], '+')
                        if all(pos is not None for pos in conserved_positions):
                            protein_length = cds_length // 3
                            abs_conserved_positions = tuple(pos + start_pos if pos is not None else None 
                                                     for pos in conserved_positions)
                            results.append((start_pos, actual_stop_pos + 3, protein_length, "+", abs_conserved_positions))

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
                    if cds_length >= self.params['min_cds_distance']:
                        conserved_positions = self._find_conserved_codon_positions_dynamic(rc_seq[start_pos:actual_stop_pos+3], '-')
                        if all(pos is not None for pos in conserved_positions):
                            protein_length = cds_length // 3
                            orig_start = len(seq) - (actual_stop_pos + 3)
                            orig_end = len(seq) - start_pos
                            abs_conserved_positions = tuple(len(seq) - (pos + start_pos) if pos is not None else None 
                                                     for pos in conserved_positions[::-1])
                            results.append((orig_start, orig_end, protein_length, "-", abs_conserved_positions))
        
        return results

    def _find_conserved_codon_positions(self, seq: str, strand: str, target_pos: int, conserve_site_range: int, amino_acid: str) -> Optional[int]:
        """
        Search for a conserved codon corresponding to the specified amino acid within the target range.
    
        Args:
            seq: Input DNA sequence
            strand: Strand orientation ('+' or '-')
            target_pos: Target position in the sequence
            conserve_site_range: Range around the target position to search
            amino_acid: Single-letter amino acid code
    
        Returns:
            Position of the conserved codon if found, else None
        """
        codon_set = get_codon_set(amino_acid)
    
        if strand == '-':
            seq = get_reverse_complement(seq)
    
        search_start = max(0, target_pos - conserve_site_range)
        search_end = min(len(seq), target_pos + conserve_site_range)
    
        for p in range(search_start, search_end - 2, 3):
            codon = seq[p:p+3]
            if codon in codon_set:
                return p
    
        return None

    def _find_conserved_codon_positions_dynamic(self, seq: str, strand: str) -> Tuple[Optional[int], ...]:
        """
        Find positions of multiple conserved codons in a sequence.

        Args:
            seq: Input DNA sequence
            strand: Strand orientation ('+' or '-')

        Returns:
            Tuple of positions for each conserved codon (None if not found)
        """
        positions = []
        for site in self.reserve_sites:
            target_pos = (site['position'] * 3) - 3
            pos = self._find_conserved_codon_positions(
                seq=seq,
                strand=strand,
                target_pos=target_pos,
                conserve_site_range=self.params['conserve_site_range'],
                amino_acid=site['amino_acid']
            )
            positions.append(pos)
        return tuple(positions)

    def _calculate_pattern_match(self, tsd_pattern: str, motif: str, distances: List[Optional[int]], left_percentage: float, right_percentage: float) -> float:
        """
        Calculate a pattern_match score based on how well the TSD pattern matches the user-defined motif.

        Args:
            tsd_pattern: The TSD sequence pattern
            motif: User-defined motif
            distances: List of distances to conserved sites
            left_percentage: Percentage of left subterminal enrichment
            right_percentage: Percentage of right subterminal enrichment

        Returns:
            A float score representing the match quality
        """
        match_score = 1.0 if motif in tsd_pattern else 0.0

        final_score = self._calculate_pattern_match_score(
            subterminal_matches=match_score,
            distances=distances,
            left_percentage=left_percentage,
            right_percentage=right_percentage
        )

        return final_score

    def _calculate_pattern_match_score(self, subterminal_matches: float, distances: List[Optional[int]], 
                                   left_percentage: float, right_percentage: float) -> float:
        """
        Calculate pattern matching scores considering thresholds.
    
        Args:
            subterminal_matches: Indicator of motif presence in TSD pattern (1.0 if present, else 0.0)
            distances: List of distances to conserved reserve sites
            left_percentage: Percentage of motif enrichment in the left subterminal region
            right_percentage: Percentage of motif enrichment in the right subterminal region
    
        Returns:
            A float score representing the match quality, ranging from 0 to 100
        """
        # Enrichment score (0-50 points)
        enrichment_score = 0.5 * ((left_percentage + right_percentage) / 2)

        # Positional Distance Score (0-50 points)
        if not distances or not any(d is not None for d in distances):
            distance_score = 0
        else:
            valid_distances = [d for d in distances if d is not None]
            max_distance = self.params['conserve_site_range']
            avg_distance = sum(valid_distances) / len(valid_distances)
            distance_score = 50 * (1 - avg_distance / max_distance)

        final_score = int(enrichment_score + distance_score)
        return final_score