#!/usr/bin/env python3
"""
TED Family Consensus - Domain boundary consensus from multiple prediction methods

This script takes a Pfam SEED alignment and fetches domain predictions from the
three TED methods (chainsaw, merizo, unidoc-ndr) for each sequence. It then
calculates family-wide consensus domain boundaries using multiple scoring strategies.

Usage:
    python ted_family_consensus.py SEED [options]

Author: Claude Code for Pfam Curation
"""

import argparse
import sys
import re
import requests
import time
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Set
import json


# ============================================================================
# Data Classes
# ============================================================================

@dataclass
class Segment:
    """A single continuous segment of a domain."""
    start: int
    end: int

    def __repr__(self):
        return f"{self.start}-{self.end}"

    def length(self) -> int:
        return self.end - self.start + 1

    def overlaps(self, other: 'Segment', min_overlap: int = 1) -> bool:
        """Check if this segment overlaps with another."""
        overlap_start = max(self.start, other.start)
        overlap_end = min(self.end, other.end)
        return (overlap_end - overlap_start + 1) >= min_overlap


@dataclass
class Domain:
    """A domain prediction, potentially with multiple segments."""
    segments: List[Segment]
    method: str
    score: float = 0.0
    domain_idx: int = 0  # Which domain in the protein (0-indexed)

    @property
    def start(self) -> int:
        """First residue of the domain (across all segments)."""
        return min(s.start for s in self.segments)

    @property
    def end(self) -> int:
        """Last residue of the domain (across all segments)."""
        return max(s.end for s in self.segments)

    @property
    def span(self) -> Tuple[int, int]:
        """The overall span (first, last) of the domain."""
        return (self.start, self.end)

    def total_length(self) -> int:
        """Total residues covered by all segments."""
        return sum(s.length() for s in self.segments)

    def __repr__(self):
        seg_str = "_".join(str(s) for s in self.segments)
        return f"Domain({seg_str}, {self.method})"


@dataclass
class SequenceEntry:
    """A sequence from the SEED alignment."""
    accession: str          # Full accession with version (e.g., A0A6N4V552.1)
    base_accession: str     # Without version (e.g., A0A6N4V552)
    region_start: int       # Start of aligned region in full sequence
    region_end: int         # End of aligned region in full sequence
    aligned_seq: str        # The aligned sequence (with gaps)

    def __repr__(self):
        return f"{self.accession}/{self.region_start}-{self.region_end}"


@dataclass
class MethodPrediction:
    """Predictions from a single method for a protein."""
    method: str
    uniprot_acc: str
    domains: List[Domain]
    nres_chain: int = 0
    score: float = 0.0


@dataclass
class ProteinPredictions:
    """All method predictions for a single protein."""
    uniprot_acc: str
    predictions: Dict[str, MethodPrediction]  # method -> prediction

    def get_all_domains(self) -> List[Domain]:
        """Get all domains from all methods."""
        domains = []
        for pred in self.predictions.values():
            domains.extend(pred.domains)
        return domains


@dataclass
class AlignmentDomain:
    """A domain mapped to alignment coordinates."""
    seq_entry: SequenceEntry
    original_domain: Domain
    ali_start: int  # Start position in alignment (1-indexed)
    ali_end: int    # End position in alignment (1-indexed)

    def __repr__(self):
        return f"AliDomain({self.seq_entry.accession}, {self.ali_start}-{self.ali_end}, {self.original_domain.method})"


@dataclass
class ConsensusResult:
    """Result of consensus calculation for a domain region."""
    ali_start: int
    ali_end: int
    method_agreement: Dict[str, int]  # method -> count of sequences supporting
    total_sequences: int
    supporting_domains: List[AlignmentDomain]
    score: float
    consensus_type: str  # e.g., "majority", "any_two", "all_three"


# ============================================================================
# SEED Alignment Parser
# ============================================================================

def parse_seed_alignment(seed_file: str) -> List[SequenceEntry]:
    """
    Parse a Pfam SEED alignment file (MUL format).

    MUL format: each named line is ACC.version/start-end followed by whitespace
    and then the aligned sequence. Alignments may be wrapped across multiple lines.
    Continuation lines (without sequence names) are appended to the previous entry.

    Gap characters can be '.' or '-'. Amino acids can be upper or lower case.

    Returns list of SequenceEntry objects.
    """
    # First pass: collect all lines grouped by sequence
    seq_data = {}  # accession -> {'region_start', 'region_end', 'seq_parts': []}
    seq_order = []  # Preserve order of first occurrence

    # Pattern to match sequence header lines: ACC.version/start-end  SEQUENCE
    # Note: no $ anchor to allow trailing whitespace
    seq_pattern = re.compile(r'^(\S+)/(\d+)-(\d+)\s+(\S+)')

    # Pattern to match continuation lines (just sequence, possibly with leading whitespace)
    # Includes: letters (upper/lower), dots, and dashes for gaps
    cont_pattern = re.compile(r'^\s*([A-Za-z.\-]+)\s*$')

    current_acc = None

    with open(seed_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n')

            # Skip empty lines and any Stockholm-style comments (for compatibility)
            if not line.strip() or line.startswith('#') or line.startswith('//'):
                continue

            # Try to match a named sequence line
            match = seq_pattern.match(line)
            if match:
                full_acc = match.group(1)
                region_start = int(match.group(2))
                region_end = int(match.group(3))
                aligned_seq = match.group(4)

                key = f"{full_acc}/{region_start}-{region_end}"

                if key not in seq_data:
                    seq_order.append(key)
                    seq_data[key] = {
                        'accession': full_acc,
                        'region_start': region_start,
                        'region_end': region_end,
                        'seq_parts': [aligned_seq]
                    }
                else:
                    # Wrapped alignment - same sequence appearing again
                    seq_data[key]['seq_parts'].append(aligned_seq)

                current_acc = key
            else:
                # Check if this is a continuation line
                cont_match = cont_pattern.match(line)
                if cont_match and current_acc:
                    seq_data[current_acc]['seq_parts'].append(cont_match.group(1))

    # Build entries from collected data
    entries = []
    for key in seq_order:
        data = seq_data[key]
        full_acc = data['accession']

        # Extract base accession (remove version suffix)
        if '.' in full_acc:
            base_acc = full_acc.rsplit('.', 1)[0]
        else:
            base_acc = full_acc

        # Concatenate all sequence parts
        aligned_seq = ''.join(data['seq_parts'])

        entry = SequenceEntry(
            accession=full_acc,
            base_accession=base_acc,
            region_start=data['region_start'],
            region_end=data['region_end'],
            aligned_seq=aligned_seq
        )
        entries.append(entry)

    return entries


def get_alignment_length(entries: List[SequenceEntry]) -> int:
    """Get the alignment length from parsed entries."""
    if not entries:
        return 0
    return len(entries[0].aligned_seq)


# ============================================================================
# TED API Client
# ============================================================================

TED_API_BASE = "https://ted.cathdb.info/api/v1"


def parse_chopping(chopping_str: str, method: str) -> List[Domain]:
    """
    Parse a TED chopping string into Domain objects.

    Format: "24-111_266-345,113-265,360-444_602-694,447-601"
    - Comma separates different domains
    - Underscore separates segments within a domain
    """
    if not chopping_str or chopping_str == '-':
        return []

    domains = []
    domain_strs = chopping_str.split(',')

    for idx, domain_str in enumerate(domain_strs):
        segments = []
        segment_strs = domain_str.split('_')

        for seg_str in segment_strs:
            if '-' in seg_str:
                parts = seg_str.split('-')
                if len(parts) == 2:
                    try:
                        start = int(parts[0])
                        end = int(parts[1])
                        segments.append(Segment(start, end))
                    except ValueError:
                        continue

        if segments:
            domain = Domain(
                segments=segments,
                method=method,
                domain_idx=idx
            )
            domains.append(domain)

    return domains


def fetch_ted_predictions(uniprot_acc: str, max_retries: int = 3,
                          retry_delay: float = 1.0,
                          include_consensus: bool = True) -> Optional[ProteinPredictions]:
    """
    Fetch domain predictions from TED methods for a UniProt accession.

    Uses the chainparse endpoint which returns chainsaw, merizo, and unidoc-ndr predictions.
    Also fetches the TED consensus if include_consensus is True.
    """
    url = f"{TED_API_BASE}/uniprot/chainparse/{uniprot_acc}"

    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=30)

            if response.status_code == 404:
                # No predictions available for this protein
                return None

            response.raise_for_status()
            data = response.json()

            predictions = {}
            nres_chain = 0

            for item in data.get('data', []):
                method = item.get('method', '')
                chopping = item.get('chopping', '')
                score = item.get('score', 0.0)
                nres = item.get('nres_chain', 0)
                if nres > nres_chain:
                    nres_chain = nres

                domains = parse_chopping(chopping, method)
                for d in domains:
                    d.score = score

                pred = MethodPrediction(
                    method=method,
                    uniprot_acc=uniprot_acc,
                    domains=domains,
                    nres_chain=nres,
                    score=score
                )
                predictions[method] = pred

            # Fetch TED consensus from the consensus endpoint
            if include_consensus:
                consensus_pred = fetch_ted_consensus(uniprot_acc, nres_chain, max_retries, retry_delay)
                if consensus_pred:
                    predictions['ted_consensus'] = consensus_pred

            return ProteinPredictions(
                uniprot_acc=uniprot_acc,
                predictions=predictions
            )

        except requests.exceptions.RequestException as e:
            if attempt < max_retries - 1:
                time.sleep(retry_delay * (attempt + 1))
            else:
                print(f"Warning: Failed to fetch predictions for {uniprot_acc}: {e}",
                      file=sys.stderr)
                return None

    return None


def fetch_ted_consensus(uniprot_acc: str, nres_chain: int = 0,
                        max_retries: int = 3, retry_delay: float = 1.0) -> Optional[MethodPrediction]:
    """
    Fetch TED consensus domain predictions for a UniProt accession.

    Uses the /summary endpoint which returns consensus domains derived by TED
    from the three prediction methods. Each domain entry has a 'chopping' field
    and 'consensus_level' (high/medium/low).
    """
    url = f"{TED_API_BASE}/uniprot/summary/{uniprot_acc}"

    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=30)

            if response.status_code == 404:
                return None

            response.raise_for_status()
            data = response.json()

            # Parse the summary response - it contains a list of consensus domains
            domains = []
            data_list = data.get('data', [])

            if not data_list:
                return None

            for idx, item in enumerate(data_list):
                chopping = item.get('chopping', '')
                if chopping:
                    # Parse the chopping string for this domain
                    domain_domains = parse_chopping(chopping, 'ted_consensus')
                    for d in domain_domains:
                        d.domain_idx = idx
                        # Use plddt as a quality score if available
                        d.score = item.get('plddt', 0.0) / 100.0  # Normalize to 0-1
                        domains.append(d)

            if not domains:
                return None

            # Get nres from the first domain if available
            nres = nres_chain
            if data_list and 'nres_domain' in data_list[0]:
                # Sum up all domain residues as rough estimate of chain length
                total_nres = sum(item.get('nres_domain', 0) for item in data_list)
                if total_nres > nres:
                    nres = total_nres

            return MethodPrediction(
                method='ted_consensus',
                uniprot_acc=uniprot_acc,
                domains=domains,
                nres_chain=nres,
                score=sum(d.score for d in domains) / len(domains) if domains else 0.0
            )

        except requests.exceptions.RequestException as e:
            if attempt < max_retries - 1:
                time.sleep(retry_delay * (attempt + 1))
            else:
                # Silently fail - consensus may not be available for all proteins
                return None

    return None


def fetch_all_predictions(entries: List[SequenceEntry],
                          verbose: bool = False) -> Dict[str, ProteinPredictions]:
    """
    Fetch TED predictions for all sequences in the alignment.

    Returns dict mapping base_accession -> ProteinPredictions
    """
    predictions = {}
    unique_accs = set(e.base_accession for e in entries)

    total = len(unique_accs)
    for i, acc in enumerate(unique_accs):
        if verbose:
            print(f"Fetching predictions for {acc} ({i+1}/{total})...", file=sys.stderr)

        pred = fetch_ted_predictions(acc)
        if pred:
            predictions[acc] = pred

        # Small delay to be nice to the API
        if i < total - 1:
            time.sleep(0.1)

    return predictions


# ============================================================================
# Coordinate Mapping
# ============================================================================

def build_seq_to_ali_map(aligned_seq: str, region_start: int) -> Dict[int, int]:
    """
    Build a mapping from sequence positions to alignment positions.

    Args:
        aligned_seq: The aligned sequence (with gaps as '.' or '-')
        region_start: The start position of this region in the full sequence

    Returns:
        Dict mapping sequence position (1-indexed) to alignment column (1-indexed)
    """
    seq_to_ali = {}
    seq_pos = region_start

    for ali_pos, char in enumerate(aligned_seq, start=1):
        if char not in '.-':  # Not a gap
            seq_to_ali[seq_pos] = ali_pos
            seq_pos += 1

    return seq_to_ali


def build_ali_to_seq_map(aligned_seq: str, region_start: int) -> Dict[int, int]:
    """
    Build a mapping from alignment positions to sequence positions.

    Args:
        aligned_seq: The aligned sequence (with gaps as '.' or '-')
        region_start: The start position of this region in the full sequence

    Returns:
        Dict mapping alignment column (1-indexed) to sequence position (1-indexed)
        Gap positions are not included in the mapping.
    """
    ali_to_seq = {}
    seq_pos = region_start

    for ali_pos, char in enumerate(aligned_seq, start=1):
        if char not in '.-':  # Not a gap
            ali_to_seq[ali_pos] = seq_pos
            seq_pos += 1

    return ali_to_seq


def map_domain_to_alignment(domain: Domain, seq_entry: SequenceEntry) -> Optional[AlignmentDomain]:
    """
    Map a domain's coordinates from sequence space to alignment space.

    For now, we only consider the overall span (first and last coordinates)
    of multi-segment domains, as per the user's request.
    """
    seq_to_ali = build_seq_to_ali_map(seq_entry.aligned_seq, seq_entry.region_start)

    domain_start = domain.start
    domain_end = domain.end

    # Check if domain overlaps with the aligned region
    if domain_end < seq_entry.region_start or domain_start > seq_entry.region_end:
        return None  # Domain is outside the aligned region

    # Clip domain to aligned region
    clipped_start = max(domain_start, seq_entry.region_start)
    clipped_end = min(domain_end, seq_entry.region_end)

    # Find alignment positions
    # For start: find the first position >= clipped_start that has a mapping
    ali_start = None
    for pos in range(clipped_start, clipped_end + 1):
        if pos in seq_to_ali:
            ali_start = seq_to_ali[pos]
            break

    # For end: find the last position <= clipped_end that has a mapping
    ali_end = None
    for pos in range(clipped_end, clipped_start - 1, -1):
        if pos in seq_to_ali:
            ali_end = seq_to_ali[pos]
            break

    if ali_start is None or ali_end is None:
        return None

    return AlignmentDomain(
        seq_entry=seq_entry,
        original_domain=domain,
        ali_start=ali_start,
        ali_end=ali_end
    )


# ============================================================================
# Consensus Scoring Strategies
# ============================================================================

def boundary_similarity_score(b1: Tuple[int, int], b2: Tuple[int, int],
                               tolerance: int = 10) -> float:
    """
    Calculate similarity between two domain boundaries.

    Score is based on how close the start and end positions are.
    Perfect match = 1.0, positions more than 'tolerance' apart = 0.0
    """
    start_diff = abs(b1[0] - b2[0])
    end_diff = abs(b1[1] - b2[1])

    start_score = max(0, 1 - start_diff / tolerance)
    end_score = max(0, 1 - end_diff / tolerance)

    return (start_score + end_score) / 2


def cluster_domains_by_position(ali_domains: List[AlignmentDomain],
                                  tolerance: int = 10) -> List[List[AlignmentDomain]]:
    """
    Cluster alignment domains that have similar boundaries.

    Uses single-linkage clustering based on boundary similarity.
    """
    if not ali_domains:
        return []

    # Sort by start position
    sorted_domains = sorted(ali_domains, key=lambda d: (d.ali_start, d.ali_end))

    clusters = []
    current_cluster = [sorted_domains[0]]

    for domain in sorted_domains[1:]:
        # Check if this domain is close enough to any domain in current cluster
        merged = False
        for cluster_domain in current_cluster:
            score = boundary_similarity_score(
                (cluster_domain.ali_start, cluster_domain.ali_end),
                (domain.ali_start, domain.ali_end),
                tolerance
            )
            if score > 0:
                current_cluster.append(domain)
                merged = True
                break

        if not merged:
            clusters.append(current_cluster)
            current_cluster = [domain]

    if current_cluster:
        clusters.append(current_cluster)

    return clusters


def calculate_cluster_consensus(cluster: List[AlignmentDomain]) -> Tuple[int, int]:
    """
    Calculate consensus boundaries for a cluster of similar domains.

    Returns median start and end positions.
    """
    starts = sorted([d.ali_start for d in cluster])
    ends = sorted([d.ali_end for d in cluster])

    # Use median
    mid = len(starts) // 2
    if len(starts) % 2 == 0:
        consensus_start = (starts[mid - 1] + starts[mid]) // 2
        consensus_end = (ends[mid - 1] + ends[mid]) // 2
    else:
        consensus_start = starts[mid]
        consensus_end = ends[mid]

    return consensus_start, consensus_end


def analyze_method_agreement(cluster: List[AlignmentDomain]) -> Dict[str, int]:
    """
    Analyze which methods contributed to a domain cluster.

    Returns dict mapping method name to count of domains from that method.
    """
    method_counts = defaultdict(int)
    for domain in cluster:
        method_counts[domain.original_domain.method] += 1
    return dict(method_counts)


def calculate_sequence_coverage(cluster: List[AlignmentDomain],
                                 total_sequences: int) -> float:
    """
    Calculate what fraction of sequences have a domain in this cluster.
    """
    unique_seqs = set(d.seq_entry.accession for d in cluster)
    return len(unique_seqs) / total_sequences if total_sequences > 0 else 0


# ============================================================================
# Consensus Strategies
# ============================================================================

class ConsensusStrategy:
    """Base class for consensus strategies."""

    name: str = "base"

    def score_cluster(self, cluster: List[AlignmentDomain],
                      total_sequences: int,
                      method_agreement: Dict[str, int]) -> float:
        """Score a domain cluster. Higher is better."""
        raise NotImplementedError

    def filter_clusters(self, clusters: List[List[AlignmentDomain]],
                        total_sequences: int,
                        min_score: float = 0.0) -> List[ConsensusResult]:
        """Filter and score clusters, returning ConsensusResult objects."""
        results = []

        for cluster in clusters:
            method_agreement = analyze_method_agreement(cluster)
            score = self.score_cluster(cluster, total_sequences, method_agreement)

            if score >= min_score:
                consensus_start, consensus_end = calculate_cluster_consensus(cluster)
                result = ConsensusResult(
                    ali_start=consensus_start,
                    ali_end=consensus_end,
                    method_agreement=method_agreement,
                    total_sequences=total_sequences,
                    supporting_domains=cluster,
                    score=score,
                    consensus_type=self.name
                )
                results.append(result)

        # Sort by score descending
        results.sort(key=lambda r: r.score, reverse=True)
        return results


class MajorityVoteStrategy(ConsensusStrategy):
    """
    Score based on majority vote across sequences.

    Domains present in more sequences get higher scores.
    """

    name = "majority"

    def score_cluster(self, cluster: List[AlignmentDomain],
                      total_sequences: int,
                      method_agreement: Dict[str, int]) -> float:
        coverage = calculate_sequence_coverage(cluster, total_sequences)
        return coverage


class MethodAgreementStrategy(ConsensusStrategy):
    """
    Score based on agreement between prediction methods.

    Domains predicted by more methods get higher scores.
    """

    name = "method_agreement"

    def score_cluster(self, cluster: List[AlignmentDomain],
                      total_sequences: int,
                      method_agreement: Dict[str, int]) -> float:
        num_methods = len(method_agreement)
        return num_methods / 3.0  # Normalize to 0-1


class CombinedStrategy(ConsensusStrategy):
    """
    Combined score using both sequence coverage and method agreement.

    score = (method_weight * method_score) + (coverage_weight * coverage_score)
    """

    name = "combined"

    def __init__(self, method_weight: float = 0.5, coverage_weight: float = 0.5):
        self.method_weight = method_weight
        self.coverage_weight = coverage_weight

    def score_cluster(self, cluster: List[AlignmentDomain],
                      total_sequences: int,
                      method_agreement: Dict[str, int]) -> float:
        coverage = calculate_sequence_coverage(cluster, total_sequences)
        method_score = len(method_agreement) / 3.0

        return (self.method_weight * method_score +
                self.coverage_weight * coverage)


class AnyTwoMethodsStrategy(ConsensusStrategy):
    """
    Require at least 2 methods to agree for a domain to be considered.
    Score by sequence coverage among agreeing predictions.
    """

    name = "any_two_methods"

    def score_cluster(self, cluster: List[AlignmentDomain],
                      total_sequences: int,
                      method_agreement: Dict[str, int]) -> float:
        if len(method_agreement) < 2:
            return 0.0
        coverage = calculate_sequence_coverage(cluster, total_sequences)
        return coverage


class AllThreeMethodsStrategy(ConsensusStrategy):
    """
    Require all 3 methods to agree for a domain to be considered.
    Score by sequence coverage among agreeing predictions.
    """

    name = "all_three_methods"

    def score_cluster(self, cluster: List[AlignmentDomain],
                      total_sequences: int,
                      method_agreement: Dict[str, int]) -> float:
        if len(method_agreement) < 3:
            return 0.0
        coverage = calculate_sequence_coverage(cluster, total_sequences)
        return coverage


class WeightedMethodStrategy(ConsensusStrategy):
    """
    Weight methods differently based on their reliability.

    Default weights give equal weight to all methods.
    """

    name = "weighted_method"

    def __init__(self, weights: Dict[str, float] = None):
        self.weights = weights or {
            'chainsaw': 1.0,
            'merizo': 1.0,
            'unidoc-ndr': 1.0
        }

    def score_cluster(self, cluster: List[AlignmentDomain],
                      total_sequences: int,
                      method_agreement: Dict[str, int]) -> float:
        weighted_sum = 0
        total_weight = sum(self.weights.values())

        for method, count in method_agreement.items():
            weight = self.weights.get(method, 1.0)
            weighted_sum += weight * (count / total_sequences)

        return weighted_sum / total_weight if total_weight > 0 else 0


# ============================================================================
# Method Consistency Scoring
# ============================================================================

@dataclass
class MethodConsistencyScore:
    """Consistency score for a single method across the family."""
    method: str
    consistency_score: float  # 0-1, higher = more consistent
    sequences_with_predictions: int
    sequences_in_main_cluster: int
    total_domains_predicted: int
    dominant_domain_count: int  # Number of domains predicted by most sequences
    boundary_std_start: float  # Std dev of domain start positions
    boundary_std_end: float    # Std dev of domain end positions


def calculate_method_consistency(entries: List[SequenceEntry],
                                  predictions: Dict[str, ProteinPredictions],
                                  tolerance: int = 10) -> Dict[str, MethodConsistencyScore]:
    """
    Calculate consistency scores for each TED method across the family.

    Consistency is measured by:
    1. How many sequences have the same number of domains
    2. How well domain boundaries cluster together
    3. What fraction of sequences agree on domain positions

    A high score means the method makes uniform predictions across the family.

    Returns:
        Dict mapping method name to MethodConsistencyScore
    """
    methods = ['chainsaw', 'merizo', 'unidoc-ndr', 'ted_consensus']
    results = {}

    for method in methods:
        # Collect alignment-mapped domains for this method
        ali_domains = []
        sequences_with_predictions = 0
        domain_counts = []  # Number of domains per sequence

        for entry in entries:
            pred = predictions.get(entry.base_accession)
            if not pred or method not in pred.predictions:
                continue

            method_pred = pred.predictions[method]
            if method_pred.domains:
                sequences_with_predictions += 1
                domain_counts.append(len(method_pred.domains))

                for domain in method_pred.domains:
                    ali_domain = map_domain_to_alignment(domain, entry)
                    if ali_domain:
                        ali_domains.append(ali_domain)

        if sequences_with_predictions == 0:
            results[method] = MethodConsistencyScore(
                method=method,
                consistency_score=0.0,
                sequences_with_predictions=0,
                sequences_in_main_cluster=0,
                total_domains_predicted=0,
                dominant_domain_count=0,
                boundary_std_start=0.0,
                boundary_std_end=0.0
            )
            continue

        # Find dominant domain count (mode)
        from collections import Counter
        count_distribution = Counter(domain_counts)
        dominant_count = count_distribution.most_common(1)[0][0]
        seqs_with_dominant_count = count_distribution[dominant_count]

        # Cluster domains by position
        clusters = cluster_domains_by_position(ali_domains, tolerance)

        # Find the largest cluster
        main_cluster = max(clusters, key=len) if clusters else []
        seqs_in_main = len(set(d.seq_entry.accession for d in main_cluster))

        # Calculate boundary standard deviations for the main cluster
        if main_cluster:
            starts = [d.ali_start for d in main_cluster]
            ends = [d.ali_end for d in main_cluster]
            import statistics
            std_start = statistics.stdev(starts) if len(starts) > 1 else 0.0
            std_end = statistics.stdev(ends) if len(ends) > 1 else 0.0
        else:
            std_start = 0.0
            std_end = 0.0

        # Calculate consistency score components
        # Component 1: Fraction of sequences in main cluster
        cluster_consistency = seqs_in_main / sequences_with_predictions if sequences_with_predictions > 0 else 0

        # Component 2: Domain count uniformity
        count_uniformity = seqs_with_dominant_count / sequences_with_predictions if sequences_with_predictions > 0 else 0

        # Component 3: Boundary tightness (inverse of standard deviation, normalized)
        max_std = 50  # Consider std > 50 as very inconsistent
        boundary_tightness = 1 - min((std_start + std_end) / 2, max_std) / max_std

        # Combined consistency score (weighted average)
        consistency_score = (
            0.5 * cluster_consistency +
            0.25 * count_uniformity +
            0.25 * boundary_tightness
        )

        results[method] = MethodConsistencyScore(
            method=method,
            consistency_score=consistency_score,
            sequences_with_predictions=sequences_with_predictions,
            sequences_in_main_cluster=seqs_in_main,
            total_domains_predicted=len(ali_domains),
            dominant_domain_count=dominant_count,
            boundary_std_start=std_start,
            boundary_std_end=std_end
        )

    return results


def get_per_method_consensus_domains(entries: List[SequenceEntry],
                                      predictions: Dict[str, ProteinPredictions],
                                      method: str,
                                      tolerance: int = 10) -> List[Tuple[int, int, float]]:
    """
    Calculate consensus domain boundaries for a single method.

    Clusters all domain predictions from this method and returns
    the median boundaries for each cluster.

    Returns:
        List of (ali_start, ali_end, coverage) tuples for consensus domains
    """
    # Collect alignment-mapped domains for this method
    ali_domains = []
    sequences_with_predictions = set()

    for entry in entries:
        pred = predictions.get(entry.base_accession)
        if not pred or method not in pred.predictions:
            continue

        method_pred = pred.predictions[method]
        if method_pred.domains:
            sequences_with_predictions.add(entry.accession)
            for domain in method_pred.domains:
                ali_domain = map_domain_to_alignment(domain, entry)
                if ali_domain:
                    ali_domains.append(ali_domain)

    if not ali_domains:
        return []

    # Cluster domains by position
    clusters = cluster_domains_by_position(ali_domains, tolerance)

    # For each cluster, calculate median boundaries and coverage
    consensus_domains = []
    total_seqs = len(sequences_with_predictions)

    for cluster in clusters:
        if not cluster:
            continue

        # Calculate median boundaries
        starts = sorted([d.ali_start for d in cluster])
        ends = sorted([d.ali_end for d in cluster])
        mid = len(starts) // 2
        if len(starts) % 2 == 0:
            consensus_start = (starts[mid - 1] + starts[mid]) // 2
            consensus_end = (ends[mid - 1] + ends[mid]) // 2
        else:
            consensus_start = starts[mid]
            consensus_end = ends[mid]

        # Calculate coverage
        seqs_in_cluster = len(set(d.seq_entry.accession for d in cluster))
        coverage = seqs_in_cluster / total_seqs if total_seqs > 0 else 0

        consensus_domains.append((consensus_start, consensus_end, coverage))

    # Sort by start position
    consensus_domains.sort(key=lambda x: x[0])
    return consensus_domains


# ============================================================================
# Main Analysis Functions
# ============================================================================

def collect_all_alignment_domains(entries: List[SequenceEntry],
                                   predictions: Dict[str, ProteinPredictions]) -> List[AlignmentDomain]:
    """
    Collect all domain predictions and map them to alignment coordinates.
    """
    ali_domains = []

    for entry in entries:
        pred = predictions.get(entry.base_accession)
        if not pred:
            continue

        for method_pred in pred.predictions.values():
            for domain in method_pred.domains:
                ali_domain = map_domain_to_alignment(domain, entry)
                if ali_domain:
                    ali_domains.append(ali_domain)

    return ali_domains


def run_consensus_analysis(entries: List[SequenceEntry],
                            predictions: Dict[str, ProteinPredictions],
                            strategies: List[ConsensusStrategy],
                            tolerance: int = 10,
                            min_score: float = 0.1) -> Dict[str, List[ConsensusResult]]:
    """
    Run consensus analysis using multiple strategies.

    Returns dict mapping strategy name to list of ConsensusResult objects.
    """
    # Collect all alignment domains
    ali_domains = collect_all_alignment_domains(entries, predictions)

    if not ali_domains:
        return {s.name: [] for s in strategies}

    # Cluster domains by position
    clusters = cluster_domains_by_position(ali_domains, tolerance)

    total_sequences = len(entries)

    # Apply each strategy
    results = {}
    for strategy in strategies:
        strategy_results = strategy.filter_clusters(clusters, total_sequences, min_score)
        results[strategy.name] = strategy_results

    return results


def get_per_sequence_consensus(entries: List[SequenceEntry],
                                predictions: Dict[str, ProteinPredictions],
                                consensus_results: List[ConsensusResult]) -> Dict[str, List[Tuple[int, int, float]]]:
    """
    For each sequence, determine which consensus domains apply.

    Returns dict mapping accession to list of (seq_start, seq_end, score) tuples.
    """
    per_seq = {}

    for entry in entries:
        ali_to_seq = build_ali_to_seq_map(entry.aligned_seq, entry.region_start)
        seq_domains = []

        for result in consensus_results:
            # Find sequence positions that correspond to consensus alignment positions
            # Find the sequence position closest to ali_start
            seq_start = None
            for ali_pos in range(result.ali_start, result.ali_end + 1):
                if ali_pos in ali_to_seq:
                    seq_start = ali_to_seq[ali_pos]
                    break

            # Find the sequence position closest to ali_end
            seq_end = None
            for ali_pos in range(result.ali_end, result.ali_start - 1, -1):
                if ali_pos in ali_to_seq:
                    seq_end = ali_to_seq[ali_pos]
                    break

            if seq_start and seq_end and seq_start <= seq_end:
                seq_domains.append((seq_start, seq_end, result.score))

        per_seq[entry.accession] = seq_domains

    return per_seq


# ============================================================================
# Report Generation
# ============================================================================

def generate_summary_report(entries: List[SequenceEntry],
                             predictions: Dict[str, ProteinPredictions],
                             results: Dict[str, List[ConsensusResult]],
                             output_file: str = None,
                             consistency_scores: Dict[str, MethodConsistencyScore] = None,
                             tolerance: int = 10) -> str:
    """
    Generate a comprehensive summary report.
    """
    lines = []
    lines.append("=" * 80)
    lines.append("TED FAMILY CONSENSUS ANALYSIS REPORT")
    lines.append("=" * 80)
    lines.append("")

    # Summary statistics
    total_seqs = len(entries)
    seqs_with_predictions = sum(1 for e in entries if e.base_accession in predictions)

    lines.append("SUMMARY")
    lines.append("-" * 40)
    lines.append(f"Total sequences in alignment: {total_seqs}")
    lines.append(f"Sequences with TED predictions: {seqs_with_predictions}")
    lines.append(f"Sequences without predictions: {total_seqs - seqs_with_predictions}")
    lines.append("")

    # Method coverage
    method_counts = defaultdict(int)
    total_domains_per_method = defaultdict(int)

    for pred in predictions.values():
        for method, method_pred in pred.predictions.items():
            if method_pred.domains:
                method_counts[method] += 1
                total_domains_per_method[method] += len(method_pred.domains)

    lines.append("METHOD COVERAGE")
    lines.append("-" * 40)
    for method in ['chainsaw', 'merizo', 'unidoc-ndr', 'ted_consensus']:
        count = method_counts.get(method, 0)
        domains = total_domains_per_method.get(method, 0)
        if count > 0 or method != 'ted_consensus':  # Always show the 3 main methods
            lines.append(f"  {method:15s}: {count:4d} proteins, {domains:4d} total domains")
    lines.append("")

    # Method Consistency Scores
    if consistency_scores is None:
        consistency_scores = calculate_method_consistency(entries, predictions, tolerance)

    lines.append("METHOD CONSISTENCY SCORES")
    lines.append("-" * 40)
    lines.append("  (Higher score = more consistent predictions across the family)")
    lines.append("")
    lines.append(f"  {'Method':<15s}  {'Score':>8s}  {'Seqs':>5s}  {'Cluster':>7s}  {'Domains':>7s}  {'StdDev':>10s}")
    lines.append("  " + "-" * 60)

    # Sort by consistency score descending
    sorted_methods = sorted(
        [(m, s) for m, s in consistency_scores.items() if s.sequences_with_predictions > 0],
        key=lambda x: x[1].consistency_score,
        reverse=True
    )

    for method, score in sorted_methods:
        std_str = f"{score.boundary_std_start:.1f}/{score.boundary_std_end:.1f}"
        lines.append(
            f"  {method:<15s}  {score.consistency_score:8.3f}  "
            f"{score.sequences_with_predictions:5d}  {score.sequences_in_main_cluster:7d}  "
            f"{score.total_domains_predicted:7d}  {std_str:>10s}"
        )

    lines.append("")
    lines.append("  Legend: Seqs = sequences with predictions, Cluster = sequences in main cluster,")
    lines.append("          Domains = total domains predicted, StdDev = boundary std dev (start/end)")
    lines.append("")

    # Best method recommendation
    if sorted_methods:
        best_method, best_score = sorted_methods[0]
        lines.append(f"  ** Best method for this family: {best_method} (score: {best_score.consistency_score:.3f}) **")
        lines.append("")

    # Results by strategy
    lines.append("CONSENSUS RESULTS BY STRATEGY")
    lines.append("-" * 40)

    for strategy_name, strategy_results in results.items():
        lines.append("")
        lines.append(f"Strategy: {strategy_name}")
        lines.append(f"  Consensus domains found: {len(strategy_results)}")

        if strategy_results:
            lines.append("")
            lines.append("  Rank  Ali.Start  Ali.End  Score   Methods                 Seq.Coverage")
            lines.append("  " + "-" * 75)

            for i, result in enumerate(strategy_results[:10], 1):  # Top 10
                methods_str = ", ".join(f"{m}:{c}" for m, c in sorted(result.method_agreement.items()))
                coverage = calculate_sequence_coverage(result.supporting_domains, result.total_sequences)
                lines.append(f"  {i:4d}  {result.ali_start:9d}  {result.ali_end:7d}  {result.score:.4f}  {methods_str:22s}  {coverage:.1%}")

    lines.append("")
    lines.append("=" * 80)

    report = "\n".join(lines)

    if output_file:
        with open(output_file, 'w') as f:
            f.write(report)

    return report


def generate_per_sequence_report(entries: List[SequenceEntry],
                                  per_seq_consensus: Dict[str, List[Tuple[int, int, float]]],
                                  output_file: str = None) -> str:
    """
    Generate a per-sequence consensus domain report.
    """
    lines = []
    lines.append("# Per-Sequence Consensus Domains")
    lines.append("# Accession\tDomain_Num\tSeq_Start\tSeq_End\tScore")

    for entry in entries:
        domains = per_seq_consensus.get(entry.accession, [])
        if domains:
            for i, (start, end, score) in enumerate(domains, 1):
                lines.append(f"{entry.accession}\t{i}\t{start}\t{end}\t{score:.4f}")
        else:
            lines.append(f"{entry.accession}\t-\t-\t-\t-")

    report = "\n".join(lines)

    if output_file:
        with open(output_file, 'w') as f:
            f.write(report)

    return report


def generate_detailed_json(entries: List[SequenceEntry],
                            predictions: Dict[str, ProteinPredictions],
                            results: Dict[str, List[ConsensusResult]],
                            output_file: str = None) -> dict:
    """
    Generate detailed JSON output for programmatic use.
    """
    output = {
        "summary": {
            "total_sequences": len(entries),
            "sequences_with_predictions": sum(1 for e in entries if e.base_accession in predictions),
        },
        "method_predictions": {},
        "consensus_results": {}
    }

    # Add method predictions
    for acc, pred in predictions.items():
        output["method_predictions"][acc] = {}
        for method, method_pred in pred.predictions.items():
            output["method_predictions"][acc][method] = {
                "domains": [
                    {"start": d.start, "end": d.end, "segments": [{"start": s.start, "end": s.end} for s in d.segments]}
                    for d in method_pred.domains
                ],
                "score": method_pred.score
            }

    # Add consensus results
    for strategy_name, strategy_results in results.items():
        output["consensus_results"][strategy_name] = [
            {
                "ali_start": r.ali_start,
                "ali_end": r.ali_end,
                "score": r.score,
                "method_agreement": r.method_agreement,
                "sequence_coverage": calculate_sequence_coverage(r.supporting_domains, r.total_sequences)
            }
            for r in strategy_results
        ]

    if output_file:
        with open(output_file, 'w') as f:
            json.dump(output, f, indent=2)

    return output


# ============================================================================
# Per-Method Reports and Visualization
# ============================================================================

def generate_per_method_reports(entries: List[SequenceEntry],
                                 predictions: Dict[str, ProteinPredictions],
                                 output_prefix: str) -> Dict[str, str]:
    """
    Generate separate per-sequence domain reports for each TED method.

    Creates one TSV file per method (chainsaw, merizo, unidoc-ndr, ted_consensus) with
    domain boundaries for each sequence.

    Args:
        entries: List of sequence entries from alignment
        predictions: Dict of protein predictions
        output_prefix: Prefix for output files (e.g., "family" -> "family_chainsaw.tsv")

    Returns:
        Dict mapping method name to output file path
    """
    methods = ['chainsaw', 'merizo', 'unidoc-ndr', 'ted_consensus']
    output_files = {}

    for method in methods:
        lines = []
        lines.append(f"# TED {method} domain predictions")
        lines.append("# Accession\tRegion\tDomain_Num\tSeq_Start\tSeq_End\tDomain_Chopping")

        for entry in entries:
            pred = predictions.get(entry.base_accession)
            region_str = f"{entry.region_start}-{entry.region_end}"

            if pred and method in pred.predictions:
                method_pred = pred.predictions[method]
                if method_pred.domains:
                    for i, domain in enumerate(method_pred.domains, 1):
                        # Build chopping string for all segments
                        chopping = "_".join(f"{s.start}-{s.end}" for s in domain.segments)
                        lines.append(f"{entry.accession}\t{region_str}\t{i}\t{domain.start}\t{domain.end}\t{chopping}")
                else:
                    lines.append(f"{entry.accession}\t{region_str}\t-\t-\t-\t-")
            else:
                lines.append(f"{entry.accession}\t{region_str}\t-\t-\t-\tno_prediction")

        output_file = f"{output_prefix}_{method}.tsv"
        with open(output_file, 'w') as f:
            f.write("\n".join(lines))
        output_files[method] = output_file

    return output_files


def create_method_visualization(entries: List[SequenceEntry],
                                 predictions: Dict[str, ProteinPredictions],
                                 method: str,
                                 output_file: str,
                                 title: str = None,
                                 view_mode: str = 'sequence'):
    """
    Create a visualization of domain predictions from a single method.

    Two view modes available:
    - 'sequence': Shows each protein at its actual length with domains at sequence
                  coordinates and SEED region as black outline box (like add_image.py)
    - 'alignment': Shows all sequences at alignment length with domains mapped to
                   alignment coordinates (useful for comparing domain positions)

    Args:
        entries: List of sequence entries from alignment
        predictions: Dict of protein predictions
        method: Method name ('chainsaw', 'merizo', or 'unidoc-ndr')
        output_file: Output PNG file path
        title: Optional title for the figure
        view_mode: 'sequence' (default) or 'alignment'
    """
    try:
        import matplotlib
        matplotlib.use('Agg')  # Non-interactive backend
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
    except ImportError:
        print(f"Warning: matplotlib not available, skipping {method} visualization", file=sys.stderr)
        return

    # Color scheme for domains (from add_image.py)
    domain_colors = ['#4A79A7', '#F28E2C', '#E15759', '#76B7B2', '#59A14F',
                     '#EDC949', '#AF7AA1', '#FF9DA7', '#9C755F', '#BAB0AB']

    # Calculate dimensions based on view mode
    row_height = 25
    margin_left = 180
    margin_right = 50
    margin_top = 50
    margin_bottom = 30

    if view_mode == 'alignment':
        # Alignment view: all sequences shown at alignment length
        ali_length = len(entries[0].aligned_seq) if entries else 0
        pixels_per_unit = 2
        max_length = ali_length
        axis_label = "Alignment Position"
    else:
        # Sequence view: proteins shown at actual length
        max_protein_length = 0
        for entry in entries:
            pred = predictions.get(entry.base_accession)
            if pred:
                for method_pred in pred.predictions.values():
                    if method_pred.nres_chain > 0:
                        max_protein_length = max(max_protein_length, method_pred.nres_chain)
                        break
            max_protein_length = max(max_protein_length, entry.region_end)
        if max_protein_length == 0:
            max_protein_length = 500
        pixels_per_unit = max(1, 1500 // max_protein_length)
        max_length = max_protein_length
        axis_label = "Sequence Position"

    fig_width = margin_left + (max_length * pixels_per_unit) + margin_right
    fig_height = margin_top + (len(entries) * row_height) + margin_bottom

    # Create figure
    fig, ax = plt.subplots(figsize=(fig_width / 100, fig_height / 100), dpi=100)
    ax.set_xlim(0, fig_width)
    ax.set_ylim(0, fig_height)
    ax.axis('off')

    # Add title
    if title:
        ax.text(fig_width / 2, fig_height - 15, title,
                ha='center', va='top', fontsize=11, weight='bold')

    # Add axis labels
    ax.text(margin_left, fig_height - margin_top + 10, "1",
            ha='left', va='bottom', fontsize=8)
    ax.text(margin_left + max_length * pixels_per_unit, fig_height - margin_top + 10,
            str(max_length), ha='right', va='bottom', fontsize=8)
    ax.text(margin_left + (max_length * pixels_per_unit) / 2, fig_height - margin_top + 10,
            axis_label, ha='center', va='bottom', fontsize=9)

    # Draw each sequence
    current_y = fig_height - margin_top - row_height / 2

    for entry in entries:
        pred = predictions.get(entry.base_accession)

        # Draw label
        label = f"{entry.accession}"
        ax.text(margin_left - 5, current_y, label,
                ha='right', va='center', fontsize=7, family='monospace')

        if view_mode == 'alignment':
            # Alignment view: backbone spans alignment length
            backbone_length = len(entry.aligned_seq)
        else:
            # Sequence view: backbone spans actual protein length
            backbone_length = entry.region_end
            if pred:
                for method_pred in pred.predictions.values():
                    if method_pred.nres_chain > 0:
                        backbone_length = method_pred.nres_chain
                        break

        # Draw backbone (grey line)
        backbone_height = 6
        backbone = patches.Rectangle(
            (margin_left, current_y - backbone_height / 2),
            backbone_length * pixels_per_unit, backbone_height,
            linewidth=0, facecolor='#CCCCCC'
        )
        ax.add_patch(backbone)

        # Draw domains for this method
        if pred and method in pred.predictions:
            method_pred = pred.predictions[method]

            for domain_idx, domain in enumerate(method_pred.domains):
                color = domain_colors[domain_idx % len(domain_colors)]

                if view_mode == 'alignment':
                    # Map domain to alignment coordinates
                    ali_domain = map_domain_to_alignment(domain, entry)
                    if ali_domain:
                        domain_x = margin_left + (ali_domain.ali_start - 1) * pixels_per_unit
                        domain_width = (ali_domain.ali_end - ali_domain.ali_start + 1) * pixels_per_unit
                    else:
                        continue
                else:
                    # Use sequence coordinates directly
                    domain_x = margin_left + (domain.start - 1) * pixels_per_unit
                    domain_width = (domain.end - domain.start + 1) * pixels_per_unit

                domain_height = 16
                domain_rect = patches.Rectangle(
                    (domain_x, current_y - domain_height / 2),
                    domain_width, domain_height,
                    linewidth=1, edgecolor='black', facecolor=color, alpha=0.8
                )
                ax.add_patch(domain_rect)

                # Add domain label if space permits
                if domain_width > 30:
                    label_text = f"D{domain_idx + 1}"
                    ax.text(domain_x + domain_width / 2, current_y,
                            label_text, ha='center', va='center',
                            fontsize=6, weight='bold', color='white')

        # Draw SEED region (black outline box)
        if view_mode == 'alignment':
            # In alignment view, SEED region spans the non-gap portion
            # For simplicity, show as full alignment width
            pass  # Skip SEED box in alignment view
        else:
            # In sequence view, show SEED region
            seed_x = margin_left + (entry.region_start - 1) * pixels_per_unit
            seed_width = (entry.region_end - entry.region_start + 1) * pixels_per_unit
            seed_height = 22
            seed_rect = patches.Rectangle(
                (seed_x, current_y - seed_height / 2),
                seed_width, seed_height,
                linewidth=2, edgecolor='black', facecolor='none'
            )
            ax.add_patch(seed_rect)

        current_y -= row_height

    # Save figure
    plt.tight_layout()
    plt.savefig(output_file, dpi=100, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Created {method} visualization: {output_file}", file=sys.stderr)


def create_all_method_visualizations(entries: List[SequenceEntry],
                                      predictions: Dict[str, ProteinPredictions],
                                      output_prefix: str,
                                      view_mode: str = 'sequence') -> Dict[str, str]:
    """
    Create visualization images for all TED methods (including ted_consensus if available).

    Args:
        entries: List of sequence entries from alignment
        predictions: Dict of protein predictions
        output_prefix: Prefix for output files (e.g., "family" -> "family_chainsaw.png")
        view_mode: 'sequence' (default) or 'alignment'

    Returns:
        Dict mapping method name to output file path
    """
    methods = ['chainsaw', 'merizo', 'unidoc-ndr', 'ted_consensus']
    output_files = {}

    for method in methods:
        # Check if this method has any predictions
        has_predictions = any(
            method in pred.predictions and pred.predictions[method].domains
            for pred in predictions.values()
        )

        if has_predictions:
            output_file = f"{output_prefix}_{method}.png"
            title = f"TED {method.replace('_', ' ').title()} Domain Predictions"
            create_method_visualization(entries, predictions, method, output_file, title, view_mode)
            output_files[method] = output_file

    return output_files


def create_consensus_visualization(entries: List[SequenceEntry],
                                    predictions: Dict[str, ProteinPredictions],
                                    output_file: str,
                                    title: str = "TED Consensus Domain Predictions",
                                    view_mode: str = 'sequence'):
    """
    Create a visualization showing all three methods side-by-side for comparison.

    Each sequence gets three rows (one per method) to easily compare predictions.

    Args:
        entries: List of sequence entries from alignment
        predictions: Dict of protein predictions
        output_file: Output PNG file path
        title: Title for the figure
        view_mode: 'sequence' (default) or 'alignment'
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
    except ImportError:
        print(f"Warning: matplotlib not available, skipping consensus visualization", file=sys.stderr)
        return

    methods = ['chainsaw', 'merizo', 'unidoc-ndr']
    method_colors = {
        'chainsaw': '#4A79A7',
        'merizo': '#F28E2C',
        'unidoc-ndr': '#E15759'
    }

    # Calculate dimensions based on view mode
    method_row_height = 12
    seq_spacing = 8
    seq_block_height = len(methods) * method_row_height + seq_spacing
    margin_left = 180
    margin_right = 80
    margin_top = 60
    margin_bottom = 30

    if view_mode == 'alignment':
        ali_length = len(entries[0].aligned_seq) if entries else 0
        pixels_per_unit = 2
        max_length = ali_length
        axis_label = f"Alignment Position (1-{ali_length})"
    else:
        max_protein_length = 0
        for entry in entries:
            pred = predictions.get(entry.base_accession)
            if pred:
                for method_pred in pred.predictions.values():
                    if method_pred.nres_chain > 0:
                        max_protein_length = max(max_protein_length, method_pred.nres_chain)
                        break
            max_protein_length = max(max_protein_length, entry.region_end)
        if max_protein_length == 0:
            max_protein_length = 500
        pixels_per_unit = max(1, 1500 // max_protein_length)
        max_length = max_protein_length
        axis_label = f"Sequence Position (1-{max_length})"

    fig_width = margin_left + (max_length * pixels_per_unit) + margin_right
    fig_height = margin_top + (len(entries) * seq_block_height) + margin_bottom

    # Create figure
    fig, ax = plt.subplots(figsize=(fig_width / 100, fig_height / 100), dpi=100)
    ax.set_xlim(0, fig_width)
    ax.set_ylim(0, fig_height)
    ax.axis('off')

    # Add title
    ax.text(fig_width / 2, fig_height - 15, title,
            ha='center', va='top', fontsize=11, weight='bold')

    # Add legend
    legend_y = fig_height - 35
    legend_x = margin_left
    for i, method in enumerate(methods):
        x = legend_x + i * 150
        legend_box = patches.Rectangle(
            (x, legend_y - 4), 12, 8,
            linewidth=1, edgecolor='black', facecolor=method_colors[method], alpha=0.8
        )
        ax.add_patch(legend_box)
        ax.text(x + 16, legend_y, method, fontsize=8, va='center')

    # Add axis labels
    ax.text(margin_left + (max_length * pixels_per_unit) / 2, fig_height - margin_top + 5,
            axis_label, ha='center', va='bottom', fontsize=9)

    # Draw each sequence
    current_y = fig_height - margin_top - seq_block_height / 2

    for entry in entries:
        pred = predictions.get(entry.base_accession)

        # Get protein length for sequence view
        protein_length = entry.region_end
        if pred:
            for method_pred in pred.predictions.values():
                if method_pred.nres_chain > 0:
                    protein_length = method_pred.nres_chain
                    break

        # Draw label (centered on sequence block)
        label = f"{entry.accession}"
        ax.text(margin_left - 5, current_y + method_row_height,
                label, ha='right', va='center', fontsize=7, family='monospace')

        # Draw each method row
        for method_idx, method in enumerate(methods):
            row_y = current_y + (len(methods) - 1 - method_idx) * method_row_height

            # Draw backbone
            backbone_height = 4
            if view_mode == 'alignment':
                backbone_length = len(entry.aligned_seq)
            else:
                backbone_length = protein_length

            backbone = patches.Rectangle(
                (margin_left, row_y - backbone_height / 2),
                backbone_length * pixels_per_unit, backbone_height,
                linewidth=0, facecolor='#EEEEEE'
            )
            ax.add_patch(backbone)

            # Draw domains
            if pred and method in pred.predictions:
                method_pred = pred.predictions[method]

                for domain in method_pred.domains:
                    if view_mode == 'alignment':
                        ali_domain = map_domain_to_alignment(domain, entry)
                        if ali_domain:
                            domain_x = margin_left + (ali_domain.ali_start - 1) * pixels_per_unit
                            domain_width = (ali_domain.ali_end - ali_domain.ali_start + 1) * pixels_per_unit
                        else:
                            continue
                    else:
                        domain_x = margin_left + (domain.start - 1) * pixels_per_unit
                        domain_width = (domain.end - domain.start + 1) * pixels_per_unit

                    domain_height = 10
                    domain_rect = patches.Rectangle(
                        (domain_x, row_y - domain_height / 2),
                        domain_width, domain_height,
                        linewidth=0.5, edgecolor='black',
                        facecolor=method_colors[method], alpha=0.8
                    )
                    ax.add_patch(domain_rect)

        # Draw SEED region in sequence view
        if view_mode != 'alignment':
            seed_x = margin_left + (entry.region_start - 1) * pixels_per_unit
            seed_width = (entry.region_end - entry.region_start + 1) * pixels_per_unit
            seed_height = seq_block_height - 4
            seed_rect = patches.Rectangle(
                (seed_x, current_y - seed_height / 2 + method_row_height / 2),
                seed_width, seed_height,
                linewidth=1.5, edgecolor='black', facecolor='none'
            )
            ax.add_patch(seed_rect)

        current_y -= seq_block_height

    # Save figure
    plt.tight_layout()
    plt.savefig(output_file, dpi=100, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Created consensus visualization: {output_file}", file=sys.stderr)


def create_computed_consensus_visualization(entries: List[SequenceEntry],
                                             predictions: Dict[str, ProteinPredictions],
                                             method: str,
                                             output_file: str,
                                             tolerance: int = 10,
                                             title: str = None):
    """
    Create a visualization showing computed consensus domains for a single method.

    Shows each sequence with:
    - Grey backbone representing the protein length
    - Colored consensus domains (computed from clustering all predictions)
    - Black outline box for SEED region

    The consensus domains are computed by clustering all domain predictions
    from this method across all sequences and taking the median boundaries.

    Args:
        entries: List of sequence entries from alignment
        predictions: Dict of protein predictions
        method: Method name ('chainsaw', 'merizo', 'unidoc-ndr', 'ted_consensus')
        output_file: Output PNG file path
        tolerance: Tolerance for boundary clustering
        title: Optional title for the figure
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
    except ImportError:
        print(f"Warning: matplotlib not available, skipping {method} consensus visualization", file=sys.stderr)
        return

    # Get computed consensus domains for this method (in alignment coordinates)
    consensus_domains = get_per_method_consensus_domains(entries, predictions, method, tolerance)

    if not consensus_domains:
        print(f"Warning: No consensus domains found for {method}", file=sys.stderr)
        return

    # Color scheme for consensus domains
    domain_colors = ['#4A79A7', '#F28E2C', '#E15759', '#76B7B2', '#59A14F',
                     '#EDC949', '#AF7AA1', '#FF9DA7', '#9C755F', '#BAB0AB']

    # Calculate dimensions
    row_height = 25
    margin_left = 180
    margin_right = 100
    margin_top = 70
    margin_bottom = 30

    # Use alignment coordinates
    ali_length = len(entries[0].aligned_seq) if entries else 0
    pixels_per_unit = 2
    max_length = ali_length

    fig_width = margin_left + (max_length * pixels_per_unit) + margin_right
    fig_height = margin_top + (len(entries) * row_height) + margin_bottom

    # Create figure
    fig, ax = plt.subplots(figsize=(fig_width / 100, fig_height / 100), dpi=100)
    ax.set_xlim(0, fig_width)
    ax.set_ylim(0, fig_height)
    ax.axis('off')

    # Add title
    if title is None:
        title = f"Computed Consensus Domains - {method}"
    ax.text(fig_width / 2, fig_height - 12, title,
            ha='center', va='top', fontsize=11, weight='bold')

    # Add legend for consensus domains
    legend_y = fig_height - 35
    legend_x = margin_left
    ax.text(legend_x - 10, legend_y, "Consensus domains:", fontsize=8, va='center')
    for i, (ali_start, ali_end, coverage) in enumerate(consensus_domains):
        x = legend_x + 100 + i * 120
        color = domain_colors[i % len(domain_colors)]
        legend_box = patches.Rectangle(
            (x, legend_y - 4), 12, 8,
            linewidth=1, edgecolor='black', facecolor=color, alpha=0.8
        )
        ax.add_patch(legend_box)
        ax.text(x + 16, legend_y, f"D{i+1} ({coverage:.0%})", fontsize=8, va='center')

    # Add axis labels
    ax.text(margin_left, fig_height - margin_top + 10, "1",
            ha='left', va='bottom', fontsize=8)
    ax.text(margin_left + max_length * pixels_per_unit, fig_height - margin_top + 10,
            str(ali_length), ha='right', va='bottom', fontsize=8)
    ax.text(margin_left + (max_length * pixels_per_unit) / 2, fig_height - margin_top + 10,
            "Alignment Position", ha='center', va='bottom', fontsize=9)

    # Draw each sequence
    current_y = fig_height - margin_top - row_height / 2

    for entry in entries:
        # Draw label
        label = f"{entry.accession}"
        ax.text(margin_left - 5, current_y, label,
                ha='right', va='center', fontsize=7, family='monospace')

        # Draw backbone (alignment length)
        backbone_length = len(entry.aligned_seq)
        backbone_height = 6
        backbone = patches.Rectangle(
            (margin_left, current_y - backbone_height / 2),
            backbone_length * pixels_per_unit, backbone_height,
            linewidth=0, facecolor='#CCCCCC'
        )
        ax.add_patch(backbone)

        # Draw consensus domains
        for i, (ali_start, ali_end, coverage) in enumerate(consensus_domains):
            color = domain_colors[i % len(domain_colors)]
            domain_x = margin_left + (ali_start - 1) * pixels_per_unit
            domain_width = (ali_end - ali_start + 1) * pixels_per_unit

            domain_height = 16
            domain_rect = patches.Rectangle(
                (domain_x, current_y - domain_height / 2),
                domain_width, domain_height,
                linewidth=1, edgecolor='black', facecolor=color, alpha=0.8
            )
            ax.add_patch(domain_rect)

            # Add domain label if space permits
            if domain_width > 30:
                label_text = f"D{i + 1}"
                ax.text(domain_x + domain_width / 2, current_y,
                        label_text, ha='center', va='center',
                        fontsize=6, weight='bold', color='white')

        # Draw SEED region (mapped to alignment coordinates)
        seq_to_ali = build_seq_to_ali_map(entry.aligned_seq, entry.region_start)
        if entry.region_start in seq_to_ali and entry.region_end in seq_to_ali:
            seed_ali_start = seq_to_ali[entry.region_start]
            seed_ali_end = seq_to_ali[entry.region_end]
            seed_x = margin_left + (seed_ali_start - 1) * pixels_per_unit
            seed_width = (seed_ali_end - seed_ali_start + 1) * pixels_per_unit
            seed_height = 22
            seed_rect = patches.Rectangle(
                (seed_x, current_y - seed_height / 2),
                seed_width, seed_height,
                linewidth=2, edgecolor='black', facecolor='none'
            )
            ax.add_patch(seed_rect)

        current_y -= row_height

    # Add summary info
    info_text = f"Consensus domains: {len(consensus_domains)}, Tolerance: {tolerance} residues"
    ax.text(fig_width / 2, 10, info_text, ha='center', va='bottom', fontsize=8, style='italic')

    # Save figure
    plt.tight_layout()
    plt.savefig(output_file, dpi=100, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Created {method} computed consensus visualization: {output_file}", file=sys.stderr)


def create_all_computed_consensus_visualizations(entries: List[SequenceEntry],
                                                  predictions: Dict[str, ProteinPredictions],
                                                  output_prefix: str,
                                                  tolerance: int = 10) -> Dict[str, str]:
    """
    Create computed consensus domain visualizations for all methods.

    These images show the consensus domains computed by clustering all predictions
    from each method, rather than the raw individual predictions.

    Args:
        entries: List of sequence entries from alignment
        predictions: Dict of protein predictions
        output_prefix: Prefix for output files
        tolerance: Tolerance for boundary clustering

    Returns:
        Dict mapping method name to output file path
    """
    methods = ['chainsaw', 'merizo', 'unidoc-ndr', 'ted_consensus']
    output_files = {}

    for method in methods:
        # Check if this method has any predictions
        has_predictions = any(
            method in pred.predictions and pred.predictions[method].domains
            for pred in predictions.values()
        )

        if has_predictions:
            output_file = f"{output_prefix}_{method}_computed_consensus.png"
            title = f"Computed Consensus Domains - {method}"
            create_computed_consensus_visualization(
                entries, predictions, method, output_file, tolerance, title
            )
            output_files[method] = output_file

    return output_files


def create_cross_method_consensus_visualization(entries: List[SequenceEntry],
                                                 consensus_results: List['ConsensusResult'],
                                                 output_file: str,
                                                 title: str = "Cross-Method Consensus Domains",
                                                 min_coverage: float = 0.0):
    """
    Create a visualization showing the cross-method consensus domains.

    These are the domains computed by clustering predictions from ALL methods
    together, as shown in the "CONSENSUS RESULTS BY STRATEGY" table.

    Args:
        entries: List of sequence entries from alignment
        consensus_results: List of ConsensusResult objects from run_consensus_analysis
        output_file: Output PNG file path
        title: Title for the figure
        min_coverage: Minimum sequence coverage to display a domain (0-1)
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
    except ImportError:
        print(f"Warning: matplotlib not available, skipping cross-method consensus visualization", file=sys.stderr)
        return

    if not consensus_results:
        print("Warning: No consensus results to visualize", file=sys.stderr)
        return

    # Filter by minimum coverage
    filtered_results = [r for r in consensus_results
                        if calculate_sequence_coverage(r.supporting_domains, r.total_sequences) >= min_coverage]

    if not filtered_results:
        print(f"Warning: No consensus domains meet minimum coverage threshold ({min_coverage:.0%})", file=sys.stderr)
        return

    # Sort by alignment start position
    filtered_results = sorted(filtered_results, key=lambda r: r.ali_start)

    # Color scheme - color by coverage level
    def get_coverage_color(coverage):
        """Return color based on coverage: green (high) -> yellow (medium) -> red (low)"""
        if coverage >= 0.7:
            return '#2E7D32'  # Dark green
        elif coverage >= 0.5:
            return '#4CAF50'  # Green
        elif coverage >= 0.3:
            return '#FFC107'  # Amber
        elif coverage >= 0.1:
            return '#FF9800'  # Orange
        else:
            return '#F44336'  # Red

    # Calculate dimensions
    row_height = 25
    margin_left = 180
    margin_right = 150
    margin_top = 90
    margin_bottom = 30

    # Use alignment coordinates
    ali_length = len(entries[0].aligned_seq) if entries else 0
    pixels_per_unit = 2
    max_length = ali_length

    fig_width = margin_left + (max_length * pixels_per_unit) + margin_right
    fig_height = margin_top + (len(entries) * row_height) + margin_bottom

    # Create figure
    fig, ax = plt.subplots(figsize=(fig_width / 100, fig_height / 100), dpi=100)
    ax.set_xlim(0, fig_width)
    ax.set_ylim(0, fig_height)
    ax.axis('off')

    # Add title
    ax.text(fig_width / 2, fig_height - 12, title,
            ha='center', va='top', fontsize=11, weight='bold')

    # Add legend for consensus domains
    legend_y = fig_height - 35
    legend_x = margin_left
    ax.text(legend_x - 10, legend_y, "Consensus domains (colored by coverage):", fontsize=8, va='center')

    for i, result in enumerate(filtered_results[:8]):  # Show up to 8 in legend
        coverage = calculate_sequence_coverage(result.supporting_domains, result.total_sequences)
        color = get_coverage_color(coverage)
        x = legend_x + 220 + i * 100
        legend_box = patches.Rectangle(
            (x, legend_y - 4), 12, 8,
            linewidth=1, edgecolor='black', facecolor=color, alpha=0.8
        )
        ax.add_patch(legend_box)
        ax.text(x + 16, legend_y, f"D{i+1} ({coverage:.0%})", fontsize=7, va='center')

    # Add second legend row for coverage scale
    legend_y2 = fig_height - 50
    ax.text(legend_x - 10, legend_y2, "Coverage scale:", fontsize=7, va='center')
    coverage_levels = [(0.7, '≥70%'), (0.5, '50-70%'), (0.3, '30-50%'), (0.1, '10-30%'), (0.05, '<10%')]
    for i, (cov, label) in enumerate(coverage_levels):
        x = legend_x + 80 + i * 90
        color = get_coverage_color(cov)
        legend_box = patches.Rectangle(
            (x, legend_y2 - 4), 12, 8,
            linewidth=1, edgecolor='black', facecolor=color, alpha=0.8
        )
        ax.add_patch(legend_box)
        ax.text(x + 16, legend_y2, label, fontsize=7, va='center')

    # Add axis labels
    ax.text(margin_left, fig_height - margin_top + 10, "1",
            ha='left', va='bottom', fontsize=8)
    ax.text(margin_left + max_length * pixels_per_unit, fig_height - margin_top + 10,
            str(ali_length), ha='right', va='bottom', fontsize=8)
    ax.text(margin_left + (max_length * pixels_per_unit) / 2, fig_height - margin_top + 10,
            "Alignment Position", ha='center', va='bottom', fontsize=9)

    # Draw each sequence
    current_y = fig_height - margin_top - row_height / 2

    for entry in entries:
        # Draw label
        label = f"{entry.accession}"
        ax.text(margin_left - 5, current_y, label,
                ha='right', va='center', fontsize=7, family='monospace')

        # Draw backbone (alignment length)
        backbone_length = len(entry.aligned_seq)
        backbone_height = 6
        backbone = patches.Rectangle(
            (margin_left, current_y - backbone_height / 2),
            backbone_length * pixels_per_unit, backbone_height,
            linewidth=0, facecolor='#CCCCCC'
        )
        ax.add_patch(backbone)

        # Draw consensus domains
        for i, result in enumerate(filtered_results):
            coverage = calculate_sequence_coverage(result.supporting_domains, result.total_sequences)
            color = get_coverage_color(coverage)

            domain_x = margin_left + (result.ali_start - 1) * pixels_per_unit
            domain_width = (result.ali_end - result.ali_start + 1) * pixels_per_unit

            domain_height = 16
            domain_rect = patches.Rectangle(
                (domain_x, current_y - domain_height / 2),
                domain_width, domain_height,
                linewidth=1, edgecolor='black', facecolor=color, alpha=0.8
            )
            ax.add_patch(domain_rect)

            # Add domain label if space permits
            if domain_width > 30:
                label_text = f"D{i + 1}"
                ax.text(domain_x + domain_width / 2, current_y,
                        label_text, ha='center', va='center',
                        fontsize=6, weight='bold', color='white')

        # Draw SEED region (mapped to alignment coordinates)
        seq_to_ali = build_seq_to_ali_map(entry.aligned_seq, entry.region_start)
        if entry.region_start in seq_to_ali and entry.region_end in seq_to_ali:
            seed_ali_start = seq_to_ali[entry.region_start]
            seed_ali_end = seq_to_ali[entry.region_end]
            seed_x = margin_left + (seed_ali_start - 1) * pixels_per_unit
            seed_width = (seed_ali_end - seed_ali_start + 1) * pixels_per_unit
            seed_height = 22
            seed_rect = patches.Rectangle(
                (seed_x, current_y - seed_height / 2),
                seed_width, seed_height,
                linewidth=2, edgecolor='black', facecolor='none'
            )
            ax.add_patch(seed_rect)

        current_y -= row_height

    # Add summary info
    info_text = f"Cross-method consensus: {len(filtered_results)} domains shown"
    ax.text(fig_width / 2, 10, info_text, ha='center', va='bottom', fontsize=8, style='italic')

    # Save figure
    plt.tight_layout()
    plt.savefig(output_file, dpi=100, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Created cross-method consensus visualization: {output_file}", file=sys.stderr)


# ============================================================================
# Mock Data for Testing
# ============================================================================

def generate_mock_predictions(entries: List[SequenceEntry]) -> Dict[str, ProteinPredictions]:
    """
    Generate mock TED predictions for testing without network access.

    Creates realistic-looking domain predictions based on sequence regions.
    Includes all four methods: chainsaw, merizo, unidoc-ndr, and ted_consensus.
    """
    import random
    random.seed(42)  # Reproducible mock data

    predictions = {}
    methods = ['chainsaw', 'merizo', 'unidoc-ndr', 'ted_consensus']

    for entry in entries:
        # Create predictions for each method with slight variations
        method_predictions = {}

        # Estimate a domain that covers most of the aligned region
        base_start = entry.region_start
        base_end = entry.region_end
        region_len = base_end - base_start + 1

        for method in methods:
            # Add random variation to boundaries (±5 residues)
            # ted_consensus has less variation (more consistent)
            if method == 'ted_consensus':
                variation = random.randint(-2, 2)
            else:
                variation = random.randint(-5, 5)

            start = max(1, base_start + variation)
            end = base_end + random.randint(-3, 3)

            if end > start:
                segments = [Segment(start, end)]
                domain = Domain(
                    segments=segments,
                    method=method,
                    score=random.uniform(0.5, 1.0),
                    domain_idx=0
                )

                method_predictions[method] = MethodPrediction(
                    method=method,
                    uniprot_acc=entry.base_accession,
                    domains=[domain],
                    nres_chain=region_len + 50,  # Assume protein is longer
                    score=random.uniform(0.5, 1.0)
                )

        if method_predictions:
            predictions[entry.base_accession] = ProteinPredictions(
                uniprot_acc=entry.base_accession,
                predictions=method_predictions
            )

    return predictions


# ============================================================================
# Main Entry Point
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="TED Family Consensus - Calculate consensus domain boundaries from multiple TED methods",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python ted_family_consensus.py SEED
    python ted_family_consensus.py SEED -o report.txt -j results.json
    python ted_family_consensus.py SEED --strategy combined --tolerance 15
    python ted_family_consensus.py SEED --all-strategies
    python ted_family_consensus.py SEED --per-method domains  # Creates domains_chainsaw.tsv, etc.
    python ted_family_consensus.py SEED --images family       # Raw prediction images
    python ted_family_consensus.py SEED --cross-method-image cross_consensus.png  # Cross-method consensus

Output Images:
    --images PREFIX creates raw prediction images showing domains as predicted by each method
    --consensus-images PREFIX creates per-method consensus images (clustering within each method)
    --cross-method-image FILE creates the cross-method consensus image showing domains
                              computed by clustering predictions from ALL methods together
                              (this visualizes the "CONSENSUS RESULTS BY STRATEGY" table)

Consensus Strategies:
    majority         - Score by fraction of sequences with the domain
    method_agreement - Score by number of methods agreeing
    combined         - Weighted combination of coverage and method agreement
    any_two_methods  - Require at least 2 methods to agree
    all_three_methods - Require all 3 methods to agree
    weighted_method  - Weight methods by reliability (customizable)

Method Consistency Scores:
    The report includes a consistency score (0-1) for each method measuring how
    uniformly it predicts domains across all sequences in the family. Higher
    scores indicate the method is more reliable for this particular family.
        """
    )

    parser.add_argument("seed_file", help="Path to SEED alignment file (MUL format)")
    parser.add_argument("-o", "--output", help="Output file for summary report")
    parser.add_argument("-s", "--seq-output", help="Output file for per-sequence consensus domains")
    parser.add_argument("-j", "--json-output", help="Output file for detailed JSON results")
    parser.add_argument("-t", "--tolerance", type=int, default=10,
                        help="Tolerance for boundary matching (default: 10 residues)")
    parser.add_argument("-m", "--min-score", type=float, default=0.1,
                        help="Minimum score threshold for consensus domains (default: 0.1)")
    parser.add_argument("--strategy", choices=['majority', 'method_agreement', 'combined',
                                                'any_two_methods', 'all_three_methods', 'weighted_method'],
                        default='combined',
                        help="Consensus strategy to use (default: combined)")
    parser.add_argument("--all-strategies", action="store_true",
                        help="Run all consensus strategies and compare results")
    parser.add_argument("--per-method", metavar="PREFIX",
                        help="Generate per-method TSV files (PREFIX_chainsaw.tsv, PREFIX_merizo.tsv, PREFIX_unidoc-ndr.tsv)")
    parser.add_argument("--images", metavar="PREFIX",
                        help="Generate visualization images (PREFIX_chainsaw.png, PREFIX_merizo.png, PREFIX_unidoc-ndr.png, PREFIX_consensus.png)")
    parser.add_argument("--consensus-images", metavar="PREFIX",
                        help="Generate computed consensus domain images for each method (PREFIX_chainsaw_computed_consensus.png, etc.)")
    parser.add_argument("--cross-method-image", metavar="FILE",
                        help="Generate cross-method consensus image showing domains from the combined strategy (alignment view)")
    parser.add_argument("--view-mode", choices=['sequence', 'alignment'], default='sequence',
                        help="Image view mode: 'sequence' (proteins at actual length with SEED box) or 'alignment' (domains mapped to alignment coordinates)")
    parser.add_argument("--mock", action="store_true",
                        help="Use mock data for testing (no API calls)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")

    args = parser.parse_args()

    # Parse SEED alignment
    if args.verbose:
        print(f"Parsing SEED alignment: {args.seed_file}", file=sys.stderr)

    entries = parse_seed_alignment(args.seed_file)

    if not entries:
        print("Error: No sequences found in SEED alignment", file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        print(f"Found {len(entries)} sequences in alignment", file=sys.stderr)

    # Fetch TED predictions (or use mock data)
    if args.mock:
        if args.verbose:
            print("Using mock data (--mock specified)...", file=sys.stderr)
        predictions = generate_mock_predictions(entries)
    else:
        if args.verbose:
            print("Fetching TED predictions from API...", file=sys.stderr)
        predictions = fetch_all_predictions(entries, verbose=args.verbose)

    if not predictions:
        print("Error: No TED predictions found for any sequence", file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        print(f"Retrieved predictions for {len(predictions)} sequences", file=sys.stderr)

    # Set up strategies
    if args.all_strategies:
        strategies = [
            MajorityVoteStrategy(),
            MethodAgreementStrategy(),
            CombinedStrategy(),
            AnyTwoMethodsStrategy(),
            AllThreeMethodsStrategy(),
            WeightedMethodStrategy()
        ]
    else:
        strategy_map = {
            'majority': MajorityVoteStrategy(),
            'method_agreement': MethodAgreementStrategy(),
            'combined': CombinedStrategy(),
            'any_two_methods': AnyTwoMethodsStrategy(),
            'all_three_methods': AllThreeMethodsStrategy(),
            'weighted_method': WeightedMethodStrategy()
        }
        strategies = [strategy_map[args.strategy]]

    # Run consensus analysis
    if args.verbose:
        print("Running consensus analysis...", file=sys.stderr)

    results = run_consensus_analysis(
        entries, predictions, strategies,
        tolerance=args.tolerance,
        min_score=args.min_score
    )

    # Calculate method consistency scores
    if args.verbose:
        print("Calculating method consistency scores...", file=sys.stderr)
    consistency_scores = calculate_method_consistency(entries, predictions, args.tolerance)

    # Generate reports
    summary_report = generate_summary_report(
        entries, predictions, results, args.output,
        consistency_scores=consistency_scores, tolerance=args.tolerance
    )

    if not args.output:
        print(summary_report)

    # Per-sequence output
    if args.seq_output:
        # Use the primary strategy's results
        primary_strategy = strategies[0].name
        primary_results = results.get(primary_strategy, [])
        per_seq = get_per_sequence_consensus(entries, predictions, primary_results)
        generate_per_sequence_report(entries, per_seq, args.seq_output)
        if args.verbose:
            print(f"Per-sequence report written to: {args.seq_output}", file=sys.stderr)

    # JSON output
    if args.json_output:
        generate_detailed_json(entries, predictions, results, args.json_output)
        if args.verbose:
            print(f"JSON results written to: {args.json_output}", file=sys.stderr)

    # Per-method TSV files
    if args.per_method:
        if args.verbose:
            print("Generating per-method domain reports...", file=sys.stderr)
        method_files = generate_per_method_reports(entries, predictions, args.per_method)
        for method, filepath in method_files.items():
            if args.verbose:
                print(f"  {method}: {filepath}", file=sys.stderr)

    # Visualization images
    if args.images:
        if args.verbose:
            print(f"Generating visualization images (view mode: {args.view_mode})...", file=sys.stderr)

        # Create per-method images
        image_files = create_all_method_visualizations(entries, predictions, args.images, args.view_mode)
        for method, filepath in image_files.items():
            if args.verbose:
                print(f"  {method}: {filepath}", file=sys.stderr)

        # Create consensus comparison image
        consensus_file = f"{args.images}_consensus.png"
        create_consensus_visualization(entries, predictions, consensus_file, view_mode=args.view_mode)
        if args.verbose:
            print(f"  consensus: {consensus_file}", file=sys.stderr)

    # Computed consensus images (per-method)
    if args.consensus_images:
        if args.verbose:
            print(f"Generating computed consensus domain images...", file=sys.stderr)

        consensus_image_files = create_all_computed_consensus_visualizations(
            entries, predictions, args.consensus_images, tolerance=args.tolerance
        )
        for method, filepath in consensus_image_files.items():
            if args.verbose:
                print(f"  {method}: {filepath}", file=sys.stderr)

    # Cross-method consensus image
    if args.cross_method_image:
        if args.verbose:
            print(f"Generating cross-method consensus image...", file=sys.stderr)

        # Get the primary strategy results (usually 'combined')
        primary_strategy = strategies[0].name
        primary_results = results.get(primary_strategy, [])

        create_cross_method_consensus_visualization(
            entries, primary_results, args.cross_method_image,
            title=f"Cross-Method Consensus Domains ({primary_strategy} strategy)"
        )

    if args.verbose:
        print("Analysis complete!", file=sys.stderr)


if __name__ == "__main__":
    main()
