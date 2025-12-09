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

    Returns list of SequenceEntry objects.
    """
    # First pass: collect all lines grouped by sequence
    seq_data = {}  # accession -> {'region_start', 'region_end', 'seq_parts': []}
    seq_order = []  # Preserve order of first occurrence

    # Pattern to match sequence header lines: ACC.version/start-end  SEQUENCE
    seq_pattern = re.compile(r'^(\S+)/(\d+)-(\d+)\s+(\S+)$')

    # Pattern to match continuation lines (just sequence, possibly with leading whitespace)
    cont_pattern = re.compile(r'^\s*([A-Za-z.]+)$')

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
                          retry_delay: float = 1.0) -> Optional[ProteinPredictions]:
    """
    Fetch domain predictions from all three TED methods for a UniProt accession.

    Uses the chainparse endpoint which returns chainsaw, merizo, and unidoc-ndr predictions.
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

            for item in data.get('data', []):
                method = item.get('method', '')
                chopping = item.get('chopping', '')
                score = item.get('score', 0.0)
                nres = item.get('nres_chain', 0)

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
                             output_file: str = None) -> str:
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
    for method in ['chainsaw', 'merizo', 'unidoc-ndr']:
        count = method_counts.get(method, 0)
        domains = total_domains_per_method.get(method, 0)
        lines.append(f"  {method:15s}: {count:4d} proteins, {domains:4d} total domains")
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
# Mock Data for Testing
# ============================================================================

def generate_mock_predictions(entries: List[SequenceEntry]) -> Dict[str, ProteinPredictions]:
    """
    Generate mock TED predictions for testing without network access.

    Creates realistic-looking domain predictions based on sequence regions.
    """
    import random
    random.seed(42)  # Reproducible mock data

    predictions = {}
    methods = ['chainsaw', 'merizo', 'unidoc-ndr']

    for entry in entries:
        # Create predictions for each method with slight variations
        method_predictions = {}

        # Estimate a domain that covers most of the aligned region
        base_start = entry.region_start
        base_end = entry.region_end
        region_len = base_end - base_start + 1

        for method in methods:
            # Add random variation to boundaries (±5 residues)
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

Consensus Strategies:
    majority         - Score by fraction of sequences with the domain
    method_agreement - Score by number of methods agreeing
    combined         - Weighted combination of coverage and method agreement
    any_two_methods  - Require at least 2 methods to agree
    all_three_methods - Require all 3 methods to agree
    weighted_method  - Weight methods by reliability (customizable)
        """
    )

    parser.add_argument("seed_file", help="Path to SEED alignment file (MUL format)")
    parser.add_argument("-o", "--output", help="Output file for summary report")
    parser.add_argument("-s", "--seq-output", help="Output file for per-sequence domains")
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

    # Generate reports
    summary_report = generate_summary_report(entries, predictions, results, args.output)

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

    if args.verbose:
        print("Analysis complete!", file=sys.stderr)


if __name__ == "__main__":
    main()
