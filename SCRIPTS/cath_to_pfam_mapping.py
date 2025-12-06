#!/usr/bin/env python3
"""
Calculate CATH to Pfam clan mapping based on domain overlaps.

This script takes sorted TED domain data (with CATH labels) and sorted Pfam
domain regions, finds overlapping domains on the same UniProt sequences,
and produces a mapping of CATH classifications to Pfam clans.

Prerequisites:
    1. ted_cath_high_sorted.tsv.gz - TED domains with CATH labels, sorted by UniProt
       Format: uniprot_acc, ted_suffix, start, end, cath_label, cath_level

    2. pfam_clan_regions_sorted.tsv.gz - Pfam domain regions for clan families, sorted by UniProt
       Format: pfamseq_acc, seq_start, seq_end, pfamA_acc

    3. pfam_clan_membership.tsv - Family to clan mapping
       Format: pfamA_acc, clan_acc, clan_id

Usage:
    python cath_to_pfam_mapping.py [options]

Output:
    - cath_to_pfam_clan_mapping.tsv: Main mapping file with confidence scores
    - cath_to_pfam_family_counts.tsv: Detailed per-family counts
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict
from typing import Dict, List, Tuple, Optional


# Default paths
DEFAULT_TED_FILE = '/nfs/production/agb/pfam/data/TED/ted_cath_high_sorted.tsv.gz'
DEFAULT_PFAM_FILE = '/nfs/production/agb/pfam/data/TED/pfam_clan_regions_sorted.tsv.gz'
DEFAULT_CLAN_FILE = '/nfs/production/agb/pfam/data/TED/pfam_clan_membership.tsv'
DEFAULT_OUTPUT_DIR = '/nfs/production/agb/pfam/data/TED'


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Calculate CATH to Pfam clan mapping based on domain overlaps'
    )
    parser.add_argument(
        '--ted-file',
        default=DEFAULT_TED_FILE,
        help=f'Path to sorted TED CATH file (default: {DEFAULT_TED_FILE})'
    )
    parser.add_argument(
        '--pfam-file',
        default=DEFAULT_PFAM_FILE,
        help=f'Path to sorted Pfam regions file (default: {DEFAULT_PFAM_FILE})'
    )
    parser.add_argument(
        '--clan-file',
        default=DEFAULT_CLAN_FILE,
        help=f'Path to clan membership file (default: {DEFAULT_CLAN_FILE})'
    )
    parser.add_argument(
        '--output-dir',
        default=DEFAULT_OUTPUT_DIR,
        help=f'Output directory (default: {DEFAULT_OUTPUT_DIR})'
    )
    parser.add_argument(
        '--min-overlap',
        type=float,
        default=0.5,
        help='Minimum overlap fraction in both directions (default: 0.5)'
    )
    parser.add_argument(
        '--progress-interval',
        type=int,
        default=1000000,
        help='Print progress every N accessions (default: 1000000)'
    )
    return parser.parse_args()


def load_clan_membership(clan_file: str) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Load clan membership mapping.

    Returns:
        family_to_clan: Dict mapping pfamA_acc -> clan_acc
        clan_to_id: Dict mapping clan_acc -> clan_id (name)
    """
    family_to_clan = {}
    clan_to_id = {}

    with open(clan_file, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                pfamA_acc, clan_acc, clan_id = parts[0], parts[1], parts[2]
                family_to_clan[pfamA_acc] = clan_acc
                clan_to_id[clan_acc] = clan_id

    print(f"Loaded {len(family_to_clan):,} family->clan mappings")
    print(f"Loaded {len(clan_to_id):,} unique clans")

    return family_to_clan, clan_to_id


def open_file(filepath: str):
    """Open a file, handling gzip if needed."""
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt')
    return open(filepath, 'r')


def parse_ted_line(line: str) -> Optional[Tuple[str, str, int, int, str, str]]:
    """
    Parse a TED domain line.

    Returns: (uniprot_acc, ted_suffix, start, end, cath_label, cath_level) or None
    """
    parts = line.rstrip('\n').split('\t')
    if len(parts) < 6:
        return None
    try:
        return (
            parts[0],           # uniprot_acc
            parts[1],           # ted_suffix
            int(parts[2]),      # start
            int(parts[3]),      # end
            parts[4],           # cath_label
            parts[5]            # cath_level
        )
    except (ValueError, IndexError):
        return None


def parse_pfam_line(line: str) -> Optional[Tuple[str, int, int, str]]:
    """
    Parse a Pfam domain line.

    Returns: (uniprot_acc, start, end, pfamA_acc) or None
    """
    parts = line.rstrip('\n').split('\t')
    if len(parts) < 4:
        return None
    try:
        return (
            parts[0],           # uniprot_acc
            int(parts[1]),      # start
            int(parts[2]),      # end
            parts[3]            # pfamA_acc
        )
    except (ValueError, IndexError):
        return None


def calculate_overlap(start1: int, end1: int, start2: int, end2: int) -> Tuple[float, float]:
    """
    Calculate reciprocal overlap fractions between two regions.

    Returns: (overlap_fraction_1, overlap_fraction_2)
        overlap_fraction_1: fraction of region 1 covered by overlap
        overlap_fraction_2: fraction of region 2 covered by overlap
    """
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)

    if overlap_start > overlap_end:
        return 0.0, 0.0

    overlap_length = overlap_end - overlap_start + 1
    length1 = end1 - start1 + 1
    length2 = end2 - start2 + 1

    return overlap_length / length1, overlap_length / length2


def find_overlapping_pairs(
    ted_domains: List[Tuple[int, int, str, str]],
    pfam_domains: List[Tuple[int, int, str]],
    min_overlap: float
) -> List[Tuple[str, str, str]]:
    """
    Find all pairs of overlapping TED and Pfam domains.

    Args:
        ted_domains: List of (start, end, cath_label, cath_level)
        pfam_domains: List of (start, end, pfamA_acc)
        min_overlap: Minimum overlap fraction required in both directions

    Returns:
        List of (cath_label, cath_level, pfamA_acc) for overlapping pairs
    """
    pairs = []

    for ted_start, ted_end, cath_label, cath_level in ted_domains:
        for pfam_start, pfam_end, pfamA_acc in pfam_domains:
            overlap_ted, overlap_pfam = calculate_overlap(
                ted_start, ted_end, pfam_start, pfam_end
            )

            if overlap_ted >= min_overlap and overlap_pfam >= min_overlap:
                pairs.append((cath_label, cath_level, pfamA_acc))

    return pairs


def process_files(
    ted_file: str,
    pfam_file: str,
    min_overlap: float,
    progress_interval: int
) -> Dict[Tuple[str, str, str], int]:
    """
    Process both sorted files and count overlapping domain pairs.

    Uses merge-join approach to stream both files with minimal memory.

    Returns:
        counts: Dict mapping (cath_label, cath_level, pfamA_acc) -> count
    """
    counts = defaultdict(int)

    ted_f = open_file(ted_file)
    pfam_f = open_file(pfam_file)

    # Current lines
    ted_line = ted_f.readline()
    pfam_line = pfam_f.readline()

    # Statistics
    accessions_processed = 0
    accessions_with_overlaps = 0
    total_overlaps = 0

    def get_acc(line: str) -> Optional[str]:
        """Extract accession from line."""
        if not line:
            return None
        return line.split('\t')[0]

    while ted_line or pfam_line:
        ted_acc = get_acc(ted_line)
        pfam_acc = get_acc(pfam_line)

        # Determine which accession to process
        if ted_acc is None:
            current_acc = pfam_acc
        elif pfam_acc is None:
            current_acc = ted_acc
        else:
            current_acc = min(ted_acc, pfam_acc)

        if current_acc is None:
            break

        # Collect all TED domains for current accession
        ted_domains = []
        while ted_line and get_acc(ted_line) == current_acc:
            parsed = parse_ted_line(ted_line)
            if parsed:
                _, _, start, end, cath_label, cath_level = parsed
                ted_domains.append((start, end, cath_label, cath_level))
            ted_line = ted_f.readline()

        # Collect all Pfam domains for current accession
        pfam_domains = []
        while pfam_line and get_acc(pfam_line) == current_acc:
            parsed = parse_pfam_line(pfam_line)
            if parsed:
                _, start, end, pfamA_acc = parsed
                pfam_domains.append((start, end, pfamA_acc))
            pfam_line = pfam_f.readline()

        # Find overlapping pairs if we have both
        if ted_domains and pfam_domains:
            pairs = find_overlapping_pairs(ted_domains, pfam_domains, min_overlap)

            if pairs:
                accessions_with_overlaps += 1
                total_overlaps += len(pairs)

                for cath_label, cath_level, pfamA_acc in pairs:
                    counts[(cath_label, cath_level, pfamA_acc)] += 1

        accessions_processed += 1

        if accessions_processed % progress_interval == 0:
            print(f"  Processed {accessions_processed:,} accessions, "
                  f"{accessions_with_overlaps:,} with overlaps, "
                  f"{total_overlaps:,} total overlaps, "
                  f"{len(counts):,} unique pairs...")

    ted_f.close()
    pfam_f.close()

    print(f"\nProcessing complete:")
    print(f"  Accessions processed: {accessions_processed:,}")
    print(f"  Accessions with overlaps: {accessions_with_overlaps:,}")
    print(f"  Total overlapping domain pairs: {total_overlaps:,}")
    print(f"  Unique (CATH, Pfam family) pairs: {len(counts):,}")

    return counts


def aggregate_to_clans(
    counts: Dict[Tuple[str, str, str], int],
    family_to_clan: Dict[str, str],
    clan_to_id: Dict[str, str]
) -> Tuple[Dict[Tuple[str, str, str], int], Dict[Tuple[str, str, str], List[Tuple[str, int]]]]:
    """
    Aggregate family-level counts to clan level.

    Returns:
        clan_counts: Dict mapping (cath_label, cath_level, clan_acc) -> total_count
        clan_families: Dict mapping (cath_label, cath_level, clan_acc) -> [(pfamA_acc, count), ...]
    """
    clan_counts = defaultdict(int)
    clan_families = defaultdict(list)

    for (cath_label, cath_level, pfamA_acc), count in counts.items():
        clan_acc = family_to_clan.get(pfamA_acc)
        if clan_acc:
            key = (cath_label, cath_level, clan_acc)
            clan_counts[key] += count
            clan_families[key].append((pfamA_acc, count))

    # Sort families by count within each clan
    for key in clan_families:
        clan_families[key].sort(key=lambda x: -x[1])

    return clan_counts, clan_families


def calculate_confidence(
    clan_counts: Dict[Tuple[str, str, str], int]
) -> Dict[Tuple[str, str], int]:
    """
    Calculate total counts per CATH label for confidence calculation.

    Returns:
        cath_totals: Dict mapping (cath_label, cath_level) -> total_count
    """
    cath_totals = defaultdict(int)

    for (cath_label, cath_level, clan_acc), count in clan_counts.items():
        cath_totals[(cath_label, cath_level)] += count

    return cath_totals


def write_outputs(
    counts: Dict[Tuple[str, str, str], int],
    clan_counts: Dict[Tuple[str, str, str], int],
    clan_families: Dict[Tuple[str, str, str], List[Tuple[str, int]]],
    cath_totals: Dict[Tuple[str, str], int],
    clan_to_id: Dict[str, str],
    output_dir: str
):
    """Write output files."""

    # Main mapping file
    mapping_file = os.path.join(output_dir, 'cath_to_pfam_clan_mapping.tsv')
    print(f"\nWriting clan mapping to {mapping_file}")

    with open(mapping_file, 'w') as f:
        # Header
        f.write("cath_label\tcath_level\tclan_acc\tclan_id\toverlap_count\tconfidence\ttop_families\n")

        # Sort by CATH label, then by count descending
        sorted_items = sorted(
            clan_counts.items(),
            key=lambda x: (x[0][0], x[0][1], -x[1])
        )

        for (cath_label, cath_level, clan_acc), count in sorted_items:
            clan_id = clan_to_id.get(clan_acc, '')
            total = cath_totals.get((cath_label, cath_level), count)
            confidence = count / total if total > 0 else 0.0

            # Get top 3 families
            families = clan_families.get((cath_label, cath_level, clan_acc), [])
            top_families = ','.join(fam for fam, _ in families[:3])

            f.write(f"{cath_label}\t{cath_level}\t{clan_acc}\t{clan_id}\t"
                    f"{count}\t{confidence:.4f}\t{top_families}\n")

    # Detailed family counts file
    family_file = os.path.join(output_dir, 'cath_to_pfam_family_counts.tsv')
    print(f"Writing family counts to {family_file}")

    with open(family_file, 'w') as f:
        # Header
        f.write("cath_label\tcath_level\tpfamA_acc\tcount\n")

        # Sort by CATH label, then by count descending
        sorted_items = sorted(
            counts.items(),
            key=lambda x: (x[0][0], x[0][1], -x[1])
        )

        for (cath_label, cath_level, pfamA_acc), count in sorted_items:
            f.write(f"{cath_label}\t{cath_level}\t{pfamA_acc}\t{count}\n")

    print(f"\nOutput files written:")
    print(f"  {mapping_file}")
    print(f"  {family_file}")


def main():
    args = parse_arguments()

    print("=" * 70)
    print("CATH to Pfam Clan Mapping")
    print("=" * 70)
    print(f"\nInput files:")
    print(f"  TED file:  {args.ted_file}")
    print(f"  Pfam file: {args.pfam_file}")
    print(f"  Clan file: {args.clan_file}")
    print(f"\nParameters:")
    print(f"  Min overlap: {args.min_overlap:.0%} (both directions)")
    print(f"  Output dir:  {args.output_dir}")

    # Check input files exist
    for filepath in [args.ted_file, args.pfam_file, args.clan_file]:
        if not os.path.exists(filepath):
            print(f"\nError: File not found: {filepath}")
            sys.exit(1)

    # Load clan membership
    print("\n" + "-" * 70)
    print("Loading clan membership...")
    family_to_clan, clan_to_id = load_clan_membership(args.clan_file)

    # Process files and count overlaps
    print("\n" + "-" * 70)
    print("Processing domain overlaps...")
    counts = process_files(
        args.ted_file,
        args.pfam_file,
        args.min_overlap,
        args.progress_interval
    )

    # Aggregate to clans
    print("\n" + "-" * 70)
    print("Aggregating to clan level...")
    clan_counts, clan_families = aggregate_to_clans(counts, family_to_clan, clan_to_id)
    cath_totals = calculate_confidence(clan_counts)

    print(f"  Unique (CATH, Clan) pairs: {len(clan_counts):,}")
    print(f"  Unique CATH labels with clan mappings: {len(cath_totals):,}")

    # Write outputs
    print("\n" + "-" * 70)
    write_outputs(
        counts,
        clan_counts,
        clan_families,
        cath_totals,
        clan_to_id,
        args.output_dir
    )

    print("\n" + "=" * 70)
    print("Done!")
    print("=" * 70)


if __name__ == '__main__':
    main()
