#!/usr/bin/env python3
"""
Identify Pfam families that map to multiple CATH H-level superfamilies.

Uses only high-confidence CATH->Clan mappings (confidence > 0.8, count >= 10)
to identify families that genuinely span multiple CATH superfamilies.

These families may need review as they could:
- Be too broad and need splitting
- Contain multiple structural domains
- Have curation issues

Usage:
    python pfam_multi_cath_analysis.py [options]

Input:
    cath_to_pfam_clan_mapping.tsv - Clan-level mappings with confidence scores
    cath_to_pfam_family_counts.tsv - Detailed family counts

Output:
    pfam_multi_cath_families.tsv - Families mapping to multiple CATH superfamilies
"""

import argparse
import os
import sys
from collections import defaultdict
from typing import Dict, List, Set, Tuple


# Default paths
DEFAULT_CLAN_MAPPING = '/nfs/production/agb/pfam/data/TED/cath_to_pfam_clan_mapping.tsv'
DEFAULT_FAMILY_COUNTS = '/nfs/production/agb/pfam/data/TED/cath_to_pfam_family_counts.tsv'
DEFAULT_OUTPUT = '/nfs/production/agb/pfam/data/TED/pfam_multi_cath_families.tsv'
DEFAULT_CLAN_FILE = '/nfs/production/agb/pfam/data/TED/pfam_clan_membership.tsv'


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Identify Pfam families mapping to multiple CATH superfamilies'
    )
    parser.add_argument(
        '--clan-mapping',
        default=DEFAULT_CLAN_MAPPING,
        help=f'Path to clan mapping file (default: {DEFAULT_CLAN_MAPPING})'
    )
    parser.add_argument(
        '--family-counts',
        default=DEFAULT_FAMILY_COUNTS,
        help=f'Path to family counts file (default: {DEFAULT_FAMILY_COUNTS})'
    )
    parser.add_argument(
        '--output',
        default=DEFAULT_OUTPUT,
        help=f'Output file path (default: {DEFAULT_OUTPUT})'
    )
    parser.add_argument(
        '--clan-file',
        default=DEFAULT_CLAN_FILE,
        help=f'Path to clan membership file (default: {DEFAULT_CLAN_FILE})'
    )
    parser.add_argument(
        '--min-confidence',
        type=float,
        default=0.8,
        help='Minimum confidence for CATH->Clan mapping (default: 0.8)'
    )
    parser.add_argument(
        '--min-count',
        type=int,
        default=10,
        help='Minimum overlap count for CATH->Clan mapping (default: 10)'
    )
    parser.add_argument(
        '--min-cath-groups',
        type=int,
        default=2,
        help='Minimum number of CATH H-groups to report (default: 2)'
    )
    return parser.parse_args()


def load_clan_membership(clan_file: str) -> Dict[str, Tuple[str, str]]:
    """
    Load clan membership mapping.

    Returns:
        Dict mapping pfamA_acc -> (clan_acc, clan_id)
    """
    family_to_clan = {}

    if not os.path.exists(clan_file):
        print(f"Warning: Clan file not found: {clan_file}")
        return family_to_clan

    with open(clan_file, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                pfamA_acc, clan_acc, clan_id = parts[0], parts[1], parts[2]
                family_to_clan[pfamA_acc] = (clan_acc, clan_id)

    return family_to_clan


def load_confident_cath_mappings(
    clan_mapping_file: str,
    min_confidence: float,
    min_count: int
) -> Set[Tuple[str, str]]:
    """
    Load high-confidence CATH -> Clan mappings.

    Returns:
        Set of (cath_label, clan_acc) tuples that meet confidence/count thresholds
    """
    confident_mappings = set()

    with open(clan_mapping_file, 'r') as f:
        # Skip header
        header = f.readline()

        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6:
                continue

            cath_label = parts[0]
            cath_level = parts[1]
            clan_acc = parts[2]
            overlap_count = int(parts[4])
            confidence = float(parts[5])

            # Only H-level (superfamily) mappings
            if cath_level != 'H':
                continue

            # Filter by confidence and count
            if confidence >= min_confidence and overlap_count >= min_count:
                confident_mappings.add((cath_label, clan_acc))

    return confident_mappings


def load_family_counts_for_confident_mappings(
    family_counts_file: str,
    confident_mappings: Set[Tuple[str, str]],
    family_to_clan: Dict[str, Tuple[str, str]]
) -> Dict[str, Dict[str, int]]:
    """
    Load family counts, but only for CATH labels that have confident clan mappings.

    Returns:
        Dict mapping pfamA_acc -> {cath_label: count}
    """
    # Build set of confident CATH labels per clan
    confident_cath_by_clan = defaultdict(set)
    for cath_label, clan_acc in confident_mappings:
        confident_cath_by_clan[clan_acc].add(cath_label)

    family_cath = defaultdict(dict)

    with open(family_counts_file, 'r') as f:
        # Skip header
        header = f.readline()

        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue

            cath_label, cath_level, pfamA_acc, count = parts[0], parts[1], parts[2], int(parts[3])

            # Only H-level
            if cath_level != 'H':
                continue

            # Get clan for this family
            clan_info = family_to_clan.get(pfamA_acc)
            if not clan_info:
                continue

            clan_acc = clan_info[0]

            # Check if this CATH label has a confident mapping to this family's clan
            if (cath_label, clan_acc) in confident_mappings:
                family_cath[pfamA_acc][cath_label] = count

    return family_cath


def analyze_multi_cath_families(
    family_cath: Dict[str, Dict[str, int]],
    family_to_clan: Dict[str, Tuple[str, str]],
    min_cath_groups: int
) -> List[Dict]:
    """
    Identify families mapping to multiple CATH superfamilies.

    Returns:
        List of dicts with family analysis results, sorted by num_cath_groups descending
    """
    results = []

    for pfamA_acc, cath_counts in family_cath.items():
        num_cath = len(cath_counts)

        if num_cath < min_cath_groups:
            continue

        # Sort CATH groups by count
        sorted_cath = sorted(cath_counts.items(), key=lambda x: -x[1])
        total_count = sum(cath_counts.values())

        # Get clan info
        clan_info = family_to_clan.get(pfamA_acc, ('', ''))

        # Calculate dominance ratio
        top_count = sorted_cath[0][1]
        dominance_ratio = top_count / total_count if total_count > 0 else 0

        results.append({
            'pfamA_acc': pfamA_acc,
            'clan_acc': clan_info[0],
            'clan_id': clan_info[1],
            'num_cath_groups': num_cath,
            'total_overlaps': total_count,
            'dominance_ratio': dominance_ratio,
            'cath_groups': sorted_cath
        })

    # Sort by number of CATH groups (descending), then by total overlaps
    results.sort(key=lambda x: (-x['num_cath_groups'], -x['total_overlaps']))

    return results


def write_output(results: List[Dict], output_file: str):
    """Write analysis results to file."""

    with open(output_file, 'w') as f:
        # Header
        f.write("pfamA_acc\tclan_acc\tclan_id\tnum_cath_groups\ttotal_overlaps\t"
                "dominance_ratio\ttop_cath\ttop_count\tall_cath_groups\n")

        for r in results:
            top_cath, top_count = r['cath_groups'][0]

            # Format all CATH groups as "label:count,label:count,..."
            all_cath = ','.join(f"{label}:{count}" for label, count in r['cath_groups'])

            f.write(f"{r['pfamA_acc']}\t{r['clan_acc']}\t{r['clan_id']}\t"
                    f"{r['num_cath_groups']}\t{r['total_overlaps']}\t"
                    f"{r['dominance_ratio']:.3f}\t{top_cath}\t{top_count}\t{all_cath}\n")


def print_summary(results: List[Dict]):
    """Print summary statistics."""

    if not results:
        print("\nNo families found mapping to multiple CATH superfamilies.")
        return

    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"Total families with multiple CATH mappings: {len(results)}")

    # Distribution of CATH group counts
    cath_dist = defaultdict(int)
    for r in results:
        cath_dist[r['num_cath_groups']] += 1

    print(f"\nDistribution by number of CATH groups:")
    for num_cath in sorted(cath_dist.keys()):
        print(f"  {num_cath} CATH groups: {cath_dist[num_cath]} families")

    # Top 20 families by CATH group count
    print(f"\nTop 20 families with most CATH superfamily mappings:")
    print(f"{'Family':<12} {'Clan':<12} {'#CATH':>6} {'Total':>8} {'Dom.':>6}  Top CATH groups")
    print("-" * 90)

    for r in results[:20]:
        top_3 = ', '.join(f"{label}" for label, count in r['cath_groups'][:3])
        print(f"{r['pfamA_acc']:<12} {r['clan_acc']:<12} {r['num_cath_groups']:>6} "
              f"{r['total_overlaps']:>8} {r['dominance_ratio']:>6.2f}  {top_3}")

    # Families with low dominance (evenly spread across CATH groups)
    low_dom = [r for r in results if r['dominance_ratio'] < 0.5 and r['num_cath_groups'] >= 2]
    if low_dom:
        print(f"\nFamilies with spread mappings (dominance < 0.5): {len(low_dom)}")
        print("These are most likely to need review/splitting:")
        print(f"{'Family':<12} {'Clan':<12} {'#CATH':>6} {'Dom.':>6}  CATH groups")
        print("-" * 80)
        for r in low_dom[:15]:
            cath_str = ', '.join(f"{label}:{count}" for label, count in r['cath_groups'][:4])
            print(f"{r['pfamA_acc']:<12} {r['clan_acc']:<12} {r['num_cath_groups']:>6} "
                  f"{r['dominance_ratio']:>6.2f}  {cath_str}")


def main():
    args = parse_arguments()

    print("=" * 70)
    print("Pfam Multi-CATH Superfamily Analysis")
    print("=" * 70)
    print(f"\nInput files:")
    print(f"  Clan mapping:   {args.clan_mapping}")
    print(f"  Family counts:  {args.family_counts}")
    print(f"  Clan membership: {args.clan_file}")
    print(f"\nOutput file: {args.output}")
    print(f"\nFilters:")
    print(f"  Min confidence: {args.min_confidence}")
    print(f"  Min count: {args.min_count}")
    print(f"  Min CATH groups: {args.min_cath_groups}")

    # Check input files exist
    for filepath in [args.clan_mapping, args.family_counts]:
        if not os.path.exists(filepath):
            print(f"\nError: Input file not found: {filepath}")
            sys.exit(1)

    # Load clan membership
    print("\n" + "-" * 70)
    print("Loading clan membership...")
    family_to_clan = load_clan_membership(args.clan_file)
    print(f"  Loaded {len(family_to_clan)} family->clan mappings")

    # Load confident CATH->Clan mappings
    print("\nLoading confident CATH->Clan mappings...")
    confident_mappings = load_confident_cath_mappings(
        args.clan_mapping, args.min_confidence, args.min_count
    )
    print(f"  Found {len(confident_mappings)} confident (CATH, Clan) pairs")

    # Load family counts for confident mappings only
    print("\nLoading family counts for confident mappings...")
    family_cath = load_family_counts_for_confident_mappings(
        args.family_counts, confident_mappings, family_to_clan
    )
    print(f"  Loaded {len(family_cath)} families with confident CATH mappings")

    # Analyze
    print("\nAnalyzing multi-CATH families...")
    results = analyze_multi_cath_families(family_cath, family_to_clan, args.min_cath_groups)

    # Output
    write_output(results, args.output)
    print(f"\nResults written to: {args.output}")

    # Print summary
    print_summary(results)

    print(f"\n{'='*70}")
    print("Done!")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()
