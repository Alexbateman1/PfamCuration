#!/usr/bin/env python3
"""
Identify Pfam families that map to multiple CATH H-level superfamilies.

These families may need review as they could:
- Be too broad and need splitting
- Contain multiple structural domains
- Have curation issues

Usage:
    python pfam_multi_cath_analysis.py [options]

Input:
    cath_to_pfam_family_counts.tsv - Output from cath_to_pfam_mapping.py

Output:
    pfam_multi_cath_families.tsv - Families mapping to multiple CATH superfamilies
"""

import argparse
import os
import sys
from collections import defaultdict
from typing import Dict, List, Tuple


# Default paths
DEFAULT_INPUT = '/nfs/production/agb/pfam/data/TED/cath_to_pfam_family_counts.tsv'
DEFAULT_OUTPUT = '/nfs/production/agb/pfam/data/TED/pfam_multi_cath_families.tsv'
DEFAULT_CLAN_FILE = '/nfs/production/agb/pfam/data/TED/pfam_clan_membership.tsv'


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Identify Pfam families mapping to multiple CATH superfamilies'
    )
    parser.add_argument(
        '--input',
        default=DEFAULT_INPUT,
        help=f'Path to family counts file (default: {DEFAULT_INPUT})'
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
        '--min-cath-groups',
        type=int,
        default=2,
        help='Minimum number of CATH H-groups to report (default: 2)'
    )
    parser.add_argument(
        '--min-count',
        type=int,
        default=10,
        help='Minimum overlap count per CATH group to consider (default: 10)'
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


def load_family_counts(input_file: str, min_count: int) -> Dict[str, Dict[str, int]]:
    """
    Load family to CATH mapping counts.

    Only considers H-level (superfamily) assignments and filters by min_count.

    Returns:
        Dict mapping pfamA_acc -> {cath_label: count}
    """
    family_cath = defaultdict(dict)

    with open(input_file, 'r') as f:
        # Skip header
        header = f.readline()

        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue

            cath_label, cath_level, pfamA_acc, count = parts[0], parts[1], parts[2], int(parts[3])

            # Only consider H-level (homologous superfamily) assignments
            if cath_level != 'H':
                continue

            # Filter by minimum count
            if count < min_count:
                continue

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

        # Calculate entropy-like measure of how spread out the mappings are
        # If one CATH dominates, ratio will be high; if spread evenly, ratio will be low
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

    # Top 10 families by CATH group count
    print(f"\nTop 10 families with most CATH superfamily mappings:")
    print(f"{'Family':<12} {'Clan':<10} {'#CATH':>6} {'Total':>8} {'Dom.':>6}  Top CATH groups")
    print("-" * 80)

    for r in results[:10]:
        top_3 = ', '.join(f"{label}" for label, count in r['cath_groups'][:3])
        print(f"{r['pfamA_acc']:<12} {r['clan_acc']:<10} {r['num_cath_groups']:>6} "
              f"{r['total_overlaps']:>8} {r['dominance_ratio']:>6.2f}  {top_3}")

    # Families with low dominance (evenly spread across CATH groups)
    low_dom = [r for r in results if r['dominance_ratio'] < 0.5 and r['num_cath_groups'] >= 3]
    if low_dom:
        print(f"\nFamilies with spread mappings (dominance < 0.5, ≥3 CATH groups): {len(low_dom)}")
        print("These may be most likely to need splitting:")
        print(f"{'Family':<12} {'Clan':<10} {'#CATH':>6} {'Dom.':>6}  CATH groups")
        print("-" * 70)
        for r in low_dom[:10]:
            cath_str = ', '.join(f"{label}:{count}" for label, count in r['cath_groups'][:4])
            print(f"{r['pfamA_acc']:<12} {r['clan_acc']:<10} {r['num_cath_groups']:>6} "
                  f"{r['dominance_ratio']:>6.2f}  {cath_str}")


def main():
    args = parse_arguments()

    print("=" * 70)
    print("Pfam Multi-CATH Superfamily Analysis")
    print("=" * 70)
    print(f"\nInput file:  {args.input}")
    print(f"Clan file:   {args.clan_file}")
    print(f"Output file: {args.output}")
    print(f"\nFilters:")
    print(f"  Min CATH groups: {args.min_cath_groups}")
    print(f"  Min count per CATH: {args.min_count}")

    # Check input file exists
    if not os.path.exists(args.input):
        print(f"\nError: Input file not found: {args.input}")
        sys.exit(1)

    # Load data
    print("\nLoading clan membership...")
    family_to_clan = load_clan_membership(args.clan_file)
    print(f"  Loaded {len(family_to_clan)} family->clan mappings")

    print("\nLoading family counts (H-level only)...")
    family_cath = load_family_counts(args.input, args.min_count)
    print(f"  Loaded {len(family_cath)} families with CATH mappings")

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
