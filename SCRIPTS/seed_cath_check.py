#!/usr/bin/env python3
"""
Check CATH classifications for sequences in a Pfam SEED alignment.

Takes a SEED alignment file and reports which CATH T-groups and H-groups
the sequences map to, based on TED domain data.

Uses grep -f for efficient lookup against the large TED data file.

Usage:
    python seed_cath_check.py SEED
    python seed_cath_check.py /path/to/SEED --ted-file /path/to/ted_cath_high_sorted.tsv.gz

Output:
    Summary of CATH groups found in the SEED sequences
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from typing import Dict, List, Set, Tuple


# Default paths
DEFAULT_TED_FILE = '/nfs/production/agb/pfam/data/TED/ted_cath_high_sorted.tsv.gz'


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Check CATH classifications for sequences in a SEED alignment'
    )
    parser.add_argument(
        'seed_file',
        help='Path to SEED alignment file (Stockholm format)'
    )
    parser.add_argument(
        '--ted-file',
        default=DEFAULT_TED_FILE,
        help=f'Path to TED CATH data file (default: {DEFAULT_TED_FILE})'
    )
    parser.add_argument(
        '--show-sequences',
        action='store_true',
        help='Show individual sequence CATH assignments'
    )
    parser.add_argument(
        '--output',
        help='Output file for detailed results (optional)'
    )
    return parser.parse_args()


def parse_seed_accessions(seed_file: str) -> List[Tuple[str, int, int]]:
    """
    Parse SEED file and extract UniProt accessions with coordinates.

    Returns:
        List of (uniprot_acc, start, end) tuples
    """
    accessions = []

    # Pattern for Stockholm format sequence lines: ACC/start-end
    # e.g., A0A009H7F2/4-116
    seq_pattern = re.compile(r'^([A-Z0-9]+)/(\d+)-(\d+)\s+')

    with open(seed_file, 'r') as f:
        for line in f:
            line = line.rstrip()

            # Skip comment lines, blank lines, and Stockholm headers
            if not line or line.startswith('#') or line.startswith('//'):
                continue

            match = seq_pattern.match(line)
            if match:
                acc = match.group(1)
                start = int(match.group(2))
                end = int(match.group(3))
                accessions.append((acc, start, end))

    return accessions


def lookup_cath_assignments(
    accessions: List[Tuple[str, int, int]],
    ted_file: str
) -> Dict[str, List[Dict]]:
    """
    Look up CATH assignments for accessions using grep -f.

    Returns:
        Dict mapping uniprot_acc -> list of TED domain dicts
    """
    if not accessions:
        return {}

    # Create temp file with unique accessions for grep -f
    unique_accs = sorted(set(acc for acc, _, _ in accessions))

    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp:
        # Write patterns that match start of line
        for acc in unique_accs:
            tmp.write(f"^{acc}\t\n")
        tmp_path = tmp.name

    try:
        # Use zgrep -f for efficient lookup
        # The pattern file has ^ACC\t to match accession at start of line
        cmd = f"zgrep -E -f {tmp_path} {ted_file}"

        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True
        )

        # Parse results
        cath_by_acc = defaultdict(list)

        for line in result.stdout.strip().split('\n'):
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) >= 6:
                acc = parts[0]
                ted_suffix = parts[1]
                start = int(parts[2])
                end = int(parts[3])
                cath_label = parts[4]
                cath_level = parts[5]

                cath_by_acc[acc].append({
                    'ted_suffix': ted_suffix,
                    'start': start,
                    'end': end,
                    'cath_label': cath_label,
                    'cath_level': cath_level
                })

        return cath_by_acc

    finally:
        os.unlink(tmp_path)


def calculate_overlap(start1: int, end1: int, start2: int, end2: int) -> float:
    """Calculate overlap fraction (of first region)."""
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)

    if overlap_start > overlap_end:
        return 0.0

    overlap_length = overlap_end - overlap_start + 1
    length1 = end1 - start1 + 1

    return overlap_length / length1


def match_seed_to_cath(
    seed_accessions: List[Tuple[str, int, int]],
    cath_by_acc: Dict[str, List[Dict]],
    min_overlap: float = 0.5
) -> List[Dict]:
    """
    Match SEED sequences to overlapping CATH domains.

    Returns:
        List of match results with SEED coords and CATH info
    """
    results = []

    for acc, seed_start, seed_end in seed_accessions:
        ted_domains = cath_by_acc.get(acc, [])

        matched_cath = []
        for ted in ted_domains:
            overlap = calculate_overlap(
                seed_start, seed_end,
                ted['start'], ted['end']
            )

            if overlap >= min_overlap:
                matched_cath.append({
                    'cath_label': ted['cath_label'],
                    'cath_level': ted['cath_level'],
                    'ted_suffix': ted['ted_suffix'],
                    'ted_start': ted['start'],
                    'ted_end': ted['end'],
                    'overlap': overlap
                })

        results.append({
            'acc': acc,
            'seed_start': seed_start,
            'seed_end': seed_end,
            'cath_matches': matched_cath
        })

    return results


def summarize_results(results: List[Dict]) -> Dict:
    """Generate summary statistics from results."""

    # Count by CATH H-group (full code)
    h_group_counts = defaultdict(int)
    # Count by CATH T-group (first 3 numbers)
    t_group_counts = defaultdict(int)
    # Count by Class.Architecture (first 2 numbers)
    ca_counts = defaultdict(int)

    # Track sequences
    seqs_with_cath = 0
    seqs_without_cath = 0
    seqs_multi_cath = 0

    # Track which sequences go to which groups
    h_group_seqs = defaultdict(list)

    for r in results:
        if not r['cath_matches']:
            seqs_without_cath += 1
            continue

        seqs_with_cath += 1

        # Get unique CATH labels for this sequence
        unique_h = set()
        for m in r['cath_matches']:
            cath = m['cath_label']
            level = m['cath_level']

            if level == 'H':
                unique_h.add(cath)
                h_group_counts[cath] += 1
                h_group_seqs[cath].append(f"{r['acc']}/{r['seed_start']}-{r['seed_end']}")

                # Extract T-group (first 3 numbers)
                parts = cath.split('.')
                if len(parts) >= 3:
                    t_group = '.'.join(parts[:3])
                    t_group_counts[t_group] += 1

                # Extract Class.Architecture
                if len(parts) >= 2:
                    ca = '.'.join(parts[:2])
                    ca_counts[ca] += 1

            elif level == 'T':
                # T-level assignment (only 3 numbers)
                t_group_counts[cath] += 1
                parts = cath.split('.')
                if len(parts) >= 2:
                    ca = '.'.join(parts[:2])
                    ca_counts[ca] += 1

        if len(unique_h) > 1:
            seqs_multi_cath += 1

    return {
        'total_sequences': len(results),
        'seqs_with_cath': seqs_with_cath,
        'seqs_without_cath': seqs_without_cath,
        'seqs_multi_cath': seqs_multi_cath,
        'h_group_counts': dict(h_group_counts),
        't_group_counts': dict(t_group_counts),
        'ca_counts': dict(ca_counts),
        'h_group_seqs': dict(h_group_seqs)
    }


def print_summary(summary: Dict, show_sequences: bool = False):
    """Print summary to stdout."""

    print(f"\n{'='*70}")
    print("CATH Classification Summary for SEED")
    print(f"{'='*70}")

    print(f"\nSequence Statistics:")
    print(f"  Total sequences in SEED: {summary['total_sequences']}")
    print(f"  Sequences with CATH assignment: {summary['seqs_with_cath']}")
    print(f"  Sequences without CATH assignment: {summary['seqs_without_cath']}")
    print(f"  Sequences with multiple H-groups: {summary['seqs_multi_cath']}")

    # Class.Architecture breakdown
    print(f"\n{'-'*70}")
    print("Class.Architecture (C.A) Distribution:")
    print(f"{'C.A':<15} {'Count':>8}")
    print("-" * 25)
    for ca, count in sorted(summary['ca_counts'].items(), key=lambda x: -x[1]):
        print(f"{ca:<15} {count:>8}")

    # T-group breakdown
    print(f"\n{'-'*70}")
    print("Topology (C.A.T) Distribution:")
    print(f"{'T-group':<20} {'Count':>8}")
    print("-" * 30)
    for t_group, count in sorted(summary['t_group_counts'].items(), key=lambda x: -x[1]):
        print(f"{t_group:<20} {count:>8}")

    # H-group breakdown
    print(f"\n{'-'*70}")
    print("Homologous Superfamily (C.A.T.H) Distribution:")
    print(f"{'H-group':<20} {'Count':>8}")
    print("-" * 30)
    for h_group, count in sorted(summary['h_group_counts'].items(), key=lambda x: -x[1]):
        print(f"{h_group:<20} {count:>8}")

    # Show sequences per H-group if requested
    if show_sequences and summary['h_group_seqs']:
        print(f"\n{'-'*70}")
        print("Sequences per H-group:")
        for h_group, seqs in sorted(summary['h_group_seqs'].items()):
            print(f"\n{h_group} ({len(seqs)} sequences):")
            for seq in seqs[:10]:  # Show first 10
                print(f"  {seq}")
            if len(seqs) > 10:
                print(f"  ... and {len(seqs) - 10} more")

    print(f"\n{'='*70}")


def write_detailed_output(results: List[Dict], output_file: str):
    """Write detailed results to file."""

    with open(output_file, 'w') as f:
        f.write("acc\tseed_start\tseed_end\tcath_label\tcath_level\tted_suffix\tted_start\tted_end\toverlap\n")

        for r in results:
            if r['cath_matches']:
                for m in r['cath_matches']:
                    f.write(f"{r['acc']}\t{r['seed_start']}\t{r['seed_end']}\t"
                            f"{m['cath_label']}\t{m['cath_level']}\t{m['ted_suffix']}\t"
                            f"{m['ted_start']}\t{m['ted_end']}\t{m['overlap']:.3f}\n")
            else:
                f.write(f"{r['acc']}\t{r['seed_start']}\t{r['seed_end']}\t-\t-\t-\t-\t-\t-\n")


def main():
    args = parse_arguments()

    # Check files exist
    if not os.path.exists(args.seed_file):
        print(f"Error: SEED file not found: {args.seed_file}")
        sys.exit(1)

    if not os.path.exists(args.ted_file):
        print(f"Error: TED file not found: {args.ted_file}")
        sys.exit(1)

    print(f"Analyzing SEED: {args.seed_file}")
    print(f"Using TED data: {args.ted_file}")

    # Parse SEED
    print("\nParsing SEED file...")
    seed_accessions = parse_seed_accessions(args.seed_file)
    print(f"  Found {len(seed_accessions)} sequences")

    if not seed_accessions:
        print("Error: No sequences found in SEED file")
        sys.exit(1)

    # Lookup CATH assignments
    print("\nLooking up CATH assignments (this may take a moment)...")
    cath_by_acc = lookup_cath_assignments(seed_accessions, args.ted_file)
    print(f"  Found TED data for {len(cath_by_acc)} accessions")

    # Match SEED to CATH
    print("\nMatching SEED coordinates to CATH domains...")
    results = match_seed_to_cath(seed_accessions, cath_by_acc)

    # Summarize
    summary = summarize_results(results)
    print_summary(summary, args.show_sequences)

    # Write detailed output if requested
    if args.output:
        write_detailed_output(results, args.output)
        print(f"\nDetailed results written to: {args.output}")

    print("\nDone!")


if __name__ == '__main__':
    main()
