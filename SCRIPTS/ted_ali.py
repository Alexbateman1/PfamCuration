#!/usr/bin/env python3
"""
Generate a report of TED domains for proteins in a Pfam alignment.

This script replaces ted_ali.pl with a Python implementation that uses
local TED data instead of API calls.

Usage:
    ted_ali.py <alignment_file> [options]

    ted_ali.py SEED
    ted_ali.py SEED --verbose
    ted_ali.py SEED --force

Output files:
    TED      - Tab-separated report of TED domains with overlap information
    TED_NUM  - Mean TED domains per protein (single number)

The TED file contains one line per TED domain found, with columns:
    Protein_Accession  TED_Domains  Domain_Count  Domain_IDs  Domain_Ranges  Overlap_with_Pfam
"""

import argparse
import os
import re
import sys
from typing import Dict, List, Tuple, Optional

# Import our local TED module
from ted_local import TEDLocal


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Generate a report of TED domains for proteins in a Pfam alignment.'
    )
    parser.add_argument(
        'alignment_file',
        help='Path to the alignment file (e.g., SEED)'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Show detailed progress information'
    )
    parser.add_argument(
        '--force', '-f',
        action='store_true',
        help='Overwrite existing output files'
    )
    parser.add_argument(
        '--output', '-o',
        default='TED',
        help='Output file name (default: TED)'
    )
    parser.add_argument(
        '--mean-file', '-m',
        default='TED_NUM',
        help='Mean domains output file name (default: TED_NUM)'
    )
    parser.add_argument(
        '--ted-data',
        default=None,
        help='Path to local TED data file (uses default if not specified)'
    )

    return parser.parse_args()


def strip_version(accession: str) -> str:
    """
    Strip version number from UniProt accession.

    Args:
        accession: UniProt accession possibly with version (e.g., 'O34341.2')

    Returns:
        Accession without version (e.g., 'O34341')
    """
    return accession.split('.')[0]


def parse_alignment_file(alignment_file: str, verbose: bool = False) -> Dict[str, Dict]:
    """
    Parse alignment file to extract sequence IDs and regions.

    Args:
        alignment_file: Path to alignment file
        verbose: Print progress information

    Returns:
        Dict mapping sequence IDs to their info (uniprot, pfam_start, pfam_end)
    """
    if verbose:
        print("Parsing alignment file...", file=sys.stderr)

    sequences = {}

    with open(alignment_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Match lines like: O34341.2/15-100 or P12345/1-200
            match = re.match(r'^(\S+)/(\d+)-(\d+)', line)
            if match:
                seq_id = match.group(0).split()[0]  # Get the full ID with region
                uniprot_with_version = match.group(1)
                pfam_start = int(match.group(2))
                pfam_end = int(match.group(3))

                sequences[seq_id] = {
                    'uniprot': uniprot_with_version,
                    'uniprot_no_version': strip_version(uniprot_with_version),
                    'pfam_start': pfam_start,
                    'pfam_end': pfam_end
                }

    if verbose:
        print(f"Found {len(sequences)} sequences in alignment.", file=sys.stderr)

    return sequences


def calculate_overlap_percentage(
    ted_start: int,
    ted_end: int,
    pfam_start: int,
    pfam_end: int
) -> float:
    """
    Calculate the percentage of the Pfam region covered by the TED domain.

    Args:
        ted_start: TED domain start position
        ted_end: TED domain end position
        pfam_start: Pfam region start position
        pfam_end: Pfam region end position

    Returns:
        Overlap percentage (0-100) relative to Pfam region length
    """
    # Check if there's any overlap
    if ted_start > pfam_end or ted_end < pfam_start:
        return 0.0

    # Calculate overlap
    overlap_start = max(ted_start, pfam_start)
    overlap_end = min(ted_end, pfam_end)
    overlap_length = overlap_end - overlap_start + 1

    # Calculate percentage relative to Pfam region
    pfam_length = pfam_end - pfam_start + 1
    overlap_percent = (overlap_length / pfam_length) * 100

    return overlap_percent


def main():
    args = parse_arguments()

    # Check input file exists
    if not os.path.exists(args.alignment_file):
        print(f"Error: Alignment file not found: {args.alignment_file}", file=sys.stderr)
        sys.exit(1)

    # Check if output files already exist
    if (os.path.exists(args.output) or os.path.exists(args.mean_file)) and not args.force:
        print(f"Error: Output file '{args.output}' or '{args.mean_file}' already exists.",
              file=sys.stderr)
        print("This indicates the script has already been run in this directory.",
              file=sys.stderr)
        print("Use --force to overwrite the existing files if needed.", file=sys.stderr)
        sys.exit(1)

    # Parse alignment file
    sequences = parse_alignment_file(args.alignment_file, args.verbose)

    if not sequences:
        print("Error: No sequences found in alignment file.", file=sys.stderr)
        sys.exit(1)

    # Initialize TED local data access
    if args.verbose:
        print("Loading TED domain data...", file=sys.stderr)

    try:
        ted = TEDLocal(args.ted_data)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # Process sequences and collect results
    seq_count = len(sequences)
    processed = 0
    with_ted = 0
    total_domains = 0
    errors = 0

    results = []

    for seq_id in sorted(sequences.keys()):
        processed += 1
        seq_info = sequences[seq_id]
        query_acc = seq_info['uniprot_no_version']

        if args.verbose:
            print(f"Fetching TED domains for {seq_id} (using {query_acc})...",
                  file=sys.stderr)
            print(f"Progress: {processed}/{seq_count} sequences processed",
                  file=sys.stderr)

        # Get TED domains from local data
        try:
            domains = ted.get_domains(query_acc)
        except Exception as e:
            if args.verbose:
                print(f"Warning: Error getting domains for {query_acc}: {e}",
                      file=sys.stderr)
            errors += 1
            continue

        domain_count = len(domains)
        if domain_count == 0:
            # Skip proteins with no TED domains
            continue

        with_ted += 1
        total_domains += domain_count

        # Process each domain separately
        pfam_start = seq_info['pfam_start']
        pfam_end = seq_info['pfam_end']

        for domain in domains:
            ted_id = domain['ted_id']
            ted_suffix = domain['ted_suffix']
            chopping = domain['chopping']
            segments = domain['segments']

            # Calculate overall domain boundaries
            if segments:
                ted_start = min(seg['start'] for seg in segments)
                ted_end = max(seg['end'] for seg in segments)
            else:
                continue

            # Calculate overlap percentage
            overlap_percent = calculate_overlap_percentage(
                ted_start, ted_end, pfam_start, pfam_end
            )

            # Format overlap as percentage string
            if overlap_percent > 0:
                overlap_str = f"{overlap_percent:.1f}%"
            else:
                overlap_str = "0"

            results.append({
                'seq_id': seq_id,
                'ted_info': f"{ted_suffix}:{chopping}",
                'domain_count': 1,  # Each line reports one domain
                'ted_id': ted_suffix,
                'domain_range': f"{ted_start}-{ted_end}",
                'overlap': overlap_str
            })

    # Write output file
    with open(args.output, 'w') as f:
        # Header
        f.write("Protein_Accession\tTED_Domains\tDomain_Count\tDomain_IDs\t"
                "Domain_Ranges\tOverlap_with_Pfam\n")

        # Data
        for r in results:
            f.write(f"{r['seq_id']}\t{r['ted_info']}\t{r['domain_count']}\t"
                    f"{r['ted_id']}\t{r['domain_range']}\t{r['overlap']}\n")

    # Calculate mean domains per protein
    mean_domains = total_domains / with_ted if with_ted > 0 else 0

    # Write mean file
    with open(args.mean_file, 'w') as f:
        f.write(f"{mean_domains:.2f}\n")

    # Print summary
    print(f"\n=== TED DOMAIN REPORT SUMMARY ===", file=sys.stderr)
    print(f"Total sequences analyzed: {seq_count}", file=sys.stderr)
    print(f"Sequences with TED domains: {with_ted}", file=sys.stderr)
    print(f"Total TED domains found: {total_domains}", file=sys.stderr)
    print(f"Mean domains per protein: {mean_domains:.2f}", file=sys.stderr)
    print(f"Errors: {errors}", file=sys.stderr)
    print(f"Results written to file: {args.output}", file=sys.stderr)
    print(f"Mean domains per protein written to: {args.mean_file}", file=sys.stderr)


if __name__ == '__main__':
    main()
