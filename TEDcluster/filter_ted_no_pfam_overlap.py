#!/usr/bin/env python3
"""
Filter TED domains that don't overlap with existing Pfam domains.

This script:
1. Reads high-confidence TED domains from the consensus TSV file
2. Fetches all Pfam domain regions from the pfamlive database
3. Filters out TED domains that overlap with any Pfam domain (even by 1 residue)
4. Outputs the non-overlapping TED domains

Usage:
    python filter_ted_no_pfam_overlap.py [options]

Output format (TSV):
    ted_id<TAB>chopping<TAB>consensus_level
"""

import argparse
import gzip
import os
import re
import subprocess
import sys
from collections import defaultdict
from typing import Dict, List, Set, Tuple

# Default paths
DEFAULT_TED_TSV = '/nfs/production/agb/pfam/data/TED/ted_365m_domain_boundaries_consensus_level.tsv.gz'
DEFAULT_OUTPUT = 'ted_high_no_pfam_overlap.tsv.gz'
DEFAULT_PFAM_CACHE = 'pfam_domains_cache.tsv'


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Filter TED domains that do not overlap with Pfam domains'
    )
    parser.add_argument(
        '--ted-file',
        default=DEFAULT_TED_TSV,
        help=f'Path to TED TSV file (default: {DEFAULT_TED_TSV})'
    )
    parser.add_argument(
        '--output',
        default=DEFAULT_OUTPUT,
        help=f'Output file path (default: {DEFAULT_OUTPUT})'
    )
    parser.add_argument(
        '--pfam-cache',
        default=DEFAULT_PFAM_CACHE,
        help=f'Pfam domains cache file (default: {DEFAULT_PFAM_CACHE})'
    )
    parser.add_argument(
        '--mysql-config',
        default='~/.my.cnf',
        help='Path to MySQL config file (default: ~/.my.cnf)'
    )
    parser.add_argument(
        '--skip-pfam-fetch',
        action='store_true',
        help='Skip fetching Pfam domains (use existing cache)'
    )
    parser.add_argument(
        '--progress-interval',
        type=int,
        default=1000000,
        help='Print progress every N TED domains (default: 1000000)'
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Force re-fetch Pfam domains even if cache exists (non-interactive)'
    )
    parser.add_argument(
        '--use-cache',
        action='store_true',
        help='Use existing Pfam cache without prompting (non-interactive)'
    )
    return parser.parse_args()


def fetch_pfam_domains(config_file: str, output_file: str) -> int:
    """
    Fetch all Pfam domain regions from pfamlive database.

    Args:
        config_file: Path to MySQL config file
        output_file: Path to save the domain data

    Returns:
        Number of domain records fetched
    """
    print("Fetching Pfam domain regions from pfamlive database...")
    print("This may take several minutes...")

    config_file = os.path.expanduser(config_file)

    # Query to get all significant Pfam domain regions
    # Columns: pfamseq_acc, seq_start, seq_end
    query = """
        SELECT pfamseq_acc, seq_start, seq_end
        FROM pfamA_reg_full_significant
        WHERE in_full = 1
    """

    cmd = [
        'mysql',
        f'--defaults-file={config_file}',
        'pfam_live',
        '--quick',
        '--silent',
        '--skip-column-names',
        '-e',
        query
    ]

    try:
        print(f"Running MySQL query...")
        with open(output_file, 'w') as out_f:
            result = subprocess.run(
                cmd,
                stdout=out_f,
                stderr=subprocess.PIPE,
                text=True,
                timeout=3600  # 1 hour timeout for large query
            )

        if result.returncode != 0:
            print(f"ERROR: MySQL query failed with code {result.returncode}")
            if result.stderr:
                print(f"  {result.stderr}")
            sys.exit(1)

        # Count lines
        with open(output_file, 'r') as f:
            count = sum(1 for _ in f)

        print(f"Fetched {count:,} Pfam domain regions to {output_file}")
        return count

    except subprocess.TimeoutExpired:
        print("ERROR: MySQL query timed out after 1 hour")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to fetch Pfam domains: {e}")
        sys.exit(1)


def load_pfam_domains(cache_file: str) -> Dict[str, List[Tuple[int, int]]]:
    """
    Load Pfam domains from cache file into memory.

    Args:
        cache_file: Path to the cached Pfam domains file

    Returns:
        Dictionary mapping UniProt accession to list of (start, end) tuples
    """
    print(f"Loading Pfam domains from {cache_file}...")

    pfam_domains = defaultdict(list)
    count = 0

    with open(cache_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n\r')
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) >= 3:
                uniprot_acc = parts[0]
                try:
                    seq_start = int(parts[1])
                    seq_end = int(parts[2])
                    pfam_domains[uniprot_acc].append((seq_start, seq_end))
                    count += 1
                except ValueError:
                    continue

            if count % 5000000 == 0:
                print(f"  Loaded {count:,} Pfam domains...")

    print(f"Loaded {count:,} Pfam domains for {len(pfam_domains):,} proteins")
    return dict(pfam_domains)


def parse_chopping(chopping: str) -> List[Tuple[int, int]]:
    """
    Parse TED chopping string into list of (start, end) tuples.

    Args:
        chopping: Chopping string like "11-41_290-389"

    Returns:
        List of (start, end) tuples
    """
    segments = []
    for segment in chopping.split('_'):
        if '-' in segment:
            try:
                start, end = segment.split('-')
                segments.append((int(start), int(end)))
            except ValueError:
                continue
    return segments


def segments_overlap(ted_segments: List[Tuple[int, int]],
                     pfam_regions: List[Tuple[int, int]]) -> bool:
    """
    Check if any TED segment overlaps with any Pfam region.

    Overlap occurs if there is any shared residue position.

    Args:
        ted_segments: List of (start, end) tuples for TED domain
        pfam_regions: List of (start, end) tuples for Pfam domains on this protein

    Returns:
        True if any overlap exists, False otherwise
    """
    for ted_start, ted_end in ted_segments:
        for pfam_start, pfam_end in pfam_regions:
            # Check for overlap: segments overlap if neither is completely before or after the other
            # Overlap exists if: max(start1, start2) <= min(end1, end2)
            if max(ted_start, pfam_start) <= min(ted_end, pfam_end):
                return True
    return False


def extract_uniprot_acc(ted_id: str) -> str:
    """
    Extract UniProt accession from TED ID.

    Args:
        ted_id: TED ID like "AF-A0A000-F1-model_v4_TED01"

    Returns:
        UniProt accession like "A0A000"
    """
    # Pattern: AF-{UniProtAcc}-F{n}-model_v{version}_TED{nn}
    match = re.match(r'^AF-([A-Z0-9]+)-F\d+-', ted_id)
    if match:
        return match.group(1)
    return None


def filter_ted_domains(ted_file: str, pfam_domains: Dict[str, List[Tuple[int, int]]],
                       output_file: str, progress_interval: int):
    """
    Filter TED domains that don't overlap with Pfam domains.

    Args:
        ted_file: Path to TED TSV file (gzipped)
        pfam_domains: Dictionary of Pfam domains by UniProt accession
        output_file: Path to output file
        progress_interval: Print progress every N records
    """
    print(f"\nFiltering TED domains from {ted_file}...")
    print("Only keeping high-confidence domains that don't overlap Pfam...")

    # Counters
    total_count = 0
    high_count = 0
    kept_count = 0
    overlap_count = 0
    no_acc_count = 0

    # Open files
    if ted_file.endswith('.gz'):
        in_f = gzip.open(ted_file, 'rt')
    else:
        in_f = open(ted_file, 'r')

    if output_file.endswith('.gz'):
        out_f = gzip.open(output_file, 'wt')
    else:
        out_f = open(output_file, 'w')

    try:
        for line in in_f:
            line = line.rstrip('\n\r')
            if not line:
                continue

            total_count += 1

            parts = line.split('\t')
            if len(parts) < 3:
                continue

            ted_id = parts[0]
            chopping = parts[1]
            consensus_level = parts[2]

            # Only process high confidence domains
            if consensus_level != 'high':
                continue

            high_count += 1

            # Extract UniProt accession
            uniprot_acc = extract_uniprot_acc(ted_id)
            if not uniprot_acc:
                no_acc_count += 1
                continue

            # Parse TED domain segments
            ted_segments = parse_chopping(chopping)
            if not ted_segments:
                continue

            # Check for overlap with Pfam domains
            pfam_regions = pfam_domains.get(uniprot_acc, [])

            if pfam_regions and segments_overlap(ted_segments, pfam_regions):
                overlap_count += 1
            else:
                # No overlap - keep this domain
                out_f.write(f"{ted_id}\t{chopping}\t{consensus_level}\n")
                kept_count += 1

            # Progress update
            if total_count % progress_interval == 0:
                print(f"  Processed {total_count:,} total, {high_count:,} high-confidence, "
                      f"kept {kept_count:,}, filtered {overlap_count:,}")

    finally:
        in_f.close()
        out_f.close()

    # Final statistics
    print(f"\n{'='*60}")
    print("FILTERING COMPLETE")
    print(f"{'='*60}")
    print(f"Total TED domains processed:       {total_count:,}")
    print(f"High-confidence domains:           {high_count:,}")
    print(f"Domains with Pfam overlap:         {overlap_count:,}")
    print(f"Domains without UniProt acc:       {no_acc_count:,}")
    print(f"Domains kept (no Pfam overlap):    {kept_count:,}")
    print(f"{'='*60}")
    print(f"Output written to: {output_file}")


def main():
    args = parse_arguments()

    print(f"\n{'='*60}")
    print("TED Domain Filter - Remove Pfam Overlapping Domains")
    print(f"{'='*60}")

    # Step 1: Fetch or load Pfam domains
    pfam_cache = args.pfam_cache

    if not args.skip_pfam_fetch:
        if os.path.exists(pfam_cache):
            if args.force:
                # Force re-fetch
                print(f"\nForce mode: Re-fetching Pfam domains...")
                fetch_pfam_domains(args.mysql_config, pfam_cache)
            elif args.use_cache:
                # Use existing cache without prompting
                print(f"\nUsing existing Pfam cache: {pfam_cache}")
            else:
                # Interactive mode
                print(f"\nPfam cache file exists: {pfam_cache}")
                response = input("Use existing cache? [Y/n]: ").strip().lower()
                if response not in ('', 'y', 'yes'):
                    fetch_pfam_domains(args.mysql_config, pfam_cache)
        else:
            fetch_pfam_domains(args.mysql_config, pfam_cache)
    else:
        if not os.path.exists(pfam_cache):
            print(f"ERROR: Pfam cache file not found: {pfam_cache}")
            print("Run without --skip-pfam-fetch to create it")
            sys.exit(1)
        print(f"\nUsing existing Pfam cache: {pfam_cache}")

    # Step 2: Load Pfam domains into memory
    pfam_domains = load_pfam_domains(pfam_cache)

    # Step 3: Filter TED domains
    filter_ted_domains(
        args.ted_file,
        pfam_domains,
        args.output,
        args.progress_interval
    )

    print("\nDone!")


if __name__ == '__main__':
    main()
