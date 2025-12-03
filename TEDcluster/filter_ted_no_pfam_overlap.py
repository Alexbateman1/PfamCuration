#!/usr/bin/env python3
"""
Filter TED domains that don't overlap with existing Pfam domains.

This script uses an external sort approach to keep memory usage low:
1. Extracts high-confidence TED domains to a temp file
2. Combines Pfam and TED domains into one file with source markers
3. Sorts by UniProt accession using Unix sort (disk-based, memory efficient)
4. Processes one protein at a time, comparing TED vs Pfam domains
5. Outputs non-overlapping TED domains

Usage:
    python filter_ted_no_pfam_overlap.py [options]

Output format (TSV):
    uniprot_acc<TAB>start<TAB>end<TAB>ted_suffix

For discontinuous domains, start is the first coordinate and end is the last.
"""

import argparse
import gzip
import os
import re
import subprocess
import sys
import tempfile
from typing import List, Tuple

# Default paths
DEFAULT_TED_TSV = '/nfs/production/agb/pfam/data/TED/ted_365m_domain_boundaries_consensus_level.tsv.gz'
DEFAULT_OUTPUT = 'ted_high_no_pfam_overlap.tsv'
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
        '--progress-interval',
        type=int,
        default=1000000,
        help='Print progress every N TED domains (default: 1000000)'
    )
    parser.add_argument(
        '--sort-memory',
        default='4G',
        help='Memory limit for Unix sort (default: 4G)'
    )
    parser.add_argument(
        '--sort-parallel',
        type=int,
        default=4,
        help='Number of parallel sort threads (default: 4)'
    )
    parser.add_argument(
        '--temp-dir',
        default=None,
        help='Temporary directory for sort (default: system temp)'
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
        print("Running MySQL query...")
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


def extract_uniprot_acc(ted_id: str) -> str:
    """
    Extract UniProt accession from TED ID.

    Args:
        ted_id: TED ID like "AF-A0A000-F1-model_v4_TED01"

    Returns:
        UniProt accession like "A0A000"
    """
    match = re.match(r'^AF-([A-Z0-9]+)-F\d+-', ted_id)
    if match:
        return match.group(1)
    return None


def extract_ted_suffix(ted_id: str) -> str:
    """
    Extract TED suffix from TED ID.

    Args:
        ted_id: TED ID like "AF-A0A000-F1-model_v4_TED01"

    Returns:
        TED suffix like "TED01"
    """
    if '_' in ted_id:
        return ted_id.split('_')[-1]
    return None


def extract_ted_domains(ted_file: str, output_file: str, progress_interval: int) -> Tuple[int, int]:
    """
    Extract high-confidence TED domains to a file.

    Output format: uniprot_acc<TAB>start<TAB>end<TAB>ted_suffix

    Args:
        ted_file: Path to TED TSV file (gzipped)
        output_file: Path to output file
        progress_interval: Print progress every N records

    Returns:
        Tuple of (total_count, high_count)
    """
    print(f"\nExtracting high-confidence TED domains from {ted_file}...")

    total_count = 0
    high_count = 0
    no_acc_count = 0

    if ted_file.endswith('.gz'):
        in_f = gzip.open(ted_file, 'rt')
    else:
        in_f = open(ted_file, 'r')

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

            # Extract UniProt accession and TED suffix
            uniprot_acc = extract_uniprot_acc(ted_id)
            if not uniprot_acc:
                no_acc_count += 1
                continue

            ted_suffix = extract_ted_suffix(ted_id)

            # Parse TED domain segments
            ted_segments = parse_chopping(chopping)
            if not ted_segments:
                continue

            # Get overall start (first coord) and end (last coord)
            overall_start = ted_segments[0][0]
            overall_end = ted_segments[-1][1]

            # Write with all segment info for overlap checking
            # Format: uniprot_acc, start, end, ted_suffix, all_segments
            segments_str = chopping
            out_f.write(f"{uniprot_acc}\t{overall_start}\t{overall_end}\t{ted_suffix}\t{segments_str}\n")
            high_count += 1

            if total_count % progress_interval == 0:
                print(f"  Processed {total_count:,} total, {high_count:,} high-confidence...")

    finally:
        in_f.close()
        out_f.close()

    print(f"Extracted {high_count:,} high-confidence TED domains")
    print(f"  (from {total_count:,} total, {no_acc_count:,} without valid accession)")

    return total_count, high_count


def create_combined_file(pfam_file: str, ted_file: str, combined_file: str):
    """
    Create a combined file with both Pfam and TED domains.

    Format: uniprot_acc<TAB>start<TAB>end<TAB>source<TAB>extra
    - source: 'P' for Pfam, 'T' for TED
    - extra: for TED, contains ted_suffix and segments_str

    Args:
        pfam_file: Path to Pfam domains file
        ted_file: Path to TED domains file
        combined_file: Path to output combined file
    """
    print(f"\nCreating combined domains file...")

    pfam_count = 0
    ted_count = 0

    with open(combined_file, 'w') as out_f:
        # Add Pfam domains (format: acc, start, end)
        print(f"  Adding Pfam domains from {pfam_file}...")
        with open(pfam_file, 'r') as f:
            for line in f:
                line = line.rstrip('\n\r')
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= 3:
                    # Format: acc, start, end, source='P', extra=''
                    out_f.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\tP\t\n")
                    pfam_count += 1

        # Add TED domains (format: acc, start, end, ted_suffix, segments)
        print(f"  Adding TED domains from {ted_file}...")
        with open(ted_file, 'r') as f:
            for line in f:
                line = line.rstrip('\n\r')
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= 5:
                    # Format: acc, start, end, source='T', extra=ted_suffix:segments
                    extra = f"{parts[3]}:{parts[4]}"
                    out_f.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\tT\t{extra}\n")
                    ted_count += 1

    print(f"Combined file created: {pfam_count:,} Pfam + {ted_count:,} TED = {pfam_count + ted_count:,} total")


def sort_file(input_file: str, output_file: str, memory: str, parallel: int, temp_dir: str = None):
    """
    Sort file by first column using Unix sort.

    Args:
        input_file: Path to input file
        output_file: Path to output sorted file
        memory: Memory limit (e.g., '4G')
        parallel: Number of parallel threads
        temp_dir: Temporary directory for sort
    """
    print(f"\nSorting combined file by UniProt accession...")
    print(f"  Memory limit: {memory}, Parallel threads: {parallel}")

    cmd = [
        'sort',
        '-t', '\t',
        '-k1,1',  # Sort by first column
        f'-S{memory}',
        f'--parallel={parallel}',
        '-o', output_file,
        input_file
    ]

    if temp_dir:
        cmd.insert(-2, f'-T{temp_dir}')

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
        if result.returncode != 0:
            print(f"ERROR: Sort failed with code {result.returncode}")
            if result.stderr:
                print(f"  {result.stderr}")
            sys.exit(1)
        print("Sort completed successfully")
    except subprocess.TimeoutExpired:
        print("ERROR: Sort timed out after 2 hours")
        sys.exit(1)


def segments_overlap(ted_segments: List[Tuple[int, int]],
                     pfam_regions: List[Tuple[int, int]]) -> bool:
    """
    Check if any TED segment overlaps with any Pfam region.

    Args:
        ted_segments: List of (start, end) tuples for TED domain
        pfam_regions: List of (start, end) tuples for Pfam domains

    Returns:
        True if any overlap exists, False otherwise
    """
    for ted_start, ted_end in ted_segments:
        for pfam_start, pfam_end in pfam_regions:
            if max(ted_start, pfam_start) <= min(ted_end, pfam_end):
                return True
    return False


def process_sorted_file(sorted_file: str, output_file: str, progress_interval: int):
    """
    Process sorted file one protein at a time.

    Args:
        sorted_file: Path to sorted combined file
        output_file: Path to output file
        progress_interval: Print progress every N proteins
    """
    print(f"\nProcessing sorted file to find non-overlapping TED domains...")

    proteins_processed = 0
    ted_total = 0
    ted_kept = 0
    ted_overlap = 0

    current_acc = None
    current_pfam_regions = []
    current_ted_domains = []

    def process_protein():
        """Process all domains for current protein."""
        nonlocal ted_total, ted_kept, ted_overlap

        if not current_acc:
            return

        for ted_start, ted_end, ted_suffix, segments_str in current_ted_domains:
            ted_total += 1

            # Parse segments for overlap checking
            ted_segments = parse_chopping(segments_str)
            if not ted_segments:
                ted_segments = [(ted_start, ted_end)]

            if current_pfam_regions and segments_overlap(ted_segments, current_pfam_regions):
                ted_overlap += 1
            else:
                # No overlap - output this domain
                out_f.write(f"{current_acc}\t{ted_start}\t{ted_end}\t{ted_suffix}\n")
                ted_kept += 1

    with open(sorted_file, 'r') as in_f, open(output_file, 'w') as out_f:
        for line in in_f:
            line = line.rstrip('\n\r')
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) < 4:
                continue

            acc = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            source = parts[3]
            extra = parts[4] if len(parts) > 4 else ''

            # Check if we've moved to a new protein
            if acc != current_acc:
                # Process previous protein
                process_protein()

                # Reset for new protein
                current_acc = acc
                current_pfam_regions = []
                current_ted_domains = []
                proteins_processed += 1

                if proteins_processed % progress_interval == 0:
                    print(f"  Processed {proteins_processed:,} proteins, "
                          f"kept {ted_kept:,}/{ted_total:,} TED domains...")

            # Add domain to current protein
            if source == 'P':
                current_pfam_regions.append((start, end))
            elif source == 'T':
                # Parse extra field: ted_suffix:segments_str
                if ':' in extra:
                    ted_suffix, segments_str = extra.split(':', 1)
                else:
                    ted_suffix = extra
                    segments_str = f"{start}-{end}"
                current_ted_domains.append((start, end, ted_suffix, segments_str))

        # Process last protein
        process_protein()

    print(f"\n{'='*60}")
    print("FILTERING COMPLETE")
    print(f"{'='*60}")
    print(f"Proteins processed:                {proteins_processed:,}")
    print(f"TED domains checked:               {ted_total:,}")
    print(f"TED domains with Pfam overlap:     {ted_overlap:,}")
    print(f"TED domains kept (no overlap):     {ted_kept:,}")
    print(f"{'='*60}")
    print(f"Output written to: {output_file}")


def main():
    args = parse_arguments()

    print(f"\n{'='*60}")
    print("TED Domain Filter - Remove Pfam Overlapping Domains")
    print("(Memory-efficient external sort approach)")
    print(f"{'='*60}")

    # Step 1: Fetch Pfam domains if cache doesn't exist
    if os.path.exists(args.pfam_cache):
        print(f"\nPfam cache already exists: {args.pfam_cache}")
        print("Skipping database fetch (delete file to re-fetch)")
    else:
        fetch_pfam_domains(args.mysql_config, args.pfam_cache)

    # Create temp directory for intermediate files
    temp_dir = args.temp_dir or tempfile.gettempdir()
    ted_extracted_file = os.path.join(temp_dir, 'ted_high_extracted.tsv')
    combined_file = os.path.join(temp_dir, 'combined_domains.tsv')
    sorted_file = os.path.join(temp_dir, 'combined_domains_sorted.tsv')

    try:
        # Step 2: Extract high-confidence TED domains
        extract_ted_domains(args.ted_file, ted_extracted_file, args.progress_interval)

        # Step 3: Create combined file
        create_combined_file(args.pfam_cache, ted_extracted_file, combined_file)

        # Step 4: Sort by UniProt accession
        sort_file(combined_file, sorted_file, args.sort_memory, args.sort_parallel, temp_dir)

        # Step 5: Process sorted file
        process_sorted_file(sorted_file, args.output, args.progress_interval)

    finally:
        # Clean up temp files
        for f in [ted_extracted_file, combined_file, sorted_file]:
            if os.path.exists(f):
                try:
                    os.remove(f)
                    print(f"Cleaned up: {f}")
                except:
                    pass

    print("\nDone!")


if __name__ == '__main__':
    main()
