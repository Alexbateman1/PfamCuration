#!/usr/bin/env python3
"""
Filter TED domains that don't overlap with existing Pfam domains.

Uses fast shell commands and merge-sort approach to minimize memory.

Usage:
    python filter_ted_no_pfam_overlap.py [options]

Output format (TSV):
    uniprot_acc<TAB>start<TAB>end<TAB>ted_suffix
"""

import argparse
import os
import subprocess
import sys

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
        '--sort-memory',
        default='1G',
        help='Memory limit for Unix sort (default: 1G)'
    )
    parser.add_argument(
        '--work-dir',
        default='.',
        help='Working directory for intermediate files (default: current dir)'
    )
    parser.add_argument(
        '--cleanup',
        action='store_true',
        help='Clean up intermediate files after successful completion'
    )
    return parser.parse_args()


def run_cmd(cmd, description):
    """Run a shell command and handle errors."""
    print(f"  Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: {description} failed with code {result.returncode}")
        if result.stderr:
            print(f"  stderr: {result.stderr}")
        sys.exit(1)
    return result


def count_lines(filepath):
    """Count lines in a file."""
    result = subprocess.run(f"wc -l < {filepath}", shell=True, capture_output=True, text=True)
    return int(result.stdout.strip())


def fetch_pfam_domains(config_file: str, output_file: str):
    """Fetch all Pfam domain regions from pfamlive database."""
    print("\nStep 1a: Fetching Pfam domain regions from pfamlive...")
    config_file = os.path.expanduser(config_file)
    query = "SELECT pfamseq_acc, seq_start, seq_end FROM pfamA_reg_full_significant WHERE in_full = 1"
    cmd = f"mysql --defaults-file={config_file} pfam_live --quick --silent --skip-column-names -e \"{query}\" > {output_file}"
    run_cmd(cmd, "MySQL query")
    print(f"  Fetched {count_lines(output_file):,} Pfam domain regions")


def fetch_pfamseq_accessions(config_file: str, output_file: str):
    """Fetch all pfamseq accessions from pfamlive database."""
    print("\nStep 1b: Fetching pfamseq accessions from pfamlive...")
    config_file = os.path.expanduser(config_file)
    query = "SELECT pfamseq_acc FROM pfamseq"
    cmd = f"mysql --defaults-file={config_file} pfam_live --quick --silent --skip-column-names -e \"{query}\" > {output_file}"
    run_cmd(cmd, "MySQL query")
    print(f"  Fetched {count_lines(output_file):,} pfamseq accessions")


def extract_ted_high(ted_file: str, output_file: str):
    """Extract high-confidence TED domains using zgrep and awk."""
    print("\nStep 2: Extracting high-confidence TED domains...")
    # Output: acc, start, end, suffix, segments
    cmd = f"""zgrep $'\\thigh$' {ted_file} | awk -F'\\t' '{{
        split($1, a, "-"); acc = a[2];
        n = split($1, b, "_"); ted_suffix = b[n];
        split($2, segs, "_");
        split(segs[1], first, "-"); start = first[1];
        n_segs = split($2, segs, "_");
        split(segs[n_segs], last, "-"); end = last[2];
        print acc "\\t" start "\\t" end "\\t" ted_suffix "\\t" $2
    }}' > {output_file}"""
    run_cmd(cmd, "TED extraction")
    print(f"  Extracted {count_lines(output_file):,} high-confidence TED domains")


def sort_file(input_file: str, output_file: str, memory: str, description: str):
    """Sort file by first column."""
    print(f"\n{description}")
    cmd = f"sort -S {memory} -k1 {input_file} > {output_file}"
    run_cmd(cmd, "Sort")
    print("  Sort completed")


def filter_ted_by_pfamseq(ted_sorted: str, pfamseq_sorted: str, output_file: str):
    """Filter TED domains to only those with accessions in pfamseq."""
    print("\nStep 4: Filtering TED domains to pfamseq proteins only...")
    # Use join to keep only TED entries where accession is in pfamseq
    # join -t $'\t' matches on first field by default
    cmd = f"join -t $'\\t' {pfamseq_sorted} {ted_sorted} > {output_file}"
    run_cmd(cmd, "Filter by pfamseq")
    print(f"  Filtered to {count_lines(output_file):,} TED domains in pfamseq")


def parse_chopping(chopping):
    """Parse chopping string into list of (start, end) tuples."""
    segments = []
    for segment in chopping.split('_'):
        if '-' in segment:
            try:
                start, end = segment.split('-')
                segments.append((int(start), int(end)))
            except ValueError:
                continue
    return segments


def segments_overlap(ted_segments, pfam_regions):
    """Check if any TED segment overlaps with any Pfam region."""
    for ted_start, ted_end in ted_segments:
        for pfam_start, pfam_end in pfam_regions:
            if max(ted_start, pfam_start) <= min(ted_end, pfam_end):
                return True
    return False


def merge_and_filter(pfam_sorted: str, ted_sorted: str, output_file: str):
    """
    Merge two sorted files and filter TED domains that don't overlap Pfam.

    Reads both files in parallel, processing one protein at a time.
    """
    print("\nStep 5: Comparing TED vs Pfam domains (one protein at a time)...")

    proteins_processed = 0
    ted_total = 0
    ted_kept = 0
    ted_overlap = 0

    pfam_f = open(pfam_sorted, 'r')
    ted_f = open(ted_sorted, 'r')
    out_f = open(output_file, 'w')

    # Current lines from each file
    pfam_line = pfam_f.readline()
    ted_line = ted_f.readline()

    def get_acc(line, source):
        if not line:
            return None
        return line.split('\t')[0]

    def parse_pfam(line):
        parts = line.rstrip('\n').split('\t')
        return parts[0], int(parts[1]), int(parts[2])

    def parse_ted(line):
        parts = line.rstrip('\n').split('\t')
        return parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4]

    while pfam_line or ted_line:
        pfam_acc = get_acc(pfam_line, 'P')
        ted_acc = get_acc(ted_line, 'T')

        # Determine which protein to process next
        if pfam_acc is None:
            current_acc = ted_acc
        elif ted_acc is None:
            current_acc = pfam_acc
        else:
            current_acc = min(pfam_acc, ted_acc)

        if current_acc is None:
            break

        # Collect all Pfam domains for current protein
        pfam_regions = []
        while pfam_line and get_acc(pfam_line, 'P') == current_acc:
            _, start, end = parse_pfam(pfam_line)
            pfam_regions.append((start, end))
            pfam_line = pfam_f.readline()

        # Collect and process all TED domains for current protein
        while ted_line and get_acc(ted_line, 'T') == current_acc:
            acc, start, end, suffix, segments_str = parse_ted(ted_line)
            ted_total += 1

            ted_segments = parse_chopping(segments_str)
            if not ted_segments:
                ted_segments = [(start, end)]

            if pfam_regions and segments_overlap(ted_segments, pfam_regions):
                ted_overlap += 1
            else:
                out_f.write(f"{acc}\t{start}\t{end}\t{suffix}\n")
                ted_kept += 1

            ted_line = ted_f.readline()

        proteins_processed += 1
        if proteins_processed % 10000000 == 0:
            print(f"  {proteins_processed:,} proteins, kept {ted_kept:,}/{ted_total:,} TED...")

    pfam_f.close()
    ted_f.close()
    out_f.close()

    print(f"\n{'='*60}")
    print("COMPLETE")
    print(f"{'='*60}")
    print(f"Proteins processed:            {proteins_processed:,}")
    print(f"TED domains checked:           {ted_total:,}")
    print(f"TED domains with overlap:      {ted_overlap:,}")
    print(f"TED domains kept:              {ted_kept:,}")
    print(f"{'='*60}")
    print(f"Output: {output_file}")


def main():
    args = parse_arguments()

    print(f"\n{'='*60}")
    print("TED Domain Filter - Remove Pfam Overlapping Domains")
    print(f"{'='*60}")

    work_dir = os.path.abspath(args.work_dir)

    # File paths
    pfam_cache = os.path.join(work_dir, args.pfam_cache)
    pfam_sorted = os.path.join(work_dir, 'pfam_domains_sorted.tsv')
    pfamseq_cache = os.path.join(work_dir, 'pfamseq_accessions.tsv')
    pfamseq_sorted = os.path.join(work_dir, 'pfamseq_accessions_sorted.tsv')
    ted_extracted = os.path.join(work_dir, 'ted_high_extracted.tsv')
    ted_sorted = os.path.join(work_dir, 'ted_high_sorted.tsv')
    ted_filtered = os.path.join(work_dir, 'ted_high_in_pfamseq.tsv')
    output_file = os.path.join(work_dir, args.output)

    # Step 1a: Pfam domains
    if os.path.exists(pfam_cache):
        print(f"\nStep 1a: Using existing Pfam cache ({count_lines(pfam_cache):,} domains)")
    else:
        fetch_pfam_domains(args.mysql_config, pfam_cache)

    # Step 1b: Pfamseq accessions
    if os.path.exists(pfamseq_cache):
        print(f"\nStep 1b: Using existing pfamseq cache ({count_lines(pfamseq_cache):,} accessions)")
    else:
        fetch_pfamseq_accessions(args.mysql_config, pfamseq_cache)

    # Step 2a: Extract TED high-confidence
    if os.path.exists(ted_extracted):
        print(f"\nStep 2a: Using existing TED extraction ({count_lines(ted_extracted):,} domains)")
    else:
        extract_ted_high(args.ted_file, ted_extracted)

    # Step 2b: Sort pfamseq accessions
    if os.path.exists(pfamseq_sorted):
        print(f"\nStep 2b: Using existing sorted pfamseq ({count_lines(pfamseq_sorted):,} accessions)")
    else:
        sort_file(pfamseq_cache, pfamseq_sorted, args.sort_memory, "Step 2b: Sorting pfamseq accessions...")

    # Step 3a: Sort Pfam file
    if os.path.exists(pfam_sorted):
        print(f"\nStep 3a: Using existing sorted Pfam file ({count_lines(pfam_sorted):,} domains)")
    else:
        sort_file(pfam_cache, pfam_sorted, args.sort_memory, "Step 3a: Sorting Pfam domains...")

    # Step 3b: Sort TED file
    if os.path.exists(ted_sorted):
        print(f"\nStep 3b: Using existing sorted TED file ({count_lines(ted_sorted):,} domains)")
    else:
        sort_file(ted_extracted, ted_sorted, args.sort_memory, "Step 3b: Sorting TED domains...")

    # Step 4: Filter TED to pfamseq proteins only
    if os.path.exists(ted_filtered):
        print(f"\nStep 4: Using existing filtered TED file ({count_lines(ted_filtered):,} domains)")
    else:
        filter_ted_by_pfamseq(ted_sorted, pfamseq_sorted, ted_filtered)

    # Step 5: Merge and filter
    merge_and_filter(pfam_sorted, ted_filtered, output_file)

    # Cleanup only if requested
    if args.cleanup:
        for f in [ted_extracted, pfam_sorted, ted_sorted, pfamseq_cache, pfamseq_sorted, ted_filtered]:
            if os.path.exists(f):
                os.remove(f)
                print(f"Cleaned up: {f}")

    print("\nDone!")


if __name__ == '__main__':
    main()
