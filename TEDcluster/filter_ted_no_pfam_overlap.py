#!/usr/bin/env python3
"""
Filter TED domains that don't overlap with existing Pfam domains.

Uses fast shell commands (zgrep, awk, sort) for processing large files.

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
        default='8G',
        help='Memory limit for Unix sort (default: 8G)'
    )
    parser.add_argument(
        '--sort-parallel',
        type=int,
        default=8,
        help='Number of parallel sort threads (default: 8)'
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


def run_cmd(cmd, description, shell=False):
    """Run a command and handle errors."""
    print(f"  Running: {cmd if shell else ' '.join(cmd)}")
    result = subprocess.run(cmd, shell=shell, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: {description} failed with code {result.returncode}")
        if result.stderr:
            print(f"  stderr: {result.stderr}")
        sys.exit(1)
    return result


def fetch_pfam_domains(config_file: str, output_file: str):
    """Fetch all Pfam domain regions from pfamlive database."""
    print("\nStep 1: Fetching Pfam domain regions from pfamlive...")

    config_file = os.path.expanduser(config_file)

    query = "SELECT pfamseq_acc, seq_start, seq_end FROM pfamA_reg_full_significant WHERE in_full = 1"

    cmd = f"mysql --defaults-file={config_file} pfam_live --quick --silent --skip-column-names -e \"{query}\" > {output_file}"

    run_cmd(cmd, "MySQL query", shell=True)

    # Count lines
    result = subprocess.run(f"wc -l < {output_file}", shell=True, capture_output=True, text=True)
    count = int(result.stdout.strip())
    print(f"  Fetched {count:,} Pfam domain regions")


def extract_ted_high(ted_file: str, output_file: str):
    """Extract high-confidence TED domains using zgrep and awk."""
    print("\nStep 2: Extracting high-confidence TED domains...")

    # Use zgrep to filter 'high' lines, awk to parse and reformat
    # Input: AF-A0A000-F1-model_v4_TED01	11-41_290-389	high
    # Output: A0A000	11	389	TED01	11-41_290-389
    cmd = f"""zgrep $'\\thigh$' {ted_file} | awk -F'\\t' '{{
        # Extract UniProt acc from ted_id (AF-XXX-F1-model_v4_TEDNN)
        split($1, a, "-");
        acc = a[2];

        # Extract TED suffix (last part after underscore)
        n = split($1, b, "_");
        ted_suffix = b[n];

        # Parse chopping to get first start and last end
        split($2, segs, "_");
        split(segs[1], first, "-");
        start = first[1];
        n_segs = split($2, segs, "_");
        split(segs[n_segs], last, "-");
        end = last[2];

        print acc "\\t" start "\\t" end "\\t" ted_suffix "\\t" $2
    }}' > {output_file}"""

    run_cmd(cmd, "TED extraction", shell=True)

    result = subprocess.run(f"wc -l < {output_file}", shell=True, capture_output=True, text=True)
    count = int(result.stdout.strip())
    print(f"  Extracted {count:,} high-confidence TED domains")


def create_combined_file(pfam_file: str, ted_file: str, combined_file: str):
    """Combine Pfam and TED domains into single file with source markers."""
    print("\nStep 3: Creating combined domains file...")

    # Pfam format: acc, start, end -> add P marker
    # TED format: acc, start, end, suffix, segments -> add T marker
    cmd = f"""awk -F'\\t' '{{print $1 "\\t" $2 "\\t" $3 "\\tP\\t"}}' {pfam_file} > {combined_file} && \
awk -F'\\t' '{{print $1 "\\t" $2 "\\t" $3 "\\tT\\t" $4 ":" $5}}' {ted_file} >> {combined_file}"""

    run_cmd(cmd, "Combine files", shell=True)

    result = subprocess.run(f"wc -l < {combined_file}", shell=True, capture_output=True, text=True)
    count = int(result.stdout.strip())
    print(f"  Combined file: {count:,} total domains")


def sort_combined_file(input_file: str, output_file: str, memory: str, parallel: int, work_dir: str):
    """Sort combined file by UniProt accession."""
    print("\nStep 4: Sorting by UniProt accession...")
    print(f"  Memory: {memory}, Parallel: {parallel}")

    cmd = f"sort -t $'\\t' -k1,1 -S {memory} --parallel={parallel} -T {work_dir} {input_file} -o {output_file}"

    run_cmd(cmd, "Sort", shell=True)
    print("  Sort completed")


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


def process_sorted_file(sorted_file: str, output_file: str):
    """Process sorted file one protein at a time."""
    print("\nStep 5: Finding non-overlapping TED domains...")

    proteins_processed = 0
    ted_total = 0
    ted_kept = 0
    ted_overlap = 0

    current_acc = None
    current_pfam_regions = []
    current_ted_domains = []

    out_f = open(output_file, 'w')

    def process_protein():
        nonlocal ted_total, ted_kept, ted_overlap

        if not current_acc:
            return

        for ted_start, ted_end, ted_suffix, segments_str in current_ted_domains:
            ted_total += 1

            ted_segments = parse_chopping(segments_str)
            if not ted_segments:
                ted_segments = [(ted_start, ted_end)]

            if current_pfam_regions and segments_overlap(ted_segments, current_pfam_regions):
                ted_overlap += 1
            else:
                out_f.write(f"{current_acc}\t{ted_start}\t{ted_end}\t{ted_suffix}\n")
                ted_kept += 1

    with open(sorted_file, 'r') as in_f:
        for line in in_f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue

            acc = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            source = parts[3]
            extra = parts[4] if len(parts) > 4 else ''

            if acc != current_acc:
                process_protein()
                current_acc = acc
                current_pfam_regions = []
                current_ted_domains = []
                proteins_processed += 1

                if proteins_processed % 10000000 == 0:
                    print(f"  {proteins_processed:,} proteins, kept {ted_kept:,}/{ted_total:,} TED...")

            if source == 'P':
                current_pfam_regions.append((start, end))
            elif source == 'T':
                if ':' in extra:
                    ted_suffix, segments_str = extra.split(':', 1)
                else:
                    ted_suffix = extra
                    segments_str = f"{start}-{end}"
                current_ted_domains.append((start, end, ted_suffix, segments_str))

        process_protein()

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
    ted_extracted = os.path.join(work_dir, 'ted_high_extracted.tsv')
    combined_file = os.path.join(work_dir, 'combined_domains.tsv')
    sorted_file = os.path.join(work_dir, 'combined_domains_sorted.tsv')
    output_file = os.path.join(work_dir, args.output)

    # Step 1: Pfam domains
    if os.path.exists(pfam_cache):
        result = subprocess.run(f"wc -l < {pfam_cache}", shell=True, capture_output=True, text=True)
        count = int(result.stdout.strip())
        print(f"\nStep 1: Using existing Pfam cache ({count:,} domains)")
    else:
        fetch_pfam_domains(args.mysql_config, pfam_cache)

    # Step 2: Extract TED high-confidence
    if os.path.exists(ted_extracted):
        result = subprocess.run(f"wc -l < {ted_extracted}", shell=True, capture_output=True, text=True)
        count = int(result.stdout.strip())
        print(f"\nStep 2: Using existing TED extraction ({count:,} domains)")
    else:
        extract_ted_high(args.ted_file, ted_extracted)

    # Step 3: Combine files
    if os.path.exists(combined_file):
        result = subprocess.run(f"wc -l < {combined_file}", shell=True, capture_output=True, text=True)
        count = int(result.stdout.strip())
        print(f"\nStep 3: Using existing combined file ({count:,} domains)")
    else:
        create_combined_file(pfam_cache, ted_extracted, combined_file)

    # Step 4: Sort
    if os.path.exists(sorted_file):
        result = subprocess.run(f"wc -l < {sorted_file}", shell=True, capture_output=True, text=True)
        count = int(result.stdout.strip())
        print(f"\nStep 4: Using existing sorted file ({count:,} domains)")
    else:
        sort_combined_file(combined_file, sorted_file, args.sort_memory, args.sort_parallel, work_dir)

    # Step 5: Process
    process_sorted_file(sorted_file, output_file)

    # Cleanup only if requested and successful
    if args.cleanup:
        for f in [ted_extracted, combined_file, sorted_file]:
            if os.path.exists(f):
                os.remove(f)
                print(f"Cleaned up: {f}")

    print("\nDone!")


if __name__ == '__main__':
    main()
