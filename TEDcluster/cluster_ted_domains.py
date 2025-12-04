#!/usr/bin/env python3
"""
Cluster TED domains using mmseqs2.

This script:
1. Extracts domain sequences directly from pfamseq FASTA file (single pass)
2. Creates mmseqs2 database and clusters
3. Generates MSAs for each cluster

Usage:
    python cluster_ted_domains.py [options]

Input: ted_high_no_pfam_overlap.tsv (from filter script)
    Format: uniprot_acc<TAB>start<TAB>end<TAB>ted_suffix

Output:
    - ted_domains.fasta: All domain sequences
    - ted_msas/: Directory with MSA files per cluster
"""

import argparse
import os
import subprocess
import sys
from collections import defaultdict

# Default paths
DEFAULT_INPUT = 'ted_high_no_pfam_overlap.tsv'
DEFAULT_FASTA = 'ted_domains.fasta'
DEFAULT_MSA_DIR = 'ted_msas'
DEFAULT_PFAMSEQ = '/nfs/production/agb/pfam/data/pfamseq/pfamseq'


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Cluster TED domains using mmseqs2'
    )
    parser.add_argument(
        '--input',
        default=DEFAULT_INPUT,
        help=f'Input TSV file (default: {DEFAULT_INPUT})'
    )
    parser.add_argument(
        '--fasta',
        default=DEFAULT_FASTA,
        help=f'Output FASTA file (default: {DEFAULT_FASTA})'
    )
    parser.add_argument(
        '--pfamseq',
        default=DEFAULT_PFAMSEQ,
        help=f'Path to pfamseq FASTA file (default: {DEFAULT_PFAMSEQ})'
    )
    parser.add_argument(
        '--msa-dir',
        default=DEFAULT_MSA_DIR,
        help=f'Output MSA directory (default: {DEFAULT_MSA_DIR})'
    )
    parser.add_argument(
        '--min-seq-id',
        type=float,
        default=0.3,
        help='Minimum sequence identity for clustering (default: 0.3)'
    )
    parser.add_argument(
        '--coverage',
        type=float,
        default=0.8,
        help='Minimum coverage for clustering (default: 0.8)'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=8,
        help='Number of threads for mmseqs (default: 8)'
    )
    parser.add_argument(
        '--tmp-dir',
        default='mmseqs_tmp',
        help='Temporary directory for mmseqs (default: mmseqs_tmp)'
    )
    parser.add_argument(
        '--skip-fetch',
        action='store_true',
        help='Skip sequence fetching (use existing FASTA)'
    )
    parser.add_argument(
        '--skip-cluster',
        action='store_true',
        help='Skip clustering (use existing clusters)'
    )
    parser.add_argument(
        '--work-dir',
        default='.',
        help='Working directory (default: current dir)'
    )
    return parser.parse_args()


def run_cmd(cmd, description, check=True):
    """Run a shell command."""
    print(f"  Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if check and result.returncode != 0:
        print(f"ERROR: {description} failed with code {result.returncode}")
        if result.stderr:
            print(f"  stderr: {result.stderr[:1000]}")
        sys.exit(1)
    return result


def count_lines(filepath):
    """Count lines in a file."""
    result = subprocess.run(f"wc -l < {filepath}", shell=True, capture_output=True, text=True)
    return int(result.stdout.strip())


def load_domain_regions(input_file):
    """
    Load TED domain regions from input file.

    Returns dict: base_acc -> list of (start, end, suffix)
    """
    print(f"  Loading domain regions from {input_file}...")

    regions = defaultdict(list)
    count = 0

    with open(input_file, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 4:
                acc = parts[0]  # Base accession without version
                start = int(parts[1])
                end = int(parts[2])
                suffix = parts[3]
                regions[acc].append((start, end, suffix))
                count += 1

    print(f"  Loaded {count:,} domain regions for {len(regions):,} proteins")
    return regions


def extract_sequences_from_pfamseq(pfamseq_file, domain_regions, output_fasta):
    """
    Extract domain sequences from pfamseq FASTA file in a single pass.

    Args:
        pfamseq_file: Path to pfamseq FASTA file
        domain_regions: Dict mapping base_acc -> list of (start, end, suffix)
        output_fasta: Output FASTA file path
    """
    print(f"\n  Reading pfamseq FASTA: {pfamseq_file}")
    print(f"  This may take a while for a large file...")

    sequences_written = 0
    proteins_matched = 0
    proteins_total = 0

    current_acc_full = None  # e.g., A0A8J8BLT5.1
    current_acc_base = None  # e.g., A0A8J8BLT5
    current_seq = []

    def process_sequence():
        """Process current sequence if we have domain regions for it."""
        nonlocal sequences_written, proteins_matched

        if current_acc_base and current_acc_base in domain_regions:
            proteins_matched += 1
            full_seq = ''.join(current_seq)

            for start, end, suffix in domain_regions[current_acc_base]:
                # Extract region (1-indexed, inclusive)
                region_seq = full_seq[start - 1:end]

                if region_seq:
                    # Write with format: >ACC.version/start-end
                    out_f.write(f">{current_acc_full}/{start}-{end}\n")
                    # Write sequence in lines of 60 characters
                    for i in range(0, len(region_seq), 60):
                        out_f.write(region_seq[i:i + 60] + '\n')
                    sequences_written += 1

    with open(pfamseq_file, 'r') as in_f, open(output_fasta, 'w') as out_f:
        for line in in_f:
            if line.startswith('>'):
                # Process previous sequence
                process_sequence()

                # Parse new header: >A0A8J8BLT5.1 A0A8J8BLT5_9EURY Description
                proteins_total += 1
                header = line[1:].split()[0]  # Get first word after >
                current_acc_full = header  # e.g., A0A8J8BLT5.1
                current_acc_base = header.split('.')[0]  # e.g., A0A8J8BLT5
                current_seq = []

                if proteins_total % 10000000 == 0:
                    print(f"    Processed {proteins_total:,} proteins, "
                          f"matched {proteins_matched:,}, wrote {sequences_written:,} sequences...")
            else:
                current_seq.append(line.rstrip('\n'))

        # Process last sequence
        process_sequence()

    print(f"  Processed {proteins_total:,} proteins from pfamseq")
    print(f"  Matched {proteins_matched:,} proteins with TED domains")
    print(f"  Wrote {sequences_written:,} domain sequences")

    return sequences_written


def extract_sequences(input_file, pfamseq_file, output_fasta):
    """
    Extract all domain sequences from pfamseq FASTA file.
    """
    print(f"\nStep 1: Extracting domain sequences...")
    print(f"  Input domains: {input_file}")
    print(f"  Pfamseq file: {pfamseq_file}")

    # Load domain regions
    domain_regions = load_domain_regions(input_file)

    # Extract sequences from pfamseq
    seq_count = extract_sequences_from_pfamseq(pfamseq_file, domain_regions, output_fasta)

    print(f"  Output: {output_fasta}")

    return seq_count


def run_mmseqs_cluster(fasta_file, tmp_dir, min_seq_id, coverage, threads, work_dir):
    """Run mmseqs clustering."""
    print(f"\nStep 2: Clustering with mmseqs...")
    print(f"  Min sequence identity: {min_seq_id}")
    print(f"  Coverage: {coverage}")
    print(f"  Threads: {threads}")

    db_name = os.path.join(work_dir, 'ted_db')
    clu_name = os.path.join(work_dir, 'ted_clu')

    # Create database
    print("\n  Creating mmseqs database...")
    cmd = f"mmseqs createdb {fasta_file} {db_name}"
    run_cmd(cmd, "createdb")

    # Cluster
    print("\n  Clustering (this may take a while)...")
    cmd = f"mmseqs cluster {db_name} {clu_name} {tmp_dir} --min-seq-id {min_seq_id} -c {coverage} --threads {threads}"
    run_cmd(cmd, "cluster")

    # Get cluster statistics
    print("\n  Generating cluster TSV...")
    tsv_file = f"{clu_name}.tsv"
    cmd = f"mmseqs createtsv {db_name} {db_name} {clu_name} {tsv_file}"
    run_cmd(cmd, "createtsv")

    # Count clusters (unique representatives)
    cmd = f"cut -f1 {tsv_file} | sort -u | wc -l"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    cluster_count = int(result.stdout.strip()) if result.returncode == 0 else 0
    print(f"  Found {cluster_count:,} clusters")

    return db_name, clu_name


def generate_msas(db_name, clu_name, msa_dir, threads):
    """Generate MSAs for each cluster."""
    print(f"\nStep 3: Generating MSAs...")

    msa_name = f"{clu_name}_msa"

    # Generate MSAs
    print("  Running result2msa...")
    cmd = f"mmseqs result2msa {db_name} {db_name} {clu_name} {msa_name} --msa-format-mode 2 --threads {threads}"
    run_cmd(cmd, "result2msa")

    # Unpack to individual files
    print(f"\n  Unpacking MSAs to {msa_dir}/...")
    os.makedirs(msa_dir, exist_ok=True)
    cmd = f"mmseqs unpackdb {msa_name} {msa_dir} --unpack-suffix .afa"
    run_cmd(cmd, "unpackdb")

    # Count MSA files
    msa_count = len([f for f in os.listdir(msa_dir) if f.endswith('.afa')])
    print(f"  Generated {msa_count:,} MSA files")

    return msa_count


def main():
    args = parse_arguments()

    print(f"\n{'='*60}")
    print("TED Domain Clustering Pipeline")
    print(f"{'='*60}")

    work_dir = os.path.abspath(args.work_dir)
    input_file = os.path.join(work_dir, args.input) if not os.path.isabs(args.input) else args.input
    fasta_file = os.path.join(work_dir, args.fasta) if not os.path.isabs(args.fasta) else args.fasta
    msa_dir = os.path.join(work_dir, args.msa_dir) if not os.path.isabs(args.msa_dir) else args.msa_dir
    tmp_dir = os.path.join(work_dir, args.tmp_dir) if not os.path.isabs(args.tmp_dir) else args.tmp_dir
    pfamseq_file = args.pfamseq

    # Check input exists
    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found: {input_file}")
        sys.exit(1)

    if not args.skip_fetch and not os.path.exists(pfamseq_file):
        print(f"ERROR: Pfamseq file not found: {pfamseq_file}")
        sys.exit(1)

    # Step 1: Extract sequences
    if args.skip_fetch:
        if not os.path.exists(fasta_file):
            print(f"ERROR: FASTA file not found: {fasta_file}")
            sys.exit(1)
        seq_result = subprocess.run(f"grep -c '^>' {fasta_file}", shell=True, capture_output=True, text=True)
        seq_count = int(seq_result.stdout.strip()) if seq_result.returncode == 0 else 0
        print(f"\nStep 1: Using existing FASTA: {fasta_file} ({seq_count:,} sequences)")
    else:
        extract_sequences(input_file, pfamseq_file, fasta_file)

    # Step 2: Cluster
    if args.skip_cluster:
        print("\nStep 2: Skipping clustering (--skip-cluster)")
        db_name = os.path.join(work_dir, 'ted_db')
        clu_name = os.path.join(work_dir, 'ted_clu')
    else:
        os.makedirs(tmp_dir, exist_ok=True)
        db_name, clu_name = run_mmseqs_cluster(
            fasta_file, tmp_dir, args.min_seq_id, args.coverage, args.threads, work_dir
        )

    # Step 3: Generate MSAs
    msa_count = generate_msas(db_name, clu_name, msa_dir, args.threads)

    print(f"\n{'='*60}")
    print("COMPLETE")
    print(f"{'='*60}")
    print(f"FASTA file: {fasta_file}")
    print(f"Cluster TSV: {clu_name}.tsv")
    print(f"MSA directory: {msa_dir}/ ({msa_count:,} clusters)")
    print(f"{'='*60}")

    print("\nDone!")


if __name__ == '__main__':
    main()
