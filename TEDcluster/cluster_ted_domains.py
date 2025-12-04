#!/usr/bin/env python3
"""
Cluster TED domains using mmseqs2.

This script:
1. Extracts domain sequences from UniProt using pfetch (in batches via xargs)
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

# Default paths
DEFAULT_INPUT = 'ted_high_no_pfam_overlap.tsv'
DEFAULT_FASTA = 'ted_domains.fasta'
DEFAULT_MSA_DIR = 'ted_msas'


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
        '--parallel-pfetch',
        type=int,
        default=32,
        help='Number of parallel pfetch processes (default: 32)'
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


def extract_sequences(input_file, output_fasta, parallel, work_dir):
    """
    Extract all domain sequences using pfetch via GNU parallel.
    Uses fetch_domain_seq.sh helper script for cleaner processing.
    """
    print(f"\nStep 1: Extracting domain sequences...")
    print(f"  Input: {input_file}")
    print(f"  Parallel pfetch processes: {parallel}")

    total = count_lines(input_file)
    print(f"  Total domains: {total:,}")

    # Find helper script (same directory as this script)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    helper_script = os.path.join(script_dir, 'fetch_domain_seq.sh')

    if not os.path.exists(helper_script):
        print(f"ERROR: Helper script not found: {helper_script}")
        sys.exit(1)

    print(f"  Running pfetch in parallel (this will take a while for {total:,} sequences)...")
    print(f"  Helper script: {helper_script}")

    # Use parallel with the helper script
    # Input format: acc\tstart\tend\tsuffix
    cmd = f"parallel --colsep '\\t' -j {parallel} --progress {helper_script} {{1}} {{2}} {{3}} {{4}} < {input_file} > {output_fasta}"

    run_cmd(cmd, "pfetch extraction")

    # Count sequences
    seq_count_result = subprocess.run(
        f"grep -c '^>' {output_fasta}", shell=True, capture_output=True, text=True
    )
    seq_count = int(seq_count_result.stdout.strip()) if seq_count_result.returncode == 0 else 0

    print(f"  Extracted {seq_count:,} sequences")
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

    # Check input exists
    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found: {input_file}")
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
        extract_sequences(input_file, fasta_file, args.parallel_pfetch, work_dir)

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
