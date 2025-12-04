#!/usr/bin/env python3
"""
Cluster TED domains using mmseqs2.

This script:
1. Extracts domain sequences from pfamseq FASTA file in batches (multiple passes)
2. Creates mmseqs2 database and clusters
3. Generates MSAs for each cluster

Uses batched processing to handle large numbers of domains without OOM:
- Processes ~1M domain regions per batch
- Reads pfamseq FASTA file once per batch
- Appends extracted sequences to output FASTA

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
        '--batch-size',
        type=int,
        default=1000000,
        help='Number of domains per batch (default: 1000000)'
    )
    parser.add_argument(
        '--skip-fetch',
        action='store_true',
        help='Skip sequence fetching (use existing FASTA) - DEPRECATED: auto-detected'
    )
    parser.add_argument(
        '--force-fetch',
        action='store_true',
        help='Force re-extraction of sequences even if FASTA exists'
    )
    parser.add_argument(
        '--skip-cluster',
        action='store_true',
        help='Skip clustering (use existing clusters)'
    )
    parser.add_argument(
        '--split-memory-limit',
        default='16G',
        help='Memory limit for mmseqs split mode (default: 16G)'
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


def count_domains(input_file):
    """Count total number of domains in input file."""
    result = subprocess.run(f"wc -l < {input_file}", shell=True, capture_output=True, text=True)
    return int(result.stdout.strip())


def load_domain_regions_batch(input_file, start_line, batch_size):
    """
    Load a batch of TED domain regions from input file.

    Args:
        input_file: Path to input TSV file
        start_line: Line number to start from (0-indexed)
        batch_size: Maximum number of lines to read

    Returns:
        tuple: (dict of base_acc -> list of (start, end, suffix), count of domains loaded)
    """
    regions = defaultdict(list)
    count = 0

    with open(input_file, 'r') as f:
        # Skip to start line
        for _ in range(start_line):
            next(f, None)

        # Read batch_size lines
        for i, line in enumerate(f):
            if i >= batch_size:
                break
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 4:
                acc = parts[0]  # Base accession without version
                start = int(parts[1])
                end = int(parts[2])
                suffix = parts[3]
                regions[acc].append((start, end, suffix))
                count += 1

    return regions, count


def extract_sequences_batch(pfamseq_file, domain_regions, output_fasta, append=False):
    """
    Extract domain sequences from pfamseq FASTA file for a batch of regions.

    Args:
        pfamseq_file: Path to pfamseq FASTA file
        domain_regions: Dict mapping base_acc -> list of (start, end, suffix)
        output_fasta: Output FASTA file path
        append: If True, append to existing file; otherwise overwrite

    Returns:
        tuple: (sequences_written, proteins_matched, proteins_total)
    """
    sequences_written = 0
    proteins_matched = 0
    proteins_total = 0

    current_acc_full = None  # e.g., A0A8J8BLT5.1
    current_acc_base = None  # e.g., A0A8J8BLT5
    current_seq = []

    mode = 'a' if append else 'w'

    def process_sequence(out_f):
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

    with open(pfamseq_file, 'r') as in_f, open(output_fasta, mode) as out_f:
        for line in in_f:
            if line.startswith('>'):
                # Process previous sequence
                process_sequence(out_f)

                # Parse new header: >A0A8J8BLT5.1 A0A8J8BLT5_9EURY Description
                proteins_total += 1
                header = line[1:].split()[0]  # Get first word after >
                current_acc_full = header  # e.g., A0A8J8BLT5.1
                current_acc_base = header.split('.')[0]  # e.g., A0A8J8BLT5
                current_seq = []
            else:
                current_seq.append(line.rstrip('\n'))

        # Process last sequence
        process_sequence(out_f)

    return sequences_written, proteins_matched, proteins_total


def extract_sequences(input_file, pfamseq_file, output_fasta, batch_size):
    """
    Extract all domain sequences from pfamseq FASTA file using batched processing.

    Processes domains in batches to avoid memory issues with large domain lists.
    Reads through pfamseq file once per batch.
    """
    print(f"\nStep 1: Extracting domain sequences (batched)...")
    print(f"  Input domains: {input_file}")
    print(f"  Pfamseq file: {pfamseq_file}")
    print(f"  Batch size: {batch_size:,} domains")

    # Count total domains
    total_domains = count_domains(input_file)
    num_batches = (total_domains + batch_size - 1) // batch_size
    print(f"  Total domains: {total_domains:,}")
    print(f"  Number of batches: {num_batches}")

    total_sequences = 0
    total_proteins_matched = 0

    for batch_num in range(num_batches):
        start_line = batch_num * batch_size
        print(f"\n  Batch {batch_num + 1}/{num_batches}: domains {start_line + 1:,} to {min(start_line + batch_size, total_domains):,}")

        # Load this batch of domain regions
        print(f"    Loading domain regions...")
        domain_regions, domains_loaded = load_domain_regions_batch(input_file, start_line, batch_size)
        print(f"    Loaded {domains_loaded:,} domains for {len(domain_regions):,} proteins")

        # Extract sequences for this batch
        print(f"    Reading pfamseq FASTA...")
        append = batch_num > 0  # Append after first batch
        seqs_written, prots_matched, prots_total = extract_sequences_batch(
            pfamseq_file, domain_regions, output_fasta, append
        )
        print(f"    Wrote {seqs_written:,} sequences (matched {prots_matched:,} proteins)")

        total_sequences += seqs_written
        total_proteins_matched += prots_matched

        # Clear memory
        del domain_regions

    print(f"\n  Total sequences extracted: {total_sequences:,}")
    print(f"  Output: {output_fasta}")

    return total_sequences


def run_mmseqs_cluster(fasta_file, tmp_dir, min_seq_id, coverage, threads, work_dir, split_memory_limit):
    """Run mmseqs clustering."""
    print(f"\nStep 2: Clustering with mmseqs...")
    print(f"  Min sequence identity: {min_seq_id}")
    print(f"  Coverage: {coverage}")
    print(f"  Threads: {threads}")
    print(f"  Memory limit: {split_memory_limit}")

    db_name = os.path.join(work_dir, 'ted_db')
    clu_name = os.path.join(work_dir, 'ted_clu')

    # Create database
    print("\n  Creating mmseqs database...")
    cmd = f"mmseqs createdb {fasta_file} {db_name}"
    run_cmd(cmd, "createdb")

    # Cluster with memory limit
    print("\n  Clustering (this may take a while)...")
    cmd = f"mmseqs cluster {db_name} {clu_name} {tmp_dir} --min-seq-id {min_seq_id} -c {coverage} --threads {threads} --split-memory-limit {split_memory_limit}"
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


def generate_msas(db_name, clu_name, threads):
    """Generate MSAs for each cluster, keeping in mmseqs database format."""
    print(f"\nStep 3: Generating MSAs...")

    msa_name = f"{clu_name}_msa"

    # Generate MSAs (stays in mmseqs database format - single file)
    print("  Running result2msa...")
    cmd = f"mmseqs result2msa {db_name} {db_name} {clu_name} {msa_name} --msa-format-mode 2 --threads {threads}"
    run_cmd(cmd, "result2msa")

    # Count clusters from the index
    index_file = f"{msa_name}.index"
    if os.path.exists(index_file):
        msa_count = count_lines(index_file)
    else:
        msa_count = 0

    print(f"  Generated {msa_count:,} MSAs (stored in mmseqs database: {msa_name})")
    print(f"\n  To extract a specific MSA by cluster rep ID:")
    print(f"    mmseqs result2flat {db_name} {db_name} {msa_name} <output.afa> --use-fasta-header")
    print(f"\n  To extract all to a single file:")
    print(f"    mmseqs unpackdb {msa_name} <output_dir> --unpack-suffix .afa")

    return msa_count, msa_name


def main():
    args = parse_arguments()

    print(f"\n{'='*60}")
    print("TED Domain Clustering Pipeline")
    print(f"{'='*60}")

    work_dir = os.path.abspath(args.work_dir)
    input_file = os.path.join(work_dir, args.input) if not os.path.isabs(args.input) else args.input
    fasta_file = os.path.join(work_dir, args.fasta) if not os.path.isabs(args.fasta) else args.fasta
    tmp_dir = os.path.join(work_dir, args.tmp_dir) if not os.path.isabs(args.tmp_dir) else args.tmp_dir
    pfamseq_file = args.pfamseq

    # Check input exists
    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found: {input_file}")
        sys.exit(1)

    # Step 1: Extract sequences
    # Auto-detect existing FASTA file unless --force-fetch is specified
    fasta_exists = os.path.exists(fasta_file) and os.path.getsize(fasta_file) > 0

    if fasta_exists and not args.force_fetch:
        seq_result = subprocess.run(f"grep -c '^>' {fasta_file}", shell=True, capture_output=True, text=True)
        seq_count = int(seq_result.stdout.strip()) if seq_result.returncode == 0 else 0
        print(f"\nStep 1: Using existing FASTA: {fasta_file} ({seq_count:,} sequences)")
        print(f"  (use --force-fetch to re-extract)")
    else:
        if not os.path.exists(pfamseq_file):
            print(f"ERROR: Pfamseq file not found: {pfamseq_file}")
            sys.exit(1)
        extract_sequences(input_file, pfamseq_file, fasta_file, args.batch_size)

    # Step 2: Cluster
    if args.skip_cluster:
        print("\nStep 2: Skipping clustering (--skip-cluster)")
        db_name = os.path.join(work_dir, 'ted_db')
        clu_name = os.path.join(work_dir, 'ted_clu')
    else:
        os.makedirs(tmp_dir, exist_ok=True)
        db_name, clu_name = run_mmseqs_cluster(
            fasta_file, tmp_dir, args.min_seq_id, args.coverage, args.threads, work_dir,
            args.split_memory_limit
        )

    # Step 3: Generate MSAs
    msa_count, msa_name = generate_msas(db_name, clu_name, args.threads)

    print(f"\n{'='*60}")
    print("COMPLETE")
    print(f"{'='*60}")
    print(f"FASTA file: {fasta_file}")
    print(f"Cluster TSV: {clu_name}.tsv")
    print(f"MSA database: {msa_name} ({msa_count:,} clusters)")
    print(f"{'='*60}")

    print("\nDone!")


if __name__ == '__main__':
    main()
