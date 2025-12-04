#!/usr/bin/env python3
"""
Generate alignments for TED domain clusters with >9 sequences.

Uses the original sequences from ted_domains.fasta (correct coordinates)
and runs create_alignment.pl to build proper MSAs.

Usage:
    python generate_cluster_alignments.py [options]
"""

import argparse
import os
import subprocess
import sys
from collections import defaultdict


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Generate alignments for large TED clusters'
    )
    parser.add_argument(
        '--cluster-tsv',
        default='ted_clu.tsv',
        help='Cluster TSV file from mmseqs (default: ted_clu.tsv)'
    )
    parser.add_argument(
        '--fasta',
        default='ted_domains.fasta',
        help='Original domain sequences FASTA (default: ted_domains.fasta)'
    )
    parser.add_argument(
        '--output-dir',
        default='cluster_alignments',
        help='Output directory for alignments (default: cluster_alignments)'
    )
    parser.add_argument(
        '--min-seqs',
        type=int,
        default=10,
        help='Minimum sequences per cluster (default: 10, i.e. >9)'
    )
    parser.add_argument(
        '--work-dir',
        default='.',
        help='Working directory (default: current dir)'
    )
    return parser.parse_args()


def load_clusters(cluster_tsv, min_seqs):
    """
    Load clusters from TSV and filter by minimum sequence count.

    Returns dict: rep_id -> list of member_ids (including rep)
    """
    print(f"Loading clusters from {cluster_tsv}...")

    clusters = defaultdict(list)

    with open(cluster_tsv, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2:
                rep_id = parts[0]
                member_id = parts[1]
                clusters[rep_id].append(member_id)

    # Filter by minimum size
    large_clusters = {rep: members for rep, members in clusters.items()
                      if len(members) >= min_seqs}

    print(f"  Total clusters: {len(clusters):,}")
    print(f"  Clusters with >= {min_seqs} sequences: {len(large_clusters):,}")

    return large_clusters


def load_fasta_index(fasta_file):
    """
    Build index of sequence positions in FASTA file.

    Returns dict: seq_id -> (file_offset, seq_length)
    """
    print(f"Indexing FASTA file {fasta_file}...")

    index = {}

    with open(fasta_file, 'r') as f:
        current_id = None
        seq_start = None
        seq_lines = []

        while True:
            line = f.readline()
            if not line:
                # End of file - save last sequence
                if current_id:
                    index[current_id] = ''.join(seq_lines)
                break

            if line.startswith('>'):
                # Save previous sequence
                if current_id:
                    index[current_id] = ''.join(seq_lines)

                # Start new sequence
                current_id = line[1:].split()[0]  # Get ID without >
                seq_lines = []
            else:
                seq_lines.append(line.rstrip('\n'))

    print(f"  Indexed {len(index):,} sequences")
    return index


def extract_cluster_fasta(cluster_members, fasta_index, output_fasta):
    """
    Extract sequences for a cluster and write to FASTA file.

    Returns number of sequences written.
    """
    written = 0

    with open(output_fasta, 'w') as f:
        for member_id in cluster_members:
            if member_id in fasta_index:
                seq = fasta_index[member_id]
                f.write(f">{member_id}\n")
                # Write sequence in 60-char lines
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i+60] + '\n')
                written += 1

    return written


def run_alignment(fasta_file, output_seed):
    """
    Run create_alignment.pl to generate alignment.
    """
    cmd = f"create_alignment.pl -fasta {fasta_file} -mu > {output_seed}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        return False, result.stderr
    return True, None


def main():
    args = parse_arguments()

    work_dir = os.path.abspath(args.work_dir)
    cluster_tsv = os.path.join(work_dir, args.cluster_tsv) if not os.path.isabs(args.cluster_tsv) else args.cluster_tsv
    fasta_file = os.path.join(work_dir, args.fasta) if not os.path.isabs(args.fasta) else args.fasta
    output_dir = os.path.join(work_dir, args.output_dir) if not os.path.isabs(args.output_dir) else args.output_dir

    print(f"\n{'='*60}")
    print("TED Cluster Alignment Generator")
    print(f"{'='*60}")

    # Check inputs exist
    if not os.path.exists(cluster_tsv):
        print(f"ERROR: Cluster TSV not found: {cluster_tsv}")
        sys.exit(1)

    if not os.path.exists(fasta_file):
        print(f"ERROR: FASTA file not found: {fasta_file}")
        sys.exit(1)

    # Load clusters
    clusters = load_clusters(cluster_tsv, args.min_seqs)

    if not clusters:
        print("No clusters found with enough sequences.")
        sys.exit(0)

    # Load FASTA sequences into memory
    # For 14M sequences this might use significant RAM, but should be manageable
    fasta_index = load_fasta_index(fasta_file)

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Process each cluster
    print(f"\nGenerating alignments in {output_dir}/...")

    successful = 0
    failed = 0

    for i, (rep_id, members) in enumerate(clusters.items()):
        # Create safe filename from rep_id
        safe_name = rep_id.replace('/', '_').replace('|', '_')

        cluster_fasta = os.path.join(output_dir, f"{safe_name}.fa")
        cluster_seed = os.path.join(output_dir, f"{safe_name}.SEED")

        # Extract sequences
        seq_count = extract_cluster_fasta(members, fasta_index, cluster_fasta)

        if seq_count < args.min_seqs:
            print(f"  Warning: {rep_id} - only found {seq_count}/{len(members)} sequences in FASTA")
            continue

        # Run alignment
        success, error = run_alignment(cluster_fasta, cluster_seed)

        if success:
            successful += 1
        else:
            failed += 1
            print(f"  Failed: {rep_id} - {error}")

        # Progress update
        if (i + 1) % 100 == 0:
            print(f"  Processed {i + 1:,}/{len(clusters):,} clusters...")

    print(f"\n{'='*60}")
    print("COMPLETE")
    print(f"{'='*60}")
    print(f"Successful alignments: {successful:,}")
    print(f"Failed alignments: {failed:,}")
    print(f"Output directory: {output_dir}/")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
