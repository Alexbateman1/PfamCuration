#!/usr/bin/env python3
"""
Generate Pfam-style curation directories for TED domain clusters.

For each cluster with >N sequences, creates a directory named after
the representative (e.g., M3XIJ8_TED01) containing:
  - FA: FASTA sequences for the cluster
  - SEED: Initial multiple sequence alignment

Uses original sequences from ted_domains.fasta (correct coordinates)
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
        description='Generate Pfam-style curation directories for TED clusters'
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
        '--domain-info',
        default='ted_high_no_pfam_overlap.tsv',
        help='Domain info TSV with TED suffixes (default: ted_high_no_pfam_overlap.tsv)'
    )
    parser.add_argument(
        '--output-dir',
        default='cluster_alignments',
        help='Output directory for cluster directories (default: cluster_alignments)'
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


def load_ted_suffixes(domain_info_file):
    """
    Load TED suffixes from domain info file.

    Returns dict: (acc, start, end) -> ted_suffix
    """
    print(f"Loading TED suffixes from {domain_info_file}...")

    suffixes = {}

    with open(domain_info_file, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 4:
                acc = parts[0]
                start = parts[1]
                end = parts[2]
                suffix = parts[3]
                suffixes[(acc, start, end)] = suffix

    print(f"  Loaded {len(suffixes):,} domain entries")
    return suffixes


def parse_seq_id(seq_id):
    """
    Parse sequence ID like 'A0A8J8BLT5.1/1-100' into components.

    Returns: (acc_with_version, acc_base, start, end)
    """
    # Split on /
    if '/' in seq_id:
        acc_part, coords = seq_id.rsplit('/', 1)
        start, end = coords.split('-')
    else:
        acc_part = seq_id
        start, end = None, None

    # Get base accession without version
    acc_base = acc_part.split('.')[0]

    return acc_part, acc_base, start, end


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

    Returns dict: seq_id -> sequence
    """
    print(f"Indexing FASTA file {fasta_file}...")

    index = {}

    with open(fasta_file, 'r') as f:
        current_id = None
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


def get_cluster_name(rep_id, ted_suffixes):
    """
    Generate cluster directory name from representative ID.

    Format: ACC_TEDXX (e.g., M3XIJ8_TED01)
    """
    acc_full, acc_base, start, end = parse_seq_id(rep_id)

    # Look up TED suffix
    suffix = ted_suffixes.get((acc_base, start, end), None)

    if suffix:
        return f"{acc_base}_{suffix}"
    else:
        # Fallback: use accession and coordinates
        return f"{acc_base}_{start}_{end}"


def main():
    args = parse_arguments()

    work_dir = os.path.abspath(args.work_dir)
    cluster_tsv = os.path.join(work_dir, args.cluster_tsv) if not os.path.isabs(args.cluster_tsv) else args.cluster_tsv
    fasta_file = os.path.join(work_dir, args.fasta) if not os.path.isabs(args.fasta) else args.fasta
    domain_info = os.path.join(work_dir, args.domain_info) if not os.path.isabs(args.domain_info) else args.domain_info
    output_dir = os.path.join(work_dir, args.output_dir) if not os.path.isabs(args.output_dir) else args.output_dir

    print(f"\n{'='*60}")
    print("TED Cluster Pfam Directory Generator")
    print(f"{'='*60}")

    # Check inputs exist
    for path, name in [(cluster_tsv, 'Cluster TSV'), (fasta_file, 'FASTA'), (domain_info, 'Domain info')]:
        if not os.path.exists(path):
            print(f"ERROR: {name} not found: {path}")
            sys.exit(1)

    # Load TED suffixes
    ted_suffixes = load_ted_suffixes(domain_info)

    # Load clusters
    clusters = load_clusters(cluster_tsv, args.min_seqs)

    if not clusters:
        print("No clusters found with enough sequences.")
        sys.exit(0)

    # Load FASTA sequences into memory
    fasta_index = load_fasta_index(fasta_file)

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Process each cluster
    print(f"\nGenerating Pfam directories in {output_dir}/...")

    successful = 0
    failed = 0

    for i, (rep_id, members) in enumerate(clusters.items()):
        # Get cluster name (ACC_TEDXX format)
        cluster_name = get_cluster_name(rep_id, ted_suffixes)

        # Create cluster directory
        cluster_dir = os.path.join(output_dir, cluster_name)
        os.makedirs(cluster_dir, exist_ok=True)

        fa_file = os.path.join(cluster_dir, 'FA')
        seed_file = os.path.join(cluster_dir, 'SEED')

        # Extract sequences
        seq_count = extract_cluster_fasta(members, fasta_index, fa_file)

        if seq_count < args.min_seqs:
            print(f"  Warning: {cluster_name} - only found {seq_count}/{len(members)} sequences in FASTA")
            continue

        # Run alignment
        success, error = run_alignment(fa_file, seed_file)

        if success:
            successful += 1
        else:
            failed += 1
            print(f"  Failed: {cluster_name} - {error}")

        # Progress update
        if (i + 1) % 100 == 0:
            print(f"  Processed {i + 1:,}/{len(clusters):,} clusters...")

    print(f"\n{'='*60}")
    print("COMPLETE")
    print(f"{'='*60}")
    print(f"Successful alignments: {successful:,}")
    print(f"Failed alignments: {failed:,}")
    print(f"Output directory: {output_dir}/")
    print(f"  Each subdirectory contains FA and SEED files")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
