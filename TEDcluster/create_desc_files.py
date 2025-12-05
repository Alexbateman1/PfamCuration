#!/usr/bin/env python3
"""
Create template DESC files for TED cluster curation directories.

For each cluster directory, creates a DESC file with standard Pfam format.
The SE (source) line is set to the TED domain name (directory name).

Usage:
    python create_desc_files.py [options]
"""

import argparse
import os
import sys


DESC_TEMPLATE = """ID   ShortName
DE   Family description
AU   Bateman A;0000-0002-6982-4660
SE   TED:{ted_name}
GA   27.00 27.00;
TC   27.00 27.00;
NC   27.00 27.00;
BM   hmmbuild  -o /dev/null HMM SEED
SM   hmmsearch -Z 90746521 -E 1000 --cpu 8 HMM pfamseq
TP   Domain
"""


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Create DESC files for TED cluster directories'
    )
    parser.add_argument(
        '--input-dir',
        default='cluster_alignments',
        help='Directory containing cluster subdirectories (default: cluster_alignments)'
    )
    parser.add_argument(
        '--overwrite',
        action='store_true',
        help='Overwrite existing DESC files'
    )
    parser.add_argument(
        '--work-dir',
        default='.',
        help='Working directory (default: current dir)'
    )
    return parser.parse_args()


def main():
    args = parse_arguments()

    work_dir = os.path.abspath(args.work_dir)
    input_dir = os.path.join(work_dir, args.input_dir) if not os.path.isabs(args.input_dir) else args.input_dir

    print(f"\n{'='*60}")
    print("TED Cluster DESC File Generator")
    print(f"{'='*60}")

    if not os.path.exists(input_dir):
        print(f"ERROR: Input directory not found: {input_dir}")
        sys.exit(1)

    # Get list of cluster directories
    cluster_dirs = [d for d in os.listdir(input_dir)
                    if os.path.isdir(os.path.join(input_dir, d))]

    print(f"Found {len(cluster_dirs):,} cluster directories")

    created = 0
    skipped = 0

    for cluster_name in sorted(cluster_dirs):
        cluster_path = os.path.join(input_dir, cluster_name)
        desc_file = os.path.join(cluster_path, 'DESC')

        # Check if DESC already exists
        if os.path.exists(desc_file) and not args.overwrite:
            skipped += 1
            continue

        # Create DESC file
        desc_content = DESC_TEMPLATE.format(ted_name=cluster_name).lstrip('\n')

        with open(desc_file, 'w') as f:
            f.write(desc_content)

        created += 1

    print(f"\n{'='*60}")
    print("COMPLETE")
    print(f"{'='*60}")
    print(f"DESC files created: {created:,}")
    print(f"Skipped (already exist): {skipped:,}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
