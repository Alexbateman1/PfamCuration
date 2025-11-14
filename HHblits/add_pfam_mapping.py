#!/usr/bin/env python3
"""
Add Pfam family mapping to HHblits results.
Maps sequence IDs to their source Pfam families.
"""

import sys
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def build_sequence_to_family_map(seed_dir):
    """Build mapping from sequence IDs to Pfam families."""
    seq_to_family = {}
    seed_dir = Path(seed_dir)

    seed_files = list(seed_dir.glob("PF*_SEED"))
    logging.info(f"Building sequence-to-family map from {len(seed_files)} SEED files...")

    for i, seed_file in enumerate(seed_files, 1):
        if i % 1000 == 0:
            logging.info(f"Progress: {i}/{len(seed_files)}")

        family_id = seed_file.stem.replace('_SEED', '')

        try:
            with open(seed_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split(None, 1)
                    if len(parts) >= 2:
                        seq_id = parts[0]
                        seq_to_family[seq_id] = family_id
        except Exception as e:
            logging.warning(f"Failed to read {seed_file}: {e}")

    logging.info(f"Mapped {len(seq_to_family)} sequences to families")
    return seq_to_family

def add_pfam_mapping(input_file, output_file, seq_to_family):
    """Add Pfam family column to results file."""
    with open(input_file, 'r') as inf, open(output_file, 'w') as outf:
        header = inf.readline().strip()
        # Add target_pfam_family column after target_family
        new_header = header.replace('target_family\t', 'target_family\ttarget_pfam_family\t')
        outf.write(new_header + '\n')

        mapped = 0
        unmapped = 0

        for line in inf:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue

            query_family = parts[0]
            target_seq = parts[1]

            # Extract sequence ID (remove /start-end)
            seq_id = target_seq.split('/')[0] if '/' in target_seq else target_seq

            # Look up Pfam family
            target_pfam = seq_to_family.get(seq_id, 'UNKNOWN')
            if target_pfam != 'UNKNOWN':
                mapped += 1
            else:
                unmapped += 1

            # Insert target_pfam_family after target_family
            new_parts = [parts[0], parts[1], target_pfam] + parts[2:]
            outf.write('\t'.join(new_parts) + '\n')

    logging.info(f"Mapped {mapped} hits, {unmapped} unmapped")

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python3 add_pfam_mapping.py <seed_dir> <input_tsv> <output_tsv>")
        sys.exit(1)

    seed_dir = sys.argv[1]
    input_file = sys.argv[2]
    output_file = sys.argv[3]

    # Build mapping
    seq_to_family = build_sequence_to_family_map(seed_dir)

    # Add mapping to results
    add_pfam_mapping(input_file, output_file, seq_to_family)

    logging.info(f"Enhanced results written to {output_file}")
