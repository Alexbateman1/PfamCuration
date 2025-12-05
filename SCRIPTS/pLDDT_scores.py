#!/usr/bin/env python3
"""
Calculate mean pLDDT scores for sequences in a SEED alignment.

Reads a SEED alignment, fetches AlphaFold models for each sequence,
and calculates the mean pLDDT for the domain region.

Usage:
    python pLDDT_scores.py [SEED_file]

Output:
    pLDDT.scores - scores file with format: score accession
"""

import argparse
import os
import sys
import re
import tempfile
import urllib.request


def parse_seed_alignment(seed_file):
    """
    Parse a Pfam SEED alignment file and extract UniProt accessions with regions.

    Returns:
        list of tuples: (uniprot_acc, start, end, full_id)
    """
    sequences = []

    with open(seed_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if not line or line.startswith('#') or line.startswith('//'):
                continue

            # Match lines like: P39607.2/21-100  SEQUENCE...
            match = re.match(r'^(\S+)/(\d+)-(\d+)\s+', line)
            if match:
                full_acc = match.group(0).split()[0]  # e.g., P39607.2/21-100
                uniprot_acc = match.group(1)  # e.g., P39607.2
                start = int(match.group(2))
                end = int(match.group(3))
                sequences.append((uniprot_acc, start, end, full_acc))

    return sequences


def download_alphafold_model(uniprot_acc, output_dir):
    """
    Download AlphaFold CIF model for a given UniProt accession.

    Args:
        uniprot_acc: UniProt accession (may include version like "P12345.2")
        output_dir: directory to save the file

    Returns:
        str: path to downloaded CIF file, or None if failed
    """
    # Strip version number from accession for AlphaFold API (P12345.2 -> P12345)
    base_acc = uniprot_acc.split('.')[0]

    url = f"https://alphafold.ebi.ac.uk/files/AF-{base_acc}-F1-model_v4.cif"
    output_file = os.path.join(output_dir, f"AF-{base_acc}-F1-model_v4.cif")

    try:
        urllib.request.urlretrieve(url, output_file)
        return output_file
    except Exception:
        # Try v3 as fallback
        try:
            url = f"https://alphafold.ebi.ac.uk/files/AF-{base_acc}-F1-model_v3.cif"
            urllib.request.urlretrieve(url, output_file)
            return output_file
        except Exception:
            return None


def calculate_mean_plddt(cif_file, start, end):
    """
    Calculate mean pLDDT score for a specific region in the structure.

    Args:
        cif_file: path to CIF file
        start: start residue position (1-indexed)
        end: end residue position (1-indexed)

    Returns:
        float: mean pLDDT score, or 0 if failed
    """
    try:
        from Bio.PDB import MMCIFParser

        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure('protein', cif_file)

        chain = next(structure.get_chains())

        plddt_scores = []
        for residue in chain.get_residues():
            res_id = residue.get_id()[1]
            if start <= res_id <= end:
                # pLDDT is stored in the B-factor field
                for atom in residue.get_atoms():
                    plddt_scores.append(atom.get_bfactor())
                    break  # Only need one atom per residue

        if plddt_scores:
            return sum(plddt_scores) / len(plddt_scores)
        return 0.0

    except Exception:
        return 0.0


def main():
    parser = argparse.ArgumentParser(
        description='Calculate mean pLDDT scores for SEED alignment sequences'
    )
    parser.add_argument(
        'seed_file',
        nargs='?',
        default='SEED',
        help='Path to SEED alignment file (default: SEED)'
    )
    parser.add_argument(
        '--output',
        default='pLDDT.scores',
        help='Output scores file (default: pLDDT.scores)'
    )

    args = parser.parse_args()

    # Check SEED file exists
    if not os.path.exists(args.seed_file):
        print(f"Error: SEED file not found: {args.seed_file}", file=sys.stderr)
        sys.exit(1)

    # Parse SEED alignment
    sequences = parse_seed_alignment(args.seed_file)

    if not sequences:
        print("Error: No sequences found in SEED file", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(sequences)} sequences in SEED alignment", file=sys.stderr)

    # Process each sequence
    scores = []

    with tempfile.TemporaryDirectory() as tmp_dir:
        for i, (uniprot_acc, start, end, full_id) in enumerate(sequences, 1):
            print(f"Processing {i}/{len(sequences)}: {full_id}...", file=sys.stderr)

            # Download AlphaFold model
            cif_file = download_alphafold_model(uniprot_acc, tmp_dir)

            if cif_file:
                # Calculate mean pLDDT
                mean_plddt = calculate_mean_plddt(cif_file, start, end)
            else:
                mean_plddt = 0.0

            scores.append((mean_plddt, full_id))

    # Write output scores file
    with open(args.output, 'w') as f:
        for score, acc in scores:
            f.write(f"{score:.1f} {acc}\n")

    print(f"Wrote {len(scores)} scores to {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
