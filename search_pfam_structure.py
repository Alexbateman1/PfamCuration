#!/usr/bin/env python3
"""
Search a new Pfam family against the existing library of Pfam AlphaFold structures.

This script:
1. Reads a SEED alignment from a Pfam curation directory
2. Downloads AlphaFold v6 models for sequences in the SEED
3. Verifies sequence matches and calculates pLDDT scores
4. Selects the best model (highest mean pLDDT over alignment region)
5. Chops the model to the alignment boundaries
6. Searches against the Pfam AlphaFold structure library using Foldseek
7. Filters results (60% overlap, E-value <= 1e-3)
8. Outputs significant matches
"""

import argparse
import os
import sys
import re
import urllib.request
import subprocess
import tempfile
from pathlib import Path
from Bio.PDB import MMCIFParser, MMCIFIO, Select
from collections import defaultdict


def parse_seed_alignment(seed_file):
    """
    Parse a Pfam SEED alignment file and extract UniProt accessions with regions.

    Returns:
        list of tuples: (uniprot_acc, start, end, sequence)
    """
    sequences = []

    with open(seed_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if not line or line.startswith('#') or line.startswith('//'):
                continue

            # Match lines like: P39607.2/21-100  SEQUENCE...
            match = re.match(r'^(\w+)\.?\d*/(\d+)-(\d+)\s+(.+)$', line)
            if match:
                uniprot_acc = match.group(1)
                start = int(match.group(2))
                end = int(match.group(3))
                sequence = match.group(4).replace('.', '').replace('-', '')
                sequences.append((uniprot_acc, start, end, sequence))

    return sequences


def download_alphafold_model(uniprot_acc, output_dir):
    """
    Download AlphaFold v6 CIF model for a given UniProt accession.

    Returns:
        str: path to downloaded CIF file, or None if failed
    """
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_acc}-F1-model_v6.cif"
    output_file = os.path.join(output_dir, f"AF-{uniprot_acc}-F1-model_v6.cif")

    try:
        print(f"Downloading {url}...")
        urllib.request.urlretrieve(url, output_file)
        return output_file
    except Exception as e:
        print(f"Failed to download {uniprot_acc}: {e}")
        return None


def get_sequence_from_cif(cif_file):
    """
    Extract the protein sequence from a CIF file.

    Returns:
        str: amino acid sequence
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('protein', cif_file)

    # Get the first chain
    chain = next(structure.get_chains())

    # Extract sequence
    aa_map = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }

    sequence = []
    for residue in chain.get_residues():
        if residue.get_resname() in aa_map:
            sequence.append(aa_map[residue.get_resname()])

    return ''.join(sequence)


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
    return 0


def verify_sequence_match(alignment_seq, structure_seq, start, end):
    """
    Verify that the alignment sequence matches the structure sequence.

    Returns:
        bool: True if sequences match
    """
    # Extract the relevant region from structure sequence
    structure_region = structure_seq[start-1:end]

    # Compare (allowing for some flexibility)
    return alignment_seq == structure_region


class RegionSelect(Select):
    """Select only residues within a specific region."""

    def __init__(self, start, end):
        self.start = start
        self.end = end

    def accept_residue(self, residue):
        res_id = residue.get_id()[1]
        return self.start <= res_id <= self.end


def chop_structure(input_cif, output_cif, start, end):
    """
    Chop a structure to keep only residues in the specified region.

    Args:
        input_cif: input CIF file path
        output_cif: output CIF file path
        start: start residue (1-indexed)
        end: end residue (1-indexed)
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('protein', input_cif)

    io = MMCIFIO()
    io.set_structure(structure)
    io.save(output_cif, RegionSelect(start, end))


def run_foldseek_search(query_cif, database_dir, output_file, tmp_dir):
    """
    Run foldseek easy-search against the Pfam structure database.

    Returns:
        str: path to foldseek output file
    """
    foldseek_output = os.path.join(tmp_dir, 'foldseek_raw.out')
    log_file = os.path.join(tmp_dir, 'foldseek.log')

    cmd = [
        'foldseek', 'easy-search',
        query_cif,
        database_dir,
        foldseek_output,
        tmp_dir
    ]

    print(f"Running foldseek search...")
    with open(log_file, 'w') as log:
        result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)

    if result.returncode != 0:
        raise RuntimeError(f"Foldseek search failed. Check {log_file}")

    return foldseek_output


def parse_and_filter_results(foldseek_output, curation_dir, output_file,
                             min_overlap=0.6, max_evalue=1e-3):
    """
    Parse foldseek results and filter by overlap and E-value.

    Outputs a TSV file with columns:
    - Curation directory
    - Query info
    - Target Pfam accession
    - Identity
    - Alignment length
    - Query length
    - Target length
    - Overlap (fraction)
    - E-value
    - Bitscore
    """
    curation_name = os.path.basename(os.path.normpath(curation_dir))

    results = []

    with open(foldseek_output, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue

            query = fields[0]
            target = fields[1]
            identity = float(fields[2])
            alnlen = int(fields[3])
            qlen = int(fields[4]) if len(fields) > 4 else 0
            tlen = int(fields[5]) if len(fields) > 5 else 0
            evalue = float(fields[10])
            bitscore = float(fields[11]) if len(fields) > 11 else 0

            # Calculate overlap as fraction of shorter sequence covered
            min_len = min(qlen, tlen) if qlen and tlen else alnlen
            overlap = alnlen / min_len if min_len > 0 else 0

            # Extract Pfam accession from target (format: PF00001_...cif)
            pfam_match = re.match(r'(PF\d+)', target)
            pfam_acc = pfam_match.group(1) if pfam_match else target

            # Apply filters
            if overlap >= min_overlap and evalue <= max_evalue:
                results.append({
                    'curation_dir': curation_name,
                    'query': query,
                    'pfam_acc': pfam_acc,
                    'target_full': target,
                    'identity': identity,
                    'alnlen': alnlen,
                    'qlen': qlen,
                    'tlen': tlen,
                    'overlap': overlap,
                    'evalue': evalue,
                    'bitscore': bitscore
                })

    # Sort by E-value
    results.sort(key=lambda x: x['evalue'])

    # Write output
    with open(output_file, 'w') as out:
        # Header
        out.write('\t'.join([
            'Curation_Dir', 'Query', 'Pfam_Accession', 'Target_Full',
            'Identity', 'Aln_Length', 'Query_Length', 'Target_Length',
            'Overlap', 'E-value', 'Bitscore'
        ]) + '\n')

        # Data
        for r in results:
            out.write('\t'.join([
                r['curation_dir'],
                r['query'],
                r['pfam_acc'],
                r['target_full'],
                f"{r['identity']:.2f}",
                str(r['alnlen']),
                str(r['qlen']),
                str(r['tlen']),
                f"{r['overlap']:.3f}",
                f"{r['evalue']:.2e}",
                f"{r['bitscore']:.1f}"
            ]) + '\n')

    return len(results)


def main():
    parser = argparse.ArgumentParser(
        description='Search a new Pfam family against existing Pfam structure library'
    )
    parser.add_argument(
        'curation_dir',
        help='Path to Pfam curation directory containing SEED file'
    )
    parser.add_argument(
        '--database',
        default='/nfs/production/agb/interpro/users/typhaine/interpro-pfam-curation-tools/pfam_add_clan_search/results/chopped_cif',
        help='Path to Pfam AlphaFold structure database (default: chopped_cif location)'
    )
    parser.add_argument(
        '--output',
        default='foldseek',
        help='Output file name (default: foldseek)'
    )

    args = parser.parse_args()

    # Check if output already exists
    output_file = os.path.join(args.curation_dir, args.output)
    if os.path.exists(output_file):
        print(f"Output file {output_file} already exists. Skipping to avoid duplication.")
        sys.exit(0)

    # Check SEED file exists
    seed_file = os.path.join(args.curation_dir, 'SEED')
    if not os.path.exists(seed_file):
        print(f"Error: SEED file not found at {seed_file}")
        sys.exit(1)

    print(f"Processing curation directory: {args.curation_dir}")

    # Parse SEED alignment
    print("Parsing SEED alignment...")
    sequences = parse_seed_alignment(seed_file)
    print(f"Found {len(sequences)} sequences in SEED alignment")

    if not sequences:
        print("Error: No sequences found in SEED file")
        sys.exit(1)

    # Create temporary directory for downloads
    with tempfile.TemporaryDirectory() as tmp_dir:
        best_model = None
        best_plddt = 0
        best_info = None

        # Process each sequence
        for uniprot_acc, start, end, alignment_seq in sequences:
            print(f"\nProcessing {uniprot_acc} ({start}-{end})...")

            # Download AlphaFold model
            cif_file = download_alphafold_model(uniprot_acc, tmp_dir)
            if not cif_file:
                continue

            # Verify sequence match
            structure_seq = get_sequence_from_cif(cif_file)
            if not verify_sequence_match(alignment_seq, structure_seq, start, end):
                print(f"Warning: Sequence mismatch for {uniprot_acc}. Skipping.")
                print(f"  Alignment: {alignment_seq[:50]}...")
                print(f"  Structure: {structure_seq[start-1:end][:50]}...")
                continue

            # Calculate mean pLDDT
            mean_plddt = calculate_mean_plddt(cif_file, start, end)
            print(f"  Mean pLDDT: {mean_plddt:.2f}")

            if mean_plddt > best_plddt:
                best_plddt = mean_plddt
                best_model = cif_file
                best_info = (uniprot_acc, start, end)

        if not best_model:
            print("\nError: No valid models found with matching sequences")
            sys.exit(1)

        print(f"\nBest model: {best_info[0]} ({best_info[1]}-{best_info[2]}) "
              f"with mean pLDDT {best_plddt:.2f}")

        # Chop the best model to alignment boundaries
        chopped_model = os.path.join(tmp_dir, 'query_chopped.cif')
        print(f"Chopping structure to region {best_info[1]}-{best_info[2]}...")
        chop_structure(best_model, chopped_model, best_info[1], best_info[2])

        # Run foldseek search
        print(f"\nSearching against database: {args.database}")
        foldseek_output = run_foldseek_search(
            chopped_model,
            args.database,
            output_file,
            tmp_dir
        )

        # Parse and filter results
        print("\nParsing and filtering results...")
        num_hits = parse_and_filter_results(
            foldseek_output,
            args.curation_dir,
            output_file
        )

        print(f"\nDone! Found {num_hits} significant matches.")
        print(f"Results written to: {output_file}")


if __name__ == '__main__':
    main()
