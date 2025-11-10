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
import shutil
from pathlib import Path
from Bio.PDB import MMCIFParser, MMCIFIO, Select
from collections import defaultdict
import mysql.connector
from configparser import ConfigParser


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


def connect_to_pfam_db(config_file='~/.my.cnf'):
    """
    Connect to the Pfam MySQL database using credentials from config file.

    Args:
        config_file: path to MySQL config file (default: ~/.my.cnf)

    Returns:
        mysql.connector connection object
    """
    config_file = os.path.expanduser(config_file)

    if not os.path.exists(config_file):
        raise FileNotFoundError(f"MySQL config file not found: {config_file}")

    # Parse MySQL config file
    config = ConfigParser()
    config.read(config_file)

    # Get credentials from [client] section
    if 'client' not in config:
        raise ValueError(f"No [client] section found in {config_file}")

    db_config = {
        'host': config.get('client', 'host', fallback='localhost'),
        'user': config.get('client', 'user'),
        'password': config.get('client', 'password', fallback=''),
        'database': 'pfam_live'
    }

    if config.has_option('client', 'port'):
        db_config['port'] = config.getint('client', 'port')

    try:
        connection = mysql.connector.connect(**db_config)
        return connection
    except mysql.connector.Error as e:
        raise RuntimeError(f"Failed to connect to Pfam database: {e}")


def get_pfam_info(pfam_accs, connection=None):
    """
    Query Pfam database to get Pfam ID and clan information for given accessions.

    Args:
        pfam_accs: list of Pfam accessions (e.g., ['PF00001', 'PF00002'])
        connection: mysql connection object (will create new one if None)

    Returns:
        dict: {pfam_acc: {'pfam_id': str, 'clan_acc': str or None}}
    """
    if not pfam_accs:
        return {}

    close_connection = False
    if connection is None:
        connection = connect_to_pfam_db()
        close_connection = True

    cursor = connection.cursor(dictionary=True)

    # Query for Pfam ID and clan information
    placeholders = ','.join(['%s'] * len(pfam_accs))
    query = f"""
        SELECT pfamA.pfamA_acc, pfamA.pfamA_id, clan_membership.clan_acc
        FROM pfamA
        LEFT JOIN clan_membership ON pfamA.pfamA_acc = clan_membership.pfamA_acc
        WHERE pfamA.pfamA_acc IN ({placeholders})
    """

    cursor.execute(query, pfam_accs)
    results = cursor.fetchall()

    # Build result dictionary
    pfam_info = {}
    for row in results:
        pfam_info[row['pfamA_acc']] = {
            'pfam_id': row['pfamA_id'],
            'clan_acc': row['clan_acc'] if row['clan_acc'] else 'No_clan'
        }

    cursor.close()
    if close_connection:
        connection.close()

    # Fill in missing entries
    for acc in pfam_accs:
        if acc not in pfam_info:
            pfam_info[acc] = {
                'pfam_id': 'Unknown',
                'clan_acc': 'Unknown'
            }

    return pfam_info


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


def calculate_mean_plddt(cif_file, start, end, verbose=False):
    """
    Calculate mean pLDDT score for a specific region in the structure.

    Args:
        cif_file: path to CIF file
        start: start residue position (1-indexed)
        end: end residue position (1-indexed)
        verbose: print debug information

    Returns:
        tuple: (mean pLDDT score, number of residues counted)
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('protein', cif_file)

    chain = next(structure.get_chains())

    plddt_scores = []
    residue_count = 0
    for residue in chain.get_residues():
        res_id = residue.get_id()[1]
        if start <= res_id <= end:
            # pLDDT is stored in the B-factor field
            for atom in residue.get_atoms():
                plddt_scores.append(atom.get_bfactor())
                residue_count += 1
                break  # Only need one atom per residue

    if verbose and plddt_scores:
        print(f"  pLDDT calculated over {residue_count} residues (region {start}-{end})")

    if plddt_scores:
        return sum(plddt_scores) / len(plddt_scores), residue_count
    return 0, 0


def verify_sequence_match(alignment_seq, structure_seq, start, end, verbose=False):
    """
    Verify that the alignment sequence (with gaps removed) matches the structure sequence.

    Args:
        alignment_seq: sequence from SEED alignment (gaps already removed)
        structure_seq: full sequence from AlphaFold structure
        start: start residue position (1-indexed)
        end: end residue position (1-indexed)
        verbose: print debug information

    Returns:
        bool: True if sequences match exactly
    """
    # Extract the relevant region from structure sequence (convert to 0-indexed)
    structure_region = structure_seq[start-1:end]

    # Compare sequences
    match = alignment_seq == structure_region

    if verbose and not match:
        print(f"  Sequence mismatch detected:")
        print(f"    Expected length: {len(alignment_seq)}")
        print(f"    Structure length: {len(structure_region)}")
        if len(alignment_seq) > 0 and len(structure_region) > 0:
            print(f"    Alignment (first 50): {alignment_seq[:50]}")
            print(f"    Structure (first 50): {structure_region[:50]}")

    return match


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


def create_foldseek_database(source_dir, db_path):
    """
    Create a foldseek database from a directory of CIF files.
    Only creates if database doesn't already exist.

    Args:
        source_dir: directory containing CIF files
        db_path: path where foldseek database should be created

    Returns:
        str: path to the database
    """
    # Check if database already exists (foldseek creates multiple files with this prefix)
    if os.path.exists(db_path) or os.path.exists(f"{db_path}.dbtype"):
        print(f"Foldseek database already exists at {db_path}")
        return db_path

    print(f"Creating foldseek database at {db_path}...")
    print(f"  Source: {source_dir}")
    print("  This is a one-time setup and may take several minutes...")

    # Create directory if it doesn't exist
    db_dir = os.path.dirname(db_path)
    os.makedirs(db_dir, exist_ok=True)

    # Run foldseek createdb
    cmd = ['foldseek', 'createdb', source_dir, db_path]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"Foldseek createdb failed: {result.stderr}")
        print("  Database created successfully!")
    except Exception as e:
        raise RuntimeError(f"Failed to create foldseek database: {e}")

    return db_path


def run_foldseek_search(query_cif, database_path, output_file, tmp_dir):
    """
    Run foldseek easy-search against the Pfam structure database.

    Args:
        query_cif: path to query structure
        database_path: path to indexed foldseek database
        output_file: where to write results
        tmp_dir: temporary directory for foldseek

    Returns:
        str: path to foldseek output file
    """
    foldseek_output = os.path.join(tmp_dir, 'foldseek_raw.out')
    log_file = os.path.join(tmp_dir, 'foldseek.log')

    cmd = [
        'foldseek', 'easy-search',
        query_cif,
        database_path,
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
    Queries Pfam database to get Pfam ID and clan information.

    Outputs a TSV file with columns:
    - Curation directory
    - Query info
    - Target Pfam accession
    - Pfam ID (name)
    - Clan accession
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

    # Query database for Pfam ID and clan information
    print("Querying Pfam database for family names and clan information...")
    unique_pfam_accs = list(set([r['pfam_acc'] for r in results]))

    try:
        pfam_info = get_pfam_info(unique_pfam_accs)
    except Exception as e:
        print(f"Warning: Could not query Pfam database: {e}")
        print("Continuing without Pfam ID and clan information...")
        pfam_info = {acc: {'pfam_id': 'Unknown', 'clan_acc': 'Unknown'}
                     for acc in unique_pfam_accs}

    # Add Pfam ID and clan info to results
    for r in results:
        pfam_acc = r['pfam_acc']
        if pfam_acc in pfam_info:
            r['pfam_id'] = pfam_info[pfam_acc]['pfam_id']
            r['clan_acc'] = pfam_info[pfam_acc]['clan_acc']
        else:
            r['pfam_id'] = 'Unknown'
            r['clan_acc'] = 'Unknown'

    # Write output
    with open(output_file, 'w') as out:
        # Header
        out.write('\t'.join([
            'Curation_Dir', 'Query', 'Pfam_Accession', 'Pfam_ID', 'Clan_Accession',
            'Target_Full', 'Identity', 'Aln_Length', 'Query_Length', 'Target_Length',
            'Overlap', 'E-value', 'Bitscore'
        ]) + '\n')

        # Data
        for r in results:
            out.write('\t'.join([
                r['curation_dir'],
                r['query'],
                r['pfam_acc'],
                r['pfam_id'],
                r['clan_acc'],
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
        '--source-dir',
        default='/nfs/production/agb/interpro/users/typhaine/interpro-pfam-curation-tools/pfam_add_clan_search/results/chopped_cif',
        help='Path to directory containing Pfam AlphaFold CIF files (default: chopped_cif location)'
    )
    parser.add_argument(
        '--database',
        default='/nfs/production/agb/pfam/curation/foldseek/pfam_db',
        help='Path to foldseek database (default: /nfs/production/agb/pfam/curation/foldseek/pfam_db)'
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

    # Limit to first 20 sequences for performance
    max_sequences = 20
    if len(sequences) > max_sequences:
        print(f"Limiting to first {max_sequences} sequences for performance (skipping {len(sequences) - max_sequences})")
        sequences = sequences[:max_sequences]

    # Create or verify foldseek database exists
    print("\nChecking foldseek database...")
    foldseek_db = create_foldseek_database(args.source_dir, args.database)

    # Check if chopped model already exists in curation directory
    saved_chopped_model = os.path.join(args.curation_dir, 'query_model.cif')

    # Create temporary directory for foldseek
    with tempfile.TemporaryDirectory() as tmp_dir:
        if os.path.exists(saved_chopped_model):
            print(f"\nFound existing chopped model: {saved_chopped_model}")
            print("Using cached model (delete this file to force reprocessing)")
            chopped_model = saved_chopped_model
        else:
            # Process SEED sequences to find best model
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

                # Verify sequence match (alignment_seq already has gaps removed)
                structure_seq = get_sequence_from_cif(cif_file)
                if not verify_sequence_match(alignment_seq, structure_seq, start, end, verbose=True):
                    print(f"Warning: Sequence mismatch for {uniprot_acc}. Skipping.")
                    continue

                # Calculate mean pLDDT (only for the alignment region)
                mean_plddt, num_residues = calculate_mean_plddt(cif_file, start, end, verbose=True)
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
            chopped_model_tmp = os.path.join(tmp_dir, 'query_chopped.cif')
            print(f"Chopping structure to region {best_info[1]}-{best_info[2]}...")
            chop_structure(best_model, chopped_model_tmp, best_info[1], best_info[2])

            # Save chopped model to curation directory for reuse
            print(f"Saving chopped model to {saved_chopped_model}")
            shutil.copy2(chopped_model_tmp, saved_chopped_model)
            chopped_model = saved_chopped_model

        # Run foldseek search
        print(f"\nSearching against database: {foldseek_db}")
        foldseek_output = run_foldseek_search(
            chopped_model,
            foldseek_db,
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
