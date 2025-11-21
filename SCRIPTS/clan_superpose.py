#!/usr/bin/env python3
"""
Create structural superposition of AlphaFold models for a Pfam clan.

This script:
1. Gets list of Pfam families in a clan from the database
2. For each family:
   - Reads the SEED alignment
   - Downloads AlphaFold models for sequences (up to 20)
   - Extracts domain regions and calculates mean pLDDT
   - Selects the best model (highest pLDDT)
   - Saves as PFXXXXX.pdb
3. Orders domains by mean pLDDT
4. Uses the highest scoring domain as reference
5. Superposes all other structures using TM-align
6. Creates a combined PDB file with all superposed structures
"""

import argparse
import os
import sys
import re
import tempfile
import shutil
import subprocess
from pathlib import Path


def connect_to_pfam_db(config_file='~/.my.cnf'):
    """
    Connect to the Pfam MySQL database using credentials from config file.

    Args:
        config_file: path to MySQL config file (default: ~/.my.cnf)

    Returns:
        mysql.connector connection object
    """
    import mysql.connector
    from configparser import ConfigParser

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


def get_clan_families(clan_acc, connection):
    """
    Get all Pfam families in a clan.

    Args:
        clan_acc: Clan accession (e.g., 'CL0004')
        connection: MySQL connection object

    Returns:
        list of tuples: (pfam_acc, pfam_id, description)
    """
    query = """
        SELECT pfamA.pfamA_acc, pfamA_id, pfamA.description
        FROM clan
        JOIN clan_membership ON clan.clan_acc = clan_membership.clan_acc
        JOIN pfamA ON clan_membership.pfamA_acc = pfamA.pfamA_acc
        WHERE clan.clan_acc = %s
        ORDER BY pfamA.pfamA_acc
    """

    cursor = connection.cursor()
    cursor.execute(query, (clan_acc,))
    families = cursor.fetchall()
    cursor.close()

    print(f"\nFound {len(families)} families in clan {clan_acc}:")
    for pfam_acc, pfam_id, desc in families:
        print(f"  {pfam_acc} ({pfam_id}): {desc}")

    return families


def fetch_seed_file(pfam_acc, work_dir):
    """
    Fetch SEED file for a Pfam family directly from SVN repository.
    This is much faster than pfco as it only downloads the SEED file.

    Args:
        pfam_acc: Pfam accession (e.g., 'PF00001')
        work_dir: Working directory where SEED file will be saved

    Returns:
        Path to SEED file, or None if fetch failed
    """
    family_dir = Path(work_dir) / pfam_acc
    seed_path = family_dir / 'SEED'

    # Check if already exists
    if seed_path.exists():
        print(f"  Family {pfam_acc} SEED file already exists")
        return str(seed_path)

    # Create family directory if needed
    family_dir.mkdir(parents=True, exist_ok=True)

    # Fetch SEED file directly from SVN
    svn_url = f"https://xfam-svn-hl.ebi.ac.uk/svn/pfam/trunk/Data/Families/{pfam_acc}/SEED"

    try:
        result = subprocess.run(
            ['svn', 'cat', svn_url],
            capture_output=True,
            text=True,
            timeout=30
        )

        if result.returncode != 0:
            print(f"  Warning: SVN fetch failed for {pfam_acc}")
            print(f"    stderr: {result.stderr}")
            return None

        # Write SEED file
        with open(seed_path, 'w') as f:
            f.write(result.stdout)

        print(f"  Fetched SEED file from SVN: {seed_path}")
        return str(seed_path)

    except subprocess.TimeoutExpired:
        print(f"  Warning: SVN fetch timed out for {pfam_acc}")
        return None
    except FileNotFoundError:
        print(f"  Error: svn command not found. Please ensure it's in your PATH")
        return None
    except Exception as e:
        print(f"  Warning: Error fetching SEED for {pfam_acc}: {e}")
        return None


def parse_seed_alignment(seed_file, max_sequences=20):
    """
    Parse a Pfam SEED alignment file and extract UniProt accessions with regions.

    Args:
        seed_file: Path to SEED file
        max_sequences: Maximum number of sequences to process

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
            match = re.match(r'^(\w+\.?\d*)/(\d+)-(\d+)\s+(.+)$', line)
            if match:
                uniprot_acc = match.group(1)
                start = int(match.group(2))
                end = int(match.group(3))
                sequence = match.group(4).replace('.', '').replace('-', '')
                sequences.append((uniprot_acc, start, end, sequence))

                if len(sequences) >= max_sequences:
                    break

    return sequences


def download_alphafold_model(uniprot_acc, output_dir):
    """
    Download AlphaFold v6 CIF model for a given UniProt accession.

    Args:
        uniprot_acc: UniProt accession (may include version like "P12345.2")
        output_dir: directory to save the file

    Returns:
        str: path to downloaded CIF file, or None if failed
    """
    import urllib.request

    # Strip version number from accession for AlphaFold API
    base_acc = uniprot_acc.split('.')[0]

    url = f"https://alphafold.ebi.ac.uk/files/AF-{base_acc}-F1-model_v6.cif"
    output_file = os.path.join(output_dir, f"AF-{base_acc}-F1-model_v6.cif")

    try:
        urllib.request.urlretrieve(url, output_file)
        return output_file
    except Exception as e:
        print(f"  Warning: Failed to download {uniprot_acc}: {e}")
        return None


def get_sequence_from_cif(cif_file):
    """
    Extract the protein sequence from a CIF file.

    Returns:
        str: amino acid sequence
    """
    from Bio.PDB import MMCIFParser

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
        float: mean pLDDT score
    """
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


def verify_sequence_match(alignment_seq, structure_seq, start, end):
    """
    Verify that the alignment sequence matches the structure sequence.

    Args:
        alignment_seq: sequence from SEED alignment (gaps already removed)
        structure_seq: full sequence from AlphaFold structure
        start: start residue position (1-indexed)
        end: end residue position (1-indexed)

    Returns:
        bool: True if sequences match exactly
    """
    # Extract the relevant region from structure sequence (convert to 0-indexed)
    structure_region = structure_seq[start-1:end]
    return alignment_seq == structure_region


def chop_structure_to_pdb(input_cif, output_pdb, start, end):
    """
    Chop a CIF structure to keep only residues in the specified region and save as PDB.

    Args:
        input_cif: input CIF file path
        output_pdb: output PDB file path
        start: start residue (1-indexed)
        end: end residue (1-indexed)
    """
    from Bio.PDB import MMCIFParser, PDBIO, Select

    class RegionSelect(Select):
        """Select only residues within a specific region."""

        def __init__(self, start, end):
            self.start = start
            self.end = end

        def accept_residue(self, residue):
            res_id = residue.get_id()[1]
            return self.start <= res_id <= self.end

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('protein', input_cif)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, RegionSelect(start, end))


def process_family(pfam_acc, pfam_id, work_dir, output_dir, tmp_dir):
    """
    Process a single Pfam family to get the best AlphaFold domain structure.

    Args:
        pfam_acc: Pfam accession
        pfam_id: Pfam ID (name)
        work_dir: Working directory where families will be checked out
        output_dir: Directory to save output PDB files
        tmp_dir: Temporary directory for downloads

    Returns:
        dict with 'pfam_acc', 'pfam_id', 'pdb_file', 'mean_plddt', or None if failed
    """
    print(f"\n{'='*60}")
    print(f"Processing {pfam_acc} ({pfam_id})")
    print(f"{'='*60}")

    # Check if PDB already exists
    output_pdb = os.path.join(output_dir, f"{pfam_acc}.pdb")
    plddt_file = os.path.join(output_dir, f"{pfam_acc}.plddt")

    if os.path.exists(output_pdb):
        print(f"PDB file already exists: {output_pdb}")
        print("Skipping processing (delete file to reprocess)")

        # Try to read pLDDT from saved file
        mean_plddt = 0.0
        if os.path.exists(plddt_file):
            try:
                with open(plddt_file, 'r') as f:
                    mean_plddt = float(f.read().strip())
                print(f"Read pLDDT score: {mean_plddt:.2f}")
            except:
                print("Warning: Could not read pLDDT file, using 0.0")

        return {
            'pfam_acc': pfam_acc,
            'pfam_id': pfam_id,
            'pdb_file': output_pdb,
            'mean_plddt': mean_plddt
        }

    # Fetch SEED file from SVN
    seed_file = fetch_seed_file(pfam_acc, work_dir)
    if not seed_file:
        return None

    # Parse SEED
    sequences = parse_seed_alignment(seed_file, max_sequences=20)
    print(f"Found {len(sequences)} sequences in SEED")

    if not sequences:
        print(f"No sequences found in SEED for {pfam_acc}")
        return None

    # Find best model
    best_model = None
    best_plddt = 0
    best_info = None

    for uniprot_acc, start, end, alignment_seq in sequences:
        print(f"\n  Processing {uniprot_acc}/{start}-{end}...")

        # Download AlphaFold model
        cif_file = download_alphafold_model(uniprot_acc, tmp_dir)
        if not cif_file:
            continue

        # Verify sequence match
        structure_seq = get_sequence_from_cif(cif_file)
        if not verify_sequence_match(alignment_seq, structure_seq, start, end):
            print(f"  Warning: Sequence mismatch, skipping")
            continue

        # Calculate mean pLDDT
        mean_plddt = calculate_mean_plddt(cif_file, start, end)
        print(f"  Mean pLDDT: {mean_plddt:.2f}")

        if mean_plddt > best_plddt:
            best_plddt = mean_plddt
            best_model = cif_file
            best_info = (uniprot_acc, start, end)

    if not best_model:
        print(f"No valid models found for {pfam_acc}")
        return None

    print(f"\nBest model: {best_info[0]}/{best_info[1]}-{best_info[2]} "
          f"with mean pLDDT {best_plddt:.2f}")

    # Save as PDB
    output_pdb = os.path.join(output_dir, f"{pfam_acc}.pdb")
    chop_structure_to_pdb(best_model, output_pdb, best_info[1], best_info[2])
    print(f"Saved domain structure to {output_pdb}")

    # Save pLDDT score to a file for future reference
    plddt_file = os.path.join(output_dir, f"{pfam_acc}.plddt")
    with open(plddt_file, 'w') as f:
        f.write(f"{best_plddt:.2f}\n")

    return {
        'pfam_acc': pfam_acc,
        'pfam_id': pfam_id,
        'pdb_file': output_pdb,
        'mean_plddt': best_plddt
    }


def extract_chain_from_pdb(input_pdb, output_pdb, chain_id='B'):
    """
    Extract a specific chain from a PDB file.

    Args:
        input_pdb: Input PDB file
        output_pdb: Output PDB file
        chain_id: Chain ID to extract (default: 'B')

    Returns:
        bool: True if successful, False otherwise
    """
    try:
        with open(input_pdb, 'r') as infile:
            with open(output_pdb, 'w') as outfile:
                for line in infile:
                    # Keep ATOM/HETATM lines for the specified chain
                    if line.startswith(('ATOM', 'HETATM')):
                        if len(line) > 21 and line[21] == chain_id:
                            outfile.write(line)
                    # Skip FATCAT REMARK lines which cause errors in Chimera
                    elif line.startswith('REMARK'):
                        # Filter out FATCAT-specific remarks
                        if any(keyword in line for keyword in [
                            'superimposed protein', 'twisted protein',
                            'result after optimizing', 'chain A', 'chain B',
                            'PF0', '.pdb'  # Also skip lines mentioning PDB files
                        ]):
                            continue
                        # Keep other REMARK lines
                        outfile.write(line)
                    # Keep other header lines
                    elif line.startswith(('HEADER', 'TITLE', 'COMPND', 'SOURCE')):
                        outfile.write(line)
                outfile.write("END\n")
        return True
    except Exception as e:
        print(f"  Error extracting chain {chain_id}: {e}")
        return False


def run_tmalign_superposition(reference_pdb, target_pdb, output_dir, target_name, tmalign_path='TMalign'):
    """
    Perform structural superposition using TM-align.

    Args:
        reference_pdb: Reference structure PDB file
        target_pdb: Target structure PDB file to superpose
        output_dir: Directory to save output files
        target_name: Name for output file (e.g., 'PF00651')
        tmalign_path: Path to TMalign executable

    Returns:
        Path to superposed PDB file, or None if failed
    """
    import subprocess

    output_prefix = os.path.join(output_dir, f"tmalign_{target_name}")

    # TM-align syntax: TMalign structure1.pdb structure2.pdb -o output_prefix
    # This superimposes structure1 onto structure2
    cmd = [
        tmalign_path,
        target_pdb,
        reference_pdb,
        '-o', output_prefix
    ]

    print(f"  Running: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

        # TM-align creates files: output_prefix, output_prefix_all, output_prefix_atm
        # The main file (output_prefix) contains the superposed target structure

        if os.path.exists(output_prefix):
            # Rename to our standard naming
            output_pdb = os.path.join(output_dir, f"{target_name}_superposed.pdb")
            import shutil
            shutil.move(output_prefix, output_pdb)

            # Extract TM-score from output
            for line in result.stdout.split('\n'):
                if 'TM-score' in line:
                    print(f"  {line.strip()}")
                    break

            print(f"  Saved superposed structure to: {output_pdb}")

            # Clean up extra files
            for suffix in ['_all', '_all_atm', '_all_atm_lig', '_atm', '_atm_lig']:
                extra_file = f"{output_prefix}{suffix}"
                if os.path.exists(extra_file):
                    os.remove(extra_file)

            return output_pdb
        else:
            print(f"  ERROR: TM-align failed - output file not found: {output_prefix}")
            print(f"  Return code: {result.returncode}")
            if result.stdout:
                print(f"  stdout: {result.stdout[:500]}")
            if result.stderr:
                print(f"  stderr: {result.stderr[:500]}")
            return None

    except subprocess.TimeoutExpired:
        print(f"  ERROR: TM-align timed out after 300 seconds")
        return None
    except FileNotFoundError:
        print(f"  ERROR: TM-align executable not found: {tmalign_path}")
        print(f"  Please ensure TM-align is installed and in your PATH")
        return None
    except Exception as e:
        print(f"  ERROR: TM-align exception: {e}")
        import traceback
        traceback.print_exc()
        return None


def run_rigid_superposition(reference_pdb, target_pdb, output_dir, target_name):
    """
    Perform rigid-body superposition using Bio.PDB (no flexible alignment).

    Args:
        reference_pdb: Reference structure PDB file
        target_pdb: Target structure PDB file to superpose
        output_dir: Directory to save output files
        target_name: Name for output file (e.g., 'PF00651')

    Returns:
        Path to superposed PDB file, or None if failed
    """
    from Bio.PDB import PDBParser, PDBIO, Superimposer

    try:
        parser = PDBParser(QUIET=True)

        # Load structures
        ref_structure = parser.get_structure('reference', reference_pdb)
        target_structure = parser.get_structure('target', target_pdb)

        # Get CA atoms from both structures
        ref_atoms = []
        target_atoms = []

        ref_model = ref_structure[0]
        target_model = target_structure[0]

        # Collect CA atoms
        for ref_chain in ref_model:
            for ref_res in ref_chain:
                if 'CA' in ref_res:
                    ref_atoms.append(ref_res['CA'])

        for target_chain in target_model:
            for target_res in target_chain:
                if 'CA' in target_res:
                    target_atoms.append(target_res['CA'])

        # Use the minimum number of atoms available
        num_atoms = min(len(ref_atoms), len(target_atoms))

        if num_atoms < 3:
            print(f"  ERROR: Not enough CA atoms for superposition (need >= 3, found {num_atoms})")
            return None

        print(f"  Superposing {num_atoms} CA atoms")

        # Perform superposition
        super_imposer = Superimposer()
        super_imposer.set_atoms(ref_atoms[:num_atoms], target_atoms[:num_atoms])

        # Apply transformation to entire target structure
        super_imposer.apply(target_structure.get_atoms())

        # Save superposed structure
        output_pdb = os.path.join(output_dir, f"{target_name}_superposed.pdb")
        io = PDBIO()
        io.set_structure(target_structure)
        io.save(output_pdb)

        print(f"  RMSD: {super_imposer.rms:.2f} Å")
        print(f"  Saved superposed structure to: {output_pdb}")

        return output_pdb

    except Exception as e:
        print(f"  ERROR: Superposition failed: {e}")
        import traceback
        traceback.print_exc()
        return None


def run_fatcat(reference_pdb, target_pdb, output_dir, target_name, fatcat_path='FATCAT'):
    """
    Run FATCAT to superpose target structure onto reference and extract chain B.

    Args:
        reference_pdb: Reference structure PDB file
        target_pdb: Target structure PDB file to superpose
        output_dir: Directory to save output files
        target_name: Name for output file (e.g., 'PF00651')
        fatcat_path: Path to FATCAT executable

    Returns:
        Path to superposed PDB file (chain B only), or None if failed
    """
    # Output prefix in the output directory
    output_prefix = os.path.join(output_dir, f"fatcat_{target_name}")

    cmd = [
        fatcat_path,
        '-p1', reference_pdb,
        '-p2', target_pdb,
        '-o', output_prefix,
        '-t'  # Output transformed PDB files
    ]

    print(f"  Running: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

        # Note: FATCAT may return exit code 1 even on success, so we check for output files instead
        # FATCAT creates .ini.twist.pdb (rigid-body) and .opt.twist.pdb (flexible)
        # Use .ini.twist.pdb for rigid-body superposition without flexibility
        fatcat_output = f"{output_prefix}.ini.twist.pdb"

        if not os.path.exists(fatcat_output):
            print(f"  ERROR: FATCAT failed - output file not found: {fatcat_output}")
            print(f"  Return code: {result.returncode}")
            if result.stdout:
                print(f"  stdout: {result.stdout[:1000]}")
            if result.stderr:
                print(f"  stderr: {result.stderr[:1000]}")
            return None

        print(f"  FATCAT output created: {fatcat_output}")

        # Extract only chain B (the superposed target structure)
        chain_b_output = os.path.join(output_dir, f"{target_name}_superposed.pdb")

        if extract_chain_from_pdb(fatcat_output, chain_b_output, chain_id='B'):
            print(f"  Extracted chain B to: {chain_b_output}")
            return chain_b_output
        else:
            print(f"  ERROR: Failed to extract chain B")
            return None

    except subprocess.TimeoutExpired:
        print(f"  ERROR: FATCAT timed out after 300 seconds")
        return None
    except FileNotFoundError:
        print(f"  ERROR: FATCAT executable not found: {fatcat_path}")
        print(f"  Please ensure FATCAT is installed and in your PATH")
        return None
    except Exception as e:
        print(f"  ERROR: FATCAT exception: {e}")
        return None


def combine_structures(pdb_files, output_file, labels=None):
    """
    Combine multiple PDB structures into a single file using MODEL blocks.
    Each structure is saved as a separate MODEL.

    Args:
        pdb_files: List of PDB file paths
        output_file: Output combined PDB file
        labels: Optional list of labels (e.g., Pfam accessions) for each structure
    """
    print(f"\nCombining {len(pdb_files)} structures:")
    for i, pdb_file in enumerate(pdb_files, 1):
        label = labels[i-1] if labels and i-1 < len(labels) else "Unknown"
        print(f"  {i}. {pdb_file} ({label})")

    with open(output_file, 'w') as out:
        model_num = 1

        for idx, pdb_file in enumerate(pdb_files):
            if not os.path.exists(pdb_file):
                print(f"  WARNING: File not found, skipping: {pdb_file}")
                continue

            label = labels[idx] if labels and idx < len(labels) else f"Model_{model_num}"
            out.write(f"MODEL     {model_num:4d}\n")

            with open(pdb_file, 'r') as f:
                for line in f:
                    # Skip existing MODEL/ENDMDL lines
                    if line.startswith('MODEL') or line.startswith('ENDMDL'):
                        continue
                    # Skip END line
                    if line.startswith('END') and not line.startswith('ENDMDL'):
                        continue
                    # Skip existing TITLE/HEADER/COMPND/SOURCE lines
                    if line.startswith(('TITLE', 'HEADER', 'COMPND', 'SOURCE')):
                        continue
                    # Skip FATCAT REMARK lines
                    if line.startswith('REMARK') and ('superimposed protein' in line or
                                                      'twisted protein' in line or
                                                      'result after optimizing' in line or
                                                      'chain A' in line or
                                                      'chain B' in line):
                        continue
                    out.write(line)

            out.write("ENDMDL\n")
            model_num += 1

        out.write("END\n")

    print(f"\nSuccessfully combined {model_num - 1} structures into {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Create structural superposition of AlphaFold models for a Pfam clan'
    )
    parser.add_argument(
        'clan_acc',
        help='Clan accession (e.g., CL0004)'
    )
    parser.add_argument(
        '--work-dir',
        default='.',
        help='Working directory where families will be checked out (default: current directory)'
    )
    parser.add_argument(
        '--output-dir',
        default='.',
        help='Output directory for PDB files (default: current directory)'
    )
    parser.add_argument(
        '--config',
        default='~/.my.cnf',
        help='Path to MySQL config file (default: ~/.my.cnf)'
    )

    args = parser.parse_args()

    # Create directories
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"Clan Superposition Pipeline")
    print(f"{'='*60}")
    print(f"Clan: {args.clan_acc}")
    print(f"Working directory: {work_dir}")
    print(f"Output directory: {output_dir}")

    # Connect to database
    print("\nConnecting to Pfam database...")
    connection = connect_to_pfam_db(args.config)

    # Get clan families
    families = get_clan_families(args.clan_acc, connection)

    if not families:
        print(f"No families found for clan {args.clan_acc}")
        sys.exit(1)

    # Process each family
    family_results = []

    with tempfile.TemporaryDirectory() as tmp_dir:
        for pfam_acc, pfam_id, description in families:
            result = process_family(pfam_acc, pfam_id, work_dir, output_dir, tmp_dir)
            if result:
                family_results.append(result)

        # Sort by pLDDT (highest first)
        family_results.sort(key=lambda x: x['mean_plddt'], reverse=True)

        print(f"\n{'='*60}")
        print(f"Summary: Processed {len(family_results)} families")
        print(f"{'='*60}")

        for result in family_results:
            print(f"  {result['pfam_acc']} ({result['pfam_id']}): "
                  f"pLDDT {result['mean_plddt']:.2f}")

        if len(family_results) == 0:
            print("No structures successfully processed")
            sys.exit(1)

        # Use highest pLDDT as reference
        reference = family_results[0]
        print(f"\nUsing {reference['pfam_acc']} as reference (pLDDT: {reference['mean_plddt']:.2f})")

        # Superpose all other structures
        print(f"\n{'='*60}")
        print(f"Superposing structures using TM-align")
        print(f"{'='*60}")

        superposed_files = []
        structure_labels = []

        # Add reference structure (chain A)
        superposed_files.append(reference['pdb_file'])
        structure_labels.append(f"{reference['pfam_acc']} (reference)")

        for result in family_results[1:]:
            print(f"\nSuperposing {result['pfam_acc']} onto reference...")

            superposed_pdb = run_tmalign_superposition(
                reference['pdb_file'],
                result['pdb_file'],
                str(output_dir),
                result['pfam_acc'],
                'TMalign'
            )

            if superposed_pdb:
                superposed_files.append(superposed_pdb)
                structure_labels.append(result['pfam_acc'])
                print(f"  Successfully superposed {result['pfam_acc']}")
            else:
                print(f"  WARNING: Skipping {result['pfam_acc']} due to TM-align failure")

        # Combine all structures into one file
        final_output = output_dir / f"{args.clan_acc}_superposed.pdb"
        combine_structures(superposed_files, str(final_output), labels=structure_labels)

        print(f"\n{'='*60}")
        print(f"SUCCESS!")
        print(f"{'='*60}")
        print(f"Combined superposed structures saved to: {final_output}")
        print(f"Individual domain structures saved in: {output_dir}")

    connection.close()


if __name__ == '__main__':
    main()
