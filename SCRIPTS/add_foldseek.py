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

# Import local TED module
from ted_local import TEDLocal


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
            match = re.match(r'^(\w+\.?\d*)/(\d+)-(\d+)\s+(.+)$', line)
            if match:
                uniprot_acc = match.group(1)
                start = int(match.group(2))
                end = int(match.group(3))
                sequence = match.group(4).replace('.', '').replace('-', '')
                sequences.append((uniprot_acc, start, end, sequence))

    return sequences


def parse_desc_ted_reference(desc_file):
    """
    Parse DESC file to find TED reference lines.

    Looks for lines like: SE   TED:Q9RYM9_TED03

    Returns:
        dict with 'uniprot_acc' and 'ted_id', or None if not found
    """
    if not os.path.exists(desc_file):
        return None

    with open(desc_file, 'r') as f:
        for line in f:
            # Match lines like: SE   TED:Q9RYM9_TED03
            match = re.match(r'^SE\s+TED:(\w+\.?\d*)_(\w+)', line)
            if match:
                return {
                    'uniprot_acc': match.group(1),
                    'ted_id': f"{match.group(1)}_{match.group(2)}"
                }

    return None


def get_ted_domain_coordinates(uniprot_acc, ted_id, ted_db_path=None):
    """
    Fetch coordinates for a specific TED domain from local database.

    Args:
        uniprot_acc: UniProt accession (e.g., "A0A0H3PIP6")
        ted_id: Internal TED identifier like "A0A0H3PIP6_TED01"
                (only the suffix like "TED01" is used for matching)
        ted_db_path: Optional path to TED SQLite database

    Returns:
        tuple: (start, end) coordinates, or None if not found
    """
    try:
        ted = TEDLocal(ted_db_path)
        domains = ted.get_domains(uniprot_acc)

        if not domains:
            return None

        # Extract just the TED suffix (e.g., "TED01" from "A0A0H3PIP6_TED01")
        ted_suffix = ted_id.split('_')[-1] if '_' in ted_id else ted_id

        # Loop through domains to find matching suffix
        for domain in domains:
            api_ted_id = domain.get('ted_id', '')

            # Match domains where the ted_id ends with our TED suffix
            if api_ted_id.endswith(ted_suffix):
                segments = domain.get('segments', [])
                if segments:
                    # Return first start and last end
                    return (min(s['start'] for s in segments), max(s['end'] for s in segments))

        return None

    except Exception as e:
        print(f"Warning: Failed to fetch TED domain coordinates: {e}")
        return None


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

    Args:
        uniprot_acc: UniProt accession (may include version like "P12345.2")
        output_dir: directory to save the file

    Returns:
        str: path to downloaded CIF file, or None if failed
    """
    import urllib.request

    # Strip version number from accession for AlphaFold API (P12345.2 -> P12345)
    base_acc = uniprot_acc.split('.')[0]

    url = f"https://alphafold.ebi.ac.uk/files/AF-{base_acc}-F1-model_v6.cif"
    output_file = os.path.join(output_dir, f"AF-{base_acc}-F1-model_v6.cif")

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
    from Bio.PDB import MMCIFParser

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


def chop_structure(input_cif, output_cif, start, end):
    """
    Chop a structure to keep only residues in the specified region.

    Args:
        input_cif: input CIF file path
        output_cif: output CIF file path
        start: start residue (1-indexed)
        end: end residue (1-indexed)
    """
    from Bio.PDB import MMCIFParser, MMCIFIO, Select

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
    import subprocess

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
    import subprocess

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


def fetch_ted_domains(uniprot_acc, protein_length_cache=None, ted_db_path=None):
    """
    Fetch TED domain information from local database.

    Args:
        uniprot_acc: UniProt accession (may include version)
        protein_length_cache: dict of {uniprot_acc: length} to avoid re-downloading
        ted_db_path: Optional path to TED SQLite database

    Returns:
        dict with 'protein_length' and 'domains' list, or None if failed
    """
    import tempfile

    # Strip version number
    base_acc = uniprot_acc.split('.')[0]

    try:
        ted = TEDLocal(ted_db_path)
        ted_domains = ted.get_domains(base_acc)

        if not ted_domains:
            return None

        # Local database doesn't have protein length, check cache or download
        protein_length = 0
        if protein_length_cache and uniprot_acc in protein_length_cache:
            protein_length = protein_length_cache[uniprot_acc]
        else:
            # Download AlphaFold model to get sequence length
            with tempfile.TemporaryDirectory() as tmp_dir:
                cif_file = download_alphafold_model(base_acc, tmp_dir)
                if cif_file:
                    sequence = get_sequence_from_cif(cif_file)
                    protein_length = len(sequence)

        domains = []

        for domain in ted_domains:
            ted_id = domain.get('ted_id', '')
            segments = domain.get('segments', [])

            if not segments:
                continue

            domains.append({
                'ted_id': ted_id,
                'segments': segments
            })

        return {
            'protein_length': protein_length if protein_length else 0,
            'domains': domains
        }

    except Exception as e:
        print(f"Warning: Failed to fetch TED domains for {uniprot_acc}: {e}")
        return None


def calculate_ted_consistency_score(seed_start, seed_end, ted_domains):
    """
    Calculate TED-Pfam consistency score for a SEED region.

    Uses bidirectional coverage with minimum:
    - Coverage_TED = (SEED residues covered by TED) / (SEED length)
    - Coverage_SEED = (TED residues within SEED) / (TED total length)
    - Score = min(Coverage_TED, Coverage_SEED)

    Args:
        seed_start: SEED region start position
        seed_end: SEED region end position
        ted_domains: list of TED domains with segments

    Returns:
        float: consistency score from 0 to 1, or 0 if no TED domains
    """
    if not ted_domains:
        return 0.0

    seed_length = seed_end - seed_start + 1
    best_score = 0.0

    for domain in ted_domains:
        segments = domain.get('segments', [])
        if not segments:
            continue

        # Treat multi-segment domain as continuous from first start to last end
        ted_start = min(seg['start'] for seg in segments)
        ted_end = max(seg['end'] for seg in segments)
        ted_length = ted_end - ted_start + 1

        # Calculate overlap
        overlap_start = max(seed_start, ted_start)
        overlap_end = min(seed_end, ted_end)

        if overlap_start <= overlap_end:
            intersection = overlap_end - overlap_start + 1

            # Bidirectional coverage
            coverage_ted = intersection / seed_length  # SEED covered by TED
            coverage_seed = intersection / ted_length   # TED within SEED

            # Take minimum for this domain
            score = min(coverage_ted, coverage_seed)
            best_score = max(best_score, score)

    return best_score


def run_webshot(uniprot_acc, curation_dir, output_file='ted_web.png'):
    """
    Capture a screenshot of the TED website protein page.

    Args:
        uniprot_acc: UniProt accession (may include version like "P12345.2")
        curation_dir: path to curation directory
        output_file: output filename (default: ted_web.png)
    """
    import subprocess

    # Strip version number from accession
    base_acc = uniprot_acc.split('.')[0]

    url = f"https://ted.cathdb.info/uniprot/{base_acc}"
    output_path = os.path.join(curation_dir, output_file)

    print(f"\nCapturing TED website screenshot for {base_acc}...")
    cmd = ['webshot', url, output_path]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Warning: webshot failed: {result.stderr}")
        else:
            print(f"Created TED website screenshot: {output_path}")
    except FileNotFoundError:
        print("Warning: webshot command not found, skipping TED website screenshot")
    except Exception as e:
        print(f"Warning: Failed to capture TED website screenshot: {e}")


def create_ted_visualization(plddt_scores, curation_dir, ted_cache=None, mean_consistency=None, output_file='ted.png'):
    """
    Create a visualization of TED domains for proteins in the SEED alignment.

    Args:
        plddt_scores: list of dicts with 'accession', 'mean_plddt', and 'consistency_score'
        curation_dir: path to curation directory
        ted_cache: dict of cached TED domain data (optional)
        mean_consistency: mean family consistency score (optional)
        output_file: output filename (default: ted.png)
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
    except ImportError:
        print("Warning: matplotlib not available, skipping TED visualization")
        return

    # Color scheme for TED domains (official TED database colors)
    domain_colors = ['#4A79A7', '#F28E2C', '#E15759', '#76B7B2', '#59A14F',
                     '#EDC949', '#AF7AA1', '#FF9DA7', '#9C755F', '#BAB0AB']

    # Build protein data for visualization
    proteins_data = []

    for score_entry in plddt_scores:
        accession = score_entry['accession']  # Format: P12345.1/10-145

        # Parse accession
        if '/' not in accession:
            continue

        parts = accession.split('/')
        uniprot_acc = parts[0]
        region_parts = parts[1].split('-')
        seed_start = int(region_parts[0])
        seed_end = int(region_parts[1])

        # Get TED data from cache or fetch if not cached
        if ted_cache and uniprot_acc in ted_cache:
            ted_data = ted_cache[uniprot_acc]
        else:
            ted_data = fetch_ted_domains(uniprot_acc)

        if ted_data:
            # Get consistency score from plddt_scores entry (already calculated)
            consistency_score = score_entry.get('consistency_score', 0.0)

            proteins_data.append({
                'accession': accession,
                'uniprot_acc': uniprot_acc,
                'seed_start': seed_start,
                'seed_end': seed_end,
                'protein_length': ted_data['protein_length'],
                'domains': ted_data['domains'],
                'mean_plddt': score_entry['mean_plddt'],
                'consistency_score': consistency_score
            })

    if not proteins_data:
        print("Warning: No TED domain data available for visualization")
        return

    # Calculate figure dimensions
    residues_per_line = 1500
    line_height = 60  # pixels per protein line
    margin_left = 200  # space for labels (increased for consistency score)
    margin_right = 50
    margin_top = 40
    margin_bottom = 40
    line_spacing = 10

    # Calculate total lines needed
    total_lines = 0
    for protein in proteins_data:
        lines_needed = (protein['protein_length'] + residues_per_line - 1) // residues_per_line
        protein['lines_needed'] = lines_needed
        total_lines += lines_needed

    # Create figure
    fig_width = margin_left + residues_per_line + margin_right
    fig_height = margin_top + (total_lines * (line_height + line_spacing)) + margin_bottom

    fig, ax = plt.subplots(figsize=(fig_width / 100, fig_height / 100), dpi=100)
    ax.set_xlim(0, fig_width)
    ax.set_ylim(0, fig_height)
    ax.axis('off')

    # Add mean consistency score as title if available
    if mean_consistency is not None:
        ax.text(fig_width / 2, fig_height - 15,
               f"Mean Family TED Consistency Score: {mean_consistency:.2f}",
               ha='center', va='top', fontsize=11, weight='bold')

    # Draw proteins
    current_y = fig_height - margin_top - line_height / 2

    for protein in proteins_data:
        protein_length = protein['protein_length']
        lines_needed = protein['lines_needed']

        for line_idx in range(lines_needed):
            line_start = line_idx * residues_per_line
            line_end = min(line_start + residues_per_line, protein_length)

            # Draw label on first line only
            if line_idx == 0:
                label = f"{protein['accession']} ({protein['mean_plddt']:.2f}) [{protein['consistency_score']:.2f}]"
                ax.text(margin_left - 10, current_y, label,
                       ha='right', va='center', fontsize=9, family='monospace')

            # Convert to pixel coordinates
            x_start = margin_left + (0 if line_idx == 0 else 0)
            x_end = margin_left + (line_end - line_start)

            # Draw protein backbone (grey line)
            backbone_height = 8
            backbone = patches.Rectangle(
                (x_start, current_y - backbone_height / 2),
                x_end - x_start, backbone_height,
                linewidth=0, facecolor='#CCCCCC'
            )
            ax.add_patch(backbone)

            # Draw TED domains
            for domain_idx, domain in enumerate(protein['domains']):
                color = domain_colors[domain_idx % len(domain_colors)]
                ted_id = domain.get('ted_id', '')

                # Extract TED domain number (e.g., "TED01" from "P12345_TED01")
                ted_label = ted_id.split('_')[-1] if '_' in ted_id else ted_id

                for segment in domain['segments']:
                    seg_start = segment['start']
                    seg_end = segment['end']

                    # Check if segment overlaps with current line
                    if seg_end >= line_start and seg_start <= line_end:
                        # Calculate visible portion
                        visible_start = max(seg_start, line_start)
                        visible_end = min(seg_end, line_end)

                        # Convert to pixel coordinates relative to line
                        domain_x = margin_left + (visible_start - line_start)
                        domain_width = visible_end - visible_start

                        domain_height = 20
                        domain_rect = patches.Rectangle(
                            (domain_x, current_y - domain_height / 2),
                            domain_width, domain_height,
                            linewidth=1, edgecolor='black', facecolor=color, alpha=0.8
                        )
                        ax.add_patch(domain_rect)

                        # Add TED domain label inside the box if there's enough space
                        if domain_width > 30:  # Only add text if box is wide enough
                            ax.text(domain_x + domain_width / 2, current_y,
                                   ted_label, ha='center', va='center',
                                   fontsize=7, weight='bold', color='white')

            # Draw SEED region
            seed_start = protein['seed_start']
            seed_end = protein['seed_end']

            if seed_end >= line_start and seed_start <= line_end:
                visible_start = max(seed_start, line_start)
                visible_end = min(seed_end, line_end)

                seed_x = margin_left + (visible_start - line_start)
                seed_width = visible_end - visible_start

                seed_height = 30
                seed_rect = patches.Rectangle(
                    (seed_x, current_y - seed_height / 2),
                    seed_width, seed_height,
                    linewidth=2, edgecolor='black', facecolor='none', alpha=0.6
                )
                ax.add_patch(seed_rect)

            current_y -= (line_height + line_spacing)

    # Save figure
    output_path = f"{curation_dir}/{output_file}"
    plt.tight_layout()
    plt.savefig(output_path, dpi=100, bbox_inches='tight')
    plt.close()

    print(f"\nCreated TED domain visualization: {output_path}")


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
    parser.add_argument(
        '--force',
        action='store_true',
        help='Force rerun even if output files already exist'
    )

    args = parser.parse_args()

    # Check SEED file exists
    seed_file = os.path.join(args.curation_dir, 'SEED')
    if not os.path.exists(seed_file):
        print(f"Error: SEED file not found at {seed_file}")
        sys.exit(1)

    print(f"Processing curation directory: {args.curation_dir}")

    # Import remaining modules only when actually needed
    import tempfile
    import shutil

    # Check for TED reference in DESC file
    desc_file = os.path.join(args.curation_dir, 'DESC')
    ted_reference = parse_desc_ted_reference(desc_file)

    # Parse SEED alignment
    print("Parsing SEED alignment...")
    sequences = parse_seed_alignment(seed_file)
    print(f"Found {len(sequences)} sequences in SEED alignment")

    if ted_reference:
        print(f"\nFound TED reference in DESC: {ted_reference['ted_id']}")
        coords = get_ted_domain_coordinates(ted_reference['uniprot_acc'], ted_reference['ted_id'])
        if coords:
            print(f"  TED domain coordinates: {coords[0]}-{coords[1]}")
            # Add TED reference to beginning of sequences list (with empty sequence as placeholder)
            sequences.insert(0, (ted_reference['uniprot_acc'], coords[0], coords[1], ''))
        else:
            print(f"  Warning: Could not fetch coordinates for {ted_reference['ted_id']}")

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
        if os.path.exists(saved_chopped_model) and not args.force:
            print(f"\nFound existing chopped model: {saved_chopped_model}")
            print("Using cached model (use --force to reprocess)")
            chopped_model = saved_chopped_model
        else:
            if args.force and os.path.exists(saved_chopped_model):
                print(f"\nForce flag set, reprocessing models (ignoring cached model)")

            # Process SEED sequences to find best model
            best_model = None
            best_plddt = 0
            best_info = None

            # Collect all pLDDT scores for output file
            plddt_scores = []

            # Cache protein lengths from downloaded models
            protein_length_cache = {}

            # Process each sequence
            for uniprot_acc, start, end, alignment_seq in sequences:
                print(f"\nProcessing {uniprot_acc} ({start}-{end})...")

                # Download AlphaFold model
                cif_file = download_alphafold_model(uniprot_acc, tmp_dir)
                if not cif_file:
                    continue

                # Get structure sequence
                structure_seq = get_sequence_from_cif(cif_file)

                # Verify sequence match only if alignment_seq is provided (not empty)
                # Empty alignment_seq indicates a TED reference from DESC file
                if alignment_seq:
                    if not verify_sequence_match(alignment_seq, structure_seq, start, end, verbose=True):
                        print(f"Warning: Sequence mismatch for {uniprot_acc}. Skipping.")
                        continue
                else:
                    print(f"  Processing TED reference from DESC file")

                # Cache protein length for later use
                protein_length_cache[uniprot_acc] = len(structure_seq)

                # Calculate mean pLDDT (only for the alignment region)
                mean_plddt, num_residues = calculate_mean_plddt(cif_file, start, end, verbose=True)
                print(f"  Mean pLDDT: {mean_plddt:.2f}")

                # Store pLDDT score for output file
                plddt_scores.append({
                    'accession': f"{uniprot_acc}/{start}-{end}",
                    'mean_plddt': mean_plddt
                })

                if mean_plddt > best_plddt:
                    best_plddt = mean_plddt
                    best_model = cif_file
                    best_info = (uniprot_acc, start, end)

            if not best_model:
                print("\nError: No valid models found with matching sequences")
                sys.exit(1)

            # Write pLDDT scores to file, sorted by highest score first
            if plddt_scores:
                # Fetch TED domains once and cache them
                print("\nFetching TED domain data...")
                ted_cache = {}
                for entry in plddt_scores:
                    accession = entry['accession']
                    if '/' in accession:
                        parts = accession.split('/')
                        uniprot_acc = parts[0]

                        # Fetch TED domains only once per protein (with protein length cache)
                        if uniprot_acc not in ted_cache:
                            ted_data = fetch_ted_domains(uniprot_acc, protein_length_cache)
                            ted_cache[uniprot_acc] = ted_data
                            if ted_data and ted_data.get('domains'):
                                print(f"  {uniprot_acc}: found {len(ted_data['domains'])} TED domain(s)")
                            else:
                                print(f"  {uniprot_acc}: no TED domains found")

                # Calculate consistency scores for each entry
                print("Calculating TED consistency scores...")
                consistency_scores = []
                for entry in plddt_scores:
                    accession = entry['accession']
                    if '/' in accession:
                        parts = accession.split('/')
                        uniprot_acc = parts[0]
                        region_parts = parts[1].split('-')
                        seed_start = int(region_parts[0])
                        seed_end = int(region_parts[1])

                        # Use cached TED data
                        ted_data = ted_cache.get(uniprot_acc)
                        if ted_data:
                            consistency_score = calculate_ted_consistency_score(
                                seed_start, seed_end, ted_data['domains']
                            )
                            entry['consistency_score'] = consistency_score
                            consistency_scores.append(consistency_score)
                            print(f"  {accession}: consistency score = {consistency_score:.2f}")
                        else:
                            entry['consistency_score'] = 0.0
                            print(f"  {accession}: no TED data, score = 0.00")
                    else:
                        entry['consistency_score'] = 0.0

                # Calculate mean family consistency score
                if consistency_scores:
                    mean_consistency = sum(consistency_scores) / len(consistency_scores)
                    print(f"\nMean family TED consistency score: {mean_consistency:.2f}")
                else:
                    mean_consistency = 0.0
                    print(f"\nMean family TED consistency score: 0.00 (no TED data available)")

                plddt_file = os.path.join(args.curation_dir, 'pLDDT')
                plddt_scores.sort(key=lambda x: x['mean_plddt'], reverse=True)

                with open(plddt_file, 'w') as f:
                    # Write header with mean consistency score
                    f.write(f"# Mean family TED consistency score: {mean_consistency:.2f}\n")
                    for entry in plddt_scores:
                        f.write(f"{entry['mean_plddt']:.2f}\t{entry['accession']}\t{entry['consistency_score']:.2f}\n")

                print(f"Wrote pLDDT scores for {len(plddt_scores)} sequences to {plddt_file}")

                # Create TED domain visualization (reusing cached TED data)
                print("\nCreating TED domain visualization...")
                create_ted_visualization(plddt_scores, args.curation_dir, ted_cache, mean_consistency)

                # Capture TED website screenshot for best model
                run_webshot(best_info[0], args.curation_dir)

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
    # EARLY EXIT: Quick check before importing heavy modules
    import argparse
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('curation_dir')
    parser.add_argument('--output', default='foldseek')
    parser.add_argument('--force', action='store_true')
    args, _ = parser.parse_known_args()

    output_file = os.path.join(args.curation_dir, args.output)
    plddt_file = os.path.join(args.curation_dir, 'pLDDT')
    query_model_file = os.path.join(args.curation_dir, 'query_model.cif')
    ted_png_file = os.path.join(args.curation_dir, 'ted.png')

    # Check if all output files exist
    all_outputs_exist = (os.path.exists(output_file) and
                        os.path.exists(plddt_file) and
                        os.path.exists(query_model_file) and
                        os.path.exists(ted_png_file))

    if all_outputs_exist and not args.force:
        print(f"Output files already exist in {args.curation_dir}. Skipping to avoid duplication.")
        print(f"Use --force to rerun anyway.")
        sys.exit(0)

    # Output files missing or force flag is set, proceed with full processing
    main()
