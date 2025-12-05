#!/usr/bin/env python3
"""
Generate TED domain and Pfam domain visualizations for SEED alignment proteins.

This script:
1. Reads a SEED alignment from a Pfam curation directory
2. Downloads AlphaFold models to get protein lengths and pLDDT scores
3. Fetches TED domain annotations for each protein
4. Fetches Pfam domain annotations from the database
5. Creates visualization images showing domain architecture

Usage:
    add_image.py <curation_dir> [options]

Example:
    add_image.py /path/to/PF06267 --max-proteins 30
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

# Import local TED module
from ted_local import TEDLocal


def parse_seed_alignment(seed_file: str) -> List[Tuple[str, int, int, str]]:
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


def parse_desc_ted_reference(desc_file: str) -> Optional[Dict]:
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


def get_ted_domain_coordinates(uniprot_acc: str, ted_id: str, ted_db_path: str = None) -> Optional[Tuple[int, int]]:
    """
    Fetch coordinates for a specific TED domain from local database.

    Args:
        uniprot_acc: UniProt accession (e.g., "A0A0H3PIP6")
        ted_id: Internal TED identifier like "A0A0H3PIP6_TED01"
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

        for domain in domains:
            api_ted_id = domain.get('ted_id', '')
            if api_ted_id.endswith(ted_suffix):
                segments = domain.get('segments', [])
                if segments:
                    return (min(s['start'] for s in segments), max(s['end'] for s in segments))

        return None

    except Exception as e:
        print(f"Warning: Failed to fetch TED domain coordinates: {e}")
        return None


def download_alphafold_model(uniprot_acc: str, output_dir: str) -> Optional[str]:
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
        print(f"  Downloading AlphaFold model for {base_acc}...")
        urllib.request.urlretrieve(url, output_file)
        return output_file
    except Exception as e:
        print(f"  Failed to download {uniprot_acc}: {e}")
        return None


def get_sequence_from_cif(cif_file: str) -> str:
    """
    Extract the protein sequence from a CIF file.

    Returns:
        str: amino acid sequence
    """
    from Bio.PDB import MMCIFParser

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('protein', cif_file)

    chain = next(structure.get_chains())

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


def calculate_mean_plddt(cif_file: str, start: int, end: int, verbose: bool = False) -> Tuple[float, int]:
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
            for atom in residue.get_atoms():
                plddt_scores.append(atom.get_bfactor())
                residue_count += 1
                break  # Only need one atom per residue

    if verbose and plddt_scores:
        print(f"    pLDDT calculated over {residue_count} residues (region {start}-{end})")

    if plddt_scores:
        return sum(plddt_scores) / len(plddt_scores), residue_count
    return 0, 0


def fetch_ted_domains(uniprot_acc: str, protein_length_cache: Dict = None, ted_db_path: str = None) -> Optional[Dict]:
    """
    Fetch TED domain information from local database.

    Args:
        uniprot_acc: UniProt accession (may include version)
        protein_length_cache: dict of {uniprot_acc: length} to avoid re-downloading
        ted_db_path: Optional path to TED SQLite database

    Returns:
        dict with 'protein_length' and 'domains' list, or None if failed
    """
    # Strip version number
    base_acc = uniprot_acc.split('.')[0]

    try:
        ted = TEDLocal(ted_db_path)
        ted_domains = ted.get_domains(base_acc)

        if not ted_domains:
            return None

        # Get protein length from cache or use 0 as placeholder
        protein_length = 0
        if protein_length_cache and uniprot_acc in protein_length_cache:
            protein_length = protein_length_cache[uniprot_acc]

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
            'protein_length': protein_length,
            'domains': domains
        }

    except Exception as e:
        print(f"Warning: Failed to fetch TED domains for {uniprot_acc}: {e}")
        return None


def fetch_pfam_domains(uniprot_acc: str, config_file: str = '~/.my.cnf') -> List[Dict]:
    """
    Fetch Pfam domain annotations for a protein from the database.

    Args:
        uniprot_acc: UniProt accession (may include version)
        config_file: Path to MySQL config file

    Returns:
        List of domain dicts with 'pfam_acc', 'pfam_id', 'start', 'end'
    """
    # Strip version number
    base_acc = uniprot_acc.split('.')[0]
    config_file = os.path.expanduser(config_file)

    try:
        query = f"""
            SELECT r.seq_start, r.seq_end, r.pfamA_acc, a.pfamA_id
            FROM pfamA_reg_full_significant r
            JOIN pfamA a ON r.pfamA_acc = a.pfamA_acc
            WHERE r.pfamseq_acc = '{base_acc}'
            AND r.in_full = 1
            ORDER BY r.seq_start
        """

        cmd = [
            'mysql',
            f'--defaults-file={config_file}',
            'pfam_live',
            '--quick',
            '--silent',
            '--skip-column-names',
            '-e',
            query
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)

        if result.returncode != 0:
            return []

        domains = []
        for line in result.stdout.strip().split('\n'):
            if line:
                parts = line.split('\t')
                if len(parts) >= 4:
                    domains.append({
                        'start': int(parts[0]),
                        'end': int(parts[1]),
                        'pfam_acc': parts[2],
                        'pfam_id': parts[3]
                    })

        return domains

    except Exception as e:
        print(f"Warning: Failed to fetch Pfam domains for {uniprot_acc}: {e}")
        return []


def calculate_ted_consistency_score(seed_start: int, seed_end: int, ted_domains: List[Dict]) -> float:
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

        ted_start = min(seg['start'] for seg in segments)
        ted_end = max(seg['end'] for seg in segments)
        ted_length = ted_end - ted_start + 1

        overlap_start = max(seed_start, ted_start)
        overlap_end = min(seed_end, ted_end)

        if overlap_start <= overlap_end:
            intersection = overlap_end - overlap_start + 1
            coverage_ted = intersection / seed_length
            coverage_seed = intersection / ted_length
            score = min(coverage_ted, coverage_seed)
            best_score = max(best_score, score)

    return best_score


def create_ted_visualization(plddt_scores: List[Dict], curation_dir: str,
                             ted_cache: Dict = None, mean_consistency: float = None,
                             output_file: str = 'ted.png'):
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
    line_height = 60
    margin_left = 200
    margin_right = 50
    margin_top = 40
    margin_bottom = 40
    line_spacing = 10

    total_lines = 0
    for protein in proteins_data:
        lines_needed = (protein['protein_length'] + residues_per_line - 1) // residues_per_line
        protein['lines_needed'] = max(1, lines_needed)
        total_lines += protein['lines_needed']

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

            x_start = margin_left
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
                ted_label = ted_id.split('_')[-1] if '_' in ted_id else ted_id

                for segment in domain['segments']:
                    seg_start = segment['start']
                    seg_end = segment['end']

                    if seg_end >= line_start and seg_start <= line_end:
                        visible_start = max(seg_start, line_start)
                        visible_end = min(seg_end, line_end)

                        domain_x = margin_left + (visible_start - line_start)
                        domain_width = visible_end - visible_start

                        domain_height = 20
                        domain_rect = patches.Rectangle(
                            (domain_x, current_y - domain_height / 2),
                            domain_width, domain_height,
                            linewidth=1, edgecolor='black', facecolor=color, alpha=0.8
                        )
                        ax.add_patch(domain_rect)

                        if domain_width > 30:
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
    output_path = os.path.join(curation_dir, output_file)
    plt.tight_layout()
    plt.savefig(output_path, dpi=100, bbox_inches='tight')
    plt.close()

    print(f"Created TED domain visualization: {output_path}")


def create_pfam_visualization(proteins_data: List[Dict], curation_dir: str,
                              output_file: str = 'pfam.png'):
    """
    Create a visualization of Pfam domains for proteins in the SEED alignment.

    Args:
        proteins_data: list of dicts with protein info and Pfam domains
        curation_dir: path to curation directory
        output_file: output filename (default: pfam.png)
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
    except ImportError:
        print("Warning: matplotlib not available, skipping Pfam visualization")
        return

    if not proteins_data:
        print("Warning: No Pfam domain data available for visualization")
        return

    # Color scheme for Pfam domains (distinct from TED colors)
    domain_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                     '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

    # Build color mapping for unique Pfam families
    pfam_families = set()
    for protein in proteins_data:
        for domain in protein.get('pfam_domains', []):
            pfam_families.add(domain['pfam_acc'])

    pfam_color_map = {}
    for idx, pfam_acc in enumerate(sorted(pfam_families)):
        pfam_color_map[pfam_acc] = domain_colors[idx % len(domain_colors)]

    # Calculate figure dimensions
    residues_per_line = 1500
    line_height = 60
    margin_left = 200
    margin_right = 50
    margin_top = 40
    margin_bottom = 60  # Extra space for legend
    line_spacing = 10

    total_lines = 0
    for protein in proteins_data:
        lines_needed = (protein['protein_length'] + residues_per_line - 1) // residues_per_line
        protein['lines_needed'] = max(1, lines_needed)
        total_lines += protein['lines_needed']

    # Add space for legend
    legend_lines = (len(pfam_families) + 3) // 4  # 4 items per row
    legend_height = legend_lines * 25 + 20

    fig_width = margin_left + residues_per_line + margin_right
    fig_height = margin_top + (total_lines * (line_height + line_spacing)) + margin_bottom + legend_height

    fig, ax = plt.subplots(figsize=(fig_width / 100, fig_height / 100), dpi=100)
    ax.set_xlim(0, fig_width)
    ax.set_ylim(0, fig_height)
    ax.axis('off')

    # Add title
    ax.text(fig_width / 2, fig_height - 15,
           "Pfam Domain Architecture",
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
                label = f"{protein['accession']} ({protein['mean_plddt']:.2f})"
                ax.text(margin_left - 10, current_y, label,
                       ha='right', va='center', fontsize=9, family='monospace')

            x_start = margin_left
            x_end = margin_left + (line_end - line_start)

            # Draw protein backbone (grey line)
            backbone_height = 8
            backbone = patches.Rectangle(
                (x_start, current_y - backbone_height / 2),
                x_end - x_start, backbone_height,
                linewidth=0, facecolor='#CCCCCC'
            )
            ax.add_patch(backbone)

            # Draw Pfam domains
            for domain in protein.get('pfam_domains', []):
                pfam_acc = domain['pfam_acc']
                pfam_id = domain['pfam_id']
                color = pfam_color_map.get(pfam_acc, '#888888')

                dom_start = domain['start']
                dom_end = domain['end']

                if dom_end >= line_start and dom_start <= line_end:
                    visible_start = max(dom_start, line_start)
                    visible_end = min(dom_end, line_end)

                    domain_x = margin_left + (visible_start - line_start)
                    domain_width = visible_end - visible_start

                    domain_height = 20
                    domain_rect = patches.Rectangle(
                        (domain_x, current_y - domain_height / 2),
                        domain_width, domain_height,
                        linewidth=1, edgecolor='black', facecolor=color, alpha=0.8
                    )
                    ax.add_patch(domain_rect)

                    # Add Pfam ID label if there's enough space
                    if domain_width > 50:
                        ax.text(domain_x + domain_width / 2, current_y,
                               pfam_id, ha='center', va='center',
                               fontsize=6, weight='bold', color='white')

            # Draw SEED region (black outline box)
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

    # Draw legend
    legend_y = margin_bottom + legend_height - 20
    legend_x = margin_left
    items_per_row = 4  # Fewer items per row to fit longer labels with accession
    item_width = (residues_per_line) // items_per_row

    ax.text(legend_x, legend_y + 15, "Legend:", fontsize=9, weight='bold')

    for idx, (pfam_acc, color) in enumerate(sorted(pfam_color_map.items())):
        row = idx // items_per_row
        col = idx % items_per_row

        x = legend_x + col * item_width
        y = legend_y - row * 25

        # Draw color box
        legend_box = patches.Rectangle(
            (x, y - 5), 15, 10,
            linewidth=1, edgecolor='black', facecolor=color, alpha=0.8
        )
        ax.add_patch(legend_box)

        # Get Pfam ID for this accession from any protein
        pfam_id = pfam_acc
        for protein in proteins_data:
            for domain in protein.get('pfam_domains', []):
                if domain['pfam_acc'] == pfam_acc:
                    pfam_id = domain['pfam_id']
                    break

        ax.text(x + 20, y, f"{pfam_acc} ({pfam_id})", fontsize=7, va='center')

    # Save figure
    output_path = os.path.join(curation_dir, output_file)
    plt.tight_layout()
    plt.savefig(output_path, dpi=100, bbox_inches='tight')
    plt.close()

    print(f"Created Pfam domain visualization: {output_path}")


def run_webshot(uniprot_acc: str, curation_dir: str, output_file: str = 'ted_web.png'):
    """
    Capture a screenshot of the TED website protein page.

    Args:
        uniprot_acc: UniProt accession (may include version like "P12345.2")
        curation_dir: path to curation directory
        output_file: output filename (default: ted_web.png)
    """
    # Strip version number from accession
    base_acc = uniprot_acc.split('.')[0]

    url = f"https://ted.cathdb.info/uniprot/{base_acc}"
    output_path = os.path.join(curation_dir, output_file)

    print(f"Capturing TED website screenshot for {base_acc}...")
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


def main():
    parser = argparse.ArgumentParser(
        description='Generate TED and Pfam domain visualizations for SEED proteins'
    )
    parser.add_argument(
        'curation_dir',
        help='Path to Pfam curation directory containing SEED file'
    )
    parser.add_argument(
        '--max-proteins',
        type=int,
        default=20,
        help='Maximum number of proteins to display (default: 20)'
    )
    parser.add_argument(
        '--ted-only',
        action='store_true',
        help='Only generate TED domain visualization'
    )
    parser.add_argument(
        '--pfam-only',
        action='store_true',
        help='Only generate Pfam domain visualization'
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Force regeneration even if output files exist'
    )
    parser.add_argument(
        '--skip-webshot',
        action='store_true',
        help='Skip capturing TED website screenshot'
    )

    args = parser.parse_args()

    # Check SEED file exists
    seed_file = os.path.join(args.curation_dir, 'SEED')
    if not os.path.exists(seed_file):
        print(f"Error: SEED file not found at {seed_file}")
        sys.exit(1)

    # Check if outputs already exist
    ted_png = os.path.join(args.curation_dir, 'ted.png')
    pfam_png = os.path.join(args.curation_dir, 'pfam.png')

    if not args.force:
        if os.path.exists(ted_png) and os.path.exists(pfam_png):
            print(f"Output files already exist in {args.curation_dir}")
            print("Use --force to regenerate")
            sys.exit(0)

    print(f"Processing curation directory: {args.curation_dir}")
    print(f"Max proteins to display: {args.max_proteins}")

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
            sequences.insert(0, (ted_reference['uniprot_acc'], coords[0], coords[1], ''))
        else:
            print(f"  Warning: Could not fetch coordinates for {ted_reference['ted_id']}")

    if not sequences:
        print("Error: No sequences found in SEED file")
        sys.exit(1)

    # Limit sequences
    if len(sequences) > args.max_proteins:
        print(f"Limiting to first {args.max_proteins} sequences (skipping {len(sequences) - args.max_proteins})")
        sequences = sequences[:args.max_proteins]

    # Process sequences
    print("\nProcessing sequences...")
    plddt_scores = []
    protein_length_cache = {}
    ted_cache = {}
    pfam_proteins_data = []

    with tempfile.TemporaryDirectory() as tmp_dir:
        for uniprot_acc, start, end, alignment_seq in sequences:
            print(f"\nProcessing {uniprot_acc}/{start}-{end}...")

            # Download AlphaFold model
            cif_file = download_alphafold_model(uniprot_acc, tmp_dir)
            if not cif_file:
                continue

            # Get structure sequence and length
            structure_seq = get_sequence_from_cif(cif_file)
            protein_length = len(structure_seq)
            protein_length_cache[uniprot_acc] = protein_length

            # Calculate mean pLDDT
            mean_plddt, num_residues = calculate_mean_plddt(cif_file, start, end, verbose=True)
            print(f"    Mean pLDDT: {mean_plddt:.2f}")

            # Fetch TED domains
            ted_data = fetch_ted_domains(uniprot_acc, protein_length_cache)
            if ted_data:
                ted_data['protein_length'] = protein_length
                ted_cache[uniprot_acc] = ted_data
                print(f"    Found {len(ted_data['domains'])} TED domain(s)")

            # Fetch Pfam domains
            pfam_domains = fetch_pfam_domains(uniprot_acc)
            print(f"    Found {len(pfam_domains)} Pfam domain(s)")

            # Calculate TED consistency score
            consistency_score = 0.0
            if ted_data:
                consistency_score = calculate_ted_consistency_score(start, end, ted_data['domains'])

            accession = f"{uniprot_acc}/{start}-{end}"

            plddt_scores.append({
                'accession': accession,
                'mean_plddt': mean_plddt,
                'consistency_score': consistency_score
            })

            pfam_proteins_data.append({
                'accession': accession,
                'uniprot_acc': uniprot_acc,
                'seed_start': start,
                'seed_end': end,
                'protein_length': protein_length,
                'pfam_domains': pfam_domains,
                'mean_plddt': mean_plddt
            })

    if not plddt_scores:
        print("\nError: No valid proteins processed")
        sys.exit(1)

    # Calculate mean consistency score
    consistency_scores = [s['consistency_score'] for s in plddt_scores]
    mean_consistency = sum(consistency_scores) / len(consistency_scores) if consistency_scores else 0.0
    print(f"\nMean family TED consistency score: {mean_consistency:.2f}")

    # Write pLDDT scores file
    plddt_file = os.path.join(args.curation_dir, 'pLDDT')
    plddt_scores.sort(key=lambda x: x['mean_plddt'], reverse=True)

    with open(plddt_file, 'w') as f:
        f.write(f"# Mean family TED consistency score: {mean_consistency:.2f}\n")
        for entry in plddt_scores:
            f.write(f"{entry['mean_plddt']:.2f}\t{entry['accession']}\t{entry['consistency_score']:.2f}\n")

    print(f"Wrote pLDDT scores to {plddt_file}")

    # Create TED visualization
    if not args.pfam_only:
        print("\nCreating TED domain visualization...")
        create_ted_visualization(plddt_scores, args.curation_dir, ted_cache, mean_consistency)

        # Capture TED website screenshot for first protein
        if not args.skip_webshot and plddt_scores:
            first_acc = plddt_scores[0]['accession'].split('/')[0]
            run_webshot(first_acc, args.curation_dir)

    # Create Pfam visualization
    if not args.ted_only:
        print("\nCreating Pfam domain visualization...")
        create_pfam_visualization(pfam_proteins_data, args.curation_dir)

    print("\nDone!")


if __name__ == '__main__':
    main()
