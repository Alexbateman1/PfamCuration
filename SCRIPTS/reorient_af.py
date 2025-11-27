#!/usr/bin/env python3
"""
Reorient AlphaFold structure to minimize visual overlap of TED domains.

This script:
1. Fetches the AlphaFold DB structural model in CIF format
2. Fetches domain definitions from TED
3. Calculates a rotation that minimizes visual overlap of TED domains
4. Outputs a PNG image with TED domains colored (ted3d.png)
5. Outputs the transformed structure file (ted_reoriented.cif)

Usage:
    reorient_af.py <uniprot_accession>
    reorient_af.py P12345
    reorient_af.py P12345.2 --output-cif my_structure.cif --output-png my_image.png

Dependencies:
    - biopython (>= 1.79)
    - numpy (>= 1.20.0)
    - matplotlib (>= 3.3.0)

    Install with: pip install -r reorient_af_requirements.txt

Algorithm:
    The rotation is calculated by sampling many viewing angles and selecting the
    one that minimizes 2D overlap between domain bounding boxes. Only C-alpha
    atoms are used for calculations, though all atoms are included in the output.

TED Domain Colors (from TED database):
    Domain 01: #4A79A7 (blue)
    Domain 02: #F28E2C (orange)
    Domain 03: #E15759 (red)
    Domain 04: #76B7B2 (teal)
    Domain 05: #59A14F (green)
    Domain 06: #EDC949 (yellow)
    Domain 07: #AF7AA1 (purple)
    Domain 08: #FF9DA7 (pink)
    Domain 09: #9C755F (brown)
    Domain 10: #BAB0AB (grey)
"""

import argparse
import sys
import os
import json
import urllib.request
import numpy as np
from pathlib import Path

# Import local TED module
from ted_local import TEDLocal


# TED domain color scheme (official TED database colors)
TED_COLORS = [
    '#4A79A7', '#F28E2C', '#E15759', '#76B7B2', '#59A14F',
    '#EDC949', '#AF7AA1', '#FF9DA7', '#9C755F', '#BAB0AB'
]


def fetch_alphafold_model(uniprot_acc, output_file='alphafold_model.cif'):
    """
    Download AlphaFold v6 CIF model for a given UniProt accession.

    Args:
        uniprot_acc: UniProt accession (may include version like "P12345.2")
        output_file: Path to save the CIF file

    Returns:
        str: Path to downloaded CIF file
    """
    # Strip version number from accession for AlphaFold API
    base_acc = uniprot_acc.split('.')[0]

    url = f"https://alphafold.ebi.ac.uk/files/AF-{base_acc}-F1-model_v6.cif"

    print(f"Fetching AlphaFold model from: {url}")

    try:
        urllib.request.urlretrieve(url, output_file)
        print(f"Downloaded AlphaFold model to: {output_file}")
        return output_file
    except Exception as e:
        print(f"Error: Failed to download AlphaFold model: {e}")
        sys.exit(1)


def fetch_ted_domains(uniprot_acc, ted_db_path=None):
    """
    Fetch TED domain information from local database.

    Args:
        uniprot_acc: UniProt accession (may include version)
        ted_db_path: Optional path to TED SQLite database

    Returns:
        list of dicts with 'ted_id', 'segments' (each segment has 'start', 'end')
    """
    print(f"Fetching TED domains from local database...")

    try:
        ted = TEDLocal(ted_db_path)
        domains = ted.get_domains(uniprot_acc)

        if not domains:
            print(f"Warning: No TED domain data found for {uniprot_acc}")
            return []

        print(f"Found {len(domains)} TED domains")
        for i, domain in enumerate(domains):
            print(f"  Domain {i+1}: {domain['ted_id']} with {len(domain['segments'])} segment(s)")

        return domains

    except Exception as e:
        print(f"Error: Failed to fetch TED domains: {e}")
        sys.exit(1)


def parse_cif_structure(cif_file):
    """
    Parse CIF file and extract C-alpha and all atom coordinates.

    Args:
        cif_file: Path to CIF file

    Returns:
        dict with 'ca_coords' (Nx3 numpy array of C-alphas), 'ca_residues' (residue IDs),
        'all_atoms_data' (list of all atom records), 'structure' (BioPython structure object)
    """
    from Bio.PDB import MMCIFParser

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('protein', cif_file)

    # Get the first chain
    chain = next(structure.get_chains())

    ca_coords = []
    ca_residues = []
    all_atoms_data = []

    for residue in chain.get_residues():
        res_id = residue.get_id()[1]

        # Extract C-alpha
        if 'CA' in residue:
            ca_atom = residue['CA']
            ca_coords.append(ca_atom.get_coord())
            ca_residues.append(res_id)

        # Store all atoms for output
        for atom in residue.get_atoms():
            all_atoms_data.append({
                'residue': residue,
                'atom': atom,
                'res_id': res_id,
                'coord': atom.get_coord().copy()
            })

    return {
        'ca_coords': np.array(ca_coords),
        'ca_residues': np.array(ca_residues),
        'all_atoms_data': all_atoms_data,
        'structure': structure
    }


def get_domain_ca_coords(structure_data, ted_domains):
    """
    Extract C-alpha coordinates for each TED domain.

    Args:
        structure_data: Dict from parse_cif_structure
        ted_domains: List of TED domain dicts

    Returns:
        list of numpy arrays, one per domain with C-alpha coordinates
    """
    ca_coords = structure_data['ca_coords']
    ca_residues = structure_data['ca_residues']

    domain_ca_coords = []

    for domain in ted_domains:
        domain_coords = []

        for segment in domain['segments']:
            start = segment['start']
            end = segment['end']

            # Find all C-alphas in this segment
            mask = (ca_residues >= start) & (ca_residues <= end)
            segment_coords = ca_coords[mask]

            if len(segment_coords) > 0:
                domain_coords.append(segment_coords)

        if domain_coords:
            all_domain_coords = np.vstack(domain_coords)
            domain_ca_coords.append(all_domain_coords)
        else:
            print(f"Warning: No C-alpha coordinates found for domain {domain['ted_id']}")
            domain_ca_coords.append(np.array([]))

    return domain_ca_coords


def rotation_matrix_from_angles(theta, phi, psi):
    """
    Create a rotation matrix from Euler angles.

    Args:
        theta: Rotation around X axis (radians)
        phi: Rotation around Y axis (radians)
        psi: Rotation around Z axis (radians)

    Returns:
        3x3 rotation matrix
    """
    # Rotation around X
    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(theta), -np.sin(theta)],
        [0, np.sin(theta), np.cos(theta)]
    ])

    # Rotation around Y
    Ry = np.array([
        [np.cos(phi), 0, np.sin(phi)],
        [0, 1, 0],
        [-np.sin(phi), 0, np.cos(phi)]
    ])

    # Rotation around Z
    Rz = np.array([
        [np.cos(psi), -np.sin(psi), 0],
        [np.sin(psi), np.cos(psi), 0],
        [0, 0, 1]
    ])

    return Rz @ Ry @ Rx


def calculate_2d_overlap_score(domain_ca_coords, rotation_matrix):
    """
    Calculate a score for 2D overlap after applying rotation and projecting to XY plane.
    Lower score means less overlap.

    Args:
        domain_ca_coords: List of numpy arrays (one per domain)
        rotation_matrix: 3x3 rotation matrix

    Returns:
        float: overlap score (lower is better)
    """
    if len(domain_ca_coords) < 2:
        return 0.0

    # Rotate all domain coordinates
    rotated_domains = []
    for coords in domain_ca_coords:
        if len(coords) > 0:
            rotated = coords @ rotation_matrix.T
            # Project to 2D (XY plane)
            rotated_2d = rotated[:, :2]
            rotated_domains.append(rotated_2d)

    if len(rotated_domains) < 2:
        return 0.0

    # Calculate pairwise overlap using bounding box overlap
    total_overlap = 0.0
    num_pairs = 0

    for i in range(len(rotated_domains)):
        for j in range(i + 1, len(rotated_domains)):
            coords_i = rotated_domains[i]
            coords_j = rotated_domains[j]

            # Get bounding boxes
            min_i = np.min(coords_i, axis=0)
            max_i = np.max(coords_i, axis=0)
            min_j = np.min(coords_j, axis=0)
            max_j = np.max(coords_j, axis=0)

            # Calculate overlap in each dimension
            overlap_x = max(0, min(max_i[0], max_j[0]) - max(min_i[0], min_j[0]))
            overlap_y = max(0, min(max_i[1], max_j[1]) - max(min_i[1], min_j[1]))

            # Overlap area
            overlap_area = overlap_x * overlap_y

            # Normalize by the smaller domain's bounding box area
            area_i = (max_i[0] - min_i[0]) * (max_i[1] - min_i[1])
            area_j = (max_j[0] - min_j[0]) * (max_j[1] - min_j[1])
            min_area = min(area_i, area_j)

            if min_area > 0:
                normalized_overlap = overlap_area / min_area
                total_overlap += normalized_overlap
                num_pairs += 1

    if num_pairs > 0:
        return total_overlap / num_pairs
    return 0.0


def find_optimal_rotation(domain_ca_coords, n_samples=1000):
    """
    Find rotation matrix that minimizes 2D overlap by sampling viewing angles.

    Args:
        domain_ca_coords: List of numpy arrays (one per domain)
        n_samples: Number of random rotations to try

    Returns:
        3x3 rotation matrix
    """
    if len(domain_ca_coords) < 2:
        print("Warning: Less than 2 domains, using identity rotation")
        return np.eye(3)

    print(f"Sampling {n_samples} viewing angles to find optimal orientation...")

    best_rotation = np.eye(3)
    best_score = float('inf')

    # Sample random rotations using Fibonacci sphere for good coverage
    golden_ratio = (1 + np.sqrt(5)) / 2

    for i in range(n_samples):
        # Generate rotation using Fibonacci spiral on sphere
        theta = 2 * np.pi * i / golden_ratio
        phi = np.arccos(1 - 2 * (i + 0.5) / n_samples)
        psi = np.random.uniform(0, 2 * np.pi)

        rotation = rotation_matrix_from_angles(theta, phi, psi)
        score = calculate_2d_overlap_score(domain_ca_coords, rotation)

        if score < best_score:
            best_score = score
            best_rotation = rotation

        if (i + 1) % 200 == 0:
            print(f"  Tested {i + 1}/{n_samples} rotations... (best overlap score: {best_score:.4f})")

    print(f"Optimal rotation found with overlap score: {best_score:.4f}")
    return best_rotation


def apply_rotation_to_structure(structure_data, rotation_matrix):
    """
    Apply rotation matrix to all atoms in the structure.

    Args:
        structure_data: Dict from parse_cif_structure
        rotation_matrix: 3x3 numpy array

    Returns:
        Updated structure_data with rotated coordinates
    """
    # Rotate C-alpha coordinates
    structure_data['ca_coords'] = structure_data['ca_coords'] @ rotation_matrix.T

    # Rotate all atom coordinates
    for atom_data in structure_data['all_atoms_data']:
        rotated_coord = atom_data['coord'] @ rotation_matrix.T
        atom_data['atom'].set_coord(rotated_coord)
        atom_data['coord'] = rotated_coord

    return structure_data


def save_structure_with_annotations(structure_data, ted_domains, output_file):
    """
    Save the transformed structure to CIF file with TED domain annotations.

    Args:
        structure_data: Dict with BioPython structure
        ted_domains: List of TED domain dicts
        output_file: Path to output CIF file
    """
    from Bio.PDB import MMCIFIO

    io = MMCIFIO()
    io.set_structure(structure_data['structure'])
    io.save(output_file)

    print(f"Saved transformed structure to: {output_file}")

    # Add domain annotations as comments (CIF files support comments with #)
    with open(output_file, 'r') as f:
        content = f.read()

    # Prepare domain annotation header
    annotations = ["# TED Domain Annotations"]
    for i, domain in enumerate(ted_domains):
        color = TED_COLORS[i % len(TED_COLORS)]
        annotations.append(f"# Domain {i+1}: {domain['ted_id']} Color: {color}")
        for seg in domain['segments']:
            annotations.append(f"#   Segment: {seg['start']}-{seg['end']}")

    # Insert annotations after the first line
    lines = content.split('\n')
    annotated_content = lines[0] + '\n' + '\n'.join(annotations) + '\n' + '\n'.join(lines[1:])

    with open(output_file, 'w') as f:
        f.write(annotated_content)


def create_3d_visualization(structure_data, ted_domains, output_file='ted3d.png'):
    """
    Create a 3D visualization of the C-alpha backbone with colored TED domains.

    Args:
        structure_data: Dict from parse_cif_structure (with rotated coords)
        ted_domains: List of TED domain dicts
        output_file: Path to output PNG file
    """
    try:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
    except ImportError:
        print("Error: matplotlib not available, skipping visualization")
        return

    ca_coords = structure_data['ca_coords']
    ca_residues = structure_data['ca_residues']

    # Create figure with white background
    fig = plt.figure(figsize=(14, 12), facecolor='white')
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('white')

    # Create a mapping of residue to domain index
    residue_to_domain = {}
    for domain_idx, domain in enumerate(ted_domains):
        for segment in domain['segments']:
            for res in range(segment['start'], segment['end'] + 1):
                residue_to_domain[res] = domain_idx

    # Draw the C-alpha backbone as connected lines, colored by domain
    for i in range(len(ca_coords) - 1):
        res_current = ca_residues[i]
        res_next = ca_residues[i + 1]

        # Determine color for this segment
        domain_idx = residue_to_domain.get(res_current, -1)

        if domain_idx >= 0:
            color = TED_COLORS[domain_idx % len(TED_COLORS)]
            linewidth = 3.0
            alpha = 0.9
        else:
            color = 'lightgray'
            linewidth = 1.5
            alpha = 0.5

        # Draw line segment
        ax.plot([ca_coords[i, 0], ca_coords[i + 1, 0]],
                [ca_coords[i, 1], ca_coords[i + 1, 1]],
                [ca_coords[i, 2], ca_coords[i + 1, 2]],
                color=color, linewidth=linewidth, alpha=alpha)

    # Add legend
    legend_elements = []
    for i, domain in enumerate(ted_domains):
        color = TED_COLORS[i % len(TED_COLORS)]
        ted_label = domain['ted_id'].split('_')[-1] if '_' in domain['ted_id'] else domain['ted_id']
        legend_elements.append(plt.Line2D([0], [0], color=color, linewidth=3, label=ted_label))

    if legend_elements:
        ax.legend(handles=legend_elements, loc='upper right', fontsize=11, framealpha=0.9)

    # Set title
    ax.set_title('AlphaFold Structure - C-alpha Backbone with TED Domains',
                 fontsize=14, weight='bold', pad=20)

    # Remove axis labels and ticks for cleaner look
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_zlabel('')

    # Make grid less prominent
    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('none')
    ax.yaxis.pane.set_edgecolor('none')
    ax.zaxis.pane.set_edgecolor('none')

    # Set viewing angle (looking down Z-axis since we optimized for XY projection)
    ax.view_init(elev=20, azim=45)

    # Equal aspect ratio
    max_range = np.array([ca_coords[:, 0].max()-ca_coords[:, 0].min(),
                          ca_coords[:, 1].max()-ca_coords[:, 1].min(),
                          ca_coords[:, 2].max()-ca_coords[:, 2].min()]).max() / 2.0

    mid_x = (ca_coords[:, 0].max()+ca_coords[:, 0].min()) * 0.5
    mid_y = (ca_coords[:, 1].max()+ca_coords[:, 1].min()) * 0.5
    mid_z = (ca_coords[:, 2].max()+ca_coords[:, 2].min()) * 0.5

    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    # Save figure
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Saved 3D visualization to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Reorient AlphaFold structure to minimize visual overlap of TED domains'
    )
    parser.add_argument(
        'uniprot_acc',
        help='UniProt accession (e.g., P12345 or P12345.2)'
    )
    parser.add_argument(
        '--output-cif',
        default='ted_reoriented.cif',
        help='Output CIF file name (default: ted_reoriented.cif)'
    )
    parser.add_argument(
        '--output-png',
        default='ted3d.png',
        help='Output PNG file name (default: ted3d.png)'
    )
    parser.add_argument(
        '--n-samples',
        type=int,
        default=1000,
        help='Number of viewing angles to sample (default: 1000)'
    )

    args = parser.parse_args()

    print("="*70)
    print("AlphaFold TED Domain Reorientation")
    print("="*70)
    print(f"UniProt accession: {args.uniprot_acc}")
    print()

    # Step 1: Fetch AlphaFold model
    print("[1/6] Fetching AlphaFold model...")
    cif_file = fetch_alphafold_model(args.uniprot_acc, 'alphafold_temp.cif')
    print()

    # Step 2: Fetch TED domains
    print("[2/6] Fetching TED domain definitions...")
    ted_domains = fetch_ted_domains(args.uniprot_acc)

    if not ted_domains:
        print("Error: No TED domains found. Exiting.")
        sys.exit(1)
    print()

    # Step 3: Parse structure (extract C-alphas)
    print("[3/6] Parsing structure and extracting C-alpha atoms...")
    structure_data = parse_cif_structure(cif_file)
    print(f"Extracted {len(structure_data['ca_coords'])} C-alpha atoms")
    print(f"Total atoms in structure: {len(structure_data['all_atoms_data'])}")
    print()

    # Step 4: Calculate optimal rotation using C-alphas only
    print("[4/6] Finding optimal rotation to minimize domain overlap...")
    domain_ca_coords = get_domain_ca_coords(structure_data, ted_domains)
    rotation_matrix = find_optimal_rotation(domain_ca_coords, n_samples=args.n_samples)
    print()

    # Step 5: Apply rotation to ALL atoms
    print("[5/6] Applying rotation to entire structure...")
    structure_data = apply_rotation_to_structure(structure_data, rotation_matrix)
    print()

    # Step 6: Save outputs
    print("[6/6] Saving outputs...")
    save_structure_with_annotations(structure_data, ted_domains, args.output_cif)
    create_3d_visualization(structure_data, ted_domains, args.output_png)
    print()

    # Clean up temporary file
    if os.path.exists('alphafold_temp.cif'):
        os.remove('alphafold_temp.cif')

    print("="*70)
    print("SUCCESS!")
    print("="*70)
    print(f"Output files:")
    print(f"  - Structure: {args.output_cif}")
    print(f"  - Visualization: {args.output_png}")
    print()


if __name__ == '__main__':
    main()
