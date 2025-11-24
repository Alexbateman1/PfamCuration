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
    The rotation is calculated using PCA (Principal Component Analysis) on the
    domain centers of mass. This finds the orientation where domains are maximally
    spread out in 2D space, minimizing visual occlusion when viewing the structure.

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


def fetch_ted_domains(uniprot_acc):
    """
    Fetch TED domain information from TED API.

    Args:
        uniprot_acc: UniProt accession (may include version)

    Returns:
        list of dicts with 'ted_id', 'segments' (each segment has 'start', 'end')
    """
    # Strip version number for TED API
    base_acc = uniprot_acc.split('.')[0]

    url = f"https://ted.cathdb.info/api/v1/uniprot/summary/{base_acc}?skip=0&limit=100"

    print(f"Fetching TED domains from: {url}")

    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            data = json.loads(response.read().decode())

            if not data or 'data' not in data:
                print(f"Warning: No TED domain data found for {uniprot_acc}")
                return []

            domains = []

            for domain in data.get('data', []):
                ted_id = domain.get('ted_id', '')
                chopping = domain.get('chopping', '')

                if not chopping:
                    continue

                # Parse chopping format: "1-100" or "1-100_150-200"
                segments = []
                for segment in chopping.split('_'):
                    if '-' in segment:
                        start, end = segment.split('-')
                        segments.append({
                            'start': int(start),
                            'end': int(end)
                        })

                domains.append({
                    'ted_id': ted_id,
                    'segments': segments
                })

            print(f"Found {len(domains)} TED domains")
            for i, domain in enumerate(domains):
                print(f"  Domain {i+1}: {domain['ted_id']} with {len(domain['segments'])} segment(s)")

            return domains

    except Exception as e:
        print(f"Error: Failed to fetch TED domains: {e}")
        sys.exit(1)


def parse_cif_structure(cif_file):
    """
    Parse CIF file and extract atom coordinates and residue information.

    Args:
        cif_file: Path to CIF file

    Returns:
        dict with 'coords' (Nx3 numpy array), 'residues' (list of residue numbers),
        'atoms' (list of atom records), 'structure' (BioPython structure object)
    """
    from Bio.PDB import MMCIFParser

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('protein', cif_file)

    # Get the first chain
    chain = next(structure.get_chains())

    coords = []
    residues = []
    atoms_data = []

    for residue in chain.get_residues():
        res_id = residue.get_id()[1]
        for atom in residue.get_atoms():
            coord = atom.get_coord()
            coords.append(coord)
            residues.append(res_id)
            atoms_data.append({
                'residue': residue,
                'atom': atom,
                'res_id': res_id,
                'coord': coord
            })

    return {
        'coords': np.array(coords),
        'residues': np.array(residues),
        'atoms_data': atoms_data,
        'structure': structure
    }


def calculate_domain_centers(structure_data, ted_domains):
    """
    Calculate the center of mass for each TED domain.

    Args:
        structure_data: Dict from parse_cif_structure
        ted_domains: List of TED domain dicts

    Returns:
        numpy array of shape (n_domains, 3) with domain centers
    """
    coords = structure_data['coords']
    residues = structure_data['residues']

    domain_centers = []

    for domain in ted_domains:
        domain_coords = []

        for segment in domain['segments']:
            start = segment['start']
            end = segment['end']

            # Find all atoms in this segment
            mask = (residues >= start) & (residues <= end)
            segment_coords = coords[mask]

            if len(segment_coords) > 0:
                domain_coords.append(segment_coords)

        if domain_coords:
            # Concatenate all segment coordinates for this domain
            all_domain_coords = np.vstack(domain_coords)
            # Calculate center of mass (simple average)
            center = np.mean(all_domain_coords, axis=0)
            domain_centers.append(center)
        else:
            print(f"Warning: No coordinates found for domain {domain['ted_id']}")
            domain_centers.append(np.array([0, 0, 0]))

    return np.array(domain_centers)


def find_optimal_rotation(domain_centers):
    """
    Find rotation matrix that maximizes 2D spread of domain centers.

    Uses PCA to find the principal components, then rotates to align
    the first two principal components with the XY plane for optimal viewing.

    Args:
        domain_centers: numpy array of shape (n_domains, 3)

    Returns:
        3x3 rotation matrix
    """
    if len(domain_centers) < 2:
        print("Warning: Less than 2 domains, using identity rotation")
        return np.eye(3)

    # Center the domain centers
    centered = domain_centers - np.mean(domain_centers, axis=0)

    # Perform PCA
    cov_matrix = np.cov(centered.T)
    eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)

    # Sort by eigenvalues (descending)
    idx = eigenvalues.argsort()[::-1]
    eigenvectors = eigenvectors[:, idx]

    # The rotation matrix aligns the principal components with X, Y, Z axes
    # We want the two largest spread directions (PC1 and PC2) to be in the XY plane
    rotation_matrix = eigenvectors.T

    # Ensure proper rotation matrix (det = 1, not -1)
    if np.linalg.det(rotation_matrix) < 0:
        rotation_matrix[:, 2] *= -1

    print(f"Optimal rotation found using PCA")
    print(f"Eigenvalues (variance along PCs): {eigenvalues[idx]}")

    return rotation_matrix


def apply_rotation_to_structure(structure_data, rotation_matrix):
    """
    Apply rotation matrix to all atoms in the structure.

    Args:
        structure_data: Dict from parse_cif_structure
        rotation_matrix: 3x3 numpy array

    Returns:
        Updated structure_data with rotated coordinates
    """
    # Rotate coordinates
    rotated_coords = structure_data['coords'] @ rotation_matrix.T

    # Update atom coordinates in BioPython structure
    for i, atom_data in enumerate(structure_data['atoms_data']):
        atom_data['atom'].set_coord(rotated_coords[i])
        atom_data['coord'] = rotated_coords[i]

    structure_data['coords'] = rotated_coords

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
    Create a 3D visualization of the structure with colored TED domains.

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

    coords = structure_data['coords']
    residues = structure_data['residues']

    # Create figure
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Plot backbone in grey
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
              c='lightgrey', s=1, alpha=0.3, label='Backbone')

    # Plot each TED domain with its color
    for i, domain in enumerate(ted_domains):
        color = TED_COLORS[i % len(TED_COLORS)]
        ted_label = domain['ted_id'].split('_')[-1] if '_' in domain['ted_id'] else domain['ted_id']

        domain_coords = []

        for segment in domain['segments']:
            start = segment['start']
            end = segment['end']

            # Find all atoms in this segment
            mask = (residues >= start) & (residues <= end)
            segment_coords = coords[mask]

            if len(segment_coords) > 0:
                domain_coords.append(segment_coords)

        if domain_coords:
            all_domain_coords = np.vstack(domain_coords)
            ax.scatter(all_domain_coords[:, 0], all_domain_coords[:, 1], all_domain_coords[:, 2],
                      c=color, s=20, alpha=0.8, label=f'{ted_label}')

    # Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('AlphaFold Structure with TED Domains', fontsize=14, weight='bold')

    # Add legend
    ax.legend(loc='upper right', fontsize=10)

    # Set viewing angle for best separation
    ax.view_init(elev=20, azim=45)

    # Equal aspect ratio
    max_range = np.array([coords[:, 0].max()-coords[:, 0].min(),
                          coords[:, 1].max()-coords[:, 1].min(),
                          coords[:, 2].max()-coords[:, 2].min()]).max() / 2.0

    mid_x = (coords[:, 0].max()+coords[:, 0].min()) * 0.5
    mid_y = (coords[:, 1].max()+coords[:, 1].min()) * 0.5
    mid_z = (coords[:, 2].max()+coords[:, 2].min()) * 0.5

    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    # Save figure
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
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

    # Step 3: Parse structure
    print("[3/6] Parsing structure...")
    structure_data = parse_cif_structure(cif_file)
    print(f"Parsed {len(structure_data['coords'])} atoms")
    print()

    # Step 4: Calculate optimal rotation
    print("[4/6] Calculating optimal rotation...")
    domain_centers = calculate_domain_centers(structure_data, ted_domains)
    rotation_matrix = find_optimal_rotation(domain_centers)
    print()

    # Step 5: Apply rotation
    print("[5/6] Applying rotation to structure...")
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
