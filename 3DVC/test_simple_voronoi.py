#!/usr/bin/env python3
"""
Simple Voronoi domain predictor - test script

Idea: Partition 3D space to minimize chain crossings while rewarding domain creation.

- No contacts, just CA coordinates
- Chain connectivity: sequential edges (1-2, 2-3, etc.)
- Crossing = when residue i and i+1 are in different partitions
- Score = domain_bonus * K - crossing_penalty * crossings
- Minimum domain size to avoid trivial solutions
"""

import argparse
import numpy as np
from scipy.spatial.distance import cdist
from pathlib import Path
import urllib.request
from Bio.PDB import MMCIFParser


CACHE_DIR = Path(__file__).parent / ".3dvc_cache"


def load_structure(uniprot_acc: str):
    """Load AlphaFold structure, return CA coords and pLDDT values."""
    CACHE_DIR.mkdir(exist_ok=True)
    cif_path = CACHE_DIR / f"{uniprot_acc}.cif"

    if not cif_path.exists():
        print(f"Downloading AlphaFold model for {uniprot_acc}...")
        for version in ["v4", "v3", "v2"]:
            url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_acc}-F1-model_{version}.cif"
            try:
                urllib.request.urlretrieve(url, cif_path)
                print(f"  Downloaded {version}")
                break
            except Exception as e:
                continue
        else:
            raise RuntimeError(f"Failed to download structure for {uniprot_acc}")

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("protein", cif_path)

    coords = []
    plddts = []
    resnums = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                if "CA" not in residue:
                    continue
                ca = residue["CA"]
                coords.append(ca.get_coord())
                plddts.append(ca.get_bfactor())
                resnums.append(residue.id[1])
        break

    return np.array(coords), np.array(plddts), np.array(resnums)


def filter_structured(coords, plddts, resnums, plddt_cutoff=70.0):
    """Keep only residues with pLDDT >= cutoff."""
    mask = plddts >= plddt_cutoff
    return coords[mask], plddts[mask], resnums[mask]


def assign_to_seeds(coords, seeds):
    """Assign each residue to nearest seed (Voronoi assignment)."""
    distances = cdist(coords, seeds)
    return np.argmin(distances, axis=1)


def count_crossings(assignments, resnums, max_intrusion_length=3):
    """
    Count chain crossings, distinguishing brief intrusions from real crossings.

    A brief intrusion is when the chain dips into another domain for a short
    stretch (≤ max_intrusion_length residues) then returns to the original domain.
    These should be penalized less than real crossings.

    Returns: (real_crossings, brief_intrusions)
    """
    sorted_indices = np.argsort(resnums)
    sorted_assignments = assignments[sorted_indices]
    sorted_resnums = resnums[sorted_indices]

    # Build list of sequential assignment changes
    # Find runs of same assignment
    runs = []  # List of (domain, start_idx, end_idx)
    if len(sorted_assignments) == 0:
        return 0, 0

    current_domain = sorted_assignments[0]
    run_start = 0

    for i in range(1, len(sorted_assignments)):
        # Check if sequential in original numbering
        if sorted_resnums[i] - sorted_resnums[i-1] != 1:
            # Gap in sequence - end current run, start new one
            runs.append((current_domain, run_start, i-1))
            current_domain = sorted_assignments[i]
            run_start = i
        elif sorted_assignments[i] != current_domain:
            # Domain change
            runs.append((current_domain, run_start, i-1))
            current_domain = sorted_assignments[i]
            run_start = i

    # Don't forget the last run
    runs.append((current_domain, run_start, len(sorted_assignments)-1))

    # Now analyze the runs to find brief intrusions vs real crossings
    real_crossings = 0
    brief_intrusions = 0

    for i in range(1, len(runs)):
        prev_domain, prev_start, prev_end = runs[i-1]
        curr_domain, curr_start, curr_end = runs[i]

        run_length = curr_end - curr_start + 1

        # Check if this is a brief intrusion (returns to previous domain)
        is_brief_intrusion = False
        if run_length <= max_intrusion_length and i + 1 < len(runs):
            next_domain, _, _ = runs[i + 1]
            if next_domain == prev_domain:
                # Pattern: A -> B -> A with short B segment
                is_brief_intrusion = True

        if is_brief_intrusion:
            brief_intrusions += 1
        else:
            real_crossings += 1

    return real_crossings, brief_intrusions


def score_partition(assignments, resnums, k, domain_bonus=50, crossing_penalty=50,
                    intrusion_penalty=10, min_size=30):
    """
    Score a partition.

    Score = domain_bonus * num_valid_domains
            - crossing_penalty * real_crossings
            - intrusion_penalty * brief_intrusions

    Domains smaller than min_size are invalid (don't count for bonus, add penalty).
    """
    # Count domain sizes
    domain_sizes = np.bincount(assignments, minlength=k)
    valid_domains = np.sum(domain_sizes >= min_size)
    invalid_domains = k - valid_domains

    # Count crossings (real vs brief intrusions)
    real_crossings, brief_intrusions = count_crossings(assignments, resnums)

    # Score: reward valid domains, penalize crossings and invalid tiny domains
    score = (
        domain_bonus * valid_domains
        - crossing_penalty * real_crossings
        - intrusion_penalty * brief_intrusions
        - 50 * invalid_domains  # Penalty for undersized domains
    )

    return score, real_crossings, brief_intrusions, valid_domains


def kmeans_optimize(coords, resnums, k, max_iter=100):
    """
    K-means style optimization for Voronoi seeds.
    Optimize seed positions to maximize score.
    """
    n = len(coords)

    # Initialize with k-means++
    seeds = np.zeros((k, 3))
    seeds[0] = coords[np.random.randint(n)]

    for i in range(1, k):
        distances = np.min(cdist(coords, seeds[:i]), axis=1)
        probs = distances ** 2
        probs /= probs.sum()
        idx = np.random.choice(n, p=probs)
        seeds[i] = coords[idx]

    best_score = -float('inf')
    best_seeds = seeds.copy()
    best_assignments = None

    for iteration in range(max_iter):
        # Assign residues to nearest seed
        assignments = assign_to_seeds(coords, seeds)

        # Compute score
        score, _, _, _ = score_partition(assignments, resnums, k)

        if score > best_score:
            best_score = score
            best_seeds = seeds.copy()
            best_assignments = assignments.copy()

        # Update seeds to centroids of assignments
        new_seeds = np.zeros_like(seeds)
        for c in range(k):
            mask = assignments == c
            if mask.sum() > 0:
                new_seeds[c] = coords[mask].mean(axis=0)
            else:
                # Empty cluster - reinitialize randomly
                new_seeds[c] = coords[np.random.randint(n)]

        # Check convergence
        if np.allclose(seeds, new_seeds, atol=0.1):
            break

        seeds = new_seeds

    return best_seeds, best_assignments, best_score


def clip_polygon_by_plane(vertices, plane_point, plane_normal):
    """Clip a polygon by a half-plane (Sutherland-Hodgman)."""
    if len(vertices) < 3:
        return []

    result = []
    n = len(vertices)

    for idx in range(n):
        v1 = vertices[idx]
        v2 = vertices[(idx + 1) % n]

        d1 = np.dot(v1 - plane_point, plane_normal)
        d2 = np.dot(v2 - plane_point, plane_normal)

        if d1 >= 0:
            result.append(v1)

        if (d1 >= 0 and d2 < 0) or (d1 < 0 and d2 >= 0):
            t = d1 / (d1 - d2)
            intersection = v1 + t * (v2 - v1)
            result.append(intersection)

    return result


def compute_voronoi_face(seeds, i, j, max_radius=50.0):
    """Compute the polygon boundary between Voronoi cells i and j."""
    midpoint = (seeds[i] + seeds[j]) / 2
    normal = seeds[j] - seeds[i]
    normal = normal / np.linalg.norm(normal)

    # Create orthogonal vectors in the bisector plane
    if abs(normal[0]) < 0.9:
        u = np.cross(normal, np.array([1.0, 0.0, 0.0]))
    else:
        u = np.cross(normal, np.array([0.0, 1.0, 0.0]))
    u = u / np.linalg.norm(u)
    v = np.cross(normal, u)

    # Initial polygon: large square
    vertices = [
        midpoint + max_radius * u + max_radius * v,
        midpoint - max_radius * u + max_radius * v,
        midpoint - max_radius * u - max_radius * v,
        midpoint + max_radius * u - max_radius * v,
    ]

    # Clip by each other seed's bisector plane
    for k in range(len(seeds)):
        if k == i or k == j:
            continue
        mid_ik = (seeds[i] + seeds[k]) / 2
        norm_ik = seeds[k] - seeds[i]
        norm_ik = norm_ik / np.linalg.norm(norm_ik)
        vertices = clip_polygon_by_plane(vertices, mid_ik, -norm_ik)
        if len(vertices) < 3:
            break

    return vertices


def generate_chimerax_commands(seeds, assignments, resnums, coords, uniprot_acc):
    """Generate ChimeraX commands for visualization."""
    lines = [
        f"# ChimeraX commands for {uniprot_acc} with K={len(seeds)} domains",
        "",
        "# Color domains",
    ]

    colors = ["red", "blue", "green", "yellow", "orange", "purple",
              "cyan", "magenta", "lime", "pink"]

    k = len(seeds)

    # Color each domain
    for c in range(k):
        mask = assignments == c
        if mask.sum() == 0:
            continue
        domain_resnums = resnums[mask]
        color = colors[c % len(colors)]

        # Build segments
        sorted_resnums = sorted(domain_resnums)
        segments = []
        start = sorted_resnums[0]
        end = sorted_resnums[0]
        for r in sorted_resnums[1:]:
            if r == end + 1:
                end = r
            else:
                segments.append((start, end))
                start = end = r
        segments.append((start, end))

        ranges = ",".join([f"{s}-{e}" for s, e in segments])
        lines.append(f"color #1:{ranges} {color}")

    lines.append("")
    lines.append("# Voronoi seed positions (as spheres)")

    for i, seed in enumerate(seeds):
        color = colors[i % len(colors)]
        lines.append(f"shape sphere radius 4 center {seed[0]:.2f},{seed[1]:.2f},{seed[2]:.2f} color {color}")

    if k >= 2:
        lines.append("")
        lines.append("# Voronoi boundary polygons")

        # Find which pairs are neighbors
        pairs_to_draw = []
        for i in range(k):
            for j in range(i + 1, k):
                midpoint = (seeds[i] + seeds[j]) / 2
                dist_ij = np.linalg.norm(seeds[i] - midpoint)
                is_neighbor = True
                for m in range(k):
                    if m == i or m == j:
                        continue
                    dist_m = np.linalg.norm(seeds[m] - midpoint)
                    if dist_m < dist_ij:
                        is_neighbor = False
                        break
                if is_neighbor:
                    pairs_to_draw.append((i, j))

        # Draw polygons as triangles
        for i, j in pairs_to_draw:
            vertices = compute_voronoi_face(seeds, i, j)
            if len(vertices) < 3:
                continue

            centroid = np.mean(vertices, axis=0)
            for idx in range(len(vertices)):
                v1 = vertices[idx]
                v2 = vertices[(idx + 1) % len(vertices)]
                lines.append(
                    f"shape triangle "
                    f"point {centroid[0]:.2f},{centroid[1]:.2f},{centroid[2]:.2f} "
                    f"point {v1[0]:.2f},{v1[1]:.2f},{v1[2]:.2f} "
                    f"point {v2[0]:.2f},{v2[1]:.2f},{v2[2]:.2f} "
                    f"color 128,128,128,51"
                )

    lines.append("")
    lines.append("# Display settings")
    lines.append("cartoon #1")
    lines.append("hide #1 atoms")

    return "\n".join(lines)


def predict_domains(coords, resnums, max_k=10, domain_bonus=50, crossing_penalty=50,
                    intrusion_penalty=10, min_size=30, output_prefix=None, uniprot_acc="protein"):
    """
    Predict domains by trying different K values.
    Saves ChimeraX files for each K if output_prefix is provided.
    """
    print(f"\nOptimizing with {len(coords)} residues...")
    print(f"Parameters: domain_bonus={domain_bonus}, crossing_penalty={crossing_penalty}, "
          f"intrusion_penalty={intrusion_penalty}, min_size={min_size}")

    best_overall_score = -float('inf')
    best_k = 1
    best_result = None

    for k in range(1, max_k + 1):
        seeds, assignments, score = kmeans_optimize(coords, resnums, k)

        _, real_crossings, brief_intrusions, valid_domains = score_partition(
            assignments, resnums, k, domain_bonus, crossing_penalty, intrusion_penalty, min_size
        )

        print(f"  K={k}: score={score:.1f}, real_crossings={real_crossings}, "
              f"intrusions={brief_intrusions}, valid_domains={valid_domains}")

        # Save ChimeraX file for this K
        if output_prefix:
            cxc_content = generate_chimerax_commands(seeds, assignments, resnums, coords, uniprot_acc)
            cxc_path = f"{output_prefix}_K{k}.cxc"
            with open(cxc_path, 'w') as f:
                f.write(cxc_content)
            print(f"       -> Saved {cxc_path}")

        if score > best_overall_score:
            best_overall_score = score
            best_k = k
            best_result = (seeds, assignments)

    print(f"\nBest: K={best_k} with score={best_overall_score:.1f}")
    return best_k, best_result[0], best_result[1]


def assignments_to_domains(assignments, resnums, min_size=30):
    """Convert assignments to domain segments."""
    domains = []
    k = assignments.max() + 1

    for c in range(k):
        mask = assignments == c
        if mask.sum() < min_size:
            continue

        domain_resnums = sorted(resnums[mask])

        # Convert to segments
        segments = []
        start = domain_resnums[0]
        end = domain_resnums[0]

        for r in domain_resnums[1:]:
            if r == end + 1:
                end = r
            else:
                segments.append((start, end))
                start = end = r
        segments.append((start, end))

        domains.append({
            'segments': segments,
            'size': len(domain_resnums),
            'chopping': '_'.join(f"{s}-{e}" for s, e in segments)
        })

    return domains


def test_protein(uniprot_acc, expected_domains=None, max_k=10, output_chimerax=False, **kwargs):
    """Test on a specific protein."""
    print(f"\n{'='*60}")
    print(f"Testing {uniprot_acc}")
    print('='*60)

    # Load structure
    coords, plddts, resnums = load_structure(uniprot_acc)
    print(f"Loaded {len(coords)} residues")

    # Filter by pLDDT
    coords, plddts, resnums = filter_structured(coords, plddts, resnums)
    print(f"After pLDDT filter: {len(coords)} residues")

    # Set output prefix for ChimeraX files
    output_prefix = uniprot_acc if output_chimerax else None

    # Predict
    best_k, seeds, assignments = predict_domains(
        coords, resnums, max_k=max_k,
        output_prefix=output_prefix, uniprot_acc=uniprot_acc,
        **kwargs
    )

    # Convert to domains
    domains = assignments_to_domains(assignments, resnums)

    print(f"\nPredicted {len(domains)} domains:")
    for i, d in enumerate(domains):
        print(f"  {i+1}: {d['chopping']} (size={d['size']})")

    if expected_domains:
        print(f"\nExpected: {expected_domains}")

    return domains


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Simple Voronoi domain predictor - minimize crossings, reward domains"
    )
    parser.add_argument("accession", help="UniProt accession")
    parser.add_argument("--chimerax", action="store_true",
                        help="Output ChimeraX .cxc files for each K value")
    parser.add_argument("--domain-bonus", type=float, default=50,
                        help="Reward for each valid domain (default: 50)")
    parser.add_argument("--crossing-penalty", type=float, default=50,
                        help="Penalty per real chain crossing (default: 50)")
    parser.add_argument("--intrusion-penalty", type=float, default=10,
                        help="Penalty per brief intrusion (default: 10)")
    parser.add_argument("--min-size", type=int, default=30,
                        help="Minimum domain size (default: 30)")
    parser.add_argument("--max-k", type=int, default=10,
                        help="Maximum number of domains to try (default: 10)")

    args = parser.parse_args()

    test_protein(
        args.accession,
        output_chimerax=args.chimerax,
        domain_bonus=args.domain_bonus,
        crossing_penalty=args.crossing_penalty,
        intrusion_penalty=args.intrusion_penalty,
        min_size=args.min_size,
        max_k=args.max_k,
    )
