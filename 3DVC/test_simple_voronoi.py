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


def count_crossings(assignments, resnums):
    """
    Count chain crossings.
    A crossing occurs when sequential residues (by resnum) are in different partitions.
    """
    crossings = 0
    sorted_indices = np.argsort(resnums)
    sorted_assignments = assignments[sorted_indices]
    sorted_resnums = resnums[sorted_indices]

    for i in range(len(sorted_resnums) - 1):
        # Check if residues are sequential in original numbering
        if sorted_resnums[i + 1] - sorted_resnums[i] == 1:
            # Sequential - check if crossing
            if sorted_assignments[i] != sorted_assignments[i + 1]:
                crossings += 1

    return crossings


def score_partition(assignments, resnums, k, domain_bonus=50, crossing_penalty=50, min_size=30):
    """
    Score a partition.

    Score = domain_bonus * num_valid_domains - crossing_penalty * crossings

    Domains smaller than min_size are invalid (don't count for bonus, add penalty).
    """
    # Count domain sizes
    domain_sizes = np.bincount(assignments, minlength=k)
    valid_domains = np.sum(domain_sizes >= min_size)
    invalid_domains = k - valid_domains

    # Count crossings
    crossings = count_crossings(assignments, resnums)

    # Score: reward valid domains, penalize crossings and invalid tiny domains
    score = (
        domain_bonus * valid_domains
        - crossing_penalty * crossings
        - 50 * invalid_domains  # Penalty for undersized domains
    )

    return score, crossings, valid_domains


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
        score, crossings, valid_domains = score_partition(assignments, resnums, k)

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


def predict_domains(coords, resnums, max_k=10, domain_bonus=50, crossing_penalty=50, min_size=30):
    """
    Predict domains by trying different K values.
    """
    print(f"\nOptimizing with {len(coords)} residues...")
    print(f"Parameters: domain_bonus={domain_bonus}, crossing_penalty={crossing_penalty}, min_size={min_size}")

    best_overall_score = -float('inf')
    best_k = 1
    best_result = None

    for k in range(1, max_k + 1):
        seeds, assignments, score = kmeans_optimize(coords, resnums, k)

        _, crossings, valid_domains = score_partition(assignments, resnums, k,
                                                       domain_bonus, crossing_penalty, min_size)

        print(f"  K={k}: score={score:.1f}, crossings={crossings}, valid_domains={valid_domains}")

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


def test_protein(uniprot_acc, expected_domains=None, max_k=10, **kwargs):
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

    # Predict
    best_k, seeds, assignments = predict_domains(coords, resnums, max_k=max_k, **kwargs)

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
    parser.add_argument("--domain-bonus", type=float, default=50,
                        help="Reward for each valid domain (default: 50)")
    parser.add_argument("--crossing-penalty", type=float, default=50,
                        help="Penalty per chain crossing (default: 50)")
    parser.add_argument("--min-size", type=int, default=30,
                        help="Minimum domain size (default: 30)")
    parser.add_argument("--max-k", type=int, default=10,
                        help="Maximum number of domains to try (default: 10)")

    args = parser.parse_args()

    test_protein(
        args.accession,
        domain_bonus=args.domain_bonus,
        crossing_penalty=args.crossing_penalty,
        min_size=args.min_size,
        max_k=args.max_k,
    )
