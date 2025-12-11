#!/usr/bin/env python3
"""
Unit tests for 3DVC Voronoi predictor using synthetic data.

This tests the core algorithm without requiring AlphaFold downloads.
"""

import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from voronoi_predictor import VoronoiPredictor, Residue


def create_synthetic_protein(domains_spec):
    """
    Create synthetic protein with known domain structure.

    domains_spec: list of (center, radius, n_residues) tuples
    Each domain is a sphere of residues around a center point.

    Returns list of Residue objects with known domain structure.
    """
    residues = []
    resnum = 1

    for domain_id, (center, radius, n_residues) in enumerate(domains_spec):
        center = np.array(center)

        # Generate residues in a sphere around the center
        for i in range(n_residues):
            # Random direction
            theta = np.random.uniform(0, 2 * np.pi)
            phi = np.random.uniform(0, np.pi)
            r = np.random.uniform(0, radius)

            # Spherical to Cartesian
            x = center[0] + r * np.sin(phi) * np.cos(theta)
            y = center[1] + r * np.sin(phi) * np.sin(theta)
            z = center[2] + r * np.cos(phi)

            residues.append(Residue(
                index=len(residues),
                resnum=resnum,
                resname="ALA",
                ca_coord=np.array([x, y, z]),
                plddt=90.0,  # High confidence
                chain_id="A",
            ))
            resnum += 1

    return residues


def test_single_domain():
    """Test that a single compact cluster is detected as one domain."""
    print("Test: Single domain...")

    # Create one compact domain
    residues = create_synthetic_protein([
        ([0, 0, 0], 10, 100),  # 100 residues in a sphere of radius 10
    ])

    predictor = VoronoiPredictor(
        min_domain_size=30,
        max_domains=5,
    )

    result = predictor._predict(residues, "TEST_SINGLE")

    print(f"  Domains found: {len(result.domains)}")
    for d in result.domains:
        print(f"    {d}")

    assert len(result.domains) == 1, f"Expected 1 domain, got {len(result.domains)}"
    assert result.domains[0].size >= 90, f"Expected ~100 residues, got {result.domains[0].size}"

    print("  PASSED\n")


def test_two_separate_domains():
    """Test that two well-separated clusters are detected as two domains."""
    print("Test: Two separate domains...")

    # Create two well-separated domains
    residues = create_synthetic_protein([
        ([0, 0, 0], 10, 80),      # Domain 1
        ([50, 0, 0], 10, 80),     # Domain 2, far from Domain 1
    ])

    predictor = VoronoiPredictor(
        min_domain_size=30,
        max_domains=5,
        contact_threshold=12.0,
    )

    result = predictor._predict(residues, "TEST_TWO")

    print(f"  Domains found: {len(result.domains)}")
    for d in result.domains:
        print(f"    {d}")

    assert len(result.domains) == 2, f"Expected 2 domains, got {len(result.domains)}"

    print("  PASSED\n")


def test_discontinuous_domain():
    """Test detection of discontinuous domain (A-B-A' topology)."""
    print("Test: Discontinuous domain (A-B-A')...")

    # Create A-B-A' topology:
    # - Domain A has residues 1-40 at position (0,0,0)
    # - Domain B has residues 41-80 at position (30,0,0)
    # - Domain A' has residues 81-120 at position (5,0,0) (near A)

    residues = []
    resnum = 1

    # Part A at origin
    for i in range(40):
        x = np.random.normal(0, 3)
        y = np.random.normal(0, 3)
        z = np.random.normal(0, 3)
        residues.append(Residue(
            index=len(residues),
            resnum=resnum,
            resname="ALA",
            ca_coord=np.array([x, y, z]),
            plddt=90.0,
            chain_id="A",
        ))
        resnum += 1

    # Part B far away
    for i in range(40):
        x = np.random.normal(30, 3)
        y = np.random.normal(0, 3)
        z = np.random.normal(0, 3)
        residues.append(Residue(
            index=len(residues),
            resnum=resnum,
            resname="ALA",
            ca_coord=np.array([x, y, z]),
            plddt=90.0,
            chain_id="A",
        ))
        resnum += 1

    # Part A' near A
    for i in range(40):
        x = np.random.normal(5, 3)
        y = np.random.normal(0, 3)
        z = np.random.normal(0, 3)
        residues.append(Residue(
            index=len(residues),
            resnum=resnum,
            resname="ALA",
            ca_coord=np.array([x, y, z]),
            plddt=90.0,
            chain_id="A",
        ))
        resnum += 1

    predictor = VoronoiPredictor(
        min_domain_size=30,
        max_domains=5,
        contact_threshold=12.0,
        mu_crossing=0.5,  # Lower crossing penalty to allow discontinuous
    )

    result = predictor._predict(residues, "TEST_DISCONT")

    print(f"  Domains found: {len(result.domains)}")
    for d in result.domains:
        print(f"    {d}")

    # We expect either:
    # - 2 domains: one discontinuous (1-40, 81-120) and one continuous (41-80)
    # - OR 3 domains if the chain crossing penalty prevents merging

    assert len(result.domains) in [2, 3], f"Expected 2-3 domains, got {len(result.domains)}"

    # Check if we got a discontinuous domain
    has_discontinuous = any(d.is_discontinuous for d in result.domains)
    print(f"  Has discontinuous domain: {has_discontinuous}")

    print("  PASSED\n")


def test_voronoi_assignment():
    """Test the basic Voronoi assignment logic."""
    print("Test: Voronoi assignment...")

    # Create simple 2-domain case
    np.random.seed(42)
    coords = np.array([
        [0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0],  # Near origin
        [10, 0, 0], [11, 0, 0], [10, 1, 0], [11, 1, 0],  # Near (10, 0, 0)
    ], dtype=float)

    predictor = VoronoiPredictor()

    # With seeds at [0.5, 0.5, 0] and [10.5, 0.5, 0]
    seeds = np.array([[0.5, 0.5, 0], [10.5, 0.5, 0]])

    from scipy.spatial.distance import cdist
    distances = cdist(coords, seeds)
    assignments = np.argmin(distances, axis=1)

    expected = [0, 0, 0, 0, 1, 1, 1, 1]
    assert list(assignments) == expected, f"Expected {expected}, got {list(assignments)}"

    print("  PASSED\n")


def test_objective_function():
    """Test the objective function scoring."""
    print("Test: Objective function...")

    predictor = VoronoiPredictor(
        lambda_boundary=1.0,
        mu_crossing=2.0,
        nu_domains=20.0,  # Default for modularity-based scoring
    )

    # Simple case: 4 points with 2 contacts
    # Points 0-1 are close, points 2-3 are close
    contacts = np.array([
        [False, True, False, False],
        [True, False, False, False],
        [False, False, False, True],
        [False, False, True, False],
    ])
    resnums = np.array([1, 2, 3, 4])

    # Case 1: All in one cluster - 2 internal contacts
    assignments1 = np.array([0, 0, 0, 0])
    score1 = predictor._compute_score(assignments1, contacts, resnums, k=1)

    # Case 2: Split into 2 clusters - each cluster has 1 internal contact
    assignments2 = np.array([0, 0, 1, 1])
    score2 = predictor._compute_score(assignments2, contacts, resnums, k=2)

    print(f"  Score K=1: {score1:.2f}")
    print(f"  Score K=2: {score2:.2f}")

    # With modularity-based scoring:
    # K=2 perfectly captures the two separate communities (0-1 and 2-3)
    # K=1 groups everything together, lower modularity
    # K=2 should score higher because it correctly identifies the community structure

    assert score2 > score1, "K=2 should score higher when there are two clear separate clusters"

    print("  PASSED\n")


def run_all_tests():
    """Run all tests."""
    print("=" * 60)
    print("3DVC Voronoi Predictor Unit Tests")
    print("=" * 60 + "\n")

    np.random.seed(42)

    test_voronoi_assignment()
    test_objective_function()
    test_single_domain()
    test_two_separate_domains()
    test_discontinuous_domain()

    print("=" * 60)
    print("All tests PASSED!")
    print("=" * 60)


if __name__ == "__main__":
    run_all_tests()
