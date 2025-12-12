#!/usr/bin/env python3
"""
Unit tests for the competitive domain growth algorithm.

Tests core algorithm components with synthetic data to verify
correctness without requiring AlphaFold downloads.
"""

import sys
import numpy as np
from pathlib import Path

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent))

from domain_growth import (
    NonConvexDomain,
    competitive_domain_growth,
    simple_merge,
    boundary_strength_merge,
    size_weighted_merge,
    spatial_merge,
    find_plddt_peaks,
    calculate_iou,
    BenchmarkDomain,
)


def create_synthetic_protein(n_residues: int = 100, n_domains: int = 2) -> tuple:
    """
    Create a synthetic protein with clear domain structure.

    Returns:
        (coords, pae_matrix, plddt, domain_boundaries)
    """
    np.random.seed(42)

    # Create coordinates for domains as clusters
    coords = np.zeros((n_residues, 3))
    domain_boundaries = []

    residues_per_domain = n_residues // n_domains

    for d in range(n_domains):
        start = d * residues_per_domain
        end = (d + 1) * residues_per_domain if d < n_domains - 1 else n_residues

        # Create a cluster for this domain
        center = np.array([d * 30.0, 0.0, 0.0])  # 30Å apart
        for i in range(start, end):
            coords[i] = center + np.random.randn(3) * 5.0  # 5Å spread

        domain_boundaries.append((start, end - 1))

    # Create PAE matrix - low PAE within domains, high between
    pae_matrix = np.full((n_residues, n_residues), 20.0)  # High default

    for start, end in domain_boundaries:
        # Low PAE within domain
        for i in range(start, end + 1):
            for j in range(start, end + 1):
                pae_matrix[i, j] = 2.0 + np.random.rand() * 2.0  # 2-4Å

    # Diagonal should be 0
    np.fill_diagonal(pae_matrix, 0.0)

    # pLDDT - high confidence throughout
    plddt = 90.0 + np.random.randn(n_residues) * 5.0
    plddt = np.clip(plddt, 70, 100)

    return coords, pae_matrix, plddt, domain_boundaries


def test_non_convex_domain():
    """Test NonConvexDomain class functionality."""
    print("Testing NonConvexDomain class...")

    n = 50
    coords = np.random.randn(n, 3) * 5.0  # Clustered coords

    # Create domain with seed at index 0
    domain = NonConvexDomain(seed_index=0, coords=coords, max_reach=10.0)

    # Test membership
    assert 0 in domain.members
    assert len(domain) == 1

    # Add members
    domain.add_member(1)
    domain.add_member(2)
    assert len(domain) == 3

    # Test is_in_reach
    # Points close to seed should be in reach
    near_idx = 1  # Already added
    far_idx = n - 1
    coords[far_idx] = np.array([100.0, 100.0, 100.0])  # Move far away

    in_reach_near = domain.is_in_reach(near_idx)
    in_reach_far = domain.is_in_reach(far_idx)

    print(f"  is_in_reach(near): {in_reach_near}")
    print(f"  is_in_reach(far): {in_reach_far}")

    assert in_reach_near == True
    assert in_reach_far == False

    print("  PASSED")


def test_desire_score():
    """Test desire_score calculation."""
    print("Testing desire_score...")

    n = 20
    coords = np.random.randn(n, 3) * 3.0  # Tight cluster

    # Create PAE matrix - uniform low PAE
    pae_matrix = np.full((n, n), 3.0)
    np.fill_diagonal(pae_matrix, 0.0)

    domain = NonConvexDomain(seed_index=0, coords=coords, max_reach=15.0)

    # Test desire score for nearby point
    score = domain.desire_score(1, pae_matrix)
    print(f"  desire_score for nearby point: {score:.3f}")

    # Score should be 1/(1 + mean_pae) = 1/(1+3) = 0.25
    assert 0.2 < score < 0.3

    # Test with far point
    coords[19] = np.array([100.0, 100.0, 100.0])
    score_far = domain.desire_score(19, pae_matrix)
    print(f"  desire_score for far point: {score_far}")

    assert score_far == float('-inf')

    print("  PASSED")


def test_competitive_growth():
    """Test competitive domain growth algorithm."""
    print("Testing competitive_domain_growth...")

    coords, pae_matrix, plddt, true_boundaries = create_synthetic_protein(
        n_residues=100, n_domains=2
    )

    # Seeds at domain centers
    seeds = [25, 75]  # Approximately center of each domain

    domains, unassigned = competitive_domain_growth(
        coords=coords,
        pae_matrix=pae_matrix,
        seeds=seeds,
        max_reach=15.0,
        acceptance_threshold=0.1,
        min_domain_size=10,
        max_iterations=50,
    )

    print(f"  Found {len(domains)} domains")
    print(f"  Unassigned residues: {len(unassigned)}")

    for i, d in enumerate(domains):
        print(f"  Domain {i}: {len(d)} residues")

    # Should find 2 domains
    assert len(domains) == 2, f"Expected 2 domains, got {len(domains)}"

    # Domains should have reasonable sizes (around 50 each)
    sizes = [len(d) for d in domains]
    assert all(s >= 30 for s in sizes), f"Domain sizes too small: {sizes}"

    print("  PASSED")


def test_simple_merge():
    """Test simple_merge strategy."""
    print("Testing simple_merge strategy...")

    n = 100
    coords = np.random.randn(n, 3) * 3.0  # Single cluster

    # Create PAE matrix - uniformly low (single domain)
    pae_matrix = np.full((n, n), 3.0)
    np.fill_diagonal(pae_matrix, 0.0)

    # Create 3 artificial domains from same cluster
    domains = [
        NonConvexDomain(0, coords, 20.0, 0),
        NonConvexDomain(33, coords, 20.0, 1),
        NonConvexDomain(66, coords, 20.0, 2),
    ]

    # Add members to each
    for i in range(33):
        domains[0].add_member(i)
    for i in range(33, 66):
        domains[1].add_member(i)
    for i in range(66, 100):
        domains[2].add_member(i)

    print(f"  Before merge: {len(domains)} domains")

    # They should merge (low inter-PAE)
    merged = simple_merge(domains, pae_matrix, merge_threshold=8.0)

    print(f"  After merge: {len(merged)} domains")

    # Should merge to 1 domain (all have low inter-PAE)
    assert len(merged) == 1, f"Expected 1 domain after merge, got {len(merged)}"
    assert len(merged[0]) == 100

    print("  PASSED")


def test_boundary_strength_merge():
    """Test boundary_strength_merge strategy."""
    print("Testing boundary_strength_merge strategy...")

    coords, pae_matrix, plddt, true_boundaries = create_synthetic_protein(
        n_residues=100, n_domains=2
    )

    # Create domains matching true boundaries
    domains = [
        NonConvexDomain(25, coords, 20.0, 0),
        NonConvexDomain(75, coords, 20.0, 1),
    ]

    for i in range(50):
        domains[0].add_member(i)
    for i in range(50, 100):
        domains[1].add_member(i)

    print(f"  Before merge: {len(domains)} domains")

    # Should NOT merge (strong boundary - high inter-PAE vs low intra-PAE)
    merged = boundary_strength_merge(domains, pae_matrix, strength_threshold=3.0)

    print(f"  After merge: {len(merged)} domains")

    # Should stay as 2 domains
    assert len(merged) == 2, f"Expected 2 domains, got {len(merged)}"

    print("  PASSED")


def test_spatial_merge():
    """Test spatial_merge strategy."""
    print("Testing spatial_merge strategy...")

    n = 100
    coords = np.zeros((n, 3))

    # Two distant clusters
    coords[:50] = np.random.randn(50, 3) * 3.0
    coords[50:] = np.random.randn(50, 3) * 3.0 + np.array([50.0, 0.0, 0.0])

    # But uniform low PAE (misleading)
    pae_matrix = np.full((n, n), 3.0)
    np.fill_diagonal(pae_matrix, 0.0)

    domains = [
        NonConvexDomain(0, coords, 20.0, 0),
        NonConvexDomain(75, coords, 20.0, 1),
    ]

    for i in range(50):
        domains[0].add_member(i)
    for i in range(50, 100):
        domains[1].add_member(i)

    print(f"  Before merge: {len(domains)} domains")

    # Should NOT merge (distant despite low PAE)
    merged = spatial_merge(domains, coords, pae_matrix,
                           pae_threshold=8.0, distance_threshold=10.0)

    print(f"  After merge: {len(merged)} domains")

    assert len(merged) == 2, f"Expected 2 domains (distant clusters), got {len(merged)}"

    print("  PASSED")


def test_plddt_peaks():
    """Test pLDDT peak finding."""
    print("Testing find_plddt_peaks...")

    n = 100
    coords = np.zeros((n, 3))
    for i in range(n):
        coords[i] = np.array([i * 1.0, 0.0, 0.0])  # Linear chain

    # Create pLDDT with peaks at 25 and 75
    plddt = np.ones(n) * 70.0
    plddt[20:30] = 95.0  # Peak 1
    plddt[70:80] = 95.0  # Peak 2

    seeds = find_plddt_peaks(coords, plddt, min_plddt=85.0, min_distance=20.0)

    print(f"  Found {len(seeds)} seeds at positions: {seeds}")

    assert len(seeds) >= 2, f"Expected at least 2 seeds, got {len(seeds)}"

    print("  PASSED")


def test_iou_calculation():
    """Test IoU calculation for benchmark domains."""
    print("Testing IoU calculation...")

    # Identical domains
    d1 = BenchmarkDomain(segments=[(10, 50)])
    d2 = BenchmarkDomain(segments=[(10, 50)])

    iou = calculate_iou(d1, d2)
    print(f"  Identical domains: IoU = {iou:.3f}")
    assert iou == 1.0

    # Partial overlap
    d3 = BenchmarkDomain(segments=[(30, 70)])

    iou2 = calculate_iou(d1, d3)
    print(f"  Partial overlap: IoU = {iou2:.3f}")
    assert 0 < iou2 < 1

    # No overlap
    d4 = BenchmarkDomain(segments=[(100, 150)])

    iou3 = calculate_iou(d1, d4)
    print(f"  No overlap: IoU = {iou3:.3f}")
    assert iou3 == 0.0

    # Discontinuous domain
    d5 = BenchmarkDomain(segments=[(10, 30), (60, 80)])
    print(f"  Discontinuous domain size: {d5.size}")
    assert d5.size == 42  # (30-10+1) + (80-60+1) = 21 + 21 = 42

    print("  PASSED")


def test_avg_internal_pae():
    """Test average internal PAE calculation."""
    print("Testing avg_internal_pae...")

    n = 10
    coords = np.random.randn(n, 3)

    pae_matrix = np.full((n, n), 5.0)
    np.fill_diagonal(pae_matrix, 0.0)

    domain = NonConvexDomain(0, coords, 15.0, 0)
    for i in range(5):
        domain.add_member(i)

    avg_pae = domain.avg_internal_pae(pae_matrix)
    print(f"  avg_internal_pae: {avg_pae:.3f}")

    # Should be approximately 5.0 (the uniform PAE value)
    assert 4.9 < avg_pae < 5.1

    print("  PASSED")


def test_domain_merge():
    """Test domain merging functionality."""
    print("Testing domain merge functionality...")

    n = 20
    coords = np.random.randn(n, 3)

    d1 = NonConvexDomain(0, coords, 15.0, 0)
    d2 = NonConvexDomain(10, coords, 15.0, 1)

    for i in range(5):
        d1.add_member(i)
    for i in range(10, 15):
        d2.add_member(i)

    print(f"  Before merge: d1 has {len(d1)} members, d2 has {len(d2)} members")

    d1.merge_with(d2)

    print(f"  After merge: d1 has {len(d1)} members")

    assert len(d1) == 10
    assert 0 in d1.members
    assert 10 in d1.members

    print("  PASSED")


def run_all_tests():
    """Run all unit tests."""
    print("=" * 60)
    print("Running Domain Growth Algorithm Unit Tests")
    print("=" * 60)
    print()

    tests = [
        test_non_convex_domain,
        test_desire_score,
        test_avg_internal_pae,
        test_domain_merge,
        test_competitive_growth,
        test_simple_merge,
        test_boundary_strength_merge,
        test_spatial_merge,
        test_plddt_peaks,
        test_iou_calculation,
    ]

    passed = 0
    failed = 0

    for test_func in tests:
        try:
            test_func()
            passed += 1
        except AssertionError as e:
            print(f"  FAILED: {e}")
            failed += 1
        except Exception as e:
            print(f"  ERROR: {e}")
            failed += 1
        print()

    print("=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 60)

    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
