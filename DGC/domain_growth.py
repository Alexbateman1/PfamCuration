#!/usr/bin/env python3
"""
Competitive Domain Growth Algorithm for Protein Domain Prediction

This module implements a domain growth algorithm that competes for residues
based on PAE connectivity and spatial proximity. Includes multiple merging
strategies to handle over-segmentation from over-seeding.

Algorithm Overview:
1. Place seeds using various strategies (pLDDT peaks, etc.)
2. Grow non-convex domains by competing for residues
3. Apply death phase to remove undersized domains
4. Merge over-segmented domains using one of 4 strategies
5. Evaluate against ground truth using IoU and boundary metrics

Usage:
    python domain_growth.py --test_file ecod_100.txt --merge_strategy simple_merge
"""

import argparse
import json
import logging
import os
import sys
import urllib.request
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
from scipy.spatial.distance import cdist

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# Domain color palette for visualization (10 distinct colors)
DOMAIN_COLORS = [
    "#4A79A7",  # Blue
    "#F28E2C",  # Orange
    "#E15759",  # Red
    "#76B7B2",  # Teal
    "#59A14F",  # Green
    "#EDC949",  # Yellow
    "#AF7AA1",  # Purple
    "#FF9DA7",  # Pink
    "#9C755F",  # Brown
    "#BAB0AB",  # Gray
]


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class Residue:
    """Single residue with coordinates and metadata."""
    index: int
    resnum: int
    resname: str
    ca_coord: np.ndarray
    plddt: float
    chain_id: str = "A"


@dataclass
class BenchmarkDomain:
    """A domain with potentially discontinuous segments for benchmarking."""
    segments: List[Tuple[int, int]]

    @property
    def residues(self) -> Set[int]:
        """Get all residue positions in this domain."""
        res = set()
        for start, end in self.segments:
            res.update(range(start, end + 1))
        return res

    @property
    def size(self) -> int:
        return len(self.residues)

    @property
    def boundaries(self) -> List[int]:
        """Get all boundary positions."""
        bounds = []
        for start, end in self.segments:
            bounds.extend([start, end])
        return bounds

    def __repr__(self):
        return "_".join(f"{s}-{e}" for s, e in self.segments)


@dataclass
class ProteinAnnotation:
    """Ground truth or predicted annotation for a protein."""
    uniprot_acc: str
    domains: List[BenchmarkDomain]
    exclude: bool = False

    @property
    def has_domains(self) -> bool:
        return len(self.domains) > 0


@dataclass
class ProteinResult:
    """Evaluation result for a single protein."""
    uniprot_acc: str
    true_domains: List[BenchmarkDomain]
    pred_domains: List[BenchmarkDomain]
    weighted_iou: float
    domain_count_error: int
    error_type: str
    boundary_f1: Dict[int, float] = field(default_factory=dict)


@dataclass
class BenchmarkResults:
    """Complete benchmark results."""
    protein_results: List[ProteinResult]
    boundary_f1: Dict[int, float]
    overall_iou: float
    domain_mae: float
    under_splits: int
    over_splits: int
    correct_domain_count: int
    perfect_predictions: int

    def summary(self) -> str:
        n = len(self.protein_results)
        lines = [
            "=" * 60,
            "Domain Growth Algorithm Benchmark Results",
            "=" * 60,
            "",
            f"Proteins evaluated: {n}",
            f"Overall length-weighted IoU: {self.overall_iou:.3f}",
            "",
            "Domain Count Analysis:",
            f"  Correct domain count: {self.correct_domain_count} ({100*self.correct_domain_count/n:.1f}%)",
            f"  Perfect predictions:  {self.perfect_predictions} ({100*self.perfect_predictions/n:.1f}%)",
            f"  Under-splits: {self.under_splits} ({100*self.under_splits/n:.1f}%)",
            f"  Over-splits:  {self.over_splits} ({100*self.over_splits/n:.1f}%)",
            f"  Domain count MAE: {self.domain_mae:.2f}",
            "",
            "Boundary F1 Score:",
        ]
        for tol in sorted(self.boundary_f1.keys()):
            lines.append(f"  +/-{tol:2d} residues: F1 = {self.boundary_f1[tol]:.3f}")
        return "\n".join(lines)


# =============================================================================
# NonConvexDomain Class
# =============================================================================

class NonConvexDomain:
    """
    A domain that can grow non-convexly in 3D space.

    Grows by incorporating residues within spatial reach that have
    low PAE to existing members.

    Attributes:
        members: Set of residue indices belonging to this domain
        coords: Reference to full CA coordinate array (Nx3)
        max_reach: Spatial distance threshold in Angstroms
        domain_id: Unique identifier for this domain
    """

    def __init__(
        self,
        seed_index: int,
        coords: np.ndarray,
        max_reach: float = 15.0,
        domain_id: int = 0,
    ):
        """
        Initialize domain with a seed residue.

        Parameters:
            seed_index: Index of the seed residue
            coords: Full CA coordinate array (Nx3)
            max_reach: Maximum distance to consider a residue "in reach"
            domain_id: Unique identifier for this domain
        """
        self.members: Set[int] = {seed_index}
        self.coords = coords
        self.max_reach = max_reach
        self.domain_id = domain_id

    def is_in_reach(self, point_idx: int) -> bool:
        """
        Check if a point is within reach of any domain member.

        Uses minimum distance to any member coordinate.

        Parameters:
            point_idx: Index of the point to check

        Returns:
            True if point is within max_reach of any member
        """
        point_coord = self.coords[point_idx]
        member_coords = self.coords[list(self.members)]

        # Calculate minimum distance to any member
        distances = np.linalg.norm(member_coords - point_coord, axis=1)
        min_dist = np.min(distances)

        return min_dist < self.max_reach

    def desire_score(self, point_idx: int, pae_matrix: np.ndarray) -> float:
        """
        Calculate how much this domain "wants" a residue.

        Score is based on PAE connectivity - lower mean PAE to members
        means higher desire.

        Parameters:
            point_idx: Index of the point to score
            pae_matrix: NxN PAE matrix (predicted aligned error)

        Returns:
            Desire score (0-1), or -inf if point not in reach
        """
        if not self.is_in_reach(point_idx):
            return float('-inf')

        member_list = list(self.members)

        # Get PAE values between point and all members
        pae_to_members = pae_matrix[point_idx, member_list]
        mean_pae = np.mean(pae_to_members)

        # Convert to score: lower PAE = higher score
        # PAE typically 0-30 Angstroms, with <5 being very good
        score = 1.0 / (1.0 + mean_pae)

        return score

    def avg_internal_pae(self, pae_matrix: np.ndarray) -> float:
        """
        Calculate average PAE between all member pairs.

        Parameters:
            pae_matrix: NxN PAE matrix

        Returns:
            Mean PAE within the domain
        """
        if len(self.members) < 2:
            return 0.0

        member_list = list(self.members)
        internal_pae = pae_matrix[np.ix_(member_list, member_list)]

        # Get upper triangle (excluding diagonal)
        n = len(member_list)
        mask = np.triu_indices(n, k=1)
        pae_values = internal_pae[mask]

        if len(pae_values) == 0:
            return 0.0

        return np.mean(pae_values)

    def avg_pae_to_domain(
        self,
        other_domain: 'NonConvexDomain',
        pae_matrix: np.ndarray
    ) -> float:
        """
        Calculate average PAE between this domain and another.

        Parameters:
            other_domain: Another NonConvexDomain
            pae_matrix: NxN PAE matrix

        Returns:
            Mean PAE between the two domains
        """
        self_members = list(self.members)
        other_members = list(other_domain.members)

        inter_pae = pae_matrix[np.ix_(self_members, other_members)]

        return np.mean(inter_pae)

    def add_member(self, idx: int):
        """Add a residue to the domain."""
        self.members.add(idx)

    def remove_member(self, idx: int):
        """Remove a residue from the domain."""
        self.members.discard(idx)

    def merge_with(self, other: 'NonConvexDomain'):
        """Merge another domain into this one."""
        self.members.update(other.members)

    def __len__(self):
        return len(self.members)


# =============================================================================
# Seeding Strategies
# =============================================================================

def find_plddt_peaks(
    coords: np.ndarray,
    plddt: np.ndarray,
    min_plddt: float = 85.0,
    min_distance: float = 15.0,
) -> List[int]:
    """
    Find local maxima in pLDDT scores as seed positions.

    Uses spatial filtering to ensure seeds are well-separated.

    Parameters:
        coords: CA coordinates (Nx3)
        plddt: Per-residue pLDDT scores (N,)
        min_plddt: Minimum pLDDT to consider as seed
        min_distance: Minimum distance between seeds in Angstroms

    Returns:
        List of residue indices to use as seeds
    """
    n = len(plddt)

    # Find high-confidence residues
    high_conf = np.where(plddt >= min_plddt)[0]

    if len(high_conf) == 0:
        # Fallback: use top 10% by pLDDT
        threshold = np.percentile(plddt, 90)
        high_conf = np.where(plddt >= threshold)[0]

    if len(high_conf) == 0:
        # Last resort: middle residue
        return [n // 2]

    # Find local maxima (pLDDT higher than neighbors within window)
    window = 10
    local_maxima = []

    for idx in high_conf:
        start = max(0, idx - window)
        end = min(n, idx + window + 1)

        if plddt[idx] >= np.max(plddt[start:end]) - 0.1:  # Allow small tolerance
            local_maxima.append(idx)

    if len(local_maxima) == 0:
        local_maxima = list(high_conf)

    # Filter by minimum distance
    seeds = []
    for idx in sorted(local_maxima, key=lambda i: -plddt[i]):  # Highest pLDDT first
        if len(seeds) == 0:
            seeds.append(idx)
            continue

        # Check distance to existing seeds
        seed_coords = coords[seeds]
        point_coord = coords[idx]
        distances = np.linalg.norm(seed_coords - point_coord, axis=1)

        if np.min(distances) >= min_distance:
            seeds.append(idx)

    return seeds


def find_grid_seeds(
    coords: np.ndarray,
    plddt: np.ndarray,
    grid_spacing: float = 20.0,
    min_plddt: float = 70.0,
) -> List[int]:
    """
    Place seeds on a regular 3D grid within the protein.

    Parameters:
        coords: CA coordinates (Nx3)
        plddt: Per-residue pLDDT scores (N,)
        grid_spacing: Distance between grid points in Angstroms
        min_plddt: Minimum pLDDT for a residue to be considered as seed

    Returns:
        List of residue indices closest to grid points
    """
    # Get bounding box
    min_coords = np.min(coords, axis=0)
    max_coords = np.max(coords, axis=0)

    # Create grid
    x_points = np.arange(min_coords[0], max_coords[0] + grid_spacing, grid_spacing)
    y_points = np.arange(min_coords[1], max_coords[1] + grid_spacing, grid_spacing)
    z_points = np.arange(min_coords[2], max_coords[2] + grid_spacing, grid_spacing)

    grid_points = []
    for x in x_points:
        for y in y_points:
            for z in z_points:
                grid_points.append([x, y, z])

    if len(grid_points) == 0:
        return [len(coords) // 2]

    grid_points = np.array(grid_points)

    # Find closest residue to each grid point (with good pLDDT)
    valid_indices = np.where(plddt >= min_plddt)[0]
    if len(valid_indices) == 0:
        valid_indices = np.arange(len(plddt))

    valid_coords = coords[valid_indices]

    seeds = set()
    for gp in grid_points:
        distances = np.linalg.norm(valid_coords - gp, axis=1)
        closest_valid_idx = np.argmin(distances)
        seeds.add(valid_indices[closest_valid_idx])

    return list(seeds)


def find_random_seeds(
    coords: np.ndarray,
    plddt: np.ndarray,
    n_seeds: int = 5,
    min_plddt: float = 85.0,
    min_distance: float = 15.0,
) -> List[int]:
    """
    Select random high-confidence residues as seeds.

    Parameters:
        coords: CA coordinates (Nx3)
        plddt: Per-residue pLDDT scores (N,)
        n_seeds: Number of seeds to select
        min_plddt: Minimum pLDDT to consider
        min_distance: Minimum distance between seeds

    Returns:
        List of residue indices
    """
    high_conf = np.where(plddt >= min_plddt)[0]

    if len(high_conf) == 0:
        threshold = np.percentile(plddt, 80)
        high_conf = np.where(plddt >= threshold)[0]

    if len(high_conf) < n_seeds:
        return list(high_conf)

    # Random selection with distance constraint
    np.random.shuffle(high_conf)

    seeds = []
    for idx in high_conf:
        if len(seeds) >= n_seeds:
            break

        if len(seeds) == 0:
            seeds.append(idx)
            continue

        seed_coords = coords[seeds]
        distances = np.linalg.norm(seed_coords - coords[idx], axis=1)

        if np.min(distances) >= min_distance:
            seeds.append(idx)

    return seeds


def get_seeds(
    coords: np.ndarray,
    plddt: np.ndarray,
    strategy: str = "plddt_peaks",
    **kwargs
) -> List[int]:
    """
    Get seed positions using specified strategy.

    Parameters:
        coords: CA coordinates (Nx3)
        plddt: Per-residue pLDDT scores (N,)
        strategy: One of "plddt_peaks", "grid", "random"
        **kwargs: Additional arguments for specific strategy

    Returns:
        List of seed residue indices
    """
    if strategy == "plddt_peaks":
        return find_plddt_peaks(coords, plddt, **kwargs)
    elif strategy == "grid":
        return find_grid_seeds(coords, plddt, **kwargs)
    elif strategy == "random":
        return find_random_seeds(coords, plddt, **kwargs)
    else:
        raise ValueError(f"Unknown seed strategy: {strategy}")


# =============================================================================
# Competitive Domain Growth Algorithm
# =============================================================================

def competitive_domain_growth(
    coords: np.ndarray,
    pae_matrix: np.ndarray,
    seeds: List[int],
    max_reach: float = 15.0,
    acceptance_threshold: float = 0.1,
    min_domain_size: int = 15,
    max_iterations: int = 100,
) -> Tuple[List[NonConvexDomain], Set[int]]:
    """
    Grow domains competitively from seed positions.

    Each iteration:
    1. GROWTH: Unassigned residues join the domain with highest desire_score
    2. DEATH: Domains smaller than min_domain_size are dissolved

    Parameters:
        coords: CA coordinates (Nx3)
        pae_matrix: NxN PAE matrix
        seeds: List of seed residue indices
        max_reach: Maximum distance for spatial reach
        acceptance_threshold: Minimum desire score to join a domain
        min_domain_size: Minimum residues for a domain to survive
        max_iterations: Maximum growth iterations

    Returns:
        (domains, unassigned): List of NonConvexDomain objects and set of unassigned indices
    """
    n_residues = len(coords)

    # Initialize domains from seeds
    domains = []
    for i, seed_idx in enumerate(seeds):
        domain = NonConvexDomain(
            seed_index=seed_idx,
            coords=coords,
            max_reach=max_reach,
            domain_id=i
        )
        domains.append(domain)

    # Track assigned residues
    assigned = set(seeds)
    unassigned = set(range(n_residues)) - assigned

    for iteration in range(max_iterations):
        changes = False

        # GROWTH PHASE
        to_assign = []
        for res_idx in list(unassigned):
            best_score = acceptance_threshold
            best_domain = None

            for domain in domains:
                score = domain.desire_score(res_idx, pae_matrix)
                if score > best_score:
                    best_score = score
                    best_domain = domain

            if best_domain is not None:
                to_assign.append((res_idx, best_domain))

        for res_idx, domain in to_assign:
            domain.add_member(res_idx)
            unassigned.discard(res_idx)
            assigned.add(res_idx)
            changes = True

        # DEATH PHASE
        dead_domains = []
        for domain in domains:
            if len(domain) < min_domain_size:
                # Return members to unassigned
                for member in domain.members:
                    unassigned.add(member)
                    assigned.discard(member)
                dead_domains.append(domain)
                changes = True

        for domain in dead_domains:
            domains.remove(domain)

        # Check convergence
        if not changes:
            break

    return domains, unassigned


# =============================================================================
# Merging Strategies
# =============================================================================

def simple_merge(
    domains: List[NonConvexDomain],
    pae_matrix: np.ndarray,
    merge_threshold: float = 8.0,
) -> List[NonConvexDomain]:
    """
    Strategy 1: Hierarchical Agglomerative Merging

    Iteratively merge domains with lowest inter-PAE until all pairs
    have inter-PAE above threshold.

    Rationale: Single domains have uniform low PAE everywhere and will
    merge completely. Multi-domain proteins stop merging at high-PAE
    domain boundaries.

    Parameters:
        domains: List of NonConvexDomain objects
        pae_matrix: NxN PAE matrix
        merge_threshold: Maximum inter-PAE for merging

    Returns:
        List of merged domains
    """
    if len(domains) <= 1:
        return domains

    while len(domains) > 1:
        # Find pair with minimum inter-domain PAE
        best_i, best_j = -1, -1
        min_inter_pae = float('inf')

        for i in range(len(domains)):
            for j in range(i + 1, len(domains)):
                inter_pae = domains[i].avg_pae_to_domain(domains[j], pae_matrix)
                if inter_pae < min_inter_pae:
                    min_inter_pae = inter_pae
                    best_i, best_j = i, j

        if min_inter_pae >= merge_threshold:
            break

        # Merge best pair
        logger.debug(f"Merging domains {best_i} and {best_j} (inter-PAE: {min_inter_pae:.2f})")
        domains[best_i].merge_with(domains[best_j])
        domains.pop(best_j)

    return domains


def boundary_strength_merge(
    domains: List[NonConvexDomain],
    pae_matrix: np.ndarray,
    strength_threshold: float = 3.0,
) -> List[NonConvexDomain]:
    """
    Strategy 2: Gradient-Based Boundary Strength Merging

    Merge only if there's no clear PAE discontinuity at the boundary.
    Compares inter-domain PAE with intra-domain PAE.

    Rationale: Real domain boundaries have inter_pae >> intra_pae.
    Same-domain segments have inter_pae approximately equal to intra_pae.

    Parameters:
        domains: List of NonConvexDomain objects
        pae_matrix: NxN PAE matrix
        strength_threshold: Minimum boundary strength to keep domains separate

    Returns:
        List of merged domains
    """
    if len(domains) <= 1:
        return domains

    changed = True
    while changed and len(domains) > 1:
        changed = False

        for i in range(len(domains)):
            if changed:
                break

            for j in range(i + 1, len(domains)):
                intra_pae_1 = domains[i].avg_internal_pae(pae_matrix)
                intra_pae_2 = domains[j].avg_internal_pae(pae_matrix)
                inter_pae = domains[i].avg_pae_to_domain(domains[j], pae_matrix)

                # Boundary strength: how much higher is inter-PAE than intra-PAE
                boundary_strength = inter_pae - max(intra_pae_1, intra_pae_2)

                if boundary_strength < strength_threshold:
                    # Weak boundary - merge
                    logger.debug(f"Merging domains {i} and {j} (boundary strength: {boundary_strength:.2f})")
                    domains[i].merge_with(domains[j])
                    domains.pop(j)
                    changed = True
                    break

    return domains


def size_weighted_merge(
    domains: List[NonConvexDomain],
    pae_matrix: np.ndarray,
    merge_threshold: float = 8.0,
) -> List[NonConvexDomain]:
    """
    Strategy 3: Size-Weighted Merging

    Prioritize merging small domains into larger ones. Small domains
    are more likely to be over-segmentation artifacts.

    Rationale: Large domains are more likely to be real structural units.
    Small domains adjacent to large domains with low PAE are likely
    fragmentation artifacts.

    Parameters:
        domains: List of NonConvexDomain objects
        pae_matrix: NxN PAE matrix
        merge_threshold: Maximum inter-PAE for merging

    Returns:
        List of merged domains
    """
    if len(domains) <= 1:
        return domains

    while len(domains) > 1:
        # Score all pairs: prioritize small domains with low inter-PAE
        candidates = []

        for i in range(len(domains)):
            for j in range(i + 1, len(domains)):
                min_size = min(len(domains[i]), len(domains[j]))
                inter_pae = domains[i].avg_pae_to_domain(domains[j], pae_matrix)

                # Score: lower inter-PAE and smaller size = better candidate
                score = inter_pae / (1.0 + min_size / 10.0)
                candidates.append((score, i, j, inter_pae))

        if not candidates:
            break

        # Get best candidate
        candidates.sort(key=lambda x: x[0])
        best_score, best_i, best_j, best_inter_pae = candidates[0]

        if best_inter_pae >= merge_threshold:
            break

        # Merge
        logger.debug(f"Merging domains {best_i} and {best_j} (score: {best_score:.2f}, inter-PAE: {best_inter_pae:.2f})")
        domains[best_i].merge_with(domains[best_j])
        domains.pop(best_j)

    return domains


def spatial_merge(
    domains: List[NonConvexDomain],
    coords: np.ndarray,
    pae_matrix: np.ndarray,
    pae_threshold: float = 8.0,
    distance_threshold: float = 5.0,
) -> List[NonConvexDomain]:
    """
    Strategy 4: Spatial + PAE Merging

    Require both low PAE AND spatial adjacency for merging.
    Prevents merging distant low-PAE regions.

    Rationale: Repeat domains can have low inter-PAE despite being
    spatially separate. This strategy requires both connectivity
    and proximity.

    Parameters:
        domains: List of NonConvexDomain objects
        coords: CA coordinates (Nx3)
        pae_matrix: NxN PAE matrix
        pae_threshold: Maximum inter-PAE for merging
        distance_threshold: Maximum inter-domain distance for merging

    Returns:
        List of merged domains
    """
    if len(domains) <= 1:
        return domains

    changed = True
    while changed and len(domains) > 1:
        changed = False

        for i in range(len(domains)):
            if changed:
                break

            for j in range(i + 1, len(domains)):
                inter_pae = domains[i].avg_pae_to_domain(domains[j], pae_matrix)

                if inter_pae >= pae_threshold:
                    continue

                # Compute minimum distance between domain surfaces
                coords_i = coords[list(domains[i].members)]
                coords_j = coords[list(domains[j].members)]

                dist_matrix = cdist(coords_i, coords_j)
                min_dist = np.min(dist_matrix)

                if min_dist < distance_threshold:
                    # Merge: both PAE and distance criteria met
                    logger.debug(f"Merging domains {i} and {j} (inter-PAE: {inter_pae:.2f}, dist: {min_dist:.2f})")
                    domains[i].merge_with(domains[j])
                    domains.pop(j)
                    changed = True
                    break

    return domains


def apply_merge_strategy(
    domains: List[NonConvexDomain],
    coords: np.ndarray,
    pae_matrix: np.ndarray,
    strategy: str,
    merge_threshold: float = 8.0,
    **kwargs
) -> List[NonConvexDomain]:
    """
    Apply specified merge strategy to domains.

    Parameters:
        domains: List of NonConvexDomain objects
        coords: CA coordinates (Nx3)
        pae_matrix: NxN PAE matrix
        strategy: One of "simple_merge", "boundary_strength", "size_weighted", "spatial_merge", "none"
        merge_threshold: Threshold parameter for merging
        **kwargs: Additional strategy-specific parameters

    Returns:
        List of merged domains
    """
    if strategy == "none":
        return domains
    elif strategy == "simple_merge":
        return simple_merge(domains, pae_matrix, merge_threshold)
    elif strategy == "boundary_strength":
        strength_threshold = kwargs.get("strength_threshold", 3.0)
        return boundary_strength_merge(domains, pae_matrix, strength_threshold)
    elif strategy == "size_weighted":
        return size_weighted_merge(domains, pae_matrix, merge_threshold)
    elif strategy == "spatial_merge":
        distance_threshold = kwargs.get("distance_threshold", 5.0)
        return spatial_merge(domains, coords, pae_matrix, merge_threshold, distance_threshold)
    else:
        raise ValueError(f"Unknown merge strategy: {strategy}")


# =============================================================================
# AlphaFold Data Loading
# =============================================================================

class AlphaFoldLoader:
    """Handle downloading and parsing AlphaFold structure and PAE data."""

    def __init__(self, cache_dir: Optional[str] = None):
        if cache_dir is None:
            self.cache_dir = Path.cwd() / ".dgc_cache"
        else:
            self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def get_structure(self, uniprot_acc: str) -> Tuple[np.ndarray, np.ndarray, List[int]]:
        """
        Get CA coordinates and pLDDT scores for a protein.

        Parameters:
            uniprot_acc: UniProt accession

        Returns:
            (coords, plddt, resnums): CA coordinates (Nx3), pLDDT scores (N,), residue numbers
        """
        pdb_path = self._get_pdb_file(uniprot_acc)
        return self._parse_pdb(pdb_path)

    def get_pae_matrix(self, uniprot_acc: str) -> np.ndarray:
        """
        Get PAE matrix for a protein.

        Parameters:
            uniprot_acc: UniProt accession

        Returns:
            NxN PAE matrix
        """
        pae_path = self._get_pae_file(uniprot_acc)
        return self._parse_pae(pae_path)

    def _get_pdb_file(self, uniprot_acc: str) -> Path:
        """Download PDB/CIF file if not cached."""
        # Try CIF first, then PDB
        for ext in [".cif", ".pdb"]:
            cached = self.cache_dir / f"{uniprot_acc}{ext}"
            if cached.exists():
                return cached

        # Download (try newest version first)
        versions = ["v6", "v4", "v3", "v2"]
        for version in versions:
            for ext, url_ext in [(".cif", ".cif"), (".pdb", ".pdb")]:
                url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_acc}-F1-model_{version}{url_ext}"
                dest = self.cache_dir / f"{uniprot_acc}{ext}"
                try:
                    urllib.request.urlretrieve(url, dest)
                    logger.info(f"Downloaded {uniprot_acc} structure ({version})")
                    return dest
                except Exception:
                    continue

        raise RuntimeError(f"Failed to download structure for {uniprot_acc}")

    def _get_pae_file(self, uniprot_acc: str) -> Path:
        """Download PAE JSON file if not cached."""
        cached = self.cache_dir / f"{uniprot_acc}_pae.json"
        if cached.exists():
            return cached

        versions = ["v6", "v4", "v3", "v2"]
        for version in versions:
            url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_acc}-F1-predicted_aligned_error_{version}.json"
            try:
                urllib.request.urlretrieve(url, cached)
                logger.info(f"Downloaded {uniprot_acc} PAE ({version})")
                return cached
            except Exception:
                continue

        raise RuntimeError(f"Failed to download PAE for {uniprot_acc}")

    def _parse_pdb(self, pdb_path: Path) -> Tuple[np.ndarray, np.ndarray, List[int]]:
        """Parse PDB/CIF file to extract CA coordinates and pLDDT."""
        from Bio.PDB import MMCIFParser, PDBParser

        if pdb_path.suffix.lower() in [".cif", ".mmcif"]:
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)

        structure = parser.get_structure("protein", str(pdb_path))

        coords = []
        plddt = []
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
                    plddt.append(ca.get_bfactor())
                    resnums.append(residue.id[1])
            break  # First model only

        return np.array(coords), np.array(plddt), resnums

    def _parse_pae(self, pae_path: Path) -> np.ndarray:
        """Parse PAE JSON file."""
        with open(pae_path) as f:
            data = json.load(f)

        # Handle different JSON formats
        if isinstance(data, list):
            # Format: [{"predicted_aligned_error": [[...]]}]
            pae = np.array(data[0]["predicted_aligned_error"])
        elif "predicted_aligned_error" in data:
            pae = np.array(data["predicted_aligned_error"])
        else:
            # Try pae key
            pae = np.array(data.get("pae", data))

        return pae


# =============================================================================
# Domain Quality Filtering
# =============================================================================

def build_contact_map(coords: np.ndarray, contact_threshold: float = 10.0) -> np.ndarray:
    """
    Build a contact map from CA coordinates.

    Parameters:
        coords: CA coordinates (Nx3)
        contact_threshold: Distance threshold for contact (Angstroms)

    Returns:
        NxN boolean contact map
    """
    distances = cdist(coords, coords)
    return distances < contact_threshold


def count_long_range_contacts(
    members: Set[int],
    contact_map: np.ndarray,
    min_seq_separation: int = 10,
) -> Tuple[int, float]:
    """
    Count long-range contacts within a domain.

    Long-range contacts are those between residues far apart in sequence,
    which distinguishes globular domains from extended structures like helices.

    Parameters:
        members: Set of residue indices in the domain
        contact_map: NxN boolean contact map
        min_seq_separation: Minimum sequence separation for "long-range"

    Returns:
        (total_lr_contacts, lr_contacts_per_residue)
    """
    if len(members) < 2:
        return 0, 0.0

    member_list = sorted(members)
    lr_contacts = 0

    for i, idx1 in enumerate(member_list):
        for idx2 in member_list[i+1:]:
            seq_sep = abs(idx2 - idx1)
            if seq_sep >= min_seq_separation and contact_map[idx1, idx2]:
                lr_contacts += 1

    lr_per_residue = lr_contacts / len(members)
    return lr_contacts, lr_per_residue


def trim_terminal_extensions(
    members: Set[int],
    contact_map: np.ndarray,
    min_seq_separation: int = 10,
    min_contacts_to_keep: int = 2,
) -> Set[int]:
    """
    Trim terminal extensions that lack long-range contacts to the domain core.

    Works inward from each terminus, removing residues that don't have
    sufficient long-range contacts to the remaining domain.

    Parameters:
        members: Set of residue indices in the domain
        contact_map: NxN boolean contact map
        min_seq_separation: Minimum sequence separation for long-range contact
        min_contacts_to_keep: Minimum long-range contacts to keep a residue

    Returns:
        Set of trimmed member indices
    """
    member_list = sorted(members)
    if len(member_list) < 10:
        return set(member_list)

    trimmed = set(member_list)

    # Trim from N-terminus
    for idx in member_list:
        if idx not in trimmed:
            continue

        lr_contacts = 0
        for other_idx in trimmed:
            if other_idx == idx:
                continue
            seq_sep = abs(other_idx - idx)
            if seq_sep >= min_seq_separation and contact_map[idx, other_idx]:
                lr_contacts += 1

        if lr_contacts < min_contacts_to_keep:
            trimmed.discard(idx)
        else:
            break

    # Trim from C-terminus
    for idx in reversed(member_list):
        if idx not in trimmed:
            continue

        lr_contacts = 0
        for other_idx in trimmed:
            if other_idx == idx:
                continue
            seq_sep = abs(other_idx - idx)
            if seq_sep >= min_seq_separation and contact_map[idx, other_idx]:
                lr_contacts += 1

        if lr_contacts < min_contacts_to_keep:
            trimmed.discard(idx)
        else:
            break

    return trimmed


def calculate_adaptive_threshold(
    pae_matrix: np.ndarray,
    base_threshold: float,
    plddt: np.ndarray = None,
    plddt_cutoff: float = 70.0,
    reference_pae: float = 5.0,
) -> float:
    """
    Calculate an adaptive PAE threshold based on the protein's PAE distribution.

    For proteins with very low PAE in domain regions, this uses relative PAE
    structure to find domain boundaries. Only considers high-pLDDT regions
    to avoid skewing stats with disordered regions.

    Parameters:
        pae_matrix: NxN PAE matrix
        base_threshold: The default threshold to adapt
        plddt: Per-residue pLDDT scores (optional, for masking IDRs)
        plddt_cutoff: Only consider residues with pLDDT >= this value
        reference_pae: Expected median PAE for "normal" proteins

    Returns:
        Adapted threshold
    """
    n = pae_matrix.shape[0]

    # If pLDDT provided, only look at PAE between high-confidence residues
    if plddt is not None:
        high_conf_mask = plddt >= plddt_cutoff
        n_high_conf = np.sum(high_conf_mask)

        if n_high_conf >= 50:  # Need enough residues for meaningful stats
            # Get PAE submatrix for high-confidence residues
            high_conf_idx = np.where(high_conf_mask)[0]
            pae_submatrix = pae_matrix[np.ix_(high_conf_idx, high_conf_idx)]

            # Exclude diagonal
            sub_n = pae_submatrix.shape[0]
            sub_mask = ~np.eye(sub_n, dtype=bool)
            pae_values = pae_submatrix[sub_mask]

            logger.info(f"Adaptive PAE: using {n_high_conf} high-pLDDT residues (>= {plddt_cutoff})")
        else:
            # Fall back to full matrix
            mask = ~np.eye(n, dtype=bool)
            pae_values = pae_matrix[mask]
            logger.info(f"Adaptive PAE: using all residues (only {n_high_conf} with pLDDT >= {plddt_cutoff})")
    else:
        mask = ~np.eye(n, dtype=bool)
        pae_values = pae_matrix[mask]

    pae_median = np.percentile(pae_values, 50)
    pae_min = np.min(pae_values)
    pae_max = np.max(pae_values)
    pae_75 = np.percentile(pae_values, 75)

    logger.info(f"PAE stats (domain regions): min={pae_min:.1f}, median={pae_median:.1f}, "
               f"75th={pae_75:.1f}, max={pae_max:.1f}")

    # For proteins with relatively low PAE in domain regions, use percentile-based threshold
    # Domain boundaries should have relatively higher PAE than within-domain regions
    if pae_median < 6.0:
        # High confidence protein - threshold should be above within-domain PAE
        # but below inter-domain PAE. Use a point between 75th percentile and max.
        # This handles cases where many values cluster at 75th percentile.
        adapted = pae_75 + (pae_max - pae_75) * 0.25
        logger.info(f"High-confidence domains: using threshold = {adapted:.1f}Å "
                   f"(75th + 25% towards max)")
    else:
        # Normal protein - use scaled threshold
        scale_factor = pae_median / reference_pae
        adapted = base_threshold * scale_factor
        adapted = max(base_threshold * 0.2, min(base_threshold * 2.0, adapted))
        logger.info(f"Adaptive threshold: {base_threshold:.1f} -> {adapted:.1f} (scale={scale_factor:.2f})")

    return adapted


def calculate_radius_of_gyration(coords: np.ndarray, members: Set[int]) -> float:
    """
    Calculate radius of gyration for a domain.

    Rg = sqrt(sum((r_i - r_mean)^2) / N)

    Lower values indicate more compact structures.

    Parameters:
        coords: Full CA coordinate array (Nx3)
        members: Set of residue indices in the domain

    Returns:
        Radius of gyration in Angstroms
    """
    if len(members) == 0:
        return 0.0

    member_coords = coords[list(members)]
    center = np.mean(member_coords, axis=0)
    distances_sq = np.sum((member_coords - center) ** 2, axis=1)
    return float(np.sqrt(np.mean(distances_sq)))


def filter_domains_by_quality(
    domains: List[NonConvexDomain],
    coords: np.ndarray,
    min_size: int = 30,
    max_rg_ratio: float = 2.0,
    min_lr_contacts_per_res: float = 0.5,
    trim_terminals: bool = True,
) -> Tuple[List[NonConvexDomain], List[dict]]:
    """
    Filter domains by size, compactness, and long-range contacts.

    Filtering steps:
    1. Trim terminal extensions lacking long-range contacts
    2. Filter by minimum size
    3. Filter by Rg ratio (compactness)
    4. Filter by long-range contacts per residue

    Parameters:
        domains: List of NonConvexDomain objects
        coords: CA coordinates (Nx3)
        min_size: Minimum residues for a valid domain
        max_rg_ratio: Maximum allowed Rg / expected_Rg ratio
        min_lr_contacts_per_res: Minimum long-range contacts per residue
        trim_terminals: Whether to trim terminal extensions

    Returns:
        (filtered_domains, rejected_info): Filtered domains and info about rejected ones
    """
    filtered = []
    rejected = []

    # Build contact map once for all domains
    contact_map = build_contact_map(coords, contact_threshold=10.0)

    for domain in domains:
        original_size = len(domain.members)

        # Step 1: Trim terminal extensions
        if trim_terminals and original_size >= 10:
            trimmed_members = trim_terminal_extensions(
                domain.members, contact_map,
                min_seq_separation=10, min_contacts_to_keep=2
            )
            n_trimmed = original_size - len(trimmed_members)
            if n_trimmed > 0:
                logger.info(f"Domain {domain.domain_id}: trimmed {n_trimmed} terminal residues")
                domain.members = trimmed_members

        size = len(domain.members)

        # Step 2: Filter by minimum size
        if size < min_size:
            rejected.append({
                'domain_id': domain.domain_id,
                'size': size,
                'reason': 'too_small',
            })
            logger.info(f"Rejecting domain {domain.domain_id}: too small ({size} < {min_size} residues)")
            continue

        # Step 3: Calculate Rg and filter by compactness
        rg = calculate_radius_of_gyration(coords, domain.members)
        expected_rg = 2.5 * (size ** (1/3))
        rg_ratio = rg / expected_rg if expected_rg > 0 else 0

        if rg_ratio > max_rg_ratio:
            rejected.append({
                'domain_id': domain.domain_id,
                'size': size,
                'rg': rg,
                'expected_rg': expected_rg,
                'rg_ratio': rg_ratio,
                'reason': 'too_extended',
            })
            logger.info(
                f"Rejecting domain {domain.domain_id}: too extended "
                f"(Rg={rg:.1f}Å, expected={expected_rg:.1f}Å, ratio={rg_ratio:.1f}x > {max_rg_ratio}x)"
            )
            continue

        # Step 4: Filter by long-range contacts
        _, lr_per_res = count_long_range_contacts(domain.members, contact_map)

        if lr_per_res < min_lr_contacts_per_res:
            rejected.append({
                'domain_id': domain.domain_id,
                'size': size,
                'lr_per_res': lr_per_res,
                'reason': 'low_lr_contacts',
            })
            logger.info(
                f"Rejecting domain {domain.domain_id}: insufficient long-range contacts "
                f"({lr_per_res:.2f}/res < {min_lr_contacts_per_res}/res)"
            )
            continue

        filtered.append(domain)

    return filtered, rejected


# =============================================================================
# Domain Prediction Pipeline
# =============================================================================

class DomainGrowthPredictor:
    """
    Complete domain prediction pipeline using competitive growth.

    Parameters:
        max_reach: Spatial reach for domain growth (Angstroms)
        acceptance_threshold: Minimum desire score to add residue
        min_domain_size: Minimum residues for domain survival during growth
        merge_strategy: One of "simple_merge", "boundary_strength",
                       "size_weighted", "spatial_merge", "none"
        merge_threshold: PAE threshold for merging
        seed_strategy: One of "plddt_peaks", "grid", "random"
        max_iterations: Maximum growth iterations
        min_final_size: Minimum residues for final domain (after merging)
        max_rg_ratio: Maximum Rg/expected_Rg ratio (filters extended structures)
        min_lr_contacts: Minimum long-range contacts per residue
        trim_terminals: Whether to trim terminal extensions
        adaptive_pae: Whether to adapt thresholds based on PAE distribution
    """

    def __init__(
        self,
        max_reach: float = 15.0,
        acceptance_threshold: float = 0.1,
        min_domain_size: int = 15,
        merge_strategy: str = "simple_merge",
        merge_threshold: float = 8.0,
        seed_strategy: str = "plddt_peaks",
        max_iterations: int = 100,
        cache_dir: Optional[str] = None,
        min_final_size: int = 30,
        max_rg_ratio: float = 2.0,
        min_lr_contacts: float = 0.5,
        trim_terminals: bool = True,
        adaptive_pae: bool = False,
    ):
        self.max_reach = max_reach
        self.acceptance_threshold = acceptance_threshold
        self.min_domain_size = min_domain_size
        self.merge_strategy = merge_strategy
        self.merge_threshold = merge_threshold
        self.seed_strategy = seed_strategy
        self.max_iterations = max_iterations
        self.min_final_size = min_final_size
        self.max_rg_ratio = max_rg_ratio
        self.min_lr_contacts = min_lr_contacts
        self.trim_terminals = trim_terminals
        self.adaptive_pae = adaptive_pae

        self.loader = AlphaFoldLoader(cache_dir)

    def predict(self, uniprot_acc: str) -> List[BenchmarkDomain]:
        """
        Predict domains for a protein.

        Parameters:
            uniprot_acc: UniProt accession

        Returns:
            List of BenchmarkDomain objects
        """
        # Load data
        coords, plddt, resnums = self.loader.get_structure(uniprot_acc)
        pae_matrix = self.loader.get_pae_matrix(uniprot_acc)

        # Calculate adaptive thresholds if enabled
        merge_threshold = self.merge_threshold
        if self.adaptive_pae:
            merge_threshold = calculate_adaptive_threshold(
                pae_matrix, self.merge_threshold,
                plddt=plddt, plddt_cutoff=70.0, reference_pae=5.0
            )
            logger.info(f"{uniprot_acc}: adaptive merge threshold = {merge_threshold:.1f}Å "
                       f"(base={self.merge_threshold})")

        # Get seeds
        seeds = get_seeds(coords, plddt, self.seed_strategy)
        logger.info(f"{uniprot_acc}: {len(seeds)} seeds from {self.seed_strategy}")

        if len(seeds) == 0:
            return []

        # Grow domains
        domains, unassigned = competitive_domain_growth(
            coords=coords,
            pae_matrix=pae_matrix,
            seeds=seeds,
            max_reach=self.max_reach,
            acceptance_threshold=self.acceptance_threshold,
            min_domain_size=self.min_domain_size,
            max_iterations=self.max_iterations,
        )

        logger.info(f"{uniprot_acc}: {len(domains)} domains after growth")

        # Apply merge strategy with potentially adapted threshold
        domains = apply_merge_strategy(
            domains=domains,
            coords=coords,
            pae_matrix=pae_matrix,
            strategy=self.merge_strategy,
            merge_threshold=merge_threshold,
        )

        logger.info(f"{uniprot_acc}: {len(domains)} domains after {self.merge_strategy}")

        # Filter by quality (size, compactness, and long-range contacts)
        domains, rejected = filter_domains_by_quality(
            domains=domains,
            coords=coords,
            min_size=self.min_final_size,
            max_rg_ratio=self.max_rg_ratio,
            min_lr_contacts_per_res=self.min_lr_contacts,
            trim_terminals=self.trim_terminals,
        )

        if rejected:
            logger.info(f"{uniprot_acc}: {len(rejected)} domains rejected by quality filter")
        logger.info(f"{uniprot_acc}: {len(domains)} domains after quality filtering")

        # Convert to BenchmarkDomain objects
        benchmark_domains = []
        for domain in domains:
            segments = self._members_to_segments(sorted(domain.members), resnums)
            benchmark_domains.append(BenchmarkDomain(segments=segments))

        return benchmark_domains

    def _members_to_segments(
        self,
        members: List[int],
        resnums: List[int],
        max_gap: int = 20,
    ) -> List[Tuple[int, int]]:
        """
        Convert member indices to residue number segments.

        Stitches together segments separated by gaps up to max_gap residues.
        This prevents over-fragmentation of domains.
        """
        if not members:
            return []

        segments = []
        start_resnum = resnums[members[0]]
        prev_resnum = start_resnum

        for idx in members[1:]:
            resnum = resnums[idx]

            # Check for gap > max_gap residues
            if resnum - prev_resnum > max_gap:
                segments.append((start_resnum, prev_resnum))
                start_resnum = resnum

            prev_resnum = resnum

        segments.append((start_resnum, prev_resnum))

        return segments


# =============================================================================
# Evaluation Metrics
# =============================================================================

def parse_test_file(filepath: str) -> List[ProteinAnnotation]:
    """
    Parse test file with ground truth annotations.

    Format:
        ACCESSION domain1 domain2 ...

    Where domains are:
        - "start-end" for continuous domains
        - "start1-end1_start2-end2" for discontinuous domains
        - "?" = exclude from evaluation
    """
    annotations = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith('#'):
                continue

            parts = line.split()
            acc = parts[0]
            domain_strs = parts[1:] if len(parts) > 1 else []

            exclude = False
            if domain_strs and domain_strs[0] == '?':
                exclude = True
                domain_strs = []

            domains = []
            for dom_str in domain_strs:
                segments = []
                for seg_str in dom_str.split('_'):
                    if '-' in seg_str:
                        start, end = seg_str.split('-')
                        segments.append((int(start), int(end)))
                if segments:
                    domains.append(BenchmarkDomain(segments=segments))

            annotations.append(ProteinAnnotation(
                uniprot_acc=acc,
                domains=domains,
                exclude=exclude,
            ))

    return annotations


def calculate_iou(domain1: BenchmarkDomain, domain2: BenchmarkDomain) -> float:
    """Calculate Intersection-over-Union between two domains."""
    res1 = domain1.residues
    res2 = domain2.residues

    intersection = len(res1 & res2)
    union = len(res1 | res2)

    if union == 0:
        return 0.0

    return intersection / union


def match_domains(
    true_domains: List[BenchmarkDomain],
    pred_domains: List[BenchmarkDomain],
) -> Tuple[List[Tuple[BenchmarkDomain, BenchmarkDomain, float]], List[BenchmarkDomain], List[BenchmarkDomain]]:
    """
    Match predicted domains to ground truth using greedy IoU matching.

    Returns:
        (matches, missing, false):
        - matches: List of (true_domain, pred_domain, iou) tuples
        - missing: Unmatched true domains
        - false: Unmatched predicted domains
    """
    if not true_domains and not pred_domains:
        return [], [], []

    if not true_domains:
        return [], [], pred_domains.copy()

    if not pred_domains:
        return [], true_domains.copy(), []

    # Calculate all pairwise IoUs
    iou_matrix = np.zeros((len(true_domains), len(pred_domains)))
    for i, true_dom in enumerate(true_domains):
        for j, pred_dom in enumerate(pred_domains):
            iou_matrix[i, j] = calculate_iou(true_dom, pred_dom)

    # Greedy matching
    matches = []
    matched_true = set()
    matched_pred = set()

    while True:
        best_iou = 0
        best_i, best_j = -1, -1

        for i in range(len(true_domains)):
            if i in matched_true:
                continue
            for j in range(len(pred_domains)):
                if j in matched_pred:
                    continue
                if iou_matrix[i, j] > best_iou:
                    best_iou = iou_matrix[i, j]
                    best_i, best_j = i, j

        if best_iou == 0 or best_i == -1:
            break

        matches.append((true_domains[best_i], pred_domains[best_j], best_iou))
        matched_true.add(best_i)
        matched_pred.add(best_j)

    missing = [true_domains[i] for i in range(len(true_domains)) if i not in matched_true]
    false = [pred_domains[j] for j in range(len(pred_domains)) if j not in matched_pred]

    return matches, missing, false


def calculate_weighted_iou(
    matches: List[Tuple[BenchmarkDomain, BenchmarkDomain, float]],
    true_domains: List[BenchmarkDomain],
) -> float:
    """Calculate length-weighted average IoU."""
    if not true_domains:
        return 1.0

    total_weight = sum(d.size for d in true_domains)
    if total_weight == 0:
        return 1.0

    weighted_sum = 0.0
    for true_dom, pred_dom, iou in matches:
        weighted_sum += iou * true_dom.size

    return weighted_sum / total_weight


def calculate_boundary_f1(
    protein_results: List[ProteinResult],
    tolerance: int,
) -> float:
    """Calculate boundary detection F1 score at given tolerance."""
    tp = 0
    fp = 0
    fn = 0

    for result in protein_results:
        true_bounds = set()
        for dom in result.true_domains:
            true_bounds.update(dom.boundaries)

        pred_bounds = set()
        for dom in result.pred_domains:
            pred_bounds.update(dom.boundaries)

        # Check each true boundary
        for tb in true_bounds:
            found = any(abs(tb - pb) <= tolerance for pb in pred_bounds)
            if found:
                tp += 1
            else:
                fn += 1

        # Check each predicted boundary
        for pb in pred_bounds:
            found = any(abs(pb - tb) <= tolerance for tb in true_bounds)
            if not found:
                fp += 1

    if tp + fp == 0 or tp + fn == 0:
        return 0.0

    precision = tp / (tp + fp)
    recall = tp / (tp + fn)

    if precision + recall == 0:
        return 0.0

    f1 = 2 * precision * recall / (precision + recall)
    return f1


def classify_error(
    true_domains: List[BenchmarkDomain],
    pred_domains: List[BenchmarkDomain],
    matches: List[Tuple[BenchmarkDomain, BenchmarkDomain, float]],
    iou_threshold: float = 0.8,
) -> str:
    """Classify the type of prediction error."""
    n_true = len(true_domains)
    n_pred = len(pred_domains)

    if n_true == 0 and n_pred == 0:
        return "no_domains_correct"
    if n_true == 0 and n_pred > 0:
        return "false_positive"
    if n_true > 0 and n_pred == 0:
        return "false_negative"

    good_matches = [m for m in matches if m[2] >= iou_threshold]

    if len(good_matches) == n_true == n_pred:
        return "correct"

    if n_pred < n_true:
        return "under_split"
    elif n_pred > n_true:
        return "over_split"
    else:
        return "boundary_error"


def evaluate_protein(
    acc: str,
    true_domains: List[BenchmarkDomain],
    pred_domains: List[BenchmarkDomain],
    iou_threshold: float = 0.8,
) -> ProteinResult:
    """Evaluate predictions for a single protein."""
    matches, missing, false = match_domains(true_domains, pred_domains)
    weighted_iou = calculate_weighted_iou(matches, true_domains)

    domain_count_error = len(pred_domains) - len(true_domains)
    error_type = classify_error(true_domains, pred_domains, matches, iou_threshold)

    return ProteinResult(
        uniprot_acc=acc,
        true_domains=true_domains,
        pred_domains=pred_domains,
        weighted_iou=weighted_iou,
        domain_count_error=domain_count_error,
        error_type=error_type,
    )


def run_benchmark(
    annotations: List[ProteinAnnotation],
    predictions: Dict[str, List[BenchmarkDomain]],
    iou_threshold: float = 0.8,
) -> BenchmarkResults:
    """Run full benchmark evaluation."""
    protein_results = []

    for ann in annotations:
        if ann.exclude:
            continue

        if ann.uniprot_acc not in predictions:
            logger.warning(f"No prediction for {ann.uniprot_acc}")
            continue

        result = evaluate_protein(
            ann.uniprot_acc,
            ann.domains,
            predictions[ann.uniprot_acc],
            iou_threshold,
        )
        protein_results.append(result)

    # Calculate overall metrics
    total_weight = sum(
        sum(d.size for d in r.true_domains)
        for r in protein_results
        if r.true_domains
    )

    if total_weight > 0:
        overall_iou = sum(
            r.weighted_iou * sum(d.size for d in r.true_domains)
            for r in protein_results
            if r.true_domains
        ) / total_weight
    else:
        overall_iou = 1.0

    # Domain count analysis
    under_splits = sum(1 for r in protein_results if r.error_type == "under_split")
    over_splits = sum(1 for r in protein_results if r.error_type == "over_split")

    correct_domain_count = sum(
        1 for r in protein_results
        if r.error_type in ["correct", "no_domains_correct", "boundary_error"]
    )

    perfect_predictions = sum(
        1 for r in protein_results
        if r.error_type in ["correct", "no_domains_correct"]
    )

    domain_mae = np.mean([abs(r.domain_count_error) for r in protein_results]) if protein_results else 0.0

    # Boundary F1 at multiple tolerances
    tolerances = [5, 10, 15, 20, 25, 30]
    boundary_f1 = {}
    for tol in tolerances:
        boundary_f1[tol] = calculate_boundary_f1(protein_results, tol)

    return BenchmarkResults(
        protein_results=protein_results,
        boundary_f1=boundary_f1,
        overall_iou=overall_iou,
        domain_mae=domain_mae,
        under_splits=under_splits,
        over_splits=over_splits,
        correct_domain_count=correct_domain_count,
        perfect_predictions=perfect_predictions,
    )


# =============================================================================
# Output Generation
# =============================================================================

def generate_chimerax_commands(
    domains: List[BenchmarkDomain],
    colors: List[str] = None,
) -> str:
    """Generate ChimeraX coloring commands for predicted domains."""
    if colors is None:
        colors = DOMAIN_COLORS

    lines = ["color all white"]

    for i, domain in enumerate(domains):
        color = colors[i % len(colors)]
        for start, end in domain.segments:
            lines.append(f"color :{start}-{end} {color}")

    lines.extend(["hide atoms", "show cartoons"])

    return "\n".join(lines)


def save_chimerax_file(
    uniprot_acc: str,
    domains: List[BenchmarkDomain],
    output_dir: Path,
):
    """Save ChimeraX commands to a .cxc file."""
    output_dir.mkdir(parents=True, exist_ok=True)
    commands = generate_chimerax_commands(domains)
    cxc_path = output_dir / f"{uniprot_acc}.cxc"
    with open(cxc_path, 'w') as f:
        f.write(commands)


def save_results(
    results: BenchmarkResults,
    output_dir: Path,
    strategy: str,
):
    """Save benchmark results to files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Per-protein results CSV
    csv_path = output_dir / f"results_{strategy}.csv"
    with open(csv_path, 'w') as f:
        f.write("accession,true_count,pred_count,weighted_iou,error_type,true_domains,pred_domains\n")
        for r in results.protein_results:
            true_str = ";".join(str(d) for d in r.true_domains)
            pred_str = ";".join(str(d) for d in r.pred_domains)
            f.write(f"{r.uniprot_acc},{len(r.true_domains)},{len(r.pred_domains)},"
                   f"{r.weighted_iou:.3f},{r.error_type},{true_str},{pred_str}\n")

    # Summary text
    summary_path = output_dir / f"summary_{strategy}.txt"
    with open(summary_path, 'w') as f:
        f.write(results.summary())
        f.write("\n\n")
        f.write("=" * 60 + "\n")
        f.write("Per-Protein Details (sorted by IoU)\n")
        f.write("=" * 60 + "\n\n")

        sorted_results = sorted(results.protein_results, key=lambda r: r.weighted_iou)
        for r in sorted_results:
            f.write(f"{r.uniprot_acc}: IoU={r.weighted_iou:.3f}, "
                   f"true={len(r.true_domains)}, pred={len(r.pred_domains)}, "
                   f"type={r.error_type}\n")
            f.write(f"  True: {'; '.join(str(d) for d in r.true_domains) or 'none'}\n")
            f.write(f"  Pred: {'; '.join(str(d) for d in r.pred_domains) or 'none'}\n")
            f.write("\n")

    logger.info(f"Results saved to {output_dir}/")


def compare_strategies(
    all_results: Dict[str, BenchmarkResults],
    output_dir: Path,
):
    """Generate comparison report across all strategies."""
    output_path = output_dir / "comparison_all_strategies.txt"

    with open(output_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("Comparison of Merging Strategies\n")
        f.write("=" * 80 + "\n\n")

        # Header
        f.write(f"{'Strategy':<20} {'IoU':>8} {'Correct%':>10} {'Perfect%':>10} "
               f"{'Under%':>8} {'Over%':>8} {'F1@20':>8}\n")
        f.write("-" * 80 + "\n")

        # Data rows
        for strategy, results in all_results.items():
            n = len(results.protein_results)
            f.write(f"{strategy:<20} "
                   f"{results.overall_iou:>8.3f} "
                   f"{100*results.correct_domain_count/n:>10.1f} "
                   f"{100*results.perfect_predictions/n:>10.1f} "
                   f"{100*results.under_splits/n:>8.1f} "
                   f"{100*results.over_splits/n:>8.1f} "
                   f"{results.boundary_f1.get(20, 0):>8.3f}\n")

        f.write("-" * 80 + "\n\n")

        # Comparison target
        f.write("Comparison Target (Leiden clustering):\n")
        f.write("  Overall IoU: 0.448\n")
        f.write("  Correct count: 34%\n")
        f.write("  Perfect predictions: 28%\n")
        f.write("  Under-splits: 37%\n")
        f.write("  Over-splits: 5%\n")
        f.write("  Boundary F1 at +/-20: 0.608\n")

    logger.info(f"Strategy comparison saved to {output_path}")


# =============================================================================
# Command Line Interface
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Competitive Domain Growth Algorithm for Protein Domain Prediction",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Test simple_merge strategy on 10 proteins
  python domain_growth.py --test-file ecod_100.txt --n-proteins 10 --merge-strategy simple_merge

  # Compare all strategies
  python domain_growth.py --test-file ecod_100.txt --n-proteins 10 --compare-all

  # Full evaluation with tuned parameters
  python domain_growth.py --test-file ecod_100.txt --merge-strategy simple_merge --max-reach 12.0
        """
    )

    # Required arguments
    parser.add_argument(
        "--test-file",
        required=True,
        help="Path to test file with ground truth (ECOD format)"
    )

    # Strategy selection
    parser.add_argument(
        "--merge-strategy",
        choices=["simple_merge", "boundary_strength", "size_weighted", "spatial_merge", "none"],
        default="simple_merge",
        help="Merging strategy to use (default: simple_merge)"
    )
    parser.add_argument(
        "--seed-strategy",
        choices=["plddt_peaks", "grid", "random"],
        default="plddt_peaks",
        help="Seeding strategy to use (default: plddt_peaks)"
    )

    # Algorithm parameters
    parser.add_argument(
        "--max-reach",
        type=float,
        default=15.0,
        help="Spatial reach for domain growth in Angstroms (default: 15.0)"
    )
    parser.add_argument(
        "--acceptance-threshold",
        type=float,
        default=0.1,
        help="Minimum desire score to add residue (default: 0.1)"
    )
    parser.add_argument(
        "--min-domain-size",
        type=int,
        default=15,
        help="Minimum residues for domain survival (default: 15)"
    )
    parser.add_argument(
        "--merge-threshold",
        type=float,
        default=8.0,
        help="PAE threshold for merging (default: 8.0)"
    )
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=100,
        help="Maximum growth iterations (default: 100)"
    )
    parser.add_argument(
        "--min-final-size",
        type=int,
        default=30,
        help="Minimum residues for final domain after filtering (default: 30)"
    )
    parser.add_argument(
        "--max-rg-ratio",
        type=float,
        default=2.0,
        help="Maximum Rg/expected_Rg ratio - filters extended structures (default: 2.0)"
    )
    parser.add_argument(
        "--min-lr-contacts",
        type=float,
        default=0.5,
        help="Minimum long-range contacts per residue (default: 0.5)"
    )
    parser.add_argument(
        "--no-trim-terminals",
        action="store_true",
        help="Disable trimming of terminal extensions"
    )
    parser.add_argument(
        "--adaptive-pae",
        action="store_true",
        help="Enable adaptive PAE thresholds based on protein's PAE distribution"
    )

    # Evaluation options
    parser.add_argument(
        "--n-proteins",
        type=int,
        default=None,
        help="Test on first N proteins (default: all)"
    )
    parser.add_argument(
        "--output-dir",
        default="./results/",
        help="Output directory for results (default: ./results/)"
    )
    parser.add_argument(
        "--compare-all",
        action="store_true",
        help="Compare all merging strategies"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Print detailed progress"
    )
    parser.add_argument(
        "--cache-dir",
        default=None,
        help="Directory for caching AlphaFold files"
    )

    return parser.parse_args()


def main():
    args = parse_args()

    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Parse ground truth
    logger.info(f"Loading ground truth from {args.test_file}...")
    annotations = parse_test_file(args.test_file)

    # Filter to requested number of proteins
    valid_annotations = [a for a in annotations if not a.exclude]
    if args.n_proteins is not None:
        valid_annotations = valid_annotations[:args.n_proteins]

    logger.info(f"Evaluating {len(valid_annotations)} proteins")

    # Output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.compare_all:
        # Compare all strategies
        strategies = ["simple_merge", "boundary_strength", "size_weighted", "spatial_merge", "none"]
        all_results = {}

        for strategy in strategies:
            logger.info(f"\n{'='*60}")
            logger.info(f"Testing strategy: {strategy}")
            logger.info(f"{'='*60}")

            predictor = DomainGrowthPredictor(
                max_reach=args.max_reach,
                acceptance_threshold=args.acceptance_threshold,
                min_domain_size=args.min_domain_size,
                merge_strategy=strategy,
                merge_threshold=args.merge_threshold,
                seed_strategy=args.seed_strategy,
                max_iterations=args.max_iterations,
                cache_dir=args.cache_dir,
                min_final_size=args.min_final_size,
                max_rg_ratio=args.max_rg_ratio,
                min_lr_contacts=args.min_lr_contacts,
                trim_terminals=not args.no_trim_terminals,
                adaptive_pae=args.adaptive_pae,
            )

            predictions = {}
            for ann in valid_annotations:
                try:
                    pred_domains = predictor.predict(ann.uniprot_acc)
                    predictions[ann.uniprot_acc] = pred_domains

                    # Save ChimeraX visualization file
                    chimerax_dir = output_dir / f"chimerax_{strategy}"
                    save_chimerax_file(ann.uniprot_acc, pred_domains, chimerax_dir)

                    # Print prediction summary
                    pred_str = "; ".join(str(d) for d in pred_domains) or "none"
                    true_str = "; ".join(str(d) for d in ann.domains) or "none"
                    logger.info(f"  {ann.uniprot_acc}: {len(pred_domains)} domains - {pred_str}")
                    if args.verbose:
                        logger.info(f"    Truth: {true_str}")

                except Exception as e:
                    logger.error(f"  {ann.uniprot_acc}: Error - {e}")
                    predictions[ann.uniprot_acc] = []

            results = run_benchmark(valid_annotations, predictions)
            all_results[strategy] = results

            logger.info(f"\n{results.summary()}")
            save_results(results, output_dir, strategy)

        # Generate comparison
        compare_strategies(all_results, output_dir)

    else:
        # Single strategy
        logger.info(f"Using merge strategy: {args.merge_strategy}")
        logger.info(f"Using seed strategy: {args.seed_strategy}")

        predictor = DomainGrowthPredictor(
            max_reach=args.max_reach,
            acceptance_threshold=args.acceptance_threshold,
            min_domain_size=args.min_domain_size,
            merge_strategy=args.merge_strategy,
            merge_threshold=args.merge_threshold,
            seed_strategy=args.seed_strategy,
            max_iterations=args.max_iterations,
            cache_dir=args.cache_dir,
            min_final_size=args.min_final_size,
            max_rg_ratio=args.max_rg_ratio,
            min_lr_contacts=args.min_lr_contacts,
            trim_terminals=not args.no_trim_terminals,
            adaptive_pae=args.adaptive_pae,
        )

        predictions = {}
        for ann in valid_annotations:
            try:
                pred_domains = predictor.predict(ann.uniprot_acc)
                predictions[ann.uniprot_acc] = pred_domains

                # Save ChimeraX visualization file
                chimerax_dir = output_dir / "chimerax"
                save_chimerax_file(ann.uniprot_acc, pred_domains, chimerax_dir)

                # Print prediction summary
                pred_str = "; ".join(str(d) for d in pred_domains) or "none"
                true_str = "; ".join(str(d) for d in ann.domains) or "none"
                logger.info(f"{ann.uniprot_acc}: {len(pred_domains)} domains predicted")
                for i, d in enumerate(pred_domains):
                    logger.info(f"  Domain {i+1}: {d}")
                if args.verbose:
                    logger.info(f"  Truth: {true_str}")

            except Exception as e:
                logger.error(f"{ann.uniprot_acc}: Error - {e}")
                import traceback
                if args.verbose:
                    traceback.print_exc()
                predictions[ann.uniprot_acc] = []

        # Run benchmark
        results = run_benchmark(valid_annotations, predictions)

        # Print and save results
        print("\n" + results.summary())
        save_results(results, output_dir, args.merge_strategy)

        # Print per-protein details
        print("\n" + "=" * 60)
        print("Per-Protein Results (sorted by IoU)")
        print("=" * 60)
        sorted_results = sorted(results.protein_results, key=lambda r: r.weighted_iou)
        for r in sorted_results:
            status = "GOOD" if r.weighted_iou >= 0.9 else ("OK" if r.weighted_iou >= 0.7 else "BAD")
            print(f"{r.uniprot_acc}: IoU={r.weighted_iou:.3f} [{status}] "
                 f"true={len(r.true_domains)}, pred={len(r.pred_domains)}, type={r.error_type}")


if __name__ == "__main__":
    main()
