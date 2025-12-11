"""
Domain scoring functions.

Defines what makes a region "domain-like" based on structural properties.
Higher scores = better domains.
"""

import numpy as np
from typing import List, Dict, Tuple, Set, Optional
from dataclasses import dataclass
import networkx as nx


@dataclass
class DomainMetrics:
    """Structural metrics for a potential domain."""
    internal_contacts: int
    external_contacts: int
    contact_density: float
    contact_ratio: float
    radius_of_gyration: float
    avg_plddt: float
    size: int
    compactness: float  # Size-normalized Rg


def compute_domain_metrics(
    start: int,
    end: int,
    graph: nx.Graph,
    coords: np.ndarray,
    plddt: np.ndarray,
    residue_to_idx: Dict[int, int]
) -> DomainMetrics:
    """
    Compute structural metrics for a potential domain region.

    Parameters:
    -----------
    start, end : int
        Domain boundaries (1-based residue numbers, inclusive)
    graph : nx.Graph
        Contact graph (nodes are 0-based indices)
    coords : np.ndarray
        CA coordinates (n_residues, 3)
    plddt : np.ndarray
        pLDDT scores per residue
    residue_to_idx : Dict[int, int]
        Map from residue number to 0-based index

    Returns:
    --------
    metrics : DomainMetrics
        Computed metrics
    """
    # Get indices for this domain
    domain_indices = set()
    for resnum in range(start, end + 1):
        if resnum in residue_to_idx:
            domain_indices.add(residue_to_idx[resnum])

    if len(domain_indices) == 0:
        return DomainMetrics(
            internal_contacts=0, external_contacts=0,
            contact_density=0, contact_ratio=0,
            radius_of_gyration=0, avg_plddt=0,
            size=0, compactness=0
        )

    domain_list = sorted(domain_indices)

    # Count internal and external contacts
    internal = 0
    external = 0

    for idx in domain_list:
        if idx not in graph:
            continue
        for neighbor in graph.neighbors(idx):
            weight = graph[idx][neighbor].get('weight', 1.0)
            if neighbor in domain_indices:
                if neighbor > idx:  # Count each edge once
                    internal += weight
            else:
                external += weight

    # Contact density = internal contacts per residue
    size = len(domain_list)
    contact_density = internal / size if size > 0 else 0

    # Contact ratio = internal / external
    contact_ratio = internal / max(external, 1)

    # Radius of gyration
    domain_coords = coords[domain_list]
    centroid = domain_coords.mean(axis=0)
    rg = np.sqrt(np.mean(np.sum((domain_coords - centroid) ** 2, axis=1)))

    # Average pLDDT
    avg_plddt = plddt[domain_list].mean()

    # Compactness: expected Rg for a globular domain scales as ~2.5 * N^0.4
    # Higher compactness score = more compact than expected
    expected_rg = 2.5 * (size ** 0.4)
    compactness = expected_rg / max(rg, 1.0)

    return DomainMetrics(
        internal_contacts=int(internal),
        external_contacts=int(external),
        contact_density=contact_density,
        contact_ratio=contact_ratio,
        radius_of_gyration=rg,
        avg_plddt=avg_plddt,
        size=size,
        compactness=compactness
    )


class DomainScorer:
    """
    Scores potential domains based on structural properties.

    The score combines multiple factors:
    - Contact density (more internal contacts = better)
    - Contact ratio (more internal than external = better)
    - Compactness (smaller Rg for size = better)
    - pLDDT confidence (higher = better)
    """

    def __init__(
        self,
        graph: nx.Graph,
        coords: np.ndarray,
        plddt: np.ndarray,
        residue_numbers: List[int],
        pae_matrix: Optional[np.ndarray] = None,
        weights: Dict[str, float] = None
    ):
        """
        Initialize scorer.

        Parameters:
        -----------
        graph : nx.Graph
            Contact graph
        coords : np.ndarray
            CA coordinates
        plddt : np.ndarray
            pLDDT scores
        residue_numbers : List[int]
            Residue numbers (1-based)
        pae_matrix : np.ndarray, optional
            PAE matrix for additional scoring
        weights : Dict[str, float], optional
            Weights for score components
        """
        self.graph = graph
        self.coords = coords
        self.plddt = plddt
        self.residue_numbers = residue_numbers
        self.pae_matrix = pae_matrix

        # Build residue number to index mapping
        self.residue_to_idx = {rn: i for i, rn in enumerate(residue_numbers)}

        # Default weights (tuned empirically)
        self.weights = weights or {
            'contact_density': 2.0,
            'contact_ratio': 1.0,
            'compactness': 1.0,
            'plddt': 0.5,
            'size_bonus': 0.1,  # Small bonus for larger domains
            'pae_internal': 1.0,  # If PAE available
        }

        # Cache computed metrics
        self._cache: Dict[Tuple[int, int], DomainMetrics] = {}

    def get_metrics(self, start: int, end: int) -> DomainMetrics:
        """Get metrics for a domain, using cache."""
        key = (start, end)
        if key not in self._cache:
            self._cache[key] = compute_domain_metrics(
                start, end, self.graph, self.coords,
                self.plddt, self.residue_to_idx
            )
        return self._cache[key]

    def score(self, start: int, end: int) -> float:
        """
        Score a potential domain.

        Higher score = better domain.

        Parameters:
        -----------
        start, end : int
            Domain boundaries (1-based, inclusive)

        Returns:
        --------
        score : float
            Domain quality score
        """
        metrics = self.get_metrics(start, end)

        if metrics.size == 0:
            return -np.inf

        # Component scores (normalized to ~0-1 range)
        w = self.weights

        # Contact density: log scale, higher = better
        density_score = np.log1p(metrics.contact_density)

        # Contact ratio: log scale, >1 is good
        ratio_score = np.log1p(metrics.contact_ratio)

        # Compactness: >1 means more compact than expected
        compact_score = np.clip(metrics.compactness, 0, 2) - 1

        # pLDDT: normalized to 0-1
        plddt_score = (metrics.avg_plddt - 50) / 50  # Center around 0

        # Size bonus (log scale to not overly favor large domains)
        size_score = np.log(metrics.size)

        score = (
            w['contact_density'] * density_score +
            w['contact_ratio'] * ratio_score +
            w['compactness'] * compact_score +
            w['plddt'] * plddt_score +
            w['size_bonus'] * size_score
        )

        # PAE bonus if available
        if self.pae_matrix is not None:
            pae_score = self._compute_pae_score(start, end)
            score += w['pae_internal'] * pae_score

        return score

    def _compute_pae_score(self, start: int, end: int) -> float:
        """
        Compute PAE-based score component.

        Low internal PAE = good (confident interactions within domain)
        High external PAE = good (independent from rest of protein)
        """
        if self.pae_matrix is None:
            return 0.0

        # Get indices
        domain_indices = []
        for resnum in range(start, end + 1):
            if resnum in self.residue_to_idx:
                domain_indices.append(self.residue_to_idx[resnum])

        if len(domain_indices) < 2:
            return 0.0

        domain_set = set(domain_indices)
        all_indices = set(range(len(self.residue_numbers)))
        external_indices = all_indices - domain_set

        # Internal PAE (lower is better)
        internal_pae = []
        for i in domain_indices:
            for j in domain_indices:
                if i != j:
                    internal_pae.append(self.pae_matrix[i, j])

        avg_internal = np.mean(internal_pae) if internal_pae else 15

        # External PAE (higher is better - indicates independence)
        external_pae = []
        for i in domain_indices:
            for j in external_indices:
                external_pae.append(self.pae_matrix[i, j])

        avg_external = np.mean(external_pae) if external_pae else 15

        # Score: want low internal, high external
        # Normalize: PAE typically 0-30
        internal_score = (15 - avg_internal) / 15  # High when internal PAE is low
        external_score = (avg_external - 10) / 10  # High when external PAE is high

        return internal_score + external_score

    def score_assignment(self, domains: List[Tuple[int, int]]) -> float:
        """Score a complete domain assignment."""
        return sum(self.score(start, end) for start, end in domains)


def create_score_function(
    graph: nx.Graph,
    coords: np.ndarray,
    plddt: np.ndarray,
    residue_numbers: List[int],
    pae_matrix: Optional[np.ndarray] = None
):
    """
    Create a scoring function for use with DP optimizer.

    Returns a callable that takes (start, end) and returns a score.
    """
    scorer = DomainScorer(graph, coords, plddt, residue_numbers, pae_matrix)
    return scorer.score
