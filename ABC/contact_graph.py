"""
Contact Graph Builder

Builds a weighted graph where:
- Nodes = Cα atoms (one per residue)
- Edges = spatial proximity (residues within distance threshold)
- Edge weights = f(distance, pLDDT_i, pLDDT_j, PAE_ij)
"""

from typing import TYPE_CHECKING, List, Optional

import networkx as nx
import numpy as np
from scipy.spatial.distance import cdist

if TYPE_CHECKING:
    from .abc_predictor import Residue


class ContactGraphBuilder:
    """
    Builds a contact graph from protein structure.

    Edge weight function:
        w_ij = (pLDDT_i × pLDDT_j / 10000) × exp(-d²_ij/σ²)

    If PAE available, reduce weight for high PAE pairs:
        w_ij *= exp(-PAE_ij / pae_scale)

    Parameters:
    -----------
    distance_threshold : float
        Maximum Cα-Cα distance for contact (default: 10Å)
    sigma : float
        Gaussian decay parameter for edge weights (default: 8.0)
    use_pae : bool
        Whether to use PAE for edge weighting (default: True)
    pae_scale : float
        Scale parameter for PAE penalty (default: 10.0)
    min_sequence_separation : int
        Minimum sequence separation for contacts (default: 6)
        Contacts closer than this in sequence are ignored
        (they're trivially close due to backbone connectivity)
    """

    def __init__(
        self,
        distance_threshold: float = 10.0,
        sigma: float = 8.0,
        use_pae: bool = True,
        pae_scale: float = 10.0,
        min_sequence_separation: int = 6,
    ):
        self.distance_threshold = distance_threshold
        self.sigma = sigma
        self.use_pae = use_pae
        self.pae_scale = pae_scale
        self.min_sequence_separation = min_sequence_separation

    def build(
        self,
        residues: List["Residue"],
        pae_matrix: Optional[np.ndarray] = None,
    ) -> nx.Graph:
        """
        Build contact graph from residues.

        Parameters:
        -----------
        residues : List[Residue]
            List of Residue objects with coordinates and pLDDT
        pae_matrix : np.ndarray, optional
            NxN PAE matrix from AlphaFold

        Returns:
        --------
        nx.Graph
            Weighted contact graph
        """
        n = len(residues)
        if n == 0:
            return nx.Graph()

        # Extract coordinates and pLDDT
        coords = np.array([r.ca_coord for r in residues])
        plddt = np.array([r.plddt for r in residues])

        # Calculate distance matrix
        dist_matrix = cdist(coords, coords)

        # Build graph
        graph = nx.Graph()

        # Add nodes with attributes
        for i, res in enumerate(residues):
            graph.add_node(
                i,
                resnum=res.resnum,
                resname=res.resname,
                plddt=res.plddt,
                chain=res.chain_id,
            )

        # Add edges
        sigma_sq = self.sigma ** 2

        for i in range(n):
            for j in range(i + 1, n):
                # Check sequence separation
                if abs(i - j) < self.min_sequence_separation:
                    continue

                # Check distance threshold
                d = dist_matrix[i, j]
                if d > self.distance_threshold:
                    continue

                # Calculate edge weight
                weight = self._calculate_edge_weight(
                    d, plddt[i], plddt[j], sigma_sq, pae_matrix, i, j
                )

                if weight > 0.01:  # Minimum weight threshold
                    graph.add_edge(i, j, weight=weight, distance=d)

        return graph

    def _calculate_edge_weight(
        self,
        distance: float,
        plddt_i: float,
        plddt_j: float,
        sigma_sq: float,
        pae_matrix: Optional[np.ndarray],
        i: int,
        j: int,
    ) -> float:
        """
        Calculate edge weight.

        w_ij = (pLDDT_i × pLDDT_j / 10000) × exp(-d²/σ²)

        If PAE available:
            w_ij *= exp(-PAE_ij / pae_scale)
        """
        # Normalize pLDDT to [0,1] range (they're usually 0-100)
        plddt_factor = (plddt_i * plddt_j) / 10000.0

        # Distance decay
        distance_factor = np.exp(-(distance ** 2) / sigma_sq)

        weight = plddt_factor * distance_factor

        # PAE penalty
        if self.use_pae and pae_matrix is not None:
            pae_value = pae_matrix[i, j]
            pae_factor = np.exp(-pae_value / self.pae_scale)
            weight *= pae_factor

        return weight

    def build_at_multiple_thresholds(
        self,
        residues: List["Residue"],
        pae_matrix: Optional[np.ndarray] = None,
        thresholds: List[float] = [8.0, 10.0, 12.0, 15.0],
    ) -> dict:
        """
        Build graphs at multiple distance thresholds.

        Useful for finding domains that are stable across scales.

        Parameters:
        -----------
        residues : List[Residue]
            List of Residue objects
        pae_matrix : np.ndarray, optional
            PAE matrix
        thresholds : list
            List of distance thresholds to try

        Returns:
        --------
        Dict[float, nx.Graph]
            Mapping of threshold to graph
        """
        graphs = {}
        original_threshold = self.distance_threshold

        for threshold in thresholds:
            self.distance_threshold = threshold
            graphs[threshold] = self.build(residues, pae_matrix)

        self.distance_threshold = original_threshold
        return graphs

    def get_graph_statistics(self, graph: nx.Graph) -> dict:
        """
        Calculate statistics about the contact graph.

        Returns:
        --------
        dict with keys:
            - n_nodes: number of nodes
            - n_edges: number of edges
            - avg_degree: average node degree
            - density: graph density
            - n_components: number of connected components
            - avg_clustering: average clustering coefficient
            - avg_weight: average edge weight
        """
        if graph.number_of_nodes() == 0:
            return {
                "n_nodes": 0,
                "n_edges": 0,
                "avg_degree": 0,
                "density": 0,
                "n_components": 0,
                "avg_clustering": 0,
                "avg_weight": 0,
            }

        degrees = [d for n, d in graph.degree()]
        weights = [d.get("weight", 1.0) for u, v, d in graph.edges(data=True)]

        return {
            "n_nodes": graph.number_of_nodes(),
            "n_edges": graph.number_of_edges(),
            "avg_degree": np.mean(degrees) if degrees else 0,
            "density": nx.density(graph),
            "n_components": nx.number_connected_components(graph),
            "avg_clustering": nx.average_clustering(graph, weight="weight"),
            "avg_weight": np.mean(weights) if weights else 0,
        }


class ContactDensityCalculator:
    """
    Calculates contact density metrics for domains.
    """

    @staticmethod
    def internal_contacts(
        graph: nx.Graph,
        residue_indices: List[int],
    ) -> float:
        """Count weighted contacts within a set of residues."""
        total = 0.0
        indices_set = set(residue_indices)

        for i in residue_indices:
            for j in graph.neighbors(i):
                if j in indices_set and j > i:  # Avoid double counting
                    total += graph[i][j].get("weight", 1.0)

        return total

    @staticmethod
    def external_contacts(
        graph: nx.Graph,
        residue_indices: List[int],
    ) -> float:
        """Count weighted contacts to residues outside the set."""
        total = 0.0
        indices_set = set(residue_indices)

        for i in residue_indices:
            for j in graph.neighbors(i):
                if j not in indices_set:
                    total += graph[i][j].get("weight", 1.0)

        return total

    @staticmethod
    def contact_density_ratio(
        graph: nx.Graph,
        residue_indices: List[int],
    ) -> float:
        """
        Calculate ratio of internal to external contacts.

        Higher values indicate a well-defined domain.
        """
        internal = ContactDensityCalculator.internal_contacts(graph, residue_indices)
        external = ContactDensityCalculator.external_contacts(graph, residue_indices)

        if external == 0:
            return float("inf") if internal > 0 else 0.0

        return internal / external

    @staticmethod
    def normalized_contact_density(
        graph: nx.Graph,
        residue_indices: List[int],
    ) -> float:
        """
        Calculate normalized contact density.

        contacts / (n * (n-1) / 2)

        This gives the fraction of possible contacts that exist.
        """
        n = len(residue_indices)
        if n < 2:
            return 0.0

        max_contacts = n * (n - 1) / 2
        internal = ContactDensityCalculator.internal_contacts(graph, residue_indices)

        return internal / max_contacts
