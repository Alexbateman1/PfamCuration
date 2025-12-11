"""
Spectral graph analysis for domain boundary detection.

The Fiedler vector (2nd eigenvector of the graph Laplacian) naturally
identifies the "weakest link" in a graph - where it's easiest to partition.
For proteins, this corresponds to domain boundaries.
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
from typing import List, Tuple, Optional
import networkx as nx


def compute_laplacian(graph: nx.Graph, normalized: bool = True) -> np.ndarray:
    """
    Compute the graph Laplacian matrix.

    Parameters:
    -----------
    graph : nx.Graph
        Contact graph with nodes as residue indices
    normalized : bool
        If True, compute normalized Laplacian (better for varying degree)

    Returns:
    --------
    L : np.ndarray
        Laplacian matrix (n x n)
    """
    n = graph.number_of_nodes()

    # Get adjacency matrix with weights
    A = nx.adjacency_matrix(graph, weight='weight').toarray().astype(float)

    # Degree matrix
    D = np.diag(A.sum(axis=1))

    if normalized:
        # Normalized Laplacian: L = I - D^(-1/2) A D^(-1/2)
        # More stable for graphs with varying node degrees
        D_inv_sqrt = np.diag(1.0 / np.sqrt(np.maximum(A.sum(axis=1), 1e-10)))
        L = np.eye(n) - D_inv_sqrt @ A @ D_inv_sqrt
    else:
        # Unnormalized Laplacian: L = D - A
        L = D - A

    return L


def compute_fiedler_vector(graph: nx.Graph, num_vectors: int = 5) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the Fiedler vector and additional eigenvectors.

    The Fiedler vector is the eigenvector corresponding to the second-smallest
    eigenvalue of the Laplacian. It provides the optimal 2-way partition.
    Additional eigenvectors provide finer partitioning information.

    Parameters:
    -----------
    graph : nx.Graph
        Contact graph
    num_vectors : int
        Number of eigenvectors to compute (including the trivial first one)

    Returns:
    --------
    eigenvalues : np.ndarray
        First k eigenvalues (ascending order)
    eigenvectors : np.ndarray
        Corresponding eigenvectors (n x k), columns are eigenvectors
    """
    L = compute_laplacian(graph, normalized=True)

    # Compute smallest eigenvalues/vectors
    # The first eigenvalue is always 0 (for connected graphs)
    # The second (Fiedler) gives the optimal 2-partition
    try:
        eigenvalues, eigenvectors = eigsh(
            sparse.csr_matrix(L),
            k=min(num_vectors, L.shape[0] - 1),
            which='SM',  # Smallest magnitude
            tol=1e-6
        )
        # Sort by eigenvalue
        idx = np.argsort(eigenvalues)
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
    except Exception:
        # Fallback to dense eigendecomposition for small graphs
        eigenvalues, eigenvectors = np.linalg.eigh(L)
        eigenvalues = eigenvalues[:num_vectors]
        eigenvectors = eigenvectors[:, :num_vectors]

    return eigenvalues, eigenvectors


def find_split_candidates_from_fiedler(
    fiedler: np.ndarray,
    residue_indices: List[int],
    min_domain_size: int = 30
) -> List[int]:
    """
    Find candidate domain boundaries from the Fiedler vector.

    Split candidates are positions where:
    1. The Fiedler vector crosses zero (sign change)
    2. The Fiedler vector has large jumps (gradient peaks)

    Parameters:
    -----------
    fiedler : np.ndarray
        Fiedler vector values for each residue
    residue_indices : List[int]
        Residue numbers (1-based) corresponding to Fiedler values
    min_domain_size : int
        Minimum distance between split candidates

    Returns:
    --------
    candidates : List[int]
        Residue positions that are candidate split points
    """
    n = len(fiedler)
    candidates = []

    # 1. Find zero crossings
    for i in range(1, n):
        if fiedler[i-1] * fiedler[i] < 0:  # Sign change
            candidates.append(residue_indices[i])

    # 2. Find gradient peaks (large jumps)
    gradient = np.abs(np.diff(fiedler))
    threshold = np.mean(gradient) + 2 * np.std(gradient)

    for i in range(len(gradient)):
        if gradient[i] > threshold:
            pos = residue_indices[i + 1]
            if pos not in candidates:
                candidates.append(pos)

    # Filter candidates that are too close together
    candidates.sort()
    filtered = []
    last_pos = -min_domain_size
    for pos in candidates:
        if pos - last_pos >= min_domain_size:
            filtered.append(pos)
            last_pos = pos

    return filtered


def recursive_spectral_partition(
    graph: nx.Graph,
    residue_indices: List[int],
    min_domain_size: int = 30,
    max_depth: int = 5,
    conductance_threshold: float = 0.3
) -> List[List[int]]:
    """
    Recursively partition the graph using spectral bisection.

    At each step:
    1. Compute Fiedler vector
    2. Partition into two groups (positive vs negative Fiedler values)
    3. Check if partition is "good" (low conductance)
    4. If good, recursively partition each half
    5. Stop when domains are small enough or partition is poor

    Parameters:
    -----------
    graph : nx.Graph
        Contact graph
    residue_indices : List[int]
        Residue numbers for graph nodes
    min_domain_size : int
        Stop partitioning when below this size
    max_depth : int
        Maximum recursion depth
    conductance_threshold : float
        Maximum conductance for a "good" partition (lower = cleaner cut)

    Returns:
    --------
    partitions : List[List[int]]
        List of residue index lists, each representing a domain
    """
    if len(residue_indices) < 2 * min_domain_size or max_depth == 0:
        return [residue_indices]

    # Compute Fiedler vector
    eigenvalues, eigenvectors = compute_fiedler_vector(graph, num_vectors=2)

    # Check algebraic connectivity (2nd eigenvalue)
    # High value = well-connected, hard to partition
    if len(eigenvalues) < 2:
        return [residue_indices]

    algebraic_connectivity = eigenvalues[1]
    if algebraic_connectivity > 0.5:  # Graph is too well-connected to split
        return [residue_indices]

    fiedler = eigenvectors[:, 1]

    # Partition by Fiedler sign
    nodes = list(graph.nodes())
    node_to_idx = {n: i for i, n in enumerate(nodes)}

    group_a = []
    group_b = []

    for i, node in enumerate(nodes):
        if fiedler[i] < 0:
            group_a.append(node)
        else:
            group_b.append(node)

    # Check partition quality using conductance
    conductance = compute_partition_conductance(graph, group_a, group_b)

    if conductance > conductance_threshold:
        # Partition is poor - don't split
        return [residue_indices]

    if len(group_a) < min_domain_size or len(group_b) < min_domain_size:
        return [residue_indices]

    # Recursively partition each group
    subgraph_a = graph.subgraph(group_a).copy()
    subgraph_b = graph.subgraph(group_b).copy()

    result = []
    result.extend(recursive_spectral_partition(
        subgraph_a, sorted(group_a), min_domain_size, max_depth - 1, conductance_threshold
    ))
    result.extend(recursive_spectral_partition(
        subgraph_b, sorted(group_b), min_domain_size, max_depth - 1, conductance_threshold
    ))

    return result


def compute_partition_conductance(
    graph: nx.Graph,
    group_a: List[int],
    group_b: List[int]
) -> float:
    """
    Compute conductance of a partition.

    Conductance = cut_weight / min(vol(A), vol(B))

    Lower conductance = better partition (fewer edges crossing the cut
    relative to the smaller partition's total edges).

    Parameters:
    -----------
    graph : nx.Graph
        Full graph
    group_a, group_b : List[int]
        Node lists for each partition

    Returns:
    --------
    conductance : float
        Partition conductance (0 = perfect cut, 1 = terrible cut)
    """
    set_a = set(group_a)
    set_b = set(group_b)

    # Cut weight = sum of edge weights crossing the partition
    cut_weight = 0
    for u, v, data in graph.edges(data=True):
        if (u in set_a and v in set_b) or (u in set_b and v in set_a):
            cut_weight += data.get('weight', 1.0)

    # Volume = sum of weighted degrees
    vol_a = sum(graph.degree(n, weight='weight') for n in group_a)
    vol_b = sum(graph.degree(n, weight='weight') for n in group_b)

    min_vol = min(vol_a, vol_b)

    if min_vol == 0:
        return 1.0

    return cut_weight / min_vol


def get_spectral_embedding(graph: nx.Graph, dims: int = 3) -> np.ndarray:
    """
    Get spectral embedding of graph nodes.

    Each node is represented by its values in the first k non-trivial
    eigenvectors. Nodes in the same domain should cluster together.

    Parameters:
    -----------
    graph : nx.Graph
        Contact graph
    dims : int
        Number of dimensions (eigenvectors) to use

    Returns:
    --------
    embedding : np.ndarray
        (n_nodes, dims) array of spectral coordinates
    """
    eigenvalues, eigenvectors = compute_fiedler_vector(graph, num_vectors=dims + 1)

    # Skip first eigenvector (trivial), use next k
    return eigenvectors[:, 1:dims + 1]
