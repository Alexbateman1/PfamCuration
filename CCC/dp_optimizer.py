"""
Dynamic Programming optimizer for domain boundary placement.

Uses DP along the sequence to find the optimal domain assignment,
scoring each possible domain by its "domainness" (contact density,
compactness, etc.) using 3D structural information.
"""

import numpy as np
from typing import List, Tuple, Dict, Optional, Callable
from dataclasses import dataclass
import networkx as nx


@dataclass
class DomainAssignment:
    """Represents a domain assignment (segmentation of the sequence)."""
    boundaries: List[int]  # Domain boundaries (split points)
    domains: List[Tuple[int, int]]  # List of (start, end) for each domain
    score: float  # Total score of this assignment
    domain_scores: List[float]  # Individual domain scores


def dp_optimal_segmentation(
    seq_length: int,
    score_fn: Callable[[int, int], float],
    min_domain_size: int = 30,
    max_domains: Optional[int] = None,
    domain_penalty: float = 0.0
) -> DomainAssignment:
    """
    Find optimal sequence segmentation using dynamic programming.

    Uses the recurrence:
        opt[i] = max over j < i of: opt[j] + score(j+1, i) - penalty

    Parameters:
    -----------
    seq_length : int
        Length of the sequence
    score_fn : Callable[[int, int], float]
        Function that scores a potential domain from start to end (1-based)
        Higher score = better domain
    min_domain_size : int
        Minimum domain size
    max_domains : int, optional
        Maximum number of domains allowed
    domain_penalty : float
        Penalty for each additional domain (MDL-inspired)
        Higher penalty = fewer, larger domains

    Returns:
    --------
    assignment : DomainAssignment
        Optimal domain assignment
    """
    n = seq_length

    # dp[i] = best score for segmenting residues 1..i
    # parent[i] = where the last domain starts for optimal segmentation ending at i
    dp = np.full(n + 1, -np.inf)
    parent = np.zeros(n + 1, dtype=int)
    dp[0] = 0  # Empty prefix has score 0

    for i in range(min_domain_size, n + 1):
        # Try all possible last domain start positions
        for j in range(0, i - min_domain_size + 1):
            domain_start = j + 1  # 1-based
            domain_end = i  # 1-based

            domain_score = score_fn(domain_start, domain_end)
            total_score = dp[j] + domain_score - domain_penalty

            if total_score > dp[i]:
                dp[i] = total_score
                parent[i] = j

    # Backtrack to find optimal segmentation
    boundaries = []
    domains = []
    domain_scores = []

    pos = n
    while pos > 0:
        start_pos = parent[pos]
        domain_start = start_pos + 1
        domain_end = pos

        domains.append((domain_start, domain_end))
        domain_scores.append(score_fn(domain_start, domain_end))

        if start_pos > 0:
            boundaries.append(domain_start)

        pos = start_pos

    # Reverse to get correct order
    domains.reverse()
    domain_scores.reverse()
    boundaries.reverse()

    return DomainAssignment(
        boundaries=boundaries,
        domains=domains,
        score=dp[n],
        domain_scores=domain_scores
    )


def dp_with_candidates(
    seq_length: int,
    candidates: List[int],
    score_fn: Callable[[int, int], float],
    min_domain_size: int = 30,
    domain_penalty: float = 0.0
) -> DomainAssignment:
    """
    DP segmentation restricted to candidate split points.

    More efficient than full DP when we have good candidate boundaries
    (e.g., from spectral analysis).

    Parameters:
    -----------
    seq_length : int
        Sequence length
    candidates : List[int]
        Candidate split positions (residue numbers where domains can start)
    score_fn : Callable[[int, int], float]
        Domain scoring function
    min_domain_size : int
        Minimum domain size
    domain_penalty : float
        Per-domain penalty

    Returns:
    --------
    assignment : DomainAssignment
        Optimal assignment using only candidate boundaries
    """
    # Add sequence start and end to candidates
    all_points = sorted(set([1] + candidates + [seq_length + 1]))

    # Filter points that are too close
    filtered_points = [all_points[0]]
    for p in all_points[1:]:
        if p - filtered_points[-1] >= min_domain_size or p == seq_length + 1:
            filtered_points.append(p)

    n_points = len(filtered_points)

    # dp[i] = best score for segmenting from start to filtered_points[i]-1
    dp = np.full(n_points, -np.inf)
    parent = np.full(n_points, -1, dtype=int)
    dp[0] = 0

    for i in range(1, n_points):
        end_pos = filtered_points[i] - 1  # Domain ends here

        for j in range(i):
            start_pos = filtered_points[j]  # Domain starts here

            if end_pos - start_pos + 1 < min_domain_size:
                continue

            domain_score = score_fn(start_pos, end_pos)
            total_score = dp[j] + domain_score - domain_penalty

            if total_score > dp[i]:
                dp[i] = total_score
                parent[i] = j

    # Backtrack
    boundaries = []
    domains = []
    domain_scores = []

    idx = n_points - 1
    while idx > 0:
        prev_idx = parent[idx]
        if prev_idx < 0:
            break

        domain_start = filtered_points[prev_idx]
        domain_end = filtered_points[idx] - 1

        domains.append((domain_start, domain_end))
        domain_scores.append(score_fn(domain_start, domain_end))

        if prev_idx > 0:
            boundaries.append(domain_start)

        idx = prev_idx

    domains.reverse()
    domain_scores.reverse()
    boundaries.reverse()

    return DomainAssignment(
        boundaries=boundaries,
        domains=domains,
        score=dp[-1],
        domain_scores=domain_scores
    )


def multi_scale_dp(
    seq_length: int,
    score_fn: Callable[[int, int], float],
    min_domain_size: int = 30,
    penalties: List[float] = None
) -> List[DomainAssignment]:
    """
    Run DP at multiple penalty scales to explore domain granularity.

    Lower penalty = more domains (finer split)
    Higher penalty = fewer domains (coarser split)

    Returns assignments at each scale for user to choose or ensemble.

    Parameters:
    -----------
    seq_length : int
        Sequence length
    score_fn : Callable[[int, int], float]
        Domain scoring function
    min_domain_size : int
        Minimum domain size
    penalties : List[float], optional
        Penalty values to try (default: [0, 5, 10, 20, 50])

    Returns:
    --------
    assignments : List[DomainAssignment]
        One assignment per penalty value
    """
    if penalties is None:
        penalties = [0, 5, 10, 20, 50]

    assignments = []
    for penalty in penalties:
        assignment = dp_optimal_segmentation(
            seq_length, score_fn, min_domain_size, domain_penalty=penalty
        )
        assignments.append(assignment)

    return assignments


def select_best_assignment(
    assignments: List[DomainAssignment],
    seq_length: int,
    mdl_weight: float = 1.0
) -> DomainAssignment:
    """
    Select best assignment using MDL-inspired criterion.

    MDL score = domain_score - mdl_weight * n_domains * log(seq_length)

    Parameters:
    -----------
    assignments : List[DomainAssignment]
        Candidate assignments
    seq_length : int
        Sequence length
    mdl_weight : float
        Weight for domain count penalty

    Returns:
    --------
    best : DomainAssignment
        Best assignment by MDL criterion
    """
    best = None
    best_mdl = -np.inf

    for assignment in assignments:
        n_domains = len(assignment.domains)
        mdl_score = assignment.score - mdl_weight * n_domains * np.log(seq_length)

        if mdl_score > best_mdl:
            best_mdl = mdl_score
            best = assignment

    return best
