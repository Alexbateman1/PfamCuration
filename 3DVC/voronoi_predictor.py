"""
3DVC Domain Predictor - Voronoi-based spatial partitioning

This module implements domain prediction using Voronoi tessellation in 3D space.
Seeds are placed and optimized such that residues near the same seed form domains.

Algorithm:
1. Filter low-confidence residues (pLDDT < 70)
2. For K = 1 to max_domains:
   a. Initialize K seeds (k-means++ on CA coordinates)
   b. Iterate until convergence:
      - Assign: each residue -> nearest seed
      - Update: move seeds to optimize objective
   c. Score this K configuration
3. Select best K (balancing fit vs parsimony)
4. Map Voronoi cells -> domain boundaries
"""

import json
import logging
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.spatial.distance import cdist

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


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
class Domain:
    """Predicted domain with segments and quality metrics."""
    domain_id: int
    segments: List[Tuple[int, int]]
    residue_indices: List[int]
    quality_metrics: Dict = field(default_factory=dict)

    @property
    def size(self) -> int:
        return len(self.residue_indices)

    @property
    def is_discontinuous(self) -> bool:
        return len(self.segments) > 1

    def to_chopping_string(self) -> str:
        return "_".join(f"{s}-{e}" for s, e in self.segments)

    def __repr__(self):
        return f"Domain({self.domain_id}: {self.to_chopping_string()}, size={self.size})"


@dataclass
class NDRRegion:
    """Non-Domain Region."""
    start: int
    end: int
    residue_indices: List[int]
    avg_plddt: float
    reason: str


@dataclass
class VoronoiPrediction:
    """Complete prediction result."""
    uniprot_acc: str
    sequence_length: int
    domains: List[Domain]
    ndr_regions: List[NDRRegion]
    parameters: Dict
    seed_positions: Optional[np.ndarray] = None  # Final seed positions

    def summary(self) -> str:
        lines = [
            f"3DVC Prediction for {self.uniprot_acc}",
            f"Sequence length: {self.sequence_length}",
            f"Domains found: {len(self.domains)}",
            f"NDR regions: {len(self.ndr_regions)}",
            "",
            "Domains:",
        ]
        for d in self.domains:
            lines.append(f"  {d.domain_id}: {d.to_chopping_string()} "
                        f"(size={d.size}, discontinuous={d.is_discontinuous})")
        return "\n".join(lines)

    def to_dict(self) -> Dict:
        return {
            "uniprot_acc": self.uniprot_acc,
            "sequence_length": self.sequence_length,
            "domains": [
                {
                    "domain_id": d.domain_id,
                    "segments": d.segments,
                    "size": d.size,
                    "chopping": d.to_chopping_string(),
                    "discontinuous": d.is_discontinuous,
                    "quality_metrics": d.quality_metrics,
                }
                for d in self.domains
            ],
            "ndr_regions": [
                {
                    "start": n.start,
                    "end": n.end,
                    "size": len(n.residue_indices),
                    "avg_plddt": n.avg_plddt,
                    "reason": n.reason,
                }
                for n in self.ndr_regions
            ],
            "parameters": self.parameters,
        }


class VoronoiPredictor:
    """
    3DVC (3D Voronoi Chop) Domain Predictor

    Uses Voronoi tessellation to partition protein structure into domains.
    Seeds are placed in 3D space and optimized to maximize domain quality.

    Parameters:
    -----------
    max_domains : int
        Maximum number of domains to consider (default: 10)
    min_domain_size : int
        Minimum residues for a valid domain (default: 30)
    ndr_plddt_cutoff : float
        pLDDT below which residues are filtered (default: 70)
    contact_threshold : float
        CA-CA distance for defining contacts (default: 10.0)
    lambda_boundary : float
        Weight for boundary contact penalty (default: 1.0)
    mu_crossing : float
        Weight for chain crossing penalty (default: 2.0)
    nu_domains : float
        Weight for domain count penalty (default: 50.0)
    max_iterations : int
        Maximum optimization iterations (default: 100)
    convergence_tol : float
        Convergence tolerance for seed movement (default: 0.1)
    cache_dir : str, optional
        Directory for caching AlphaFold files
    """

    def __init__(
        self,
        max_domains: int = 10,
        min_domain_size: int = 30,
        ndr_plddt_cutoff: float = 70.0,
        contact_threshold: float = 10.0,
        lambda_boundary: float = 2.0,
        mu_crossing: float = 5.0,
        nu_domains: float = 100.0,
        max_iterations: int = 100,
        convergence_tol: float = 0.1,
        cache_dir: Optional[str] = None,
    ):
        self.max_domains = max_domains
        self.min_domain_size = min_domain_size
        self.ndr_plddt_cutoff = ndr_plddt_cutoff
        self.contact_threshold = contact_threshold
        self.lambda_boundary = lambda_boundary
        self.mu_crossing = mu_crossing
        self.nu_domains = nu_domains
        self.max_iterations = max_iterations
        self.convergence_tol = convergence_tol

        if cache_dir is None:
            self.cache_dir = Path.cwd() / ".3dvc_cache"
        else:
            self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    @property
    def parameters(self) -> Dict:
        return {
            "max_domains": self.max_domains,
            "min_domain_size": self.min_domain_size,
            "ndr_plddt_cutoff": self.ndr_plddt_cutoff,
            "contact_threshold": self.contact_threshold,
            "lambda_boundary": self.lambda_boundary,
            "mu_crossing": self.mu_crossing,
            "nu_domains": self.nu_domains,
            "max_iterations": self.max_iterations,
            "convergence_tol": self.convergence_tol,
        }

    def predict_from_file(
        self,
        pdb_path: str,
        uniprot_acc: Optional[str] = None,
    ) -> VoronoiPrediction:
        """Predict domains from a PDB/CIF file."""
        path = Path(pdb_path)
        if uniprot_acc is None:
            name = path.stem
            if name.startswith("AF-"):
                parts = name.split("-")
                uniprot_acc = parts[1] if len(parts) > 1 else name
            else:
                uniprot_acc = name

        logger.info(f"Parsing structure from {pdb_path}")
        residues = self._parse_structure(pdb_path)
        logger.info(f"Parsed {len(residues)} residues")

        return self._predict(residues, uniprot_acc)

    def predict_from_uniprot(self, uniprot_acc: str) -> VoronoiPrediction:
        """Predict domains by downloading AlphaFold model."""
        cif_path = self._get_alphafold_file(uniprot_acc)
        return self.predict_from_file(str(cif_path), uniprot_acc)

    def _get_alphafold_file(self, uniprot_acc: str) -> Path:
        """Download AlphaFold structure if not cached."""
        cif_path = self.cache_dir / f"{uniprot_acc}.cif"

        if cif_path.exists():
            logger.info(f"Using cached structure for {uniprot_acc}")
            return cif_path

        logger.info(f"Downloading AlphaFold model for {uniprot_acc}")

        versions = ["v6", "v4", "v3", "v2"]
        for version in versions:
            url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_acc}-F1-model_{version}.cif"
            try:
                urllib.request.urlretrieve(url, cif_path)
                logger.info(f"Downloaded structure ({version})")
                return cif_path
            except Exception:
                continue

        raise RuntimeError(f"Failed to download structure for {uniprot_acc}")

    def _parse_structure(self, pdb_path: str) -> List[Residue]:
        """Parse PDB/CIF file and extract CA coordinates and pLDDT."""
        from Bio.PDB import MMCIFParser, PDBParser

        path = Path(pdb_path)
        if path.suffix.lower() in [".cif", ".mmcif"]:
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)

        structure = parser.get_structure("protein", pdb_path)

        residues = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] != " ":
                        continue
                    if "CA" not in residue:
                        continue

                    ca = residue["CA"]
                    residues.append(Residue(
                        index=len(residues),
                        resnum=residue.id[1],
                        resname=residue.resname,
                        ca_coord=np.array(ca.get_coord()),
                        plddt=ca.get_bfactor(),
                        chain_id=chain.id,
                    ))
            break

        return residues

    def _predict(self, residues: List[Residue], uniprot_acc: str) -> VoronoiPrediction:
        """Main prediction pipeline using Voronoi partitioning."""
        n_residues = len(residues)

        # Step 1: Filter low-pLDDT residues
        logger.info("Filtering low-confidence residues...")
        structured_indices = [
            i for i, r in enumerate(residues)
            if r.plddt >= self.ndr_plddt_cutoff
        ]
        ndr_indices = set(range(n_residues)) - set(structured_indices)

        logger.info(f"Structured residues: {len(structured_indices)} ({len(ndr_indices)} filtered)")

        if len(structured_indices) < self.min_domain_size:
            logger.warning("Too few structured residues for domain prediction")
            return VoronoiPrediction(
                uniprot_acc=uniprot_acc,
                sequence_length=n_residues,
                domains=[],
                ndr_regions=self._build_ndr_regions(residues, ndr_indices),
                parameters=self.parameters,
            )

        # Step 2: Build contact matrix for structured residues
        logger.info("Building contact matrix...")
        coords = np.array([residues[i].ca_coord for i in structured_indices])
        resnums = np.array([residues[i].resnum for i in structured_indices])

        dist_matrix = cdist(coords, coords)
        contacts = (dist_matrix < self.contact_threshold) & (dist_matrix > 0)

        # Step 3: Try different numbers of domains
        logger.info("Optimizing Voronoi partitions...")
        best_score = -np.inf
        best_k = 1
        best_seeds = None
        best_assignments = None

        # Determine max K to try based on protein size
        max_k = min(self.max_domains, len(structured_indices) // self.min_domain_size)
        max_k = max(1, max_k)

        for k in range(1, max_k + 1):
            logger.info(f"  Trying K={k} domains...")

            seeds, assignments, score = self._optimize_voronoi(
                coords, contacts, resnums, k
            )

            logger.info(f"    K={k}: score={score:.2f}")

            if score > best_score:
                best_score = score
                best_k = k
                best_seeds = seeds
                best_assignments = assignments

        logger.info(f"Best partition: K={best_k} with score={best_score:.2f}")

        # Step 4: Convert assignments to domains
        domains = self._build_domains(
            residues, structured_indices, best_assignments, coords
        )

        # Step 5: Post-process domains
        domains = self._merge_small_domains(domains, residues, contacts, structured_indices)
        domains = self._filter_domains(domains, residues, contacts, structured_indices)

        # Collect NDR indices from filtered domains
        domain_indices = set()
        for d in domains:
            domain_indices.update(d.residue_indices)
        ndr_indices = set(range(n_residues)) - domain_indices

        # Renumber domains
        for i, d in enumerate(domains):
            d.domain_id = i + 1

        # Build NDR regions
        ndr_regions = self._build_ndr_regions(residues, ndr_indices)

        logger.info(f"Final: {len(domains)} domains, {len(ndr_regions)} NDR regions")

        return VoronoiPrediction(
            uniprot_acc=uniprot_acc,
            sequence_length=n_residues,
            domains=domains,
            ndr_regions=ndr_regions,
            parameters=self.parameters,
            seed_positions=best_seeds,
        )

    def _kmeans_plusplus_init(self, coords: np.ndarray, k: int) -> np.ndarray:
        """Initialize seeds using k-means++ algorithm."""
        n = len(coords)
        seeds = np.zeros((k, 3))

        # First seed: random point
        seeds[0] = coords[np.random.randint(n)]

        for i in range(1, k):
            # Calculate distances to nearest existing seed
            distances = np.min(cdist(coords, seeds[:i]), axis=1)

            # Sample proportional to distance squared
            probs = distances ** 2
            probs /= probs.sum()

            idx = np.random.choice(n, p=probs)
            seeds[i] = coords[idx]

        return seeds

    def _optimize_voronoi(
        self,
        coords: np.ndarray,
        contacts: np.ndarray,
        resnums: np.ndarray,
        k: int,
    ) -> Tuple[np.ndarray, np.ndarray, float]:
        """
        Optimize Voronoi partition for K domains.

        Returns (seeds, assignments, score)
        """
        n = len(coords)

        # Handle single domain case
        if k == 1:
            seeds = np.array([coords.mean(axis=0)])
            assignments = np.zeros(n, dtype=int)
            score = self._compute_score(assignments, contacts, resnums, k, coords)
            return seeds, assignments, score

        # Initialize with k-means++
        best_seeds = None
        best_assignments = None
        best_score = -np.inf

        # Run multiple restarts
        n_restarts = 5
        for restart in range(n_restarts):
            seeds = self._kmeans_plusplus_init(coords, k)

            # Optimize
            for iteration in range(self.max_iterations):
                # Assignment step: each point to nearest seed
                distances = cdist(coords, seeds)
                assignments = np.argmin(distances, axis=1)

                # Check for empty clusters
                unique_assignments = np.unique(assignments)
                if len(unique_assignments) < k:
                    # Reinitialize empty clusters
                    for c in range(k):
                        if c not in unique_assignments:
                            # Place seed at furthest point from all seeds
                            min_distances = np.min(distances, axis=1)
                            seeds[c] = coords[np.argmax(min_distances)]
                    continue

                # Update step: move seeds to optimize objective
                new_seeds = np.zeros_like(seeds)
                for c in range(k):
                    mask = assignments == c
                    if mask.sum() > 0:
                        # Move toward centroid of assigned points
                        new_seeds[c] = coords[mask].mean(axis=0)
                    else:
                        new_seeds[c] = seeds[c]

                # Apply gradient from objective function
                new_seeds = self._apply_objective_gradient(
                    new_seeds, coords, assignments, contacts, resnums
                )

                # Check convergence
                movement = np.linalg.norm(new_seeds - seeds)
                seeds = new_seeds

                if movement < self.convergence_tol:
                    break

            # Final assignment
            distances = cdist(coords, seeds)
            assignments = np.argmin(distances, axis=1)

            # Score this partition
            score = self._compute_score(assignments, contacts, resnums, k, coords)

            if score > best_score:
                best_score = score
                best_seeds = seeds.copy()
                best_assignments = assignments.copy()

        return best_seeds, best_assignments, best_score

    def _apply_objective_gradient(
        self,
        seeds: np.ndarray,
        coords: np.ndarray,
        assignments: np.ndarray,
        contacts: np.ndarray,
        resnums: np.ndarray,
    ) -> np.ndarray:
        """
        Apply gradient from objective function to adjust seeds.

        Move seeds to increase internal contacts and reduce boundary contacts.
        """
        k = len(seeds)
        new_seeds = seeds.copy()

        for c in range(k):
            mask = assignments == c
            if mask.sum() == 0:
                continue

            cluster_indices = np.where(mask)[0]
            cluster_coords = coords[mask]

            # Calculate "contact pull" - move toward points with many internal contacts
            contact_weights = np.zeros(len(cluster_indices))
            for i, idx in enumerate(cluster_indices):
                # Count contacts with other cluster members
                internal_contacts = contacts[idx, mask].sum()
                # Count contacts with other clusters (boundary)
                boundary_contacts = contacts[idx, ~mask].sum()

                # Weight: prefer points with high internal, low boundary contacts
                contact_weights[i] = internal_contacts - 0.5 * boundary_contacts

            # Normalize weights (shift to positive)
            contact_weights = contact_weights - contact_weights.min() + 1
            contact_weights /= contact_weights.sum()

            # Move seed toward weighted centroid
            weighted_centroid = np.average(cluster_coords, axis=0, weights=contact_weights)

            # Blend with simple centroid
            alpha = 0.3  # Learning rate
            new_seeds[c] = (1 - alpha) * seeds[c] + alpha * weighted_centroid

        return new_seeds

    def _compute_score(
        self,
        assignments: np.ndarray,
        contacts: np.ndarray,
        resnums: np.ndarray,
        k: int,
        coords: np.ndarray = None,
    ) -> float:
        """
        Compute objective function score for a partition.

        Score = internal_contacts - lambda * boundary_contacts
                - mu * chain_crossings - nu * K - rho * void_penalty

        The void penalty penalizes domains with large spatial extent but sparse contacts.
        """
        n = len(assignments)

        # Count internal contacts and compute void penalty for each cluster
        total_internal = 0
        void_penalty = 0
        for c in range(k):
            mask = assignments == c
            cluster_size = mask.sum()
            if cluster_size == 0:
                continue

            cluster_contacts = contacts[np.ix_(mask, mask)]
            internal = cluster_contacts.sum() / 2
            total_internal += internal

            # Void penalty based on spatial spread
            # If coords provided, penalize clusters where points are far from centroid
            if coords is not None and cluster_size > 1:
                cluster_coords = coords[mask]
                centroid = cluster_coords.mean(axis=0)
                distances = np.sqrt(np.sum((cluster_coords - centroid) ** 2, axis=1))
                avg_distance = distances.mean()
                # Penalty: large average distance = scattered cluster
                void_penalty += avg_distance * cluster_size * 0.5

        # Count boundary contacts (between different clusters)
        boundary_contacts = 0
        for i in range(n):
            for j in range(i + 1, n):
                if contacts[i, j] and assignments[i] != assignments[j]:
                    boundary_contacts += 1

        # Count chain crossings
        # A crossing occurs when residues i, j, k have i < j < k in sequence
        # but i and k are in same cluster while j is in different cluster
        chain_crossings = 0
        for c in range(k):
            mask = assignments == c
            cluster_resnums = sorted(resnums[mask])

            if len(cluster_resnums) < 2:
                continue

            # For each pair of residues in cluster, count different-cluster residues between them
            for i, r1 in enumerate(cluster_resnums[:-1]):
                for r2 in cluster_resnums[i+1:]:
                    # Count residues with resnum between r1 and r2 that are NOT in this cluster
                    between_mask = (resnums > r1) & (resnums < r2) & (~mask)
                    chain_crossings += between_mask.sum()

        # Normalize chain crossings by sequence length
        chain_crossings = chain_crossings / max(n, 1)

        # Compute final score
        # Maximize internal contacts, minimize boundary contacts, domain count, and void
        score = (
            total_internal
            - self.lambda_boundary * boundary_contacts
            - self.mu_crossing * chain_crossings
            - self.nu_domains * k
            - void_penalty
        )

        return score

    def _build_domains(
        self,
        residues: List[Residue],
        structured_indices: List[int],
        assignments: np.ndarray,
        coords: np.ndarray,
    ) -> List[Domain]:
        """Convert Voronoi assignments to Domain objects."""
        k = int(assignments.max()) + 1
        domains = []

        for c in range(k):
            mask = assignments == c
            if mask.sum() == 0:
                continue

            # Get original residue indices
            cluster_orig_indices = [structured_indices[i] for i in np.where(mask)[0]]

            # Convert to segments
            segments = self._indices_to_segments(residues, cluster_orig_indices)

            domains.append(Domain(
                domain_id=c,
                segments=segments,
                residue_indices=sorted(cluster_orig_indices),
            ))

        return domains

    def _indices_to_segments(
        self,
        residues: List[Residue],
        indices: List[int],
    ) -> List[Tuple[int, int]]:
        """Convert residue indices to contiguous segments."""
        if not indices:
            return []

        indices = sorted(indices)
        segments = []

        start_idx = indices[0]
        prev_idx = indices[0]

        for idx in indices[1:]:
            # Check for gap > 3 residues
            gap = residues[idx].resnum - residues[prev_idx].resnum
            if gap > 3:
                segments.append((
                    residues[start_idx].resnum,
                    residues[prev_idx].resnum,
                ))
                start_idx = idx
            prev_idx = idx

        segments.append((
            residues[start_idx].resnum,
            residues[prev_idx].resnum,
        ))

        return segments

    def _merge_small_domains(
        self,
        domains: List[Domain],
        residues: List[Residue],
        contacts: np.ndarray,
        structured_indices: List[int],
    ) -> List[Domain]:
        """Merge small domains into their most-connected neighbor."""
        if len(domains) <= 1:
            return domains

        # Create index mapping
        idx_to_structured = {orig: i for i, orig in enumerate(structured_indices)}

        changed = True
        while changed:
            changed = False
            domains_sorted = sorted(domains, key=lambda d: d.size)

            for i, domain in enumerate(domains_sorted):
                if domain.size >= self.min_domain_size:
                    continue

                # Find most connected neighbor domain
                best_neighbor = None
                best_contacts = 0

                for other in domains_sorted:
                    if other is domain:
                        continue

                    # Count contacts between domains
                    contact_count = 0
                    for idx1 in domain.residue_indices:
                        if idx1 not in idx_to_structured:
                            continue
                        si1 = idx_to_structured[idx1]
                        for idx2 in other.residue_indices:
                            if idx2 not in idx_to_structured:
                                continue
                            si2 = idx_to_structured[idx2]
                            if contacts[si1, si2]:
                                contact_count += 1

                    if contact_count > best_contacts:
                        best_contacts = contact_count
                        best_neighbor = other

                if best_neighbor is not None:
                    # Merge into neighbor
                    new_indices = sorted(set(
                        best_neighbor.residue_indices + domain.residue_indices
                    ))
                    best_neighbor.residue_indices = new_indices
                    best_neighbor.segments = self._indices_to_segments(residues, new_indices)
                    domains.remove(domain)
                    changed = True
                    break

        return domains

    def _filter_domains(
        self,
        domains: List[Domain],
        residues: List[Residue],
        contacts: np.ndarray,
        structured_indices: List[int],
    ) -> List[Domain]:
        """Filter domains by quality metrics."""
        idx_to_structured = {orig: i for i, orig in enumerate(structured_indices)}

        filtered = []
        for domain in domains:
            # Compute quality metrics
            metrics = self._compute_domain_metrics(
                domain, residues, contacts, idx_to_structured
            )
            domain.quality_metrics = metrics

            # Filter by size
            if domain.size < self.min_domain_size:
                logger.info(f"Filtering domain {domain.segments}: too small ({domain.size})")
                continue

            # Filter by contact ratio (internal/boundary)
            contact_ratio = metrics.get('contact_ratio', 0)
            if contact_ratio < 0.5:  # More boundary than internal contacts
                logger.info(f"Filtering domain {domain.segments}: low contact ratio ({contact_ratio:.2f})")
                continue

            filtered.append(domain)

        return filtered

    def _compute_domain_metrics(
        self,
        domain: Domain,
        residues: List[Residue],
        contacts: np.ndarray,
        idx_to_structured: Dict[int, int],
    ) -> Dict:
        """Compute quality metrics for a domain."""
        indices_set = set(domain.residue_indices)

        # Count contacts
        internal = 0
        boundary = 0

        for idx1 in domain.residue_indices:
            if idx1 not in idx_to_structured:
                continue
            si1 = idx_to_structured[idx1]

            for idx2, si2 in idx_to_structured.items():
                if si2 <= si1:
                    continue
                if not contacts[si1, si2]:
                    continue

                if idx2 in indices_set:
                    internal += 1
                else:
                    boundary += 1

        contact_ratio = internal / boundary if boundary > 0 else float('inf')

        # Average pLDDT
        avg_plddt = np.mean([residues[i].plddt for i in domain.residue_indices])

        # Radius of gyration
        coords = np.array([residues[i].ca_coord for i in domain.residue_indices])
        centroid = coords.mean(axis=0)
        rg = np.sqrt(np.mean(np.sum((coords - centroid) ** 2, axis=1)))

        return {
            'internal_contacts': internal,
            'boundary_contacts': boundary,
            'contact_ratio': contact_ratio,
            'avg_plddt': avg_plddt,
            'radius_of_gyration': rg,
        }

    def _build_ndr_regions(
        self,
        residues: List[Residue],
        ndr_indices: set,
    ) -> List[NDRRegion]:
        """Build NDR regions from indices."""
        if not ndr_indices:
            return []

        ndr_sorted = sorted(ndr_indices)
        regions = []

        current = [ndr_sorted[0]]
        for idx in ndr_sorted[1:]:
            if idx - current[-1] <= 5:
                current.append(idx)
            else:
                if len(current) >= 3:
                    avg_plddt = np.mean([residues[i].plddt for i in current])
                    reason = "low_plddt" if avg_plddt < self.ndr_plddt_cutoff else "unassigned"
                    regions.append(NDRRegion(
                        start=residues[current[0]].resnum,
                        end=residues[current[-1]].resnum,
                        residue_indices=current,
                        avg_plddt=avg_plddt,
                        reason=reason,
                    ))
                current = [idx]

        if len(current) >= 3:
            avg_plddt = np.mean([residues[i].plddt for i in current])
            reason = "low_plddt" if avg_plddt < self.ndr_plddt_cutoff else "unassigned"
            regions.append(NDRRegion(
                start=residues[current[0]].resnum,
                end=residues[current[-1]].resnum,
                residue_indices=current,
                avg_plddt=avg_plddt,
                reason=reason,
            ))

        return regions
