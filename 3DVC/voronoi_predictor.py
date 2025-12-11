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
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.spatial.distance import cdist

# Set up logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Add console handler if not already present
if not logger.handlers:
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    logger.addHandler(console_handler)


def setup_file_logging(log_file: str = "3dvc_predictions.log"):
    """
    Set up file logging for prediction results.

    Parameters:
    -----------
    log_file : str
        Path to log file (default: 3dvc_predictions.log)
    """
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s\t%(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    return file_handler


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
    boundary_pairs: Optional[List[Tuple[int, int]]] = None  # Domain pairs with shared boundaries

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

    def _clip_polygon_by_plane(self, vertices: List[np.ndarray],
                                plane_point: np.ndarray, plane_normal: np.ndarray) -> List[np.ndarray]:
        """
        Clip a polygon by a half-plane.
        Keep vertices on the side where dot(v - plane_point, plane_normal) >= 0
        """
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

            # Edge crosses plane
            if (d1 >= 0 and d2 < 0) or (d1 < 0 and d2 >= 0):
                # Compute intersection
                t = d1 / (d1 - d2)
                intersection = v1 + t * (v2 - v1)
                result.append(intersection)

        return result

    def _compute_voronoi_face(self, seeds: np.ndarray, i: int, j: int,
                              max_radius: float = 50.0) -> List[np.ndarray]:
        """
        Compute the polygon boundary between Voronoi cells i and j.
        Returns list of vertices forming the polygon.
        """
        # Bisector plane between i and j
        midpoint = (seeds[i] + seeds[j]) / 2
        normal = seeds[j] - seeds[i]
        normal = normal / np.linalg.norm(normal)

        # Create two orthogonal vectors in the bisector plane
        if abs(normal[0]) < 0.9:
            u = np.cross(normal, np.array([1.0, 0.0, 0.0]))
        else:
            u = np.cross(normal, np.array([0.0, 1.0, 0.0]))
        u = u / np.linalg.norm(u)
        v = np.cross(normal, u)

        # Initial polygon: large square on the bisector plane
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

            # Bisector plane between i and k
            mid_ik = (seeds[i] + seeds[k]) / 2
            norm_ik = seeds[k] - seeds[i]
            norm_ik = norm_ik / np.linalg.norm(norm_ik)

            # Clip polygon by this half-plane (keep side closer to i)
            vertices = self._clip_polygon_by_plane(vertices, mid_ik, -norm_ik)

            if len(vertices) < 3:
                break

        return vertices

    def chimerax_commands(self, structure_model: str = "#1", plane_radius: float = 50.0) -> str:
        """
        Generate ChimeraX commands to visualize Voronoi boundaries.

        Parameters:
        -----------
        structure_model : str
            ChimeraX model ID for the structure (default: "#1")
        plane_radius : float
            Max radius for boundary polygons in Angstroms (default: 50.0)

        Returns:
        --------
        String of ChimeraX commands
        """
        if self.seed_positions is None or len(self.seed_positions) < 2:
            return "# No Voronoi boundaries to display (single domain or no seeds)"

        lines = [
            f"# ChimeraX commands for 3DVC visualization of {self.uniprot_acc}",
            f"# {len(self.domains)} domains predicted",
            "",
            "# Color domains",
        ]

        # Domain colors
        colors = ["red", "blue", "green", "yellow", "orange", "purple",
                  "cyan", "magenta", "lime", "pink"]

        # Color each domain - use comma-separated residue ranges
        for i, domain in enumerate(self.domains):
            color = colors[i % len(colors)]
            # ChimeraX syntax: color #1:10-50,60-80 red
            ranges = ",".join([f"{s}-{e}" for s, e in domain.segments])
            lines.append(f"color {structure_model}:{ranges} {color}")

        lines.append("")
        lines.append("# Voronoi seed positions (as spheres)")

        # Add seed markers
        for i, seed in enumerate(self.seed_positions):
            color = colors[i % len(colors)]
            lines.append(f"shape sphere radius 4 center {seed[0]:.2f},{seed[1]:.2f},{seed[2]:.2f} "
                        f"color {color}")

        lines.append("")
        lines.append("# Voronoi boundary polygons (as triangles)")

        # Compute which seed pairs share a Voronoi face
        n_seeds = len(self.seed_positions)

        pairs_to_draw = []
        for i in range(n_seeds):
            for j in range(i + 1, n_seeds):
                midpoint = (self.seed_positions[i] + self.seed_positions[j]) / 2
                dist_ij = np.linalg.norm(self.seed_positions[i] - midpoint)

                # Check if any other seed is closer to midpoint
                is_neighbor = True
                for k in range(n_seeds):
                    if k == i or k == j:
                        continue
                    dist_k = np.linalg.norm(self.seed_positions[k] - midpoint)
                    if dist_k < dist_ij:
                        is_neighbor = False
                        break

                if is_neighbor:
                    pairs_to_draw.append((i, j))

        # For each Voronoi face, compute polygon and render as triangles
        for i, j in pairs_to_draw:
            vertices = self._compute_voronoi_face(self.seed_positions, i, j, plane_radius)

            if len(vertices) < 3:
                continue

            # Triangulate using fan from centroid
            centroid = np.mean(vertices, axis=0)

            for idx in range(len(vertices)):
                v1 = vertices[idx]
                v2 = vertices[(idx + 1) % len(vertices)]

                # ChimeraX shape triangle using point option for each vertex
                lines.append(
                    f"shape triangle "
                    f"point {centroid[0]:.2f},{centroid[1]:.2f},{centroid[2]:.2f} "
                    f"point {v1[0]:.2f},{v1[1]:.2f},{v1[2]:.2f} "
                    f"point {v2[0]:.2f},{v2[1]:.2f},{v2[2]:.2f} "
                    f"color 128,128,128,51"
                )

        lines.append("")
        lines.append("# Display settings")
        lines.append(f"cartoon {structure_model}")
        lines.append(f"hide {structure_model} atoms")

        return "\n".join(lines)

    def save_chimerax_script(self, output_path: str, **kwargs):
        """Save ChimeraX commands to a file."""
        commands = self.chimerax_commands(**kwargs)
        with open(output_path, 'w') as f:
            f.write(commands)
        return output_path

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
        lambda_boundary: float = 2.0,  # No longer used (modularity handles this)
        mu_crossing: float = 5.0,
        nu_domains: float = 20.0,  # Reduced - modularity-based scoring needs less penalty
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
        all_scores = {}  # Track all scores for logging

        # Determine max K to try based on protein size
        n_structured = len(structured_indices)
        max_k = min(self.max_domains, n_structured // self.min_domain_size)
        max_k = max(1, max_k)

        for k in range(1, max_k + 1):
            seeds, assignments, score = self._optimize_voronoi(
                coords, contacts, resnums, k
            )

            score_per_res = score / n_structured if n_structured > 0 else 0
            all_scores[k] = (score, score_per_res)
            logger.info(f"  K={k}: score={score:.2f} ({score_per_res:.2f}/res)")

            if score > best_score:
                best_score = score
                best_k = k
                best_seeds = seeds
                best_assignments = assignments

        # Log summary line for easy parsing
        scores_str = " ".join([f"K{k}={s:.1f}" for k, (s, _) in sorted(all_scores.items())])
        best_score_per_res = best_score / n_structured if n_structured > 0 else 0
        logger.info(f"SCORES\t{uniprot_acc}\tn_struct={n_structured}\tbest_K={best_k}\t"
                   f"best_score={best_score:.2f}\tscore_per_res={best_score_per_res:.2f}\t{scores_str}")

        # Check if best score meets minimum threshold (per-residue)
        # If score is too negative, the structure doesn't form coherent domains
        n_structured = len(structured_indices)
        score_per_residue = best_score / n_structured if n_structured > 0 else 0
        min_score_per_residue = -20.0  # Threshold for acceptable domain structure

        if score_per_residue < min_score_per_residue:
            logger.info(f"Score per residue ({score_per_residue:.2f}) below threshold "
                       f"({min_score_per_residue}), predicting 0 domains")
            logger.info(f"RESULT\t{uniprot_acc}\tn_domains=0\t"
                       f"score_per_res={score_per_residue:.2f}\tdomains=NONE (below threshold)")
            return VoronoiPrediction(
                uniprot_acc=uniprot_acc,
                sequence_length=n_residues,
                domains=[],
                ndr_regions=self._build_ndr_regions(residues, set(range(n_residues))),
                parameters=self.parameters,
            )

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

        # Log final result summary (easy to grep)
        domain_strs = [d.to_chopping_string() for d in domains]
        domains_output = "; ".join(domain_strs) if domain_strs else "NONE"
        logger.info(f"RESULT\t{uniprot_acc}\tn_domains={len(domains)}\t"
                   f"score_per_res={best_score_per_res:.2f}\tdomains={domains_output}")

        # Compute boundary pairs: which domains share contacts (for ChimeraX visualization)
        boundary_pairs = []
        if len(domains) > 1:
            # Build mapping from original residue index to domain index
            res_to_domain = {}
            for di, dom in enumerate(domains):
                for ri in dom.residue_indices:
                    res_to_domain[ri] = di

            # Check contacts between domains
            for di in range(len(domains)):
                for dj in range(di + 1, len(domains)):
                    has_contact = False
                    for ri in domains[di].residue_indices:
                        if has_contact:
                            break
                        if ri not in [structured_indices[k] for k in range(len(structured_indices))]:
                            continue
                        si = structured_indices.index(ri) if ri in structured_indices else -1
                        if si < 0:
                            continue
                        for rj in domains[dj].residue_indices:
                            if rj not in structured_indices:
                                continue
                            sj = structured_indices.index(rj) if rj in structured_indices else -1
                            if sj < 0:
                                continue
                            if contacts[si, sj]:
                                has_contact = True
                                break
                    if has_contact:
                        boundary_pairs.append((di, dj))

        # Compute seed positions from final domain centroids (not from optimization)
        # This ensures seeds match the actual domains after merging/filtering
        final_seeds = []
        for dom in domains:
            domain_coords = np.array([residues[ri].ca_coord for ri in dom.residue_indices])
            centroid = np.mean(domain_coords, axis=0)
            final_seeds.append(centroid)
        final_seeds = np.array(final_seeds) if final_seeds else None

        return VoronoiPrediction(
            uniprot_acc=uniprot_acc,
            sequence_length=n_residues,
            domains=domains,
            ndr_regions=ndr_regions,
            parameters=self.parameters,
            seed_positions=final_seeds,
            boundary_pairs=boundary_pairs,
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
        Compute objective function score for a partition using modularity.

        Modularity Q measures whether there are more internal contacts than
        expected by chance given the degree distribution. This handles
        tightly-packed domains better than raw contact counting.

        Score = modularity * scale - mu * chain_crossings - nu * K
        """
        n = len(assignments)

        # Compute modularity
        # Q = (1/2m) * sum_ij [A_ij - (k_i * k_j)/(2m)] * delta(c_i, c_j)
        degrees = contacts.sum(axis=1)  # k_i = number of contacts for each residue
        m = contacts.sum() / 2  # Total edges (each contact counted once)

        if m == 0:
            return -float('inf')

        modularity = 0.0
        for c in range(k):
            mask = assignments == c
            if mask.sum() == 0:
                continue

            # Sum of A_ij for pairs in same cluster
            cluster_contacts = contacts[np.ix_(mask, mask)]
            actual_internal = cluster_contacts.sum() / 2

            # Sum of (k_i * k_j) / (2m) for pairs in same cluster (expected by null model)
            cluster_degrees = degrees[mask]
            expected_internal = (cluster_degrees.sum() ** 2) / (4 * m)

            modularity += (actual_internal - expected_internal) / m

        # Scale modularity to be comparable to other terms (modularity is typically 0-1)
        modularity_score = modularity * 1000

        # Count chain crossings
        chain_crossings = 0
        for c in range(k):
            mask = assignments == c
            cluster_resnums = sorted(resnums[mask])

            if len(cluster_resnums) < 2:
                continue

            for i, r1 in enumerate(cluster_resnums[:-1]):
                for r2 in cluster_resnums[i+1:]:
                    between_mask = (resnums > r1) & (resnums < r2) & (~mask)
                    chain_crossings += between_mask.sum()

        chain_crossings = chain_crossings / max(n, 1)

        # Compute final score
        # Modularity handles internal vs boundary naturally, so no lambda_boundary term
        score = (
            modularity_score
            - self.mu_crossing * chain_crossings
            - self.nu_domains * k
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

            # Convert to segments and stitch over short NDR gaps only
            # Pass structured_indices so we don't bridge over other domains
            segments = self._indices_to_segments(residues, cluster_orig_indices)
            structured_set = set(structured_indices)
            segments = self._stitch_segments(segments, max_gap=30, structured_indices=structured_set)

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

    def _stitch_segments(
        self,
        segments: List[Tuple[int, int]],
        max_gap: int = 30,
        structured_indices: set = None,
    ) -> List[Tuple[int, int]]:
        """
        Stitch together segments separated by short NDR gaps.

        Only bridges gaps where the intervening residues are NOT structured
        (i.e., they are low-pLDDT NDR regions, not other domains).

        Parameters:
        -----------
        segments : list of (start, end) tuples
        max_gap : int
            Maximum gap size to bridge (default: 30 residues)
        structured_indices : set
            Set of residue indices that are structured (pLDDT >= threshold).
            If provided, only bridge gaps where gap residues are NOT in this set.

        Returns:
        --------
        List of merged segments
        """
        if len(segments) <= 1:
            return segments

        # Sort by start position
        segments = sorted(segments)
        merged = [segments[0]]

        for start, end in segments[1:]:
            prev_start, prev_end = merged[-1]
            gap = start - prev_end - 1

            if gap <= max_gap:
                # Check if gap contains structured residues (other domains)
                can_bridge = True
                if structured_indices is not None:
                    gap_residues = set(range(prev_end + 1, start))
                    if gap_residues & structured_indices:
                        # Gap contains residues from other domains - don't bridge
                        can_bridge = False

                if can_bridge:
                    # Merge: extend previous segment to include this one
                    merged[-1] = (prev_start, end)
                else:
                    merged.append((start, end))
            else:
                merged.append((start, end))

        return merged

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
