"""
CCC (Claude Code Chop) Domain Predictor

Main orchestration class that combines:
1. Spectral graph analysis to find candidate boundaries
2. Dynamic programming to optimize domain assignment
3. Structural scoring to evaluate domain quality
"""

import json
import logging
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import networkx as nx

from .spectral import (
    compute_fiedler_vector,
    find_split_candidates_from_fiedler,
    recursive_spectral_partition,
    get_spectral_embedding
)
from .dp_optimizer import (
    dp_optimal_segmentation,
    dp_with_candidates,
    multi_scale_dp,
    select_best_assignment,
    DomainAssignment
)
from .scoring import DomainScorer, create_score_function

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
    """Predicted domain."""
    domain_id: int
    start: int
    end: int
    size: int
    score: float
    metrics: Dict = field(default_factory=dict)

    def __str__(self):
        return f"{self.start}-{self.end}"


@dataclass
class CCCPrediction:
    """Complete prediction result."""
    uniprot_acc: str
    sequence_length: int
    domains: List[Domain]
    method_used: str
    spectral_info: Dict = field(default_factory=dict)

    def summary(self) -> str:
        lines = [
            f"CCC Prediction for {self.uniprot_acc}",
            f"Sequence length: {self.sequence_length}",
            f"Domains found: {len(self.domains)}",
            f"Method: {self.method_used}",
            "",
            "Domains:",
        ]
        for d in self.domains:
            lines.append(f"  {d.domain_id}: {d.start}-{d.end} (size={d.size}, score={d.score:.2f})")
        return "\n".join(lines)


class CCCPredictor:
    """
    Claude Code Chop Domain Predictor

    Uses spectral graph analysis + dynamic programming to identify
    protein domains from AlphaFold structures.
    """

    def __init__(
        self,
        distance_threshold: float = 10.0,
        min_domain_size: int = 30,
        domain_penalty: float = 2.0,  # Lower penalty to allow more domains
        use_pae: bool = True,
        cache_dir: Optional[str] = None,
        plddt_threshold: float = 70.0,  # Filter disordered regions
    ):
        self.distance_threshold = distance_threshold
        self.min_domain_size = min_domain_size
        self.domain_penalty = domain_penalty
        self.use_pae = use_pae
        self.plddt_threshold = plddt_threshold

        if cache_dir is None:
            # Try to use existing ABC cache first, then fall back to CCC cache
            abc_cache = Path.cwd() / ".abc_cache"
            if abc_cache.exists():
                self.cache_dir = abc_cache
            else:
                self.cache_dir = Path.cwd() / ".ccc_cache"
        else:
            self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def predict_from_uniprot(self, uniprot_acc: str) -> CCCPrediction:
        """
        Predict domains for a UniProt accession.

        Downloads AlphaFold structure if needed.
        """
        # Get structure files
        cif_path, pae_path = self._get_alphafold_files(uniprot_acc)

        # Parse structure
        residues = self._parse_structure(cif_path)
        logger.info(f"Parsed {len(residues)} residues")

        # Load PAE if available
        pae_matrix = None
        if self.use_pae and pae_path.exists():
            pae_matrix = self._load_pae(pae_path)
            logger.info(f"Loaded PAE matrix")

        # Build contact graph
        graph = self._build_contact_graph(residues)
        logger.info(f"Built contact graph: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")

        # Run prediction
        return self._predict(uniprot_acc, residues, graph, pae_matrix)

    def _predict(
        self,
        uniprot_acc: str,
        residues: List[Residue],
        graph: nx.Graph,
        pae_matrix: Optional[np.ndarray] = None
    ) -> CCCPrediction:
        """Core prediction logic."""

        seq_length = len(residues)
        residue_numbers = [r.resnum for r in residues]
        coords = np.array([r.ca_coord for r in residues])
        plddt = np.array([r.plddt for r in residues])

        # Create scoring function (uses full graph)
        score_fn = create_score_function(
            graph, coords, plddt, residue_numbers, pae_matrix
        )

        spectral_info = {}

        # Filter out disordered regions (low pLDDT) before spectral analysis
        # This prevents Fiedler from finding disorder/structure boundaries instead of domain boundaries
        structured_mask = plddt >= self.plddt_threshold
        structured_indices = np.where(structured_mask)[0]
        logger.info(f"Structured residues (pLDDT >= {self.plddt_threshold}): {len(structured_indices)}/{len(residues)}")

        if len(structured_indices) < self.min_domain_size:
            # Mostly disordered - return no domains
            logger.warning("Protein is mostly disordered, no domains predicted")
            return CCCPrediction(
                uniprot_acc=uniprot_acc,
                sequence_length=residue_numbers[-1],
                domains=[],
                method_used="none",
                spectral_info={'reason': 'mostly_disordered'}
            )

        # Build filtered graph with only structured residues
        filtered_graph = graph.subgraph(structured_indices).copy()
        # Relabel nodes to be 0-indexed within the filtered graph
        old_to_new = {old: new for new, old in enumerate(sorted(structured_indices))}
        filtered_graph = nx.relabel_nodes(filtered_graph, old_to_new)
        filtered_residue_numbers = [residue_numbers[i] for i in structured_indices]

        # Method 1: Spectral partitioning for candidate boundaries (on filtered graph)
        logger.info("Computing spectral analysis on structured regions...")
        eigenvalues, eigenvectors = compute_fiedler_vector(filtered_graph, num_vectors=5)

        spectral_info['algebraic_connectivity'] = float(eigenvalues[1]) if len(eigenvalues) > 1 else 0
        spectral_info['num_eigenvectors'] = len(eigenvalues)

        fiedler = eigenvectors[:, 1] if eigenvectors.shape[1] > 1 else np.zeros(len(filtered_residue_numbers))

        # Find candidate split points from Fiedler vector (using filtered residue numbers)
        raw_candidates = find_split_candidates_from_fiedler(
            fiedler, filtered_residue_numbers, self.min_domain_size
        )
        logger.info(f"Found {len(raw_candidates)} raw candidate split points")

        # Filter candidates: only keep splits that create valid domains (contact_ratio >= 1)
        candidates = self._filter_candidates_by_contact_ratio(
            raw_candidates, residue_numbers, graph, min_ratio=1.0
        )
        logger.info(f"After contact-ratio filtering: {len(candidates)} candidates")
        spectral_info['raw_candidates'] = raw_candidates
        spectral_info['candidate_boundaries'] = candidates

        # Method 2: Skip recursive spectral partitioning (produces overlapping domains)
        # Just use the candidate-based methods below

        # Method 3: DP optimization with spectral candidates
        logger.info("Running DP optimization with spectral candidates...")
        if candidates:
            dp_result = dp_with_candidates(
                seq_length=residue_numbers[-1],
                candidates=candidates,
                score_fn=score_fn,
                min_domain_size=self.min_domain_size,
                domain_penalty=self.domain_penalty
            )
        else:
            # Fall back to full DP
            dp_result = dp_optimal_segmentation(
                seq_length=residue_numbers[-1],
                score_fn=score_fn,
                min_domain_size=self.min_domain_size,
                domain_penalty=self.domain_penalty
            )

        logger.info(f"DP found {len(dp_result.domains)} domains")

        # Method 4: Multi-scale DP with candidates (fast version)
        # Only try a few penalty values using candidate-restricted DP
        logger.info("Running multi-scale DP with candidates...")
        multi_results = []
        for penalty in [0, 5, 20, 50]:
            if candidates:
                result = dp_with_candidates(
                    seq_length=residue_numbers[-1],
                    candidates=candidates,
                    score_fn=score_fn,
                    min_domain_size=self.min_domain_size,
                    domain_penalty=penalty
                )
            else:
                result = dp_optimal_segmentation(
                    seq_length=residue_numbers[-1],
                    score_fn=score_fn,
                    min_domain_size=self.min_domain_size,
                    domain_penalty=penalty
                )
            multi_results.append(result)

        best_multi = select_best_assignment(multi_results, residue_numbers[-1])
        logger.info(f"Multi-scale DP selected {len(best_multi.domains)} domains")

        # Also try direct candidate-to-domain conversion
        candidate_domains = self._candidates_to_domains(candidates, residue_numbers[-1])
        logger.info(f"Direct candidate split: {len(candidate_domains)} domains")

        # Choose best result
        # Compare methods (excluding spectral - it produces overlapping domains that inflate scores)
        results = [
            ("dp_spectral", dp_result.domains),
            ("dp_multiscale", best_multi.domains),
            ("candidates", candidate_domains),
        ]

        # Score each result
        best_score = -np.inf
        best_method = "dp_spectral"
        best_domains = dp_result.domains

        scorer = DomainScorer(graph, coords, plddt, residue_numbers, pae_matrix)

        for method, domains in results:
            if not domains:
                continue
            total_score = scorer.score_assignment(domains)
            logger.info(f"  {method}: {len(domains)} domains, score={total_score:.2f}")

            if total_score > best_score:
                best_score = total_score
                best_method = method
                best_domains = domains

        # Build final prediction
        final_domains = []
        for i, (start, end) in enumerate(best_domains):
            metrics = scorer.get_metrics(start, end)
            final_domains.append(Domain(
                domain_id=i + 1,
                start=start,
                end=end,
                size=end - start + 1,
                score=scorer.score(start, end),
                metrics={
                    'contact_density': metrics.contact_density,
                    'contact_ratio': metrics.contact_ratio,
                    'radius_of_gyration': metrics.radius_of_gyration,
                    'compactness': metrics.compactness,
                    'avg_plddt': metrics.avg_plddt,
                }
            ))

        return CCCPrediction(
            uniprot_acc=uniprot_acc,
            sequence_length=residue_numbers[-1],
            domains=final_domains,
            method_used=best_method,
            spectral_info=spectral_info
        )

    def _partition_to_domains(
        self,
        partitions: List[List[int]],
        residue_numbers: List[int]
    ) -> List[Tuple[int, int]]:
        """Convert spectral partitions (index lists) to domain boundaries."""
        domains = []
        for part in partitions:
            if len(part) < self.min_domain_size:
                continue
            # Get residue numbers for this partition
            resnums = sorted([residue_numbers[i] for i in part])
            domains.append((resnums[0], resnums[-1]))

        # Remove overlapping domains - keep only non-nested ones
        domains = sorted(domains, key=lambda d: (d[0], -d[1]))  # Sort by start, then larger first
        non_overlapping = []
        for d in domains:
            # Check if this domain is contained within any existing domain
            is_nested = False
            for existing in non_overlapping:
                if d[0] >= existing[0] and d[1] <= existing[1]:
                    is_nested = True
                    break
            if not is_nested:
                non_overlapping.append(d)

        return sorted(non_overlapping)

    def _candidates_to_domains(
        self,
        candidates: List[int],
        seq_length: int
    ) -> List[Tuple[int, int]]:
        """Convert candidate split points directly to non-overlapping domains."""
        if not candidates:
            return [(1, seq_length)]

        # Add start and end
        boundaries = sorted(set([1] + candidates + [seq_length + 1]))

        domains = []
        for i in range(len(boundaries) - 1):
            start = boundaries[i]
            end = boundaries[i + 1] - 1
            if end - start + 1 >= self.min_domain_size:
                domains.append((start, end))

        return domains

    def _filter_candidates_by_contact_ratio(
        self,
        candidates: List[int],
        residue_numbers: List[int],
        graph: nx.Graph,
        min_ratio: float = 1.0
    ) -> List[int]:
        """
        Filter candidate split points by contact ratio.

        Only keep candidates where BOTH sides of the split would have
        contact_ratio >= min_ratio (more internal than external contacts).
        """
        if not candidates:
            return []

        seq_start = residue_numbers[0]
        seq_end = residue_numbers[-1]
        resnum_to_idx = {rn: i for i, rn in enumerate(residue_numbers)}

        def compute_contact_ratio(start_res: int, end_res: int) -> float:
            """Compute contact ratio for a region."""
            # Get indices for this region
            region_indices = set()
            for rn in range(start_res, end_res + 1):
                if rn in resnum_to_idx:
                    region_indices.add(resnum_to_idx[rn])

            if len(region_indices) < self.min_domain_size:
                return 0.0

            internal = 0
            external = 0
            for idx in region_indices:
                if idx not in graph:
                    continue
                for neighbor in graph.neighbors(idx):
                    weight = graph[idx][neighbor].get('weight', 1.0)
                    if neighbor in region_indices:
                        if neighbor > idx:
                            internal += weight
                    else:
                        external += weight

            return internal / max(external, 0.1)

        # Test each candidate
        valid_candidates = []
        sorted_candidates = sorted(candidates)

        for i, cand in enumerate(sorted_candidates):
            # Left region: from start (or previous candidate) to cand-1
            left_start = sorted_candidates[i-1] if i > 0 else seq_start
            left_end = cand - 1

            # Right region: from cand to end (or next candidate)
            right_start = cand
            right_end = sorted_candidates[i+1] - 1 if i < len(sorted_candidates) - 1 else seq_end

            left_ratio = compute_contact_ratio(left_start, left_end)
            right_ratio = compute_contact_ratio(right_start, right_end)

            if left_ratio >= min_ratio and right_ratio >= min_ratio:
                valid_candidates.append(cand)
                logger.debug(f"  Candidate {cand}: left_ratio={left_ratio:.2f}, right_ratio={right_ratio:.2f} - VALID")
            else:
                logger.debug(f"  Candidate {cand}: left_ratio={left_ratio:.2f}, right_ratio={right_ratio:.2f} - rejected")

        return valid_candidates

    def _get_alphafold_files(self, uniprot_acc: str) -> Tuple[Path, Path]:
        """Download or locate AlphaFold structure files."""
        cif_path = self.cache_dir / f"{uniprot_acc}.cif"
        pae_path = self.cache_dir / f"{uniprot_acc}_pae.json"

        if cif_path.exists():
            logger.info(f"Using cached structure for {uniprot_acc}")
            return cif_path, pae_path

        # Try downloading from AlphaFold DB
        logger.info(f"Downloading AlphaFold model for {uniprot_acc}")
        base_url = "https://alphafold.ebi.ac.uk/files"

        versions = ['v6', 'v4', 'v3', 'v2']
        for version in versions:
            try:
                # Download structure
                cif_url = f"{base_url}/AF-{uniprot_acc}-F1-model_{version}.cif"
                urllib.request.urlretrieve(cif_url, cif_path)
                logger.info(f"Downloaded structure ({version})")

                # Download PAE
                pae_url = f"{base_url}/AF-{uniprot_acc}-F1-predicted_aligned_error_{version}.json"
                try:
                    urllib.request.urlretrieve(pae_url, pae_path)
                    logger.info(f"Downloaded PAE ({version})")
                except Exception:
                    logger.warning(f"Could not download PAE, continuing without it")

                return cif_path, pae_path
            except Exception:
                continue

        raise RuntimeError(f"Failed to download structure for {uniprot_acc} (tried versions: {versions})")

    def _parse_structure(self, cif_path: Path) -> List[Residue]:
        """Parse mmCIF structure file."""
        residues = []

        with open(cif_path) as f:
            content = f.read()

        # Simple mmCIF parser for ATOM records
        in_atom_site = False
        columns = {}

        for line in content.split('\n'):
            if line.startswith('_atom_site.'):
                col_name = line.split('.')[1].strip()
                columns[col_name] = len(columns)
                in_atom_site = True
            elif in_atom_site and line.startswith('ATOM') or line.startswith('HETATM'):
                parts = line.split()
                if len(parts) > max(columns.values()):
                    atom_name = parts[columns.get('label_atom_id', 3)]
                    if atom_name == 'CA':
                        residues.append(Residue(
                            index=len(residues),
                            resnum=int(parts[columns.get('label_seq_id', 8)]),
                            resname=parts[columns.get('label_comp_id', 5)],
                            ca_coord=np.array([
                                float(parts[columns.get('Cartn_x', 10)]),
                                float(parts[columns.get('Cartn_y', 11)]),
                                float(parts[columns.get('Cartn_z', 12)])
                            ]),
                            plddt=float(parts[columns.get('B_iso_or_equiv', 14)]),
                            chain_id=parts[columns.get('label_asym_id', 6)]
                        ))
            elif in_atom_site and line.startswith('#'):
                in_atom_site = False

        return residues

    def _load_pae(self, pae_path: Path) -> np.ndarray:
        """Load PAE matrix from JSON."""
        with open(pae_path) as f:
            data = json.load(f)

        if isinstance(data, list) and 'predicted_aligned_error' in data[0]:
            return np.array(data[0]['predicted_aligned_error'])
        elif 'predicted_aligned_error' in data:
            return np.array(data['predicted_aligned_error'])
        else:
            raise ValueError(f"Unknown PAE format in {pae_path}")

    def _build_contact_graph(self, residues: List[Residue]) -> nx.Graph:
        """Build contact graph from residue coordinates."""
        n = len(residues)
        coords = np.array([r.ca_coord for r in residues])

        # Compute pairwise distances
        from scipy.spatial.distance import cdist
        distances = cdist(coords, coords)

        # Build graph
        graph = nx.Graph()
        graph.add_nodes_from(range(n))

        for i in range(n):
            for j in range(i + 1, n):
                if distances[i, j] <= self.distance_threshold:
                    # Weight by distance (closer = stronger)
                    weight = np.exp(-(distances[i, j] / 8.0) ** 2)
                    graph.add_edge(i, j, weight=weight, distance=distances[i, j])

        return graph
