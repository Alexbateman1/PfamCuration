"""
ABC Domain Predictor - Main orchestration class

This module provides the main ABCPredictor class that coordinates
contact graph building, clustering, quality assessment, and NDR detection.
"""

import json
import logging
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from .contact_graph import ContactGraphBuilder
from .domain_quality import DomainQualityAssessor
from .ndr_detector import NDRDetector
from .visualize import DomainVisualizer

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class Residue:
    """Single residue with coordinates and metadata."""
    index: int  # 0-based index in the structure
    resnum: int  # Residue number (1-based, from PDB)
    resname: str  # 3-letter residue name
    ca_coord: np.ndarray  # Cα coordinates
    plddt: float  # pLDDT confidence score (0-100)
    chain_id: str = "A"


@dataclass
class Domain:
    """Predicted domain with segments and quality metrics."""
    domain_id: int
    segments: List[Tuple[int, int]]  # List of (start, end) residue numbers
    residue_indices: List[int]  # Indices of all residues in domain
    quality_metrics: Dict = field(default_factory=dict)

    @property
    def size(self) -> int:
        """Total number of residues in domain."""
        return len(self.residue_indices)

    @property
    def is_discontinuous(self) -> bool:
        """Check if domain has multiple segments."""
        return len(self.segments) > 1

    def to_chopping_string(self) -> str:
        """Convert to standard chopping format (e.g., '24-111_266-345')."""
        return "_".join(f"{s}-{e}" for s, e in self.segments)

    def __repr__(self):
        return f"Domain({self.domain_id}: {self.to_chopping_string()}, size={self.size})"


@dataclass
class NDRRegion:
    """Non-Domain Region (linker, disordered, etc.)."""
    start: int
    end: int
    residue_indices: List[int]
    avg_plddt: float
    reason: str  # "low_plddt", "isolated", "unassigned"

    @property
    def size(self) -> int:
        return len(self.residue_indices)


@dataclass
class ABCPrediction:
    """Complete prediction result for a protein."""
    uniprot_acc: str
    sequence_length: int
    domains: List[Domain]
    ndr_regions: List[NDRRegion]
    parameters: Dict

    def summary(self) -> str:
        """Generate human-readable summary."""
        lines = [
            f"ABC Prediction for {self.uniprot_acc}",
            f"Sequence length: {self.sequence_length}",
            f"Domains found: {len(self.domains)}",
            f"NDR regions: {len(self.ndr_regions)}",
            "",
            "Domains:",
        ]
        for d in self.domains:
            lines.append(f"  {d.domain_id}: {d.to_chopping_string()} "
                        f"(size={d.size}, discontinuous={d.is_discontinuous})")
            if d.quality_metrics:
                rg = d.quality_metrics.get('radius_of_gyration', 'N/A')
                cd = d.quality_metrics.get('contact_density_ratio', 'N/A')
                plddt = d.quality_metrics.get('avg_plddt', 'N/A')
                lines.append(f"     Rg={rg:.1f}Å, ContactRatio={cd:.2f}, pLDDT={plddt:.1f}")

        if self.ndr_regions:
            lines.append("\nNDR Regions:")
            for ndr in self.ndr_regions:
                lines.append(f"  {ndr.start}-{ndr.end} ({ndr.reason}, pLDDT={ndr.avg_plddt:.1f})")

        return "\n".join(lines)

    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization."""
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
                    "size": n.size,
                    "avg_plddt": n.avg_plddt,
                    "reason": n.reason,
                }
                for n in self.ndr_regions
            ],
            "parameters": self.parameters,
        }


class ABCPredictor:
    """
    ABC (Alex Bateman Chop) Domain Predictor

    Uses contact graph analysis and community detection to identify
    protein domains from AlphaFold structures.

    Parameters:
    -----------
    distance_threshold : float
        Maximum Cα-Cα distance for contact (default: 10Å)
    min_domain_size : int
        Minimum residues for a valid domain (default: 30)
    ndr_plddt_cutoff : float
        pLDDT below which residues may be NDR (default: 70)
    clustering_method : str
        'leiden' or 'louvain' (default: 'leiden')
    resolution : float
        Clustering resolution parameter (default: 1.0)
    sigma : float
        Gaussian decay parameter for edge weights (default: 8.0)
    use_pae : bool
        Whether to use PAE for edge weighting if available (default: True)
    """

    def __init__(
        self,
        distance_threshold: float = 10.0,
        min_domain_size: int = 30,
        ndr_plddt_cutoff: float = 70.0,
        clustering_method: str = "leiden",
        resolution: float = 1.0,
        sigma: float = 8.0,
        use_pae: bool = True,
    ):
        self.distance_threshold = distance_threshold
        self.min_domain_size = min_domain_size
        self.ndr_plddt_cutoff = ndr_plddt_cutoff
        self.clustering_method = clustering_method
        self.resolution = resolution
        self.sigma = sigma
        self.use_pae = use_pae

        # Component modules
        self.graph_builder = ContactGraphBuilder(
            distance_threshold=distance_threshold,
            sigma=sigma,
            use_pae=use_pae,
        )
        self.quality_assessor = DomainQualityAssessor()
        self.ndr_detector = NDRDetector(plddt_cutoff=ndr_plddt_cutoff)
        self.visualizer = DomainVisualizer()

    @property
    def parameters(self) -> Dict:
        """Current parameter settings."""
        return {
            "distance_threshold": self.distance_threshold,
            "min_domain_size": self.min_domain_size,
            "ndr_plddt_cutoff": self.ndr_plddt_cutoff,
            "clustering_method": self.clustering_method,
            "resolution": self.resolution,
            "sigma": self.sigma,
            "use_pae": self.use_pae,
        }

    def predict_from_file(
        self,
        pdb_path: str,
        pae_path: Optional[str] = None,
        uniprot_acc: Optional[str] = None,
    ) -> ABCPrediction:
        """
        Predict domains from a PDB/CIF file.

        Parameters:
        -----------
        pdb_path : str
            Path to PDB or mmCIF file
        pae_path : str, optional
            Path to PAE JSON file
        uniprot_acc : str, optional
            UniProt accession (extracted from filename if not provided)

        Returns:
        --------
        ABCPrediction
        """
        path = Path(pdb_path)
        if uniprot_acc is None:
            # Try to extract from filename (e.g., AF-P12345-F1-model_v4.cif)
            name = path.stem
            if name.startswith("AF-"):
                parts = name.split("-")
                uniprot_acc = parts[1] if len(parts) > 1 else name
            else:
                uniprot_acc = name

        # Parse structure
        logger.info(f"Parsing structure from {pdb_path}")
        residues = self._parse_structure(pdb_path)
        logger.info(f"Parsed {len(residues)} residues")

        # Load PAE if available
        pae_matrix = None
        if pae_path and Path(pae_path).exists():
            logger.info(f"Loading PAE from {pae_path}")
            pae_matrix = self._load_pae(pae_path)

        return self._predict(residues, pae_matrix, uniprot_acc)

    def predict_from_uniprot(self, uniprot_acc: str) -> ABCPrediction:
        """
        Predict domains by downloading AlphaFold model from EBI.

        Parameters:
        -----------
        uniprot_acc : str
            UniProt accession (e.g., 'P12345')

        Returns:
        --------
        ABCPrediction
        """
        import tempfile

        # Download structure - try multiple AlphaFold DB versions
        logger.info(f"Downloading AlphaFold model for {uniprot_acc}")

        with tempfile.TemporaryDirectory() as tmpdir:
            cif_path = Path(tmpdir) / f"{uniprot_acc}.cif"
            pae_path = Path(tmpdir) / f"{uniprot_acc}_pae.json"

            # Try different AlphaFold DB versions (newest first)
            versions = ["v6", "v4", "v3", "v2"]
            cif_downloaded = False

            for version in versions:
                cif_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_acc}-F1-model_{version}.cif"
                try:
                    urllib.request.urlretrieve(cif_url, cif_path)
                    logger.info(f"Downloaded structure ({version})")
                    cif_downloaded = True

                    # Try to get PAE for same version
                    pae_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_acc}-F1-predicted_aligned_error_{version}.json"
                    try:
                        urllib.request.urlretrieve(pae_url, pae_path)
                        logger.info(f"Downloaded PAE ({version})")
                    except Exception:
                        logger.warning(f"Could not download PAE for {uniprot_acc}, continuing without it")
                        pae_path = None

                    break
                except Exception:
                    continue

            if not cif_downloaded:
                raise RuntimeError(f"Failed to download structure for {uniprot_acc} (tried versions: {versions})")

            return self.predict_from_file(str(cif_path), str(pae_path) if pae_path else None, uniprot_acc)

    def _parse_structure(self, pdb_path: str) -> List[Residue]:
        """Parse PDB/CIF file and extract Cα coordinates and pLDDT scores."""
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
                for i, residue in enumerate(chain):
                    # Skip heteroatoms (water, ligands)
                    if residue.id[0] != " ":
                        continue

                    # Get Cα atom
                    if "CA" not in residue:
                        continue

                    ca = residue["CA"]
                    ca_coord = np.array(ca.get_coord())

                    # pLDDT is stored in B-factor for AlphaFold models
                    plddt = ca.get_bfactor()

                    residues.append(Residue(
                        index=len(residues),
                        resnum=residue.id[1],
                        resname=residue.resname,
                        ca_coord=ca_coord,
                        plddt=plddt,
                        chain_id=chain.id,
                    ))
            break  # Only process first model

        return residues

    def _load_pae(self, pae_path: str) -> Optional[np.ndarray]:
        """Load PAE matrix from AlphaFold JSON file."""
        try:
            with open(pae_path) as f:
                data = json.load(f)

            # AlphaFold v4 format
            if isinstance(data, list) and len(data) > 0:
                data = data[0]

            if "predicted_aligned_error" in data:
                return np.array(data["predicted_aligned_error"])
            elif "pae" in data:
                return np.array(data["pae"])
            else:
                logger.warning("Could not find PAE data in JSON file")
                return None
        except Exception as e:
            logger.warning(f"Error loading PAE: {e}")
            return None

    def _predict(
        self,
        residues: List[Residue],
        pae_matrix: Optional[np.ndarray],
        uniprot_acc: str,
    ) -> ABCPrediction:
        """Main prediction pipeline."""
        n_residues = len(residues)

        # Step 1: Build contact graph
        logger.info("Building contact graph...")
        graph = self.graph_builder.build(residues, pae_matrix)
        logger.info(f"Graph has {graph.number_of_nodes()} nodes and {graph.number_of_edges()} edges")

        # Step 2: Cluster the graph
        logger.info(f"Clustering with {self.clustering_method}...")
        cluster_assignments = self._cluster_graph(graph)

        # Count clusters
        unique_clusters = set(cluster_assignments.values())
        logger.info(f"Found {len(unique_clusters)} initial clusters")

        # Step 3: Identify NDRs (low pLDDT, isolated residues)
        logger.info("Identifying NDR regions...")
        ndr_residue_indices = self.ndr_detector.detect(residues, cluster_assignments, graph)

        # Step 4: Build initial domains from clusters
        raw_domains = self._build_domains_from_clusters(residues, cluster_assignments, ndr_residue_indices)

        # Step 5: Evaluate quality and filter
        logger.info("Assessing domain quality...")
        domains = []
        for domain in raw_domains:
            metrics = self.quality_assessor.assess(domain, residues, graph)
            domain.quality_metrics = metrics

            # Filter by minimum size
            if domain.size >= self.min_domain_size:
                domains.append(domain)
            else:
                # Add small clusters to NDR
                ndr_residue_indices.update(domain.residue_indices)

        # Step 6: Refinement - merge over-split domains
        logger.info("Refining domain boundaries...")
        domains = self._refine_domains(domains, residues, graph)

        # Renumber domains
        for i, domain in enumerate(domains):
            domain.domain_id = i + 1

        # Step 7: Build NDR regions from collected indices
        ndr_regions = self._build_ndr_regions(residues, ndr_residue_indices, domains)

        logger.info(f"Final result: {len(domains)} domains, {len(ndr_regions)} NDR regions")

        return ABCPrediction(
            uniprot_acc=uniprot_acc,
            sequence_length=n_residues,
            domains=domains,
            ndr_regions=ndr_regions,
            parameters=self.parameters,
        )

    def _cluster_graph(self, graph) -> Dict[int, int]:
        """Apply community detection to the graph."""
        import networkx as nx

        if graph.number_of_nodes() == 0:
            return {}

        # Convert to undirected for clustering
        if graph.is_directed():
            graph = graph.to_undirected()

        if self.clustering_method == "leiden":
            try:
                import leidenalg as la
                import igraph as ig

                # Convert NetworkX to igraph
                mapping = {n: i for i, n in enumerate(graph.nodes())}
                reverse_mapping = {i: n for n, i in mapping.items()}

                edges = [(mapping[u], mapping[v]) for u, v in graph.edges()]
                weights = [graph[u][v].get("weight", 1.0) for u, v in graph.edges()]

                ig_graph = ig.Graph(n=len(mapping), edges=edges, directed=False)
                ig_graph.es["weight"] = weights

                # Run Leiden with Modularity (more robust than CPM for variable weights)
                # Use RBConfigurationVertexPartition for resolution parameter support
                partition = la.find_partition(
                    ig_graph,
                    la.RBConfigurationVertexPartition,
                    weights="weight",
                    resolution_parameter=self.resolution,
                )

                return {reverse_mapping[i]: partition.membership[i] for i in range(len(mapping))}

            except ImportError:
                logger.warning("leidenalg not installed, falling back to Louvain")
                self.clustering_method = "louvain"

        if self.clustering_method == "louvain":
            try:
                from networkx.algorithms.community import louvain_communities

                communities = louvain_communities(
                    graph,
                    weight="weight",
                    resolution=self.resolution,
                )

                assignments = {}
                for i, community in enumerate(communities):
                    for node in community:
                        assignments[node] = i
                return assignments

            except ImportError:
                logger.warning("Louvain not available, using greedy modularity")
                from networkx.algorithms.community import greedy_modularity_communities

                communities = greedy_modularity_communities(graph, weight="weight")
                assignments = {}
                for i, community in enumerate(communities):
                    for node in community:
                        assignments[node] = i
                return assignments

        raise ValueError(f"Unknown clustering method: {self.clustering_method}")

    def _build_domains_from_clusters(
        self,
        residues: List[Residue],
        cluster_assignments: Dict[int, int],
        ndr_indices: set,
    ) -> List[Domain]:
        """Convert cluster assignments to Domain objects."""
        # Group residues by cluster
        clusters = {}
        for idx, cluster_id in cluster_assignments.items():
            if idx in ndr_indices:
                continue
            if cluster_id not in clusters:
                clusters[cluster_id] = []
            clusters[cluster_id].append(idx)

        domains = []
        for cluster_id, indices in clusters.items():
            indices = sorted(indices)
            if not indices:
                continue

            # Convert indices to residue numbers and find segments
            segments = self._indices_to_segments(residues, indices)

            domains.append(Domain(
                domain_id=cluster_id,
                segments=segments,
                residue_indices=indices,
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
            # Check if contiguous (allowing small gaps of 1-2 residues)
            gap = idx - prev_idx
            if gap > 3:  # New segment if gap > 3 residues
                # End previous segment
                segments.append((
                    residues[start_idx].resnum,
                    residues[prev_idx].resnum,
                ))
                start_idx = idx
            prev_idx = idx

        # Add final segment
        segments.append((
            residues[start_idx].resnum,
            residues[prev_idx].resnum,
        ))

        return segments

    def _refine_domains(
        self,
        domains: List[Domain],
        residues: List[Residue],
        graph,
    ) -> List[Domain]:
        """Refine domain boundaries by merging over-split domains."""
        if len(domains) <= 1:
            return domains

        # Calculate inter-domain contact density
        n_domains = len(domains)
        contact_matrix = np.zeros((n_domains, n_domains))

        for i, dom_i in enumerate(domains):
            for j, dom_j in enumerate(domains):
                if i >= j:
                    continue

                # Count contacts between domains
                contacts = 0
                for idx_i in dom_i.residue_indices:
                    for idx_j in dom_j.residue_indices:
                        if graph.has_edge(idx_i, idx_j):
                            contacts += graph[idx_i][idx_j].get("weight", 1.0)

                # Normalize by number of boundary residues
                boundary_size = min(len(dom_i.residue_indices), len(dom_j.residue_indices))
                if boundary_size > 0:
                    contact_matrix[i, j] = contacts / boundary_size
                    contact_matrix[j, i] = contact_matrix[i, j]

        # Merge domains with high inter-domain contact density
        merge_threshold = 0.5  # Threshold for merging
        merged = set()
        final_domains = []

        for i, domain in enumerate(domains):
            if i in merged:
                continue

            # Find domains to merge with this one
            to_merge = [i]
            for j in range(i + 1, n_domains):
                if j in merged:
                    continue
                if contact_matrix[i, j] > merge_threshold:
                    to_merge.append(j)
                    merged.add(j)

            if len(to_merge) == 1:
                final_domains.append(domain)
            else:
                # Merge domains
                all_indices = []
                for idx in to_merge:
                    all_indices.extend(domains[idx].residue_indices)
                all_indices = sorted(set(all_indices))

                segments = self._indices_to_segments(residues, all_indices)
                merged_domain = Domain(
                    domain_id=domain.domain_id,
                    segments=segments,
                    residue_indices=all_indices,
                )

                # Recalculate quality metrics
                metrics = self.quality_assessor.assess(merged_domain, residues, graph)
                merged_domain.quality_metrics = metrics

                final_domains.append(merged_domain)

        return final_domains

    def _build_ndr_regions(
        self,
        residues: List[Residue],
        ndr_indices: set,
        domains: List[Domain],
    ) -> List[NDRRegion]:
        """Build NDR regions from collected indices."""
        # Get all indices covered by domains
        domain_indices = set()
        for domain in domains:
            domain_indices.update(domain.residue_indices)

        # All remaining indices are NDR
        all_ndr = set(range(len(residues))) - domain_indices
        all_ndr.update(ndr_indices)

        if not all_ndr:
            return []

        # Group into contiguous regions
        ndr_indices_sorted = sorted(all_ndr)
        regions = []

        start_idx = ndr_indices_sorted[0]
        current_region = [start_idx]

        for idx in ndr_indices_sorted[1:]:
            if idx - current_region[-1] <= 5:  # Allow gaps of up to 5 residues
                current_region.append(idx)
            else:
                # Save current region and start new one
                if len(current_region) >= 3:  # Minimum NDR size
                    avg_plddt = np.mean([residues[i].plddt for i in current_region])
                    reason = "low_plddt" if avg_plddt < self.ndr_plddt_cutoff else "unassigned"
                    regions.append(NDRRegion(
                        start=residues[current_region[0]].resnum,
                        end=residues[current_region[-1]].resnum,
                        residue_indices=current_region,
                        avg_plddt=avg_plddt,
                        reason=reason,
                    ))
                current_region = [idx]

        # Handle last region
        if len(current_region) >= 3:
            avg_plddt = np.mean([residues[i].plddt for i in current_region])
            reason = "low_plddt" if avg_plddt < self.ndr_plddt_cutoff else "unassigned"
            regions.append(NDRRegion(
                start=residues[current_region[0]].resnum,
                end=residues[current_region[-1]].resnum,
                residue_indices=current_region,
                avg_plddt=avg_plddt,
                reason=reason,
            ))

        return regions

    def visualize(
        self,
        prediction: ABCPrediction,
        output_path: Optional[str] = None,
        residues: Optional[List[Residue]] = None,
    ) -> str:
        """Generate visualization and ChimeraX commands."""
        return self.visualizer.visualize(prediction, output_path, residues)

    def multi_threshold_analysis(
        self,
        pdb_path: str,
        pae_path: Optional[str] = None,
        thresholds: List[float] = [8.0, 10.0, 12.0, 15.0],
    ) -> Dict[float, ABCPrediction]:
        """
        Run predictions at multiple distance thresholds to find stable domains.

        Parameters:
        -----------
        pdb_path : str
            Path to structure file
        pae_path : str, optional
            Path to PAE file
        thresholds : list
            Distance thresholds to test

        Returns:
        --------
        Dict mapping threshold to prediction
        """
        results = {}

        for threshold in thresholds:
            logger.info(f"Running with distance threshold = {threshold}Å")
            self.distance_threshold = threshold
            self.graph_builder.distance_threshold = threshold

            prediction = self.predict_from_file(pdb_path, pae_path)
            results[threshold] = prediction

            logger.info(f"  Threshold {threshold}Å: {len(prediction.domains)} domains")

        return results
