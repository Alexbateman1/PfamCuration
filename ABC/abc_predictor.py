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
    min_contact_ratio : float
        Minimum internal/external contact ratio for valid domain (default: 1.5)
        Domains with lower ratio are rejected as they have more external
        contacts than internal, suggesting they are not compact units.
    ndr_plddt_cutoff : float
        pLDDT below which residues may be NDR (default: 70)
    clustering_method : str
        'leiden' or 'louvain' (default: 'leiden')
    resolution : float
        Clustering resolution parameter (default: 0.5)
    sigma : float
        Gaussian decay parameter for edge weights (default: 8.0)
    use_pae : bool
        Whether to use PAE for edge weighting if available (default: True)
    min_bridge_seq_sep : int
        Minimum sequence separation for a contact to be considered a "bridge"
        for discontinuous domain detection (default: 50)
    min_bridge_count : int
        Minimum number of bridges needed between domains to merge them (default: 2)
    cache_dir : str, optional
        Directory to cache downloaded AlphaFold files (default: ./.abc_cache)
        Set to None to disable caching
    use_dssp : bool
        Whether to use DSSP hydrogen bond information (default: False)
        When enabled, long-range H-bonds (beta sheets) add/boost graph edges
    """

    def __init__(
        self,
        distance_threshold: float = 10.0,
        min_domain_size: int = 30,
        min_contact_ratio: float = 1.5,
        ndr_plddt_cutoff: float = 70.0,
        clustering_method: str = "leiden",
        resolution: float = 0.5,
        sigma: float = 8.0,
        use_pae: bool = True,
        min_bridge_seq_sep: int = 50,
        min_bridge_count: int = 2,
        cache_dir: Optional[str] = None,
        use_dssp: bool = False,
    ):
        self.distance_threshold = distance_threshold
        self.min_domain_size = min_domain_size
        self.min_contact_ratio = min_contact_ratio
        self.ndr_plddt_cutoff = ndr_plddt_cutoff
        self.clustering_method = clustering_method
        self.resolution = resolution
        self.sigma = sigma
        self.use_pae = use_pae
        self.min_bridge_seq_sep = min_bridge_seq_sep
        self.min_bridge_count = min_bridge_count
        self.use_dssp = use_dssp

        # Set up cache directory in current working directory
        if cache_dir is None:
            self.cache_dir = Path.cwd() / ".abc_cache"
        else:
            self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

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
            "min_contact_ratio": self.min_contact_ratio,
            "ndr_plddt_cutoff": self.ndr_plddt_cutoff,
            "clustering_method": self.clustering_method,
            "resolution": self.resolution,
            "sigma": self.sigma,
            "use_pae": self.use_pae,
            "min_bridge_seq_sep": self.min_bridge_seq_sep,
            "min_bridge_count": self.min_bridge_count,
            "use_dssp": self.use_dssp,
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

        # Get DSSP data if enabled
        dssp_data = None
        if self.use_dssp:
            from .dssp_parser import get_dssp_for_structure
            dssp_data = get_dssp_for_structure(pdb_path, cache_dir=self.cache_dir)

        return self._predict(residues, pae_matrix, uniprot_acc, dssp_data)

    def _get_alphafold_files(self, uniprot_acc: str) -> Tuple[Path, Optional[Path]]:
        """
        Get AlphaFold structure and PAE files, using cache if available.

        Parameters:
        -----------
        uniprot_acc : str
            UniProt accession

        Returns:
        --------
        Tuple[Path, Optional[Path]]
            Paths to CIF file and PAE file (PAE may be None)
        """
        cif_path = self.cache_dir / f"{uniprot_acc}.cif"
        pae_path = self.cache_dir / f"{uniprot_acc}_pae.json"

        # Check if already cached
        if cif_path.exists():
            logger.info(f"Using cached structure for {uniprot_acc}")
            return cif_path, pae_path if pae_path.exists() else None

        # Download from AlphaFold DB
        logger.info(f"Downloading AlphaFold model for {uniprot_acc}")

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

        return cif_path, pae_path if pae_path and pae_path.exists() else None

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
        cif_path, pae_path = self._get_alphafold_files(uniprot_acc)
        return self.predict_from_file(str(cif_path), str(pae_path) if pae_path else None, uniprot_acc)

    def predict_from_uniprot_debug(self, uniprot_acc: str, debug_region: Optional[str] = None) -> ABCPrediction:
        """
        Predict domains with debug output to understand classification decisions.

        Parameters:
        -----------
        uniprot_acc : str
            UniProt accession
        debug_region : str, optional
            Region to focus debug output on (e.g., '1014-1085')
        """
        # Parse debug region
        region_start, region_end = None, None
        if debug_region:
            parts = debug_region.split('-')
            if len(parts) == 2:
                region_start, region_end = int(parts[0]), int(parts[1])

        # Get structure files (using cache)
        cif_path, pae_path = self._get_alphafold_files(uniprot_acc)

        # Parse structure
        residues = self._parse_structure(str(cif_path))
        pae_matrix = self._load_pae(str(pae_path)) if pae_path else None

        # DEBUG: Show pLDDT for region of interest
        print("\n" + "="*60)
        print("DEBUG: pLDDT Analysis")
        print("="*60)

        if region_start and region_end:
            print(f"\nRegion of interest: {region_start}-{region_end}")
            print("\npLDDT values in region:")
            region_plddts = []
            for res in residues:
                if region_start <= res.resnum <= region_end:
                    region_plddts.append(res.plddt)
                    if res.resnum % 10 == 0 or res.resnum == region_start or res.resnum == region_end:
                        print(f"  Residue {res.resnum}: pLDDT = {res.plddt:.1f}")

            print(f"\nRegion statistics ({region_start}-{region_end}):")
            print(f"  Average pLDDT: {np.mean(region_plddts):.1f}")
            print(f"  Min pLDDT: {np.min(region_plddts):.1f}")
            print(f"  Max pLDDT: {np.max(region_plddts):.1f}")
            print(f"  Residues above 70: {sum(1 for p in region_plddts if p >= 70)}/{len(region_plddts)}")
            print(f"  Residues above 60: {sum(1 for p in region_plddts if p >= 60)}/{len(region_plddts)}")

            # Show context around the region
            print(f"\nContext around region:")
            for res in residues:
                if region_start - 20 <= res.resnum <= region_end + 20:
                    marker = " <--" if region_start <= res.resnum <= region_end else ""
                    if res.resnum % 5 == 0 or res.resnum in [region_start, region_end]:
                        print(f"  {res.resnum}: pLDDT={res.plddt:.1f}{marker}")

        # Build graph and cluster
        print("\n" + "="*60)
        print("DEBUG: Clustering Analysis")
        print("="*60)

        graph = self.graph_builder.build(residues, pae_matrix)
        cluster_assignments = self._cluster_graph(graph)

        # Show what clusters the region belongs to
        if region_start and region_end:
            print(f"\nCluster assignments in region {region_start}-{region_end}:")
            cluster_counts = {}
            for res in residues:
                if region_start <= res.resnum <= region_end:
                    idx = res.index
                    cluster = cluster_assignments.get(idx, -1)
                    cluster_counts[cluster] = cluster_counts.get(cluster, 0) + 1

            for cluster, count in sorted(cluster_counts.items(), key=lambda x: -x[1]):
                print(f"  Cluster {cluster}: {count} residues")

            # Show contacts within region
            print(f"\nContact analysis for region:")
            region_indices = [res.index for res in residues if region_start <= res.resnum <= region_end]
            internal_contacts = 0
            external_contacts = 0
            for idx in region_indices:
                if idx in graph:
                    for neighbor in graph.neighbors(idx):
                        if neighbor in region_indices:
                            internal_contacts += 1
                        else:
                            external_contacts += 1
            internal_contacts //= 2  # Each edge counted twice
            print(f"  Internal contacts: {internal_contacts}")
            print(f"  External contacts: {external_contacts}")

        # Run NDR detection and show which methods flag the region
        print("\n" + "="*60)
        print("DEBUG: NDR Detection Analysis")
        print("="*60)

        low_plddt_ndr = self.ndr_detector._detect_low_plddt(residues)
        isolated_ndr = self.ndr_detector._detect_isolated(residues, graph)
        terminal_ndr = self.ndr_detector._detect_terminal_disorder(residues)

        if region_start and region_end:
            region_indices = set(res.index for res in residues if region_start <= res.resnum <= region_end)

            low_plddt_in_region = region_indices & low_plddt_ndr
            isolated_in_region = region_indices & isolated_ndr
            terminal_in_region = region_indices & terminal_ndr

            print(f"\nNDR flags in region {region_start}-{region_end}:")
            print(f"  Flagged by low_plddt: {len(low_plddt_in_region)}/{len(region_indices)} residues")
            print(f"  Flagged by isolated: {len(isolated_in_region)}/{len(region_indices)} residues")
            print(f"  Flagged by terminal_disorder: {len(terminal_in_region)}/{len(region_indices)} residues")

            total_flagged = region_indices & (low_plddt_ndr | isolated_ndr | terminal_ndr)
            print(f"  Total flagged as NDR: {len(total_flagged)}/{len(region_indices)} residues")

        print("\n" + "="*60)
        print("Proceeding with normal prediction...")
        print("="*60 + "\n")

        # Now run normal prediction
        return self._predict(residues, pae_matrix, uniprot_acc)

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
        dssp_data=None,
    ) -> ABCPrediction:
        """Main prediction pipeline."""
        n_residues = len(residues)

        # Step 1: Build contact graph
        logger.info("Building contact graph...")
        graph = self.graph_builder.build(residues, pae_matrix, dssp_data)
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

        # Step 3b: Rescue structured regions with moderate pLDDT
        # Some domains have moderate pLDDT (50-70) but are still structured
        ndr_residue_indices = self._rescue_structured_regions(
            residues, cluster_assignments, graph, ndr_residue_indices
        )

        # Step 4: Build initial domains from clusters
        raw_domains = self._build_domains_from_clusters(residues, cluster_assignments, ndr_residue_indices)

        # Step 5: Evaluate quality and filter by SIZE only
        # (ContactRatio filtering happens AFTER merging, so repeat units can be combined first)
        logger.info("Assessing domain quality...")
        domains = []
        rejected_for_size = []
        for domain in raw_domains:
            metrics = self.quality_assessor.assess(domain, residues, graph)
            domain.quality_metrics = metrics

            # Filter by minimum size only at this stage
            if domain.size < self.min_domain_size:
                rejected_for_size.append({
                    'segments': domain.segments,
                    'size': domain.size,
                    'contact_ratio': metrics.get('contact_density_ratio', 0),
                    'plddt': metrics.get('avg_plddt', 0),
                })
                ndr_residue_indices.update(domain.residue_indices)
            else:
                domains.append(domain)

        if rejected_for_size:
            logger.info(f"Rejected {len(rejected_for_size)} domains for being too small:")
            for rej in rejected_for_size:
                logger.info(f"  {rej['segments']}: size={rej['size']}, CR={rej['contact_ratio']:.2f}, pLDDT={rej['plddt']:.1f}")

        # Step 6: Refinement - merge over-split domains
        logger.info("Refining domain boundaries...")
        domains = self._refine_domains(domains, residues, graph)

        # Step 6b: Resolve interdigitations and enforce continuity
        domains, extra_ndr = self._enforce_domain_continuity(domains, residues)
        ndr_residue_indices.update(extra_ndr)

        # Step 6c: NOW filter by ContactRatio (after merging is complete)
        # This allows repeat protein blades to be merged before CR assessment
        # Note: CR threshold is adjusted for larger domains since they naturally
        # have more surface area and thus more external contacts
        if self.min_contact_ratio > 0:
            filtered_domains = []
            rejected_for_cr = []
            for domain in domains:
                # Re-assess quality after merging
                metrics = self.quality_assessor.assess(domain, residues, graph)
                domain.quality_metrics = metrics
                contact_ratio = metrics.get('contact_density_ratio', 0)

                # Size-aware CR threshold: larger domains get a lower threshold
                # Rationale: a 200-residue domain has ~2x the surface area of a 50-residue domain
                # so it naturally has more external contacts and lower CR
                if domain.size <= 50:
                    effective_cr_threshold = self.min_contact_ratio
                elif domain.size <= 100:
                    effective_cr_threshold = self.min_contact_ratio * 0.8  # 20% lower
                elif domain.size <= 200:
                    effective_cr_threshold = self.min_contact_ratio * 0.7  # 30% lower
                else:
                    effective_cr_threshold = self.min_contact_ratio * 0.6  # 40% lower for large domains

                if contact_ratio < effective_cr_threshold and contact_ratio != float('inf'):
                    rejected_for_cr.append({
                        'segments': domain.segments,
                        'size': domain.size,
                        'contact_ratio': contact_ratio,
                        'plddt': metrics.get('avg_plddt', 0),
                    })
                    logger.info(
                        f"Rejecting domain {domain.segments}: ContactRatio={contact_ratio:.2f} "
                        f"< {effective_cr_threshold:.2f} (size={domain.size}, pLDDT={metrics.get('avg_plddt', 0):.1f})"
                    )
                    ndr_residue_indices.update(domain.residue_indices)
                else:
                    filtered_domains.append(domain)

            if rejected_for_cr:
                logger.info(f"Rejected {len(rejected_for_cr)} domains for low ContactRatio:")
                for rej in rejected_for_cr:
                    logger.info(f"  {rej['segments']}: CR={rej['contact_ratio']:.2f}, size={rej['size']}, pLDDT={rej['plddt']:.1f}")

            domains = filtered_domains

        # Step 6d: Filter by compactness (Rg vs expected for globular domain)
        # Extended structures like alpha helices have high Rg relative to their size
        # For globular proteins: Rg ≈ 2.5 * N^(1/3) Å
        # Reject domains where actual Rg > 2.0 * expected Rg
        filtered_domains = []
        rejected_for_rg = []
        for domain in domains:
            metrics = domain.quality_metrics
            rg = metrics.get('radius_of_gyration', 0)

            # Expected Rg for a globular domain of this size
            expected_rg = 2.5 * (domain.size ** (1/3))
            rg_ratio = rg / expected_rg if expected_rg > 0 else 0

            # Reject if Rg is more than 2x expected (too extended)
            if rg_ratio > 2.0:
                rejected_for_rg.append({
                    'segments': domain.segments,
                    'size': domain.size,
                    'rg': rg,
                    'expected_rg': expected_rg,
                    'rg_ratio': rg_ratio,
                    'plddt': metrics.get('avg_plddt', 0),
                })
                logger.info(
                    f"Rejecting domain {domain.segments}: too extended "
                    f"(Rg={rg:.1f}Å, expected={expected_rg:.1f}Å, ratio={rg_ratio:.1f}x)"
                )
                ndr_residue_indices.update(domain.residue_indices)
            else:
                filtered_domains.append(domain)

        if rejected_for_rg:
            logger.info(f"Rejected {len(rejected_for_rg)} domains for being too extended (high Rg):")
            for rej in rejected_for_rg:
                logger.info(f"  {rej['segments']}: Rg={rej['rg']:.1f}Å (expected {rej['expected_rg']:.1f}Å, {rej['rg_ratio']:.1f}x)")

        domains = filtered_domains

        # Step 6e: Extend domains to absorb adjacent high-pLDDT unassigned regions
        # This rescues structured regions that were clustered separately but should
        # be part of an adjacent domain (common in repeat proteins like WD40)
        domains = self._extend_domains_to_structured_regions(
            domains, residues, ndr_residue_indices, graph
        )

        # Step 6f: If no domains remain, create domains from high-pLDDT regions
        # This handles cases where clustering fragments the structure and all
        # assembled domains get rejected, but there's clearly structured content
        if not domains:
            domains = self._create_domains_from_high_plddt_regions(
                residues, ndr_residue_indices, graph
            )

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
                # Set seed for reproducibility (Leiden is stochastic by default)
                partition = la.find_partition(
                    ig_graph,
                    la.RBConfigurationVertexPartition,
                    weights="weight",
                    resolution_parameter=self.resolution,
                    seed=42,
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
                    seed=42,
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

    def _rescue_structured_regions(
        self,
        residues: List[Residue],
        cluster_assignments: Dict[int, int],
        graph,
        ndr_indices: set,
    ) -> set:
        """
        Rescue structured regions that were marked as NDR due to moderate pLDDT.

        Some domains have pLDDT in the 50-70 range but are still well-structured.
        We rescue regions that:
        1. Form a coherent cluster (same cluster assignment)
        2. Have good internal contact density
        3. Have reasonable pLDDT (>= 50)
        4. Have sufficient size (>= min_domain_size)

        Returns:
            Updated NDR indices with rescued residues removed
        """
        if not ndr_indices:
            return ndr_indices

        # Group NDR residues by their cluster assignment
        ndr_by_cluster = {}
        for idx in ndr_indices:
            cluster = cluster_assignments.get(idx, -1)
            if cluster == -1:
                continue  # Skip unassigned residues
            if cluster not in ndr_by_cluster:
                ndr_by_cluster[cluster] = []
            ndr_by_cluster[cluster].append(idx)

        rescued = set()
        moderate_plddt_threshold = 50.0  # Lower threshold for potential rescue

        for cluster_id, indices in ndr_by_cluster.items():
            if len(indices) < self.min_domain_size:
                continue

            # Check average pLDDT of the cluster
            avg_plddt = np.mean([residues[i].plddt for i in indices])
            if avg_plddt < moderate_plddt_threshold:
                continue  # Too low pLDDT, likely truly disordered

            # Check if the cluster is contiguous (not scattered)
            sorted_indices = sorted(indices)
            max_gap = 0
            for i in range(1, len(sorted_indices)):
                gap = sorted_indices[i] - sorted_indices[i-1]
                max_gap = max(max_gap, gap)

            if max_gap > 20:  # Too fragmented
                continue

            # Check internal contact density
            internal_contacts = 0
            external_contacts = 0
            for idx in indices:
                if idx not in graph:
                    continue
                for neighbor in graph.neighbors(idx):
                    if neighbor in indices:
                        internal_contacts += 1
                    else:
                        external_contacts += 1

            internal_contacts //= 2  # Each edge counted twice

            # Calculate contact density ratio
            if external_contacts > 0:
                contact_ratio = internal_contacts / external_contacts
            else:
                contact_ratio = internal_contacts if internal_contacts > 0 else 0

            # Calculate internal contact density (contacts per residue)
            contacts_per_residue = (internal_contacts * 2) / len(indices)

            # Rescue if:
            # - Good contact ratio (internal > external) OR
            # - High internal contact density (>= 3 contacts/residue on average)
            # - AND pLDDT is not terrible
            should_rescue = (
                (contact_ratio >= 1.0 or contacts_per_residue >= 3.0) and
                avg_plddt >= moderate_plddt_threshold
            )

            if should_rescue:
                logger.info(
                    f"Rescuing cluster {cluster_id} from NDR: "
                    f"{len(indices)} residues, pLDDT={avg_plddt:.1f}, "
                    f"contact_ratio={contact_ratio:.2f}, contacts/res={contacts_per_residue:.1f}"
                )
                rescued.update(indices)

        # Return NDR indices with rescued residues removed
        return ndr_indices - rescued

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
        reassigned_to_ndr = set()

        for cluster_id, indices in clusters.items():
            indices = sorted(indices)
            if not indices:
                continue

            # Convert indices to residue numbers and find segments
            segments = self._indices_to_segments(residues, indices)

            # Filter out tiny isolated segments (e.g., single residues far from core)
            filtered_segments, removed_resnums = self._filter_isolated_segments(
                segments,
                min_segment_size=5,
                max_isolation_distance=30,
            )

            # Track removed residues for NDR assignment
            for resnum in removed_resnums:
                # Convert resnum back to index
                for idx in indices:
                    if residues[idx].resnum == resnum:
                        reassigned_to_ndr.add(idx)
                        break

            # Update indices to only include kept segments
            kept_indices = [
                idx for idx in indices
                if idx not in reassigned_to_ndr
            ]

            if filtered_segments and kept_indices:
                domains.append(Domain(
                    domain_id=cluster_id,
                    segments=filtered_segments,
                    residue_indices=kept_indices,
                ))

        # Add reassigned residues to ndr_indices
        ndr_indices.update(reassigned_to_ndr)

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

    def _filter_isolated_segments(
        self,
        segments: List[Tuple[int, int]],
        min_segment_size: int = 5,
        max_isolation_distance: int = 30,
    ) -> Tuple[List[Tuple[int, int]], List[int]]:
        """
        Filter out isolated segments that are far from the domain core.

        The core is built starting from the largest segment, then iteratively
        adding segments that are close to the current core. Segments that
        remain far from the core are removed regardless of their size.

        Returns:
            (filtered_segments, removed_residue_indices)
        """
        if len(segments) <= 1:
            return segments, []

        # Calculate segment sizes and sort by size (largest first)
        seg_with_sizes = [(end - start + 1, start, end) for start, end in segments]
        seg_with_sizes.sort(reverse=True)

        # Start core with the largest segment
        largest_size, largest_start, largest_end = seg_with_sizes[0]
        core_segments = [(largest_start, largest_end)]
        remaining = [(s, e, sz) for sz, s, e in seg_with_sizes[1:]]

        # Iteratively add segments that are close to current core
        # Small segments need to be VERY close (adjacent), larger segments can be further
        changed = True
        while changed and remaining:
            changed = False
            still_remaining = []

            for start, end, size in remaining:
                # Check distance to nearest core segment
                min_dist = float('inf')
                for core_start, core_end in core_segments:
                    if end < core_start:
                        dist = core_start - end
                    elif start > core_end:
                        dist = start - core_end
                    else:
                        dist = 0  # Overlapping
                    min_dist = min(min_dist, dist)

                # Small segments (< 10 residues) must be truly adjacent (gap <= 5)
                # Larger segments can be up to max_isolation_distance away
                if size < 10:
                    allowed_distance = 5  # Tiny segments must be adjacent
                else:
                    allowed_distance = max_isolation_distance

                if min_dist <= allowed_distance:
                    # Add to core - it's close enough
                    core_segments.append((start, end))
                    changed = True
                else:
                    still_remaining.append((start, end, size))

            remaining = still_remaining

        # All remaining segments are too far from core - remove them
        removed_residues = []
        for start, end, size in remaining:
            # Only remove if it's a small segment OR very far from core
            min_dist = float('inf')
            for core_start, core_end in core_segments:
                if end < core_start:
                    dist = core_start - end
                elif start > core_end:
                    dist = start - core_end
                else:
                    dist = 0
                min_dist = min(min_dist, dist)

            # Remove if: small segment, OR far away (>50 residues)
            if size < min_segment_size or min_dist > 50:
                removed_residues.extend(range(start, end + 1))
            else:
                # Large segment that's moderately far - keep it but warn
                core_segments.append((start, end))

        # Sort segments by start position
        core_segments.sort()
        return core_segments, removed_residues

    def _refine_domains(
        self,
        domains: List[Domain],
        residues: List[Residue],
        graph,
    ) -> List[Domain]:
        """
        Refine domain boundaries by merging over-split domains.

        This is especially important for repeat proteins (e.g., WD40 beta-propellers)
        where the clustering algorithm may split each repeat into a separate domain.

        We merge domains that:
        1. Have low ContactRatio (< 4) - indicating many external contacts
        2. Have high inter-domain contact density with each other
        3. Are adjacent or nearly adjacent in sequence
        """
        if len(domains) <= 1:
            return domains

        # Sort domains by start position
        domains = sorted(domains, key=lambda d: min(s[0] for s in d.segments))

        # Calculate inter-domain contact density matrix
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

        # First pass: standard merge for high inter-domain contact density
        merge_threshold = 0.5
        merged = set()
        intermediate_domains = []

        for i, domain in enumerate(domains):
            if i in merged:
                continue

            to_merge = [i]
            for j in range(i + 1, n_domains):
                if j in merged:
                    continue
                if contact_matrix[i, j] > merge_threshold:
                    to_merge.append(j)
                    merged.add(j)

            if len(to_merge) == 1:
                intermediate_domains.append(domain)
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
                metrics = self.quality_assessor.assess(merged_domain, residues, graph)
                merged_domain.quality_metrics = metrics
                intermediate_domains.append(merged_domain)

        # Second pass: merge adjacent domains with low ContactRatio
        # This handles repeat proteins like WD40 beta-propellers
        if len(intermediate_domains) > 1:
            intermediate_domains = self._merge_repeat_domains(
                intermediate_domains, residues, graph
            )

        # Third pass: merge non-adjacent domains connected by bridges
        # This handles discontinuous domains like A-B-A' topology
        if len(intermediate_domains) > 1:
            intermediate_domains = self._merge_bridged_domains(
                intermediate_domains, residues, graph
            )

        return intermediate_domains

    def _merge_repeat_domains(
        self,
        domains: List[Domain],
        residues: List[Residue],
        graph,
    ) -> List[Domain]:
        """
        Merge domains that appear to be parts of a repeat protein.

        Repeat proteins (WD40, TPR, HEAT, etc.) have domains with:
        - Low ContactRatio (< 4) because repeats contact each other
        - High inter-domain contacts between adjacent domains
        - Similar sizes (roughly)
        """
        if len(domains) <= 1:
            return domains

        # Sort by position
        domains = sorted(domains, key=lambda d: min(s[0] for s in d.segments))

        # Threshold for "low" contact ratio - domains below this might be repeat units
        low_contact_ratio_threshold = 4.0

        # Check if we should merge adjacent low-ContactRatio domains
        # Build a list of domains to merge together
        merge_groups = []
        current_group = [0]

        for i in range(1, len(domains)):
            prev_domain = domains[i - 1]
            curr_domain = domains[i]

            # Get ContactRatio for both domains
            prev_ratio = prev_domain.quality_metrics.get('contact_density_ratio', float('inf'))
            curr_ratio = curr_domain.quality_metrics.get('contact_density_ratio', float('inf'))

            # Check if domains are adjacent in sequence
            prev_end = max(e for s, e in prev_domain.segments)
            curr_start = min(s for s, e in curr_domain.segments)
            gap = curr_start - prev_end

            # Calculate inter-domain contacts
            inter_contacts = 0
            for idx_i in prev_domain.residue_indices:
                for idx_j in curr_domain.residue_indices:
                    if graph.has_edge(idx_i, idx_j):
                        inter_contacts += 1

            # Normalize by smaller domain size
            smaller_size = min(len(prev_domain.residue_indices), len(curr_domain.residue_indices))
            normalized_contacts = inter_contacts / smaller_size if smaller_size > 0 else 0

            # Merge if:
            # 1. Both domains have low ContactRatio (likely parts of repeat)
            # 2. They are adjacent or nearly adjacent (gap < 30 residues)
            # 3. They have significant inter-domain contacts (> 0.3 per residue)
            should_merge = (
                prev_ratio < low_contact_ratio_threshold and
                curr_ratio < low_contact_ratio_threshold and
                gap < 30 and
                normalized_contacts > 0.3
            )

            if should_merge:
                current_group.append(i)
            else:
                merge_groups.append(current_group)
                current_group = [i]

        merge_groups.append(current_group)

        # Now merge domains within each group
        final_domains = []
        for group in merge_groups:
            if len(group) == 1:
                final_domains.append(domains[group[0]])
            else:
                # Merge all domains in this group
                all_indices = []
                for idx in group:
                    all_indices.extend(domains[idx].residue_indices)
                all_indices = sorted(set(all_indices))

                segments = self._indices_to_segments(residues, all_indices)
                merged_domain = Domain(
                    domain_id=domains[group[0]].domain_id,
                    segments=segments,
                    residue_indices=all_indices,
                )
                metrics = self.quality_assessor.assess(merged_domain, residues, graph)
                merged_domain.quality_metrics = metrics

                logger.info(
                    f"Merged {len(group)} repeat-like domains into one: "
                    f"size={merged_domain.size}, ContactRatio={metrics.get('contact_density_ratio', 0):.2f}"
                )

                final_domains.append(merged_domain)

        return final_domains

    def _merge_bridged_domains(
        self,
        domains: List[Domain],
        residues: List[Residue],
        graph,
    ) -> List[Domain]:
        """
        Merge non-adjacent domains that are connected by long-range contacts (bridges).

        This handles discontinuous domains with A-B-A' topology where:
        - Domain A has segments at positions 1-50 and 200-250
        - Domain B (inserted) is at positions 51-199
        - A and A' are connected by bridges (long-range contacts)

        The key insight is that bridges indicate structural continuity despite
        sequence discontinuity. If two non-adjacent domains share many bridges,
        they likely belong to the same discontinuous fold.
        """
        if len(domains) <= 1:
            return domains

        # Sort domains by start position
        domains = sorted(domains, key=lambda d: min(s[0] for s in d.segments))

        # Detect bridges between all domain pairs
        # A bridge is when residue i from domain A contacts residue j from domain B
        # where |i - j| > min_sequence_separation
        min_sequence_separation = self.min_bridge_seq_sep
        min_bridge_strength = 3.0  # Minimum total edge weight for bridge
        min_bridges_to_merge = self.min_bridge_count

        n_domains = len(domains)
        bridge_counts = np.zeros((n_domains, n_domains))
        bridge_strengths = np.zeros((n_domains, n_domains))

        for i, dom_i in enumerate(domains):
            for j, dom_j in enumerate(domains):
                if i >= j:
                    continue

                # Count bridges between domain i and domain j
                bridges = 0
                total_strength = 0.0

                for idx_i in dom_i.residue_indices:
                    for idx_j in dom_j.residue_indices:
                        if not graph.has_edge(idx_i, idx_j):
                            continue

                        # Check sequence separation
                        seq_dist = abs(residues[idx_i].resnum - residues[idx_j].resnum)
                        if seq_dist >= min_sequence_separation:
                            bridges += 1
                            total_strength += graph[idx_i][idx_j].get("weight", 1.0)

                bridge_counts[i, j] = bridges
                bridge_counts[j, i] = bridges
                bridge_strengths[i, j] = total_strength
                bridge_strengths[j, i] = total_strength

        # Find domain pairs to merge based on bridge connections
        # Use Union-Find to track merge groups
        parent = list(range(n_domains))

        def find(x):
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]

        def union(x, y):
            px, py = find(x), find(y)
            if px != py:
                parent[px] = py

        # Merge domains connected by strong bridges
        for i in range(n_domains):
            for j in range(i + 1, n_domains):
                # Check if domains are connected by enough bridges
                if (bridge_counts[i, j] >= min_bridges_to_merge and
                    bridge_strengths[i, j] >= min_bridge_strength):

                    # Additional check: domains should not overlap too much in sequence
                    # (if they overlap, they're probably already merged correctly)
                    i_range = (min(s[0] for s in domains[i].segments),
                               max(s[1] for s in domains[i].segments))
                    j_range = (min(s[0] for s in domains[j].segments),
                               max(s[1] for s in domains[j].segments))

                    # Calculate overlap
                    overlap_start = max(i_range[0], j_range[0])
                    overlap_end = min(i_range[1], j_range[1])
                    overlap = max(0, overlap_end - overlap_start)

                    # Only merge if minimal overlap (discontinuous case)
                    i_size = i_range[1] - i_range[0]
                    j_size = j_range[1] - j_range[0]
                    min_size = min(i_size, j_size)

                    if overlap < min_size * 0.3:  # Less than 30% overlap
                        logger.info(
                            f"Bridge connection: domains {i} ({domains[i].segments}) and "
                            f"{j} ({domains[j].segments}): {int(bridge_counts[i, j])} bridges, "
                            f"strength={bridge_strengths[i, j]:.1f}"
                        )
                        union(i, j)

        # Group domains by their root parent
        groups = {}
        for i in range(n_domains):
            root = find(i)
            if root not in groups:
                groups[root] = []
            groups[root].append(i)

        # Merge domains within each group
        final_domains = []
        for root, group in groups.items():
            if len(group) == 1:
                final_domains.append(domains[group[0]])
            else:
                # Merge all domains in this group
                all_indices = []
                for idx in group:
                    all_indices.extend(domains[idx].residue_indices)
                all_indices = sorted(set(all_indices))

                segments = self._indices_to_segments(residues, all_indices)
                merged_domain = Domain(
                    domain_id=domains[group[0]].domain_id,
                    segments=segments,
                    residue_indices=all_indices,
                )
                metrics = self.quality_assessor.assess(merged_domain, residues, graph)
                merged_domain.quality_metrics = metrics

                logger.info(
                    f"Merged {len(group)} bridge-connected domains into discontinuous domain: "
                    f"segments={segments}, size={merged_domain.size}, "
                    f"ContactRatio={metrics.get('contact_density_ratio', 0):.2f}"
                )

                final_domains.append(merged_domain)

        return final_domains

    def _enforce_domain_continuity(
        self,
        domains: List[Domain],
        residues: List[Residue],
    ) -> Tuple[List[Domain], set]:
        """
        Enforce domain continuity by:
        1. Resolving interdigitations (small segments of one domain within another's range)
        2. Filling small internal gaps within domains
        3. Converting boundary interdigitations to NDR

        Returns:
            (cleaned_domains, residues_to_add_to_ndr)
        """
        if not domains:
            return domains, set()

        extra_ndr = set()

        # First, identify each domain's "primary range" (min to max residue number)
        domain_ranges = []
        for domain in domains:
            all_resnums = []
            for start, end in domain.segments:
                all_resnums.extend(range(start, end + 1))
            if all_resnums:
                domain_ranges.append((min(all_resnums), max(all_resnums), domain))

        # Sort domains by their start position
        domain_ranges.sort(key=lambda x: x[0])

        # Check for interdigitations: small segments of domain A within domain B's range
        cleaned_domains = []
        for i, (start_i, end_i, domain_i) in enumerate(domain_ranges):
            new_segments = []
            removed_indices = set()

            for seg_start, seg_end in domain_i.segments:
                seg_size = seg_end - seg_start + 1
                seg_dominated = False

                # Check if this segment is "dominated" by another domain
                for j, (start_j, end_j, domain_j) in enumerate(domain_ranges):
                    if i == j:
                        continue

                    # Is this segment entirely within another domain's main range?
                    # And is that other domain much larger in this region?
                    if seg_start >= start_j and seg_end <= end_j:
                        # This segment is within domain_j's range
                        # Check if domain_j has substantial presence here
                        j_residues_in_range = sum(
                            1 for s, e in domain_j.segments
                            for r in range(s, e + 1)
                            if r >= seg_start - 20 and r <= seg_end + 20
                        )

                        # If the other domain has more residues nearby, this segment is an intrusion
                        if j_residues_in_range > seg_size * 2 and seg_size < 15:
                            seg_dominated = True
                            # Mark these residues for NDR
                            for idx in domain_i.residue_indices:
                                if seg_start <= residues[idx].resnum <= seg_end:
                                    removed_indices.add(idx)
                                    extra_ndr.add(idx)
                            break

                if not seg_dominated:
                    new_segments.append((seg_start, seg_end))

            if new_segments:
                # Merge adjacent segments (fill small gaps)
                merged_segments = self._merge_adjacent_segments(new_segments, max_gap=20)

                # Update domain
                new_indices = [idx for idx in domain_i.residue_indices if idx not in removed_indices]

                # Also add indices for filled gaps
                for seg_start, seg_end in merged_segments:
                    for idx, res in enumerate(residues):
                        if seg_start <= res.resnum <= seg_end and idx not in new_indices:
                            new_indices.append(idx)

                new_indices = sorted(set(new_indices))

                if len(new_indices) >= self.min_domain_size:
                    cleaned_domains.append(Domain(
                        domain_id=domain_i.domain_id,
                        segments=merged_segments,
                        residue_indices=new_indices,
                        quality_metrics=domain_i.quality_metrics,
                    ))
                else:
                    # Domain too small after cleaning, convert to NDR
                    extra_ndr.update(new_indices)

        return cleaned_domains, extra_ndr

    def _merge_adjacent_segments(
        self,
        segments: List[Tuple[int, int]],
        max_gap: int = 10,
    ) -> List[Tuple[int, int]]:
        """Merge segments that are separated by small gaps."""
        if not segments:
            return []

        sorted_segs = sorted(segments)
        merged = [sorted_segs[0]]

        for start, end in sorted_segs[1:]:
            prev_start, prev_end = merged[-1]

            # If gap is small enough, merge
            if start - prev_end <= max_gap:
                merged[-1] = (prev_start, end)
            else:
                merged.append((start, end))

        return merged

    def _create_domains_from_high_plddt_regions(
        self,
        residues: List[Residue],
        ndr_indices: set,
        graph,
    ) -> List[Domain]:
        """
        Create domains from high-pLDDT regions when all other domains were rejected.

        This is a fallback for cases like A6NG13 where clustering fragments the
        structure and all assembled domains get rejected, but there's clearly
        structured content (high overall pLDDT).
        """
        # Find contiguous high-pLDDT regions
        high_plddt_indices = []
        for idx, res in enumerate(residues):
            if res.plddt >= self.ndr_plddt_cutoff:
                high_plddt_indices.append(idx)

        if not high_plddt_indices:
            return []

        # Group into contiguous stretches
        stretches = []
        current_stretch = [high_plddt_indices[0]]

        for idx in high_plddt_indices[1:]:
            # Allow gaps up to 15 residues within a stretch
            if idx - current_stretch[-1] <= 15:
                current_stretch.append(idx)
            else:
                if len(current_stretch) >= self.min_domain_size:
                    stretches.append(current_stretch)
                current_stretch = [idx]

        if len(current_stretch) >= self.min_domain_size:
            stretches.append(current_stretch)

        if not stretches:
            return []

        logger.info(f"Creating {len(stretches)} domain(s) from high-pLDDT regions (fallback)")

        domains = []
        for i, stretch in enumerate(stretches):
            # Fill in any small gaps within the stretch
            filled_indices = list(range(min(stretch), max(stretch) + 1))

            segments = self._indices_to_segments(residues, filled_indices)
            avg_plddt = np.mean([residues[idx].plddt for idx in filled_indices])

            domain = Domain(
                domain_id=i + 1,
                segments=segments,
                residue_indices=filled_indices,
                quality_metrics={},
            )

            # Assess quality
            metrics = self.quality_assessor.assess(domain, residues, graph)
            domain.quality_metrics = metrics

            start = residues[filled_indices[0]].resnum
            end = residues[filled_indices[-1]].resnum
            logger.info(
                f"  Created domain {start}-{end} "
                f"(size={len(filled_indices)}, pLDDT={avg_plddt:.1f}, CR={metrics.get('contact_ratio', 0):.2f})"
            )

            domains.append(domain)

        return domains

    def _extend_domains_to_structured_regions(
        self,
        domains: List[Domain],
        residues: List[Residue],
        ndr_indices: set,
        graph,
    ) -> List[Domain]:
        """
        Extend domains to absorb adjacent high-pLDDT unassigned regions.

        This handles cases where a structured region (high pLDDT) was clustered
        separately from an adjacent domain but should be part of it. Common in
        repeat proteins like WD40 where blades may be split.

        The function iterates until no more extensions can be made:
        1. Finds unassigned regions with high average pLDDT
        2. Checks if they're adjacent to existing domains
        3. If they have contacts with the adjacent domain, extends the domain
        4. Repeats to catch regions that become adjacent after extension
        """
        if not domains:
            return domains

        max_iterations = 10  # Safety limit
        iteration = 0

        while iteration < max_iterations:
            iteration += 1
            made_extension = False

            # Get all residue indices currently in domains
            domain_indices = set()
            for domain in domains:
                domain_indices.update(domain.residue_indices)

            # Find unassigned indices (not in domains)
            all_indices = set(range(len(residues)))
            unassigned = all_indices - domain_indices

            if not unassigned:
                break

            # Group unassigned indices into contiguous stretches
            unassigned_sorted = sorted(unassigned)
            stretches = []
            current_stretch = [unassigned_sorted[0]]

            for idx in unassigned_sorted[1:]:
                if idx - current_stretch[-1] <= 3:  # Allow small gaps
                    current_stretch.append(idx)
                else:
                    if len(current_stretch) >= 10:  # Minimum size to consider
                        stretches.append(current_stretch)
                    current_stretch = [idx]

            if len(current_stretch) >= 10:
                stretches.append(current_stretch)

            # For each stretch, find high-pLDDT sub-regions that should be rescued
            high_plddt_stretches = []
            for stretch in stretches:
                # Split stretch into high and low pLDDT sub-regions
                # This handles cases like 679-806 where 679-737 is low pLDDT
                # but 738-806 is high pLDDT
                current_sub = []
                for idx in stretch:
                    if residues[idx].plddt >= self.ndr_plddt_cutoff:
                        current_sub.append(idx)
                    else:
                        if len(current_sub) >= 10:  # Minimum size for high-pLDDT sub-region
                            high_plddt_stretches.append(current_sub)
                        current_sub = []
                if len(current_sub) >= 10:
                    high_plddt_stretches.append(current_sub)

            # Now process each high-pLDDT sub-region
            for stretch in high_plddt_stretches:
                avg_plddt = np.mean([residues[i].plddt for i in stretch])
                stretch_start = residues[stretch[0]].resnum
                stretch_end = residues[stretch[-1]].resnum

                logger.info(
                    f"Found high-pLDDT unassigned region: {stretch_start}-{stretch_end} "
                    f"(pLDDT={avg_plddt:.1f}, size={len(stretch)})"
                )

                # Find adjacent domain(s)
                best_domain = None
                best_contacts = 0

                for domain in domains:
                    dom_start = min(s[0] for s in domain.segments)
                    dom_end = max(s[1] for s in domain.segments)

                    # Check if stretch is adjacent (within 50 residues)
                    if stretch_end < dom_start - 50 or stretch_start > dom_end + 50:
                        continue

                    # Count contacts between stretch and domain
                    contacts = 0
                    for idx in stretch:
                        if idx not in graph:
                            continue
                        for neighbor in graph.neighbors(idx):
                            if neighbor in domain.residue_indices:
                                contacts += 1

                    if contacts > best_contacts:
                        best_contacts = contacts
                        best_domain = domain

                # Extend the domain if sufficient contacts
                min_contacts_per_residue = 0.5  # At least 0.5 contacts per residue
                if best_domain is not None and best_contacts >= len(stretch) * min_contacts_per_residue:
                    logger.info(
                        f"  Extending domain {best_domain.segments} to include "
                        f"{stretch_start}-{stretch_end} ({best_contacts} contacts)"
                    )

                    # Add stretch residues to domain
                    new_indices = sorted(set(best_domain.residue_indices + stretch))

                    # Rebuild segments
                    new_segments = self._indices_to_segments(residues, new_indices)

                    # Update domain
                    best_domain.residue_indices = new_indices
                    best_domain.segments = new_segments

                    # Re-assess quality
                    metrics = self.quality_assessor.assess(best_domain, residues, graph)
                    best_domain.quality_metrics = metrics

                    made_extension = True

            # If no extensions were made this iteration, we're done
            if not made_extension:
                break

        if iteration > 1:
            logger.info(f"Domain extension completed after {iteration} iterations")

        # Post-process: fill internal gaps in domains
        # Fill small gaps unconditionally, larger gaps only if high pLDDT
        max_gap_unconditional = 10  # Always fill gaps up to this size
        max_gap_high_plddt = 20     # Fill larger gaps if pLDDT >= threshold
        for domain in domains:
            if len(domain.segments) > 1:
                # Check gaps between segments
                filled_indices = list(domain.residue_indices)
                for i in range(len(domain.segments) - 1):
                    seg1_end = domain.segments[i][1]
                    seg2_start = domain.segments[i + 1][0]
                    gap_size = seg2_start - seg1_end - 1

                    if gap_size <= 0:
                        continue

                    # Find gap residue indices and calculate pLDDT
                    gap_indices = []
                    for idx, res in enumerate(residues):
                        if seg1_end < res.resnum < seg2_start:
                            gap_indices.append(idx)

                    if not gap_indices:
                        continue

                    gap_plddt = np.mean([residues[i].plddt for i in gap_indices])

                    # Decide whether to fill
                    should_fill = False
                    if gap_size <= max_gap_unconditional:
                        should_fill = True
                    elif gap_size <= max_gap_high_plddt and gap_plddt >= self.ndr_plddt_cutoff:
                        should_fill = True

                    if should_fill:
                        filled_indices.extend(gap_indices)
                        logger.info(
                            f"Filling internal gap {seg1_end+1}-{seg2_start-1} "
                            f"({gap_size} residues, pLDDT={gap_plddt:.1f}) in domain"
                        )

                if len(filled_indices) > len(domain.residue_indices):
                    domain.residue_indices = sorted(set(filled_indices))
                    domain.segments = self._indices_to_segments(residues, domain.residue_indices)
                    metrics = self.quality_assessor.assess(domain, residues, graph)
                    domain.quality_metrics = metrics

        # Post-process: merge adjacent domains with sufficient contacts
        # This handles cases where separate clusters should be one domain (e.g., WD40 blades)
        merged = True
        while merged:
            merged = False
            domains_sorted = sorted(domains, key=lambda d: min(s[0] for s in d.segments))
            for i in range(len(domains_sorted) - 1):
                d1 = domains_sorted[i]
                d2 = domains_sorted[i + 1]
                d1_end = max(s[1] for s in d1.segments)
                d2_start = min(s[0] for s in d2.segments)

                # Check if adjacent (within 20 residues)
                if d2_start - d1_end > 20:
                    continue

                # Count contacts between domains
                contacts = 0
                d1_indices = set(d1.residue_indices)
                d2_indices = set(d2.residue_indices)
                for idx in d1_indices:
                    if idx not in graph:
                        continue
                    for neighbor in graph.neighbors(idx):
                        if neighbor in d2_indices:
                            contacts += 1

                # Merge if sufficient contacts (0.3 per residue of smaller domain)
                min_size = min(len(d1_indices), len(d2_indices))
                if contacts >= min_size * 0.3:
                    logger.info(
                        f"Merging adjacent domains {d1.segments} and {d2.segments} "
                        f"({contacts} contacts)"
                    )
                    # Merge into d1
                    new_indices = sorted(set(d1.residue_indices + d2.residue_indices))
                    d1.residue_indices = new_indices
                    d1.segments = self._indices_to_segments(residues, new_indices)
                    metrics = self.quality_assessor.assess(d1, residues, graph)
                    d1.quality_metrics = metrics
                    domains.remove(d2)
                    merged = True
                    break

        return domains

    def _build_ndr_regions(
        self,
        residues: List[Residue],
        ndr_indices: set,
        domains: List[Domain],
    ) -> List[NDRRegion]:
        """Build NDR regions from collected indices."""
        # Get all indices covered by domains - use SEGMENTS not residue_indices
        # This ensures continuous domains don't have internal NDR holes
        domain_indices = set()
        for domain in domains:
            for seg_start, seg_end in domain.segments:
                # Find all residue indices within this segment range
                for idx, res in enumerate(residues):
                    if seg_start <= res.resnum <= seg_end:
                        domain_indices.add(idx)

        # All remaining indices are NDR
        all_ndr = set(range(len(residues))) - domain_indices
        all_ndr.update(ndr_indices)
        # Remove any indices that are within domain segments
        all_ndr -= domain_indices

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
