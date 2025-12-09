"""
ABC Domain Predictor V2 - Bridge-Based Approach

This version uses bridge detection to identify discontinuous domains and
handle nested domain topologies correctly.

Key differences from V1:
1. Detects long-range "bridges" (contacts between sequence-distant residues)
2. Builds nesting hierarchy from bridges before clustering
3. Respects nesting boundaries during domain assembly
4. Can properly identify A-B-A' (inserted domain) topologies

The bridge-based approach naturally handles:
- Simple inserted domains: A-B-A'
- Nested insertions: A-B-C-B'-A'
- Multiple insertions: A-B-A'-C-A''
"""

import json
import logging
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

# Reuse common components from existing modules
from .contact_graph import ContactGraphBuilder, ContactDensityCalculator
from .domain_quality import DomainQualityAssessor
from .visualize import DomainVisualizer
from .bridge_detector import BridgeDetector, Bridge, DomainRegion, NestingNode

# Import data classes from v1 (these could be factored out to a shared module)
from .abc_predictor import Residue, Domain, NDRRegion, ABCPrediction

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ABCPredictorV2:
    """
    ABC Domain Predictor V2 - Bridge-Based Nested Domain Detection

    Uses bridge detection to identify discontinuous domains and properly
    handle nested domain topologies.

    Parameters:
    -----------
    distance_threshold : float
        Maximum Cα-Cα distance for contact (default: 10Å)
    min_domain_size : int
        Minimum residues for a valid domain (default: 30)
    min_contact_ratio : float
        Minimum internal/external contact ratio (default: 1.0)
    ndr_plddt_cutoff : float
        pLDDT below which residues may be NDR (default: 70)
    min_bridge_separation : int
        Minimum sequence separation to consider a bridge (default: 50)
    min_bridge_contacts : int
        Minimum contacts to support a bridge (default: 3)
    sigma : float
        Gaussian decay for edge weights (default: 8.0)
    use_pae : bool
        Whether to use PAE for edge weighting (default: True)
    """

    def __init__(
        self,
        distance_threshold: float = 10.0,
        min_domain_size: int = 30,
        min_contact_ratio: float = 1.0,
        ndr_plddt_cutoff: float = 70.0,
        min_bridge_separation: int = 50,
        min_bridge_contacts: int = 3,
        sigma: float = 8.0,
        use_pae: bool = True,
    ):
        self.distance_threshold = distance_threshold
        self.min_domain_size = min_domain_size
        self.min_contact_ratio = min_contact_ratio
        self.ndr_plddt_cutoff = ndr_plddt_cutoff
        self.min_bridge_separation = min_bridge_separation
        self.min_bridge_contacts = min_bridge_contacts
        self.sigma = sigma
        self.use_pae = use_pae

        # Reuse components from v1
        self.graph_builder = ContactGraphBuilder(
            distance_threshold=distance_threshold,
            sigma=sigma,
            use_pae=use_pae,
        )
        self.quality_assessor = DomainQualityAssessor()
        self.visualizer = DomainVisualizer()

        # New bridge detector
        self.bridge_detector = BridgeDetector(
            min_sequence_separation=min_bridge_separation,
            min_bridge_contacts=min_bridge_contacts,
        )

    @property
    def parameters(self) -> Dict:
        """Current parameter settings."""
        return {
            "distance_threshold": self.distance_threshold,
            "min_domain_size": self.min_domain_size,
            "min_contact_ratio": self.min_contact_ratio,
            "ndr_plddt_cutoff": self.ndr_plddt_cutoff,
            "min_bridge_separation": self.min_bridge_separation,
            "min_bridge_contacts": self.min_bridge_contacts,
            "sigma": self.sigma,
            "use_pae": self.use_pae,
            "version": "v2_bridge",
        }

    def predict_from_uniprot(self, uniprot_acc: str) -> ABCPrediction:
        """Predict domains by downloading AlphaFold model."""
        import tempfile

        logger.info(f"Downloading AlphaFold model for {uniprot_acc}")

        with tempfile.TemporaryDirectory() as tmpdir:
            cif_path = Path(tmpdir) / f"{uniprot_acc}.cif"
            pae_path = Path(tmpdir) / f"{uniprot_acc}_pae.json"

            versions = ["v6", "v4", "v3", "v2"]
            cif_downloaded = False

            for version in versions:
                cif_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_acc}-F1-model_{version}.cif"
                try:
                    urllib.request.urlretrieve(cif_url, cif_path)
                    logger.info(f"Downloaded structure ({version})")
                    cif_downloaded = True

                    pae_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_acc}-F1-predicted_aligned_error_{version}.json"
                    try:
                        urllib.request.urlretrieve(pae_url, pae_path)
                        logger.info(f"Downloaded PAE ({version})")
                    except Exception:
                        logger.warning(f"Could not download PAE, continuing without it")
                        pae_path = None
                    break
                except Exception:
                    continue

            if not cif_downloaded:
                raise RuntimeError(f"Failed to download structure for {uniprot_acc}")

            return self.predict_from_file(str(cif_path), str(pae_path) if pae_path else None, uniprot_acc)

    def predict_from_file(
        self,
        pdb_path: str,
        pae_path: Optional[str] = None,
        uniprot_acc: Optional[str] = None,
    ) -> ABCPrediction:
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

        pae_matrix = None
        if pae_path and Path(pae_path).exists():
            logger.info(f"Loading PAE from {pae_path}")
            pae_matrix = self._load_pae(pae_path)

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

    def _load_pae(self, pae_path: str) -> Optional[np.ndarray]:
        """Load PAE matrix from AlphaFold JSON file."""
        try:
            with open(pae_path) as f:
                data = json.load(f)

            if isinstance(data, list) and len(data) > 0:
                data = data[0]

            if "predicted_aligned_error" in data:
                return np.array(data["predicted_aligned_error"])
            elif "pae" in data:
                return np.array(data["pae"])
            else:
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
        """Main prediction pipeline using bridge-based approach."""
        n_residues = len(residues)

        # Step 1: Build contact graph
        logger.info("Building contact graph...")
        graph = self.graph_builder.build(residues, pae_matrix)
        logger.info(f"Graph has {graph.number_of_nodes()} nodes and {graph.number_of_edges()} edges")

        # Step 2: Detect bridges (long-range contacts)
        logger.info("Detecting bridges...")
        bridges = self.bridge_detector.detect_bridges(residues, graph)
        logger.info(f"Found {len(bridges)} bridges")

        for bridge in bridges:
            logger.info(f"  Bridge: {bridge.start_resnum}-{bridge.end_resnum} "
                       f"(span={bridge.sequence_distance}, contacts={bridge.n_contacts}, "
                       f"strength={bridge.contact_strength:.2f})")

        # Step 3: Build nesting hierarchy
        logger.info("Building nesting hierarchy...")
        hierarchy = self.bridge_detector.build_nesting_hierarchy(bridges, n_residues)

        # Log hierarchy
        self._log_hierarchy(hierarchy)

        # Step 4: Identify NDR regions (low pLDDT)
        logger.info("Identifying NDR regions...")
        ndr_indices = self._detect_ndr(residues, graph)

        # Step 5: Build domains from bridges and continuous regions
        logger.info("Building domains from bridge analysis...")
        domains = self._build_domains_from_bridges(
            residues, graph, bridges, hierarchy, ndr_indices
        )

        # Step 6: Handle unassigned regions (not part of any bridge)
        domains, ndr_indices = self._handle_unassigned_regions(
            residues, graph, bridges, domains, ndr_indices
        )

        # Step 7: Filter domains by quality
        logger.info("Filtering domains by quality...")
        domains, ndr_indices = self._filter_domains(domains, residues, graph, ndr_indices)

        # Renumber domains
        for i, domain in enumerate(domains):
            domain.domain_id = i + 1

        # Step 8: Build NDR regions
        ndr_regions = self._build_ndr_regions(residues, ndr_indices, domains)

        logger.info(f"Final result: {len(domains)} domains, {len(ndr_regions)} NDR regions")

        return ABCPrediction(
            uniprot_acc=uniprot_acc,
            sequence_length=n_residues,
            domains=domains,
            ndr_regions=ndr_regions,
            parameters=self.parameters,
        )

    def _log_hierarchy(self, node: NestingNode, indent: int = 0):
        """Log the nesting hierarchy for debugging."""
        prefix = "  " * indent
        if node.bridge:
            logger.info(f"{prefix}Bridge: {node.bridge.start_resnum}-{node.bridge.end_resnum} "
                       f"(level {node.level})")
        else:
            logger.info(f"{prefix}Root (level {node.level})")

        for child in node.children:
            self._log_hierarchy(child, indent + 1)

    def _detect_ndr(self, residues: List[Residue], graph) -> Set[int]:
        """Detect NDR residues based on pLDDT and isolation."""
        ndr = set()

        # Low pLDDT regions
        for i, res in enumerate(residues):
            if res.plddt < self.ndr_plddt_cutoff - 20:  # Very low threshold
                ndr.add(i)

        # Isolated residues
        for i in range(len(residues)):
            if i not in graph:
                ndr.add(i)
                continue

            nonlocal_contacts = sum(
                1 for j in graph.neighbors(i) if abs(j - i) >= 10
            )
            if nonlocal_contacts < 2:
                ndr.add(i)

        return ndr

    def _build_domains_from_bridges(
        self,
        residues: List[Residue],
        graph,
        bridges: List[Bridge],
        hierarchy: NestingNode,
        ndr_indices: Set[int],
    ) -> List[Domain]:
        """
        Build domains based on bridge analysis.

        Each bridge defines a discontinuous domain. The domain includes:
        - The region at the start of the bridge
        - The region at the end of the bridge
        - But NOT the regions inside nested bridges
        """
        domains = []
        assigned = set(ndr_indices)  # Start with NDR as assigned

        def process_bridge_domain(node: NestingNode, domain_id: int) -> Optional[Domain]:
            """Process a bridge node into a domain."""
            if node.bridge is None:
                return None

            bridge = node.bridge

            # Find child bridge spans to exclude
            child_spans = []
            for child in node.children:
                if child.bridge:
                    child_spans.append((child.bridge.start_idx, child.bridge.end_idx))

            # Sort child spans
            child_spans.sort()

            # Collect residues for this domain (excluding children)
            domain_indices = []

            # Walk from bridge start to end, skipping child spans and NDR
            current = bridge.start_idx
            for child_start, child_end in child_spans:
                # Add residues from current to just before child
                for idx in range(current, child_start):
                    if idx not in assigned and idx not in ndr_indices:
                        domain_indices.append(idx)
                        assigned.add(idx)
                current = child_end + 1

            # Add remaining residues after last child
            for idx in range(current, bridge.end_idx + 1):
                if idx not in assigned and idx not in ndr_indices:
                    domain_indices.append(idx)
                    assigned.add(idx)

            if len(domain_indices) < self.min_domain_size:
                return None

            # Convert to segments
            segments = self._indices_to_segments(residues, sorted(domain_indices))

            domain = Domain(
                domain_id=domain_id,
                segments=segments,
                residue_indices=sorted(domain_indices),
            )

            return domain

        # Process bridges from outermost to innermost
        domain_id = 1

        def process_node(node: NestingNode):
            nonlocal domain_id
            nonlocal domains

            if node.bridge:
                domain = process_bridge_domain(node, domain_id)
                if domain:
                    domains.append(domain)
                    domain_id += 1

            # Process children (inner domains)
            for child in node.children:
                process_node(child)

        process_node(hierarchy)

        return domains

    def _handle_unassigned_regions(
        self,
        residues: List[Residue],
        graph,
        bridges: List[Bridge],
        domains: List[Domain],
        ndr_indices: Set[int],
    ) -> Tuple[List[Domain], Set[int]]:
        """
        Handle regions not covered by any bridge.

        These are continuous regions that may be:
        - Independent continuous domains
        - NDR (if low quality)
        """
        n = len(residues)

        # Find assigned indices
        assigned = set(ndr_indices)
        for domain in domains:
            assigned.update(domain.residue_indices)

        # Find unassigned regions
        unassigned = [i for i in range(n) if i not in assigned]

        if not unassigned:
            return domains, ndr_indices

        # Group into contiguous regions
        regions = []
        current_region = [unassigned[0]]

        for idx in unassigned[1:]:
            if idx - current_region[-1] <= 3:
                current_region.append(idx)
            else:
                if len(current_region) >= 5:
                    regions.append(current_region)
                current_region = [idx]

        if len(current_region) >= 5:
            regions.append(current_region)

        # Evaluate each region
        new_domains = list(domains)
        new_ndr = set(ndr_indices)

        for region_indices in regions:
            # Check if this region is a valid domain
            if len(region_indices) < self.min_domain_size:
                new_ndr.update(region_indices)
                continue

            # Calculate quality metrics
            avg_plddt = np.mean([residues[i].plddt for i in region_indices])

            if avg_plddt < self.ndr_plddt_cutoff:
                new_ndr.update(region_indices)
                continue

            # Calculate contact ratio
            internal = ContactDensityCalculator.internal_contacts(graph, region_indices)
            external = ContactDensityCalculator.external_contacts(graph, region_indices)
            contact_ratio = internal / external if external > 0 else float('inf')

            if contact_ratio < self.min_contact_ratio:
                # Low contact ratio - might be part of an adjacent domain
                # For now, mark as NDR
                logger.info(f"Region {residues[region_indices[0]].resnum}-{residues[region_indices[-1]].resnum} "
                           f"has low CR={contact_ratio:.2f}, marking as unassigned")
                new_ndr.update(region_indices)
                continue

            # Create domain for this region
            segments = self._indices_to_segments(residues, region_indices)
            domain = Domain(
                domain_id=len(new_domains) + 1,
                segments=segments,
                residue_indices=region_indices,
            )
            new_domains.append(domain)
            logger.info(f"Created continuous domain: {segments} "
                       f"(CR={contact_ratio:.2f}, pLDDT={avg_plddt:.1f})")

        return new_domains, new_ndr

    def _filter_domains(
        self,
        domains: List[Domain],
        residues: List[Residue],
        graph,
        ndr_indices: Set[int],
    ) -> Tuple[List[Domain], Set[int]]:
        """Filter domains by quality metrics."""
        filtered = []
        new_ndr = set(ndr_indices)

        for domain in domains:
            metrics = self.quality_assessor.assess(domain, residues, graph)
            domain.quality_metrics = metrics

            contact_ratio = metrics.get('contact_density_ratio', 0)
            avg_plddt = metrics.get('avg_plddt', 0)

            # Filter by size
            if domain.size < self.min_domain_size:
                logger.info(f"Rejecting domain {domain.segments}: too small "
                           f"(size={domain.size})")
                new_ndr.update(domain.residue_indices)
                continue

            # Filter by contact ratio
            if contact_ratio < self.min_contact_ratio and contact_ratio != float('inf'):
                logger.info(f"Rejecting domain {domain.segments}: low CR "
                           f"(CR={contact_ratio:.2f}, size={domain.size}, pLDDT={avg_plddt:.1f})")
                new_ndr.update(domain.residue_indices)
                continue

            filtered.append(domain)
            logger.info(f"Accepted domain {domain.segments}: "
                       f"CR={contact_ratio:.2f}, size={domain.size}, pLDDT={avg_plddt:.1f}")

        return filtered, new_ndr

    def _indices_to_segments(
        self,
        residues: List[Residue],
        indices: List[int],
    ) -> List[Tuple[int, int]]:
        """Convert residue indices to segments."""
        if not indices:
            return []

        indices = sorted(indices)
        segments = []

        start_idx = indices[0]
        prev_idx = indices[0]

        for idx in indices[1:]:
            if idx - prev_idx > 3:  # New segment if gap > 3
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

    def _build_ndr_regions(
        self,
        residues: List[Residue],
        ndr_indices: Set[int],
        domains: List[Domain],
    ) -> List[NDRRegion]:
        """Build NDR regions from indices."""
        # Get domain coverage
        domain_indices = set()
        for domain in domains:
            for start, end in domain.segments:
                for idx, res in enumerate(residues):
                    if start <= res.resnum <= end:
                        domain_indices.add(idx)

        # All non-domain residues
        all_ndr = set(range(len(residues))) - domain_indices

        if not all_ndr:
            return []

        # Group into contiguous regions
        sorted_ndr = sorted(all_ndr)
        regions = []
        current = [sorted_ndr[0]]

        for idx in sorted_ndr[1:]:
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

    def visualize(
        self,
        prediction: ABCPrediction,
        output_path: Optional[str] = None,
        residues: Optional[List[Residue]] = None,
    ) -> str:
        """Generate visualization."""
        return self.visualizer.visualize(prediction, output_path, residues)
