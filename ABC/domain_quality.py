"""
Domain Quality Assessment

Evaluates domain quality based on:
- Radius of gyration (compactness)
- Contact density ratio (internal/external contacts)
- Average pLDDT
- Secondary structure content (if available)
"""

from typing import TYPE_CHECKING, Dict, List, Optional

import networkx as nx
import numpy as np

if TYPE_CHECKING:
    from .abc_predictor import Domain, Residue

from .contact_graph import ContactDensityCalculator


class DomainQualityAssessor:
    """
    Assesses quality of predicted domains.

    Quality metrics include:
    - radius_of_gyration: Compactness measure (lower = more compact)
    - contact_density_ratio: Internal/external contacts (higher = better defined)
    - normalized_contact_density: Fraction of possible internal contacts
    - avg_plddt: Average confidence score
    - min_plddt: Minimum confidence
    - size: Number of residues
    - n_segments: Number of discontinuous segments
    - segment_sizes: Size of each segment
    """

    def assess(
        self,
        domain: "Domain",
        residues: List["Residue"],
        graph: nx.Graph,
    ) -> Dict:
        """
        Calculate all quality metrics for a domain.

        Parameters:
        -----------
        domain : Domain
            Domain to assess
        residues : List[Residue]
            All residues in the protein
        graph : nx.Graph
            Contact graph

        Returns:
        --------
        Dict of quality metrics
        """
        indices = domain.residue_indices
        coords = np.array([residues[i].ca_coord for i in indices])
        plddt_values = np.array([residues[i].plddt for i in indices])

        metrics = {}

        # Compactness metrics
        metrics["radius_of_gyration"] = self._radius_of_gyration(coords)
        metrics["max_dimension"] = self._max_dimension(coords)

        # Contact density metrics
        metrics["contact_density_ratio"] = ContactDensityCalculator.contact_density_ratio(
            graph, indices
        )
        metrics["normalized_contact_density"] = ContactDensityCalculator.normalized_contact_density(
            graph, indices
        )
        metrics["internal_contacts"] = ContactDensityCalculator.internal_contacts(
            graph, indices
        )
        metrics["external_contacts"] = ContactDensityCalculator.external_contacts(
            graph, indices
        )

        # pLDDT metrics
        metrics["avg_plddt"] = float(np.mean(plddt_values))
        metrics["min_plddt"] = float(np.min(plddt_values))
        metrics["std_plddt"] = float(np.std(plddt_values))

        # Size metrics
        metrics["size"] = len(indices)
        metrics["n_segments"] = len(domain.segments)
        metrics["segment_sizes"] = [end - start + 1 for start, end in domain.segments]

        # Sequence continuity
        metrics["sequence_coverage"] = self._sequence_coverage(domain.segments)

        # Overall quality score (composite)
        metrics["quality_score"] = self._calculate_quality_score(metrics)

        return metrics

    def _radius_of_gyration(self, coords: np.ndarray) -> float:
        """
        Calculate radius of gyration.

        Rg = sqrt(sum((r_i - r_mean)^2) / N)

        Lower values indicate more compact structures.
        """
        if len(coords) == 0:
            return 0.0

        center = np.mean(coords, axis=0)
        distances_sq = np.sum((coords - center) ** 2, axis=1)
        return float(np.sqrt(np.mean(distances_sq)))

    def _max_dimension(self, coords: np.ndarray) -> float:
        """Calculate maximum dimension (diameter) of the domain."""
        if len(coords) < 2:
            return 0.0

        from scipy.spatial.distance import pdist

        distances = pdist(coords)
        return float(np.max(distances))

    def _sequence_coverage(self, segments: List[tuple]) -> float:
        """
        Calculate what fraction of the domain's sequence range is covered.

        For continuous domains this is 1.0.
        For discontinuous domains it's less.
        """
        if not segments:
            return 0.0

        total_residues = sum(end - start + 1 for start, end in segments)
        min_res = min(s[0] for s in segments)
        max_res = max(s[1] for s in segments)
        span = max_res - min_res + 1

        return total_residues / span

    def _calculate_quality_score(self, metrics: Dict) -> float:
        """
        Calculate overall quality score (0-100).

        Combines multiple metrics into a single score.
        """
        score = 0.0

        # pLDDT contribution (0-50 points)
        # Scale avg_plddt from 0-100 to 0-50
        score += min(50, metrics["avg_plddt"] / 2)

        # Contact density ratio contribution (0-25 points)
        # Higher is better, capped at ratio of 5
        cdr = metrics["contact_density_ratio"]
        if cdr == float("inf"):
            score += 25
        else:
            score += min(25, cdr * 5)

        # Compactness contribution (0-15 points)
        # Lower Rg is better, but depends on size
        # Use Rg/sqrt(N) as a normalized metric
        rg = metrics["radius_of_gyration"]
        n = metrics["size"]
        if n > 0:
            normalized_rg = rg / np.sqrt(n)
            # Expect normalized Rg around 2-4 Å for compact domains
            if normalized_rg < 3:
                score += 15
            elif normalized_rg < 5:
                score += 10
            elif normalized_rg < 7:
                score += 5

        # Continuity bonus (0-10 points)
        coverage = metrics["sequence_coverage"]
        score += coverage * 10

        return float(score)


class DomainComparer:
    """
    Compare domains from different predictions.
    """

    @staticmethod
    def boundary_similarity(
        domain1_segments: List[tuple],
        domain2_segments: List[tuple],
        tolerance: int = 10,
    ) -> float:
        """
        Calculate boundary similarity between two domains.

        Parameters:
        -----------
        domain1_segments : List[tuple]
            Segments of first domain [(start, end), ...]
        domain2_segments : List[tuple]
            Segments of second domain
        tolerance : int
            Allowed difference in boundary positions

        Returns:
        --------
        Similarity score 0-1
        """
        # Get all boundaries
        boundaries1 = []
        for start, end in domain1_segments:
            boundaries1.extend([start, end])

        boundaries2 = []
        for start, end in domain2_segments:
            boundaries2.extend([start, end])

        if not boundaries1 or not boundaries2:
            return 0.0

        # Calculate how many boundaries match within tolerance
        matches = 0
        used = set()

        for b1 in boundaries1:
            for i, b2 in enumerate(boundaries2):
                if i in used:
                    continue
                if abs(b1 - b2) <= tolerance:
                    matches += 1
                    used.add(i)
                    break

        # Normalize by total boundaries
        total = len(boundaries1) + len(boundaries2)
        return 2 * matches / total if total > 0 else 0.0

    @staticmethod
    def residue_overlap(
        indices1: List[int],
        indices2: List[int],
    ) -> float:
        """
        Calculate Jaccard index of residue overlap.

        Parameters:
        -----------
        indices1 : List[int]
            Residue indices of first domain
        indices2 : List[int]
            Residue indices of second domain

        Returns:
        --------
        Jaccard index (0-1)
        """
        set1 = set(indices1)
        set2 = set(indices2)

        intersection = len(set1 & set2)
        union = len(set1 | set2)

        return intersection / union if union > 0 else 0.0


class SecondaryStructureAnalyzer:
    """
    Analyze secondary structure content of domains.

    Requires DSSP to be installed.
    """

    @staticmethod
    def analyze_ss_content(
        pdb_path: str,
        residue_indices: List[int],
    ) -> Optional[Dict]:
        """
        Analyze secondary structure content using DSSP.

        Parameters:
        -----------
        pdb_path : str
            Path to PDB/CIF file
        residue_indices : List[int]
            Indices of residues in domain

        Returns:
        --------
        Dict with secondary structure fractions, or None if DSSP fails
        """
        try:
            from Bio.PDB import PDBParser, MMCIFParser
            from Bio.PDB.DSSP import DSSP
        except ImportError:
            return None

        try:
            from pathlib import Path
            path = Path(pdb_path)

            if path.suffix.lower() in [".cif", ".mmcif"]:
                parser = MMCIFParser(QUIET=True)
            else:
                parser = PDBParser(QUIET=True)

            structure = parser.get_structure("protein", pdb_path)
            model = structure[0]

            dssp = DSSP(model, pdb_path, dssp="mkdssp")

            # Count SS types for our residues
            ss_counts = {"H": 0, "E": 0, "C": 0}  # Helix, Sheet, Coil
            total = 0

            for key in dssp.keys():
                chain, res_id = key
                res_num = res_id[1]

                # Check if this residue is in our domain
                # (simplified - assumes indices map to residue numbers)
                if res_num - 1 in residue_indices:  # Convert to 0-based
                    ss = dssp[key][2]

                    # Map DSSP codes to simplified categories
                    if ss in ["H", "G", "I"]:  # Helix types
                        ss_counts["H"] += 1
                    elif ss in ["E", "B"]:  # Sheet types
                        ss_counts["E"] += 1
                    else:  # Everything else is coil
                        ss_counts["C"] += 1

                    total += 1

            if total == 0:
                return None

            return {
                "helix_fraction": ss_counts["H"] / total,
                "sheet_fraction": ss_counts["E"] / total,
                "coil_fraction": ss_counts["C"] / total,
                "total_assigned": total,
            }

        except Exception:
            return None
