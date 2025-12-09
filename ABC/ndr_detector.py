"""
NDR (Non-Domain Region) Detector

Identifies regions that are not part of any domain:
- Low pLDDT regions (disordered/flexible)
- Isolated residues with few contacts
- Linker regions between domains
"""

from typing import TYPE_CHECKING, Dict, List, Set

import networkx as nx
import numpy as np

if TYPE_CHECKING:
    from .abc_predictor import Residue


class NDRDetector:
    """
    Detects Non-Domain Regions (NDRs).

    NDRs include:
    1. Low confidence regions (pLDDT < cutoff)
    2. Isolated residues (few contacts)
    3. Small clusters below minimum domain size
    4. Unassigned residues after domain assignment

    Parameters:
    -----------
    plddt_cutoff : float
        pLDDT below which residues are considered NDR candidates (default: 70)
    min_contacts : int
        Minimum contacts for a residue to not be isolated (default: 3)
    isolation_window : int
        Window size for checking local isolation (default: 10)
    """

    def __init__(
        self,
        plddt_cutoff: float = 70.0,
        min_contacts: int = 3,
        isolation_window: int = 10,
    ):
        self.plddt_cutoff = plddt_cutoff
        self.min_contacts = min_contacts
        self.isolation_window = isolation_window

    def detect(
        self,
        residues: List["Residue"],
        cluster_assignments: Dict[int, int],
        graph: nx.Graph,
    ) -> Set[int]:
        """
        Detect NDR residues.

        Parameters:
        -----------
        residues : List[Residue]
            All residues in protein
        cluster_assignments : Dict[int, int]
            Mapping of residue index to cluster ID
        graph : nx.Graph
            Contact graph

        Returns:
        --------
        Set of residue indices identified as NDR
        """
        ndr_indices = set()

        # 1. Low pLDDT residues
        low_plddt = self._detect_low_plddt(residues)
        ndr_indices.update(low_plddt)

        # 2. Isolated residues
        isolated = self._detect_isolated(residues, graph)
        ndr_indices.update(isolated)

        # 3. Terminal regions with low confidence
        terminal = self._detect_terminal_disorder(residues)
        ndr_indices.update(terminal)

        return ndr_indices

    def _detect_low_plddt(self, residues: List["Residue"]) -> Set[int]:
        """
        Identify residues with low pLDDT scores.

        Uses a sliding window approach to identify regions
        (not just individual residues) with low confidence.
        """
        n = len(residues)
        ndr = set()

        # Sliding window for low pLDDT regions
        window_size = 5
        for i in range(n - window_size + 1):
            window_plddt = [residues[j].plddt for j in range(i, i + window_size)]
            avg_plddt = np.mean(window_plddt)

            if avg_plddt < self.plddt_cutoff:
                # Mark entire window as NDR
                for j in range(i, i + window_size):
                    ndr.add(j)

        # Also mark individual very low pLDDT residues
        for i, res in enumerate(residues):
            if res.plddt < self.plddt_cutoff - 10:  # Extra stringent cutoff
                ndr.add(i)

        return ndr

    def _detect_isolated(
        self,
        residues: List["Residue"],
        graph: nx.Graph,
    ) -> Set[int]:
        """
        Identify isolated residues with few long-range contacts.
        """
        n = len(residues)
        isolated = set()

        for i in range(n):
            if i not in graph:
                isolated.add(i)
                continue

            # Count non-local contacts
            nonlocal_contacts = 0
            for neighbor in graph.neighbors(i):
                # Only count contacts to residues far in sequence
                if abs(neighbor - i) >= self.isolation_window:
                    nonlocal_contacts += 1

            if nonlocal_contacts < self.min_contacts:
                isolated.add(i)

        return isolated

    def _detect_terminal_disorder(self, residues: List["Residue"]) -> Set[int]:
        """
        Identify disordered N- and C-terminal regions.

        Terminal regions often have lower pLDDT and should be marked as NDR
        if they don't fold into the protein core.
        """
        n = len(residues)
        ndr = set()

        # Check N-terminus
        for i in range(min(30, n)):
            if residues[i].plddt < self.plddt_cutoff:
                ndr.add(i)
            else:
                break  # Stop when we hit good confidence

        # Check C-terminus
        for i in range(n - 1, max(n - 31, -1), -1):
            if residues[i].plddt < self.plddt_cutoff:
                ndr.add(i)
            else:
                break

        return ndr

    def classify_ndr_type(
        self,
        residues: List["Residue"],
        ndr_indices: List[int],
        graph: nx.Graph,
    ) -> str:
        """
        Classify the type of NDR region.

        Returns:
        --------
        str: One of:
            - "disordered": Low pLDDT, likely intrinsically disordered
            - "linker": Between domains, may have some structure
            - "terminal": N or C-terminal extension
            - "isolated": Few contacts, possibly an artifact
        """
        if not ndr_indices:
            return "unknown"

        n_residues = len(residues)
        avg_plddt = np.mean([residues[i].plddt for i in ndr_indices])

        # Check if terminal
        min_idx = min(ndr_indices)
        max_idx = max(ndr_indices)

        if min_idx < 10:
            return "terminal"
        if max_idx > n_residues - 10:
            return "terminal"

        # Check pLDDT
        if avg_plddt < 50:
            return "disordered"

        # Check contacts
        total_contacts = 0
        for i in ndr_indices:
            if i in graph:
                total_contacts += graph.degree(i)

        avg_contacts = total_contacts / len(ndr_indices) if ndr_indices else 0

        if avg_contacts < 2:
            return "isolated"

        return "linker"


class LinkerAnalyzer:
    """
    Specialized analysis of linker regions between domains.
    """

    @staticmethod
    def identify_linkers(
        residues: List["Residue"],
        domain_assignments: Dict[int, int],
        min_linker_length: int = 3,
    ) -> List[Dict]:
        """
        Identify linker regions between assigned domains.

        Parameters:
        -----------
        residues : List[Residue]
            All residues
        domain_assignments : Dict[int, int]
            Mapping of residue index to domain ID (-1 for unassigned)
        min_linker_length : int
            Minimum length to be considered a linker

        Returns:
        --------
        List of linker dictionaries with start, end, length, connecting domains
        """
        n = len(residues)
        linkers = []

        # Find gaps between domains
        in_linker = False
        linker_start = None
        prev_domain = None

        for i in range(n):
            domain = domain_assignments.get(i, -1)

            if domain == -1:  # Not assigned to a domain
                if not in_linker:
                    in_linker = True
                    linker_start = i
            else:
                if in_linker:
                    # End of linker
                    linker_end = i - 1
                    length = linker_end - linker_start + 1

                    if length >= min_linker_length:
                        # Calculate properties
                        linker_residues = list(range(linker_start, linker_end + 1))
                        avg_plddt = np.mean([residues[j].plddt for j in linker_residues])

                        linkers.append({
                            "start": residues[linker_start].resnum,
                            "end": residues[linker_end].resnum,
                            "length": length,
                            "from_domain": prev_domain,
                            "to_domain": domain,
                            "avg_plddt": avg_plddt,
                            "residue_indices": linker_residues,
                        })

                    in_linker = False

                prev_domain = domain

        # Handle trailing linker
        if in_linker:
            linker_end = n - 1
            length = linker_end - linker_start + 1

            if length >= min_linker_length:
                linker_residues = list(range(linker_start, linker_end + 1))
                avg_plddt = np.mean([residues[j].plddt for j in linker_residues])

                linkers.append({
                    "start": residues[linker_start].resnum,
                    "end": residues[linker_end].resnum,
                    "length": length,
                    "from_domain": prev_domain,
                    "to_domain": None,
                    "avg_plddt": avg_plddt,
                    "residue_indices": linker_residues,
                })

        return linkers

    @staticmethod
    def linker_flexibility_score(
        residues: List["Residue"],
        linker_indices: List[int],
    ) -> float:
        """
        Estimate linker flexibility based on pLDDT and amino acid composition.

        Returns:
        --------
        Flexibility score 0-1 (higher = more flexible)
        """
        if not linker_indices:
            return 0.0

        # pLDDT-based flexibility (lower pLDDT = more flexible)
        avg_plddt = np.mean([residues[i].plddt for i in linker_indices])
        plddt_flexibility = 1.0 - (avg_plddt / 100.0)

        # Amino acid based flexibility
        # Flexible residues: G, P, S, D, E, N, Q
        flexible_aa = {"GLY", "PRO", "SER", "ASP", "GLU", "ASN", "GLN"}
        n_flexible = sum(1 for i in linker_indices if residues[i].resname in flexible_aa)
        aa_flexibility = n_flexible / len(linker_indices)

        # Combined score (weighted average)
        return 0.7 * plddt_flexibility + 0.3 * aa_flexibility
