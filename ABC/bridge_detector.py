"""
Bridge Detector for Nested Domain Detection

Detects "bridges" - long-range contacts between residues far apart in sequence
but close in space. These bridges indicate discontinuous domains and reveal
nesting relationships.

For a nested structure like A-B-C-B'-A':
- Bridge A<->A' spans the entire nested region
- Bridge B<->B' is nested within A<->A'
- C is the innermost inserted domain

This module identifies bridges and builds a nesting hierarchy that can be
used to properly assign domains without incorrectly merging nested domains.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple, TYPE_CHECKING

import networkx as nx
import numpy as np

if TYPE_CHECKING:
    from .abc_predictor import Residue


@dataclass(frozen=True)
class Bridge:
    """A long-range contact between sequence-distant residues."""
    start_idx: int  # Index of first residue (earlier in sequence)
    end_idx: int    # Index of second residue (later in sequence)
    start_resnum: int
    end_resnum: int
    sequence_distance: int  # How far apart in sequence
    contact_strength: float  # Sum of contact weights between the regions
    n_contacts: int  # Number of individual contacts supporting this bridge

    @property
    def span(self) -> Tuple[int, int]:
        """Return the sequence range spanned by this bridge."""
        return (self.start_idx, self.end_idx)

    def contains(self, other: 'Bridge') -> bool:
        """Check if this bridge contains another bridge."""
        return self.start_idx < other.start_idx and self.end_idx > other.end_idx

    def overlaps(self, other: 'Bridge') -> bool:
        """Check if this bridge overlaps with another (but neither contains the other)."""
        # Overlapping but not nested
        if self.contains(other) or other.contains(self):
            return False
        return not (self.end_idx < other.start_idx or other.end_idx < self.start_idx)


@dataclass
class DomainRegion:
    """A region identified as belonging to a domain."""
    indices: List[int]  # Residue indices
    resnums: List[int]  # Residue numbers
    region_type: str    # 'continuous', 'outer_flank', 'inserted'
    parent_bridge: Optional[Bridge] = None  # The bridge this region is part of (for outer domains)
    nesting_level: int = 0  # 0 = outermost, higher = more nested

    @property
    def start(self) -> int:
        return min(self.resnums)

    @property
    def end(self) -> int:
        return max(self.resnums)

    @property
    def size(self) -> int:
        return len(self.indices)


@dataclass
class NestingNode:
    """A node in the nesting hierarchy tree."""
    bridge: Optional[Bridge]  # None for root or continuous domains
    children: List['NestingNode'] = field(default_factory=list)
    domain_regions: List[DomainRegion] = field(default_factory=list)
    level: int = 0

    def add_child(self, child: 'NestingNode'):
        child.level = self.level + 1
        self.children.append(child)


class BridgeDetector:
    """
    Detects bridges and builds nesting hierarchy.

    Parameters:
    -----------
    min_sequence_separation : int
        Minimum sequence distance to consider a contact as a "bridge" (default: 50)
    min_bridge_contacts : int
        Minimum number of contacts to support a bridge (default: 3)
    min_bridge_strength : float
        Minimum total contact weight for a bridge (default: 1.0)
    bridge_region_size : int
        Size of regions at each end of bridge to check for contacts (default: 10)
    """

    def __init__(
        self,
        min_sequence_separation: int = 50,
        min_bridge_contacts: int = 3,
        min_bridge_strength: float = 1.0,
        bridge_region_size: int = 10,
    ):
        self.min_sequence_separation = min_sequence_separation
        self.min_bridge_contacts = min_bridge_contacts
        self.min_bridge_strength = min_bridge_strength
        self.bridge_region_size = bridge_region_size

    def detect_bridges(
        self,
        residues: List['Residue'],
        graph: nx.Graph,
    ) -> List[Bridge]:
        """
        Detect all significant bridges in the contact graph.

        A bridge is a cluster of contacts between two sequence-distant regions.
        We look for regions that have multiple strong contacts to each other
        despite being far apart in sequence.

        Parameters:
        -----------
        residues : List[Residue]
            All residues in the protein
        graph : nx.Graph
            Contact graph with weighted edges

        Returns:
        --------
        List[Bridge] sorted by span (largest first)
        """
        n = len(residues)
        bridges = []

        # Scan for potential bridge regions
        # We look at windows of residues and check for long-range contacts
        window_size = self.bridge_region_size

        # For efficiency, first identify all long-range contacts
        long_range_contacts = []
        for i, j, data in graph.edges(data=True):
            seq_dist = abs(i - j)
            if seq_dist >= self.min_sequence_separation:
                long_range_contacts.append((i, j, data.get('weight', 1.0)))

        if not long_range_contacts:
            return []

        # Group contacts into potential bridges
        # A bridge connects two regions - we cluster nearby contacts
        contact_by_region = {}

        for i, j, weight in long_range_contacts:
            # Ensure i < j
            if i > j:
                i, j = j, i

            # Round to region boundaries
            region_i = (i // window_size) * window_size
            region_j = (j // window_size) * window_size

            key = (region_i, region_j)
            if key not in contact_by_region:
                contact_by_region[key] = {'contacts': [], 'total_weight': 0}

            contact_by_region[key]['contacts'].append((i, j, weight))
            contact_by_region[key]['total_weight'] += weight

        # Filter and create Bridge objects
        for (region_i, region_j), data in contact_by_region.items():
            n_contacts = len(data['contacts'])
            total_weight = data['total_weight']

            if n_contacts >= self.min_bridge_contacts and total_weight >= self.min_bridge_strength:
                # Find the actual residue range involved
                all_i = [c[0] for c in data['contacts']]
                all_j = [c[1] for c in data['contacts']]

                # Bridge spans from the region around min(all_i) to the region around max(all_j)
                start_idx = min(all_i)
                end_idx = max(all_j)

                bridges.append(Bridge(
                    start_idx=start_idx,
                    end_idx=end_idx,
                    start_resnum=residues[start_idx].resnum,
                    end_resnum=residues[end_idx].resnum,
                    sequence_distance=end_idx - start_idx,
                    contact_strength=total_weight,
                    n_contacts=n_contacts,
                ))

        # Sort by span (largest first) - outer bridges before inner
        bridges.sort(key=lambda b: b.sequence_distance, reverse=True)

        # Merge overlapping bridges that likely represent the same domain closure
        bridges = self._merge_overlapping_bridges(bridges, residues, graph)

        return bridges

    def _merge_overlapping_bridges(
        self,
        bridges: List[Bridge],
        residues: List['Residue'],
        graph: nx.Graph,
    ) -> List[Bridge]:
        """Merge bridges that overlap and likely represent the same domain."""
        if len(bridges) <= 1:
            return bridges

        merged = []
        used = set()

        for i, bridge_i in enumerate(bridges):
            if i in used:
                continue

            # Find all bridges that should merge with this one
            to_merge = [bridge_i]

            for j, bridge_j in enumerate(bridges):
                if j <= i or j in used:
                    continue

                # Merge if one contains the other OR they significantly overlap
                if bridge_i.contains(bridge_j):
                    # bridge_j is nested inside bridge_i - don't merge, it's a real nesting
                    continue
                elif bridge_j.contains(bridge_i):
                    # Shouldn't happen since we sorted by span
                    continue
                elif bridge_i.overlaps(bridge_j):
                    # Overlapping bridges - might be same domain
                    # Check if their endpoints are close
                    start_close = abs(bridge_i.start_idx - bridge_j.start_idx) < self.bridge_region_size * 2
                    end_close = abs(bridge_i.end_idx - bridge_j.end_idx) < self.bridge_region_size * 2

                    if start_close or end_close:
                        to_merge.append(bridge_j)
                        used.add(j)

            # Create merged bridge
            if len(to_merge) == 1:
                merged.append(bridge_i)
            else:
                # Combine all the bridges
                all_starts = [b.start_idx for b in to_merge]
                all_ends = [b.end_idx for b in to_merge]
                total_strength = sum(b.contact_strength for b in to_merge)
                total_contacts = sum(b.n_contacts for b in to_merge)

                new_start = min(all_starts)
                new_end = max(all_ends)

                merged.append(Bridge(
                    start_idx=new_start,
                    end_idx=new_end,
                    start_resnum=residues[new_start].resnum,
                    end_resnum=residues[new_end].resnum,
                    sequence_distance=new_end - new_start,
                    contact_strength=total_strength,
                    n_contacts=total_contacts,
                ))

        # Re-sort by span
        merged.sort(key=lambda b: b.sequence_distance, reverse=True)
        return merged

    def build_nesting_hierarchy(
        self,
        bridges: List[Bridge],
        n_residues: int,
    ) -> NestingNode:
        """
        Build a tree structure representing the nesting hierarchy.

        Parameters:
        -----------
        bridges : List[Bridge]
            Bridges sorted by span (largest first)
        n_residues : int
            Total number of residues

        Returns:
        --------
        Root NestingNode containing the hierarchy
        """
        root = NestingNode(bridge=None, level=0)

        if not bridges:
            return root

        # Process bridges from largest to smallest
        # Each bridge either becomes a child of root or a child of a containing bridge
        nodes = {None: root}  # Map bridge to its node

        for bridge in bridges:
            # Find the smallest containing bridge (if any)
            parent_bridge = None
            for potential_parent in bridges:
                if potential_parent is bridge:
                    continue
                if potential_parent.contains(bridge):
                    # Check if this is the smallest container
                    if parent_bridge is None or potential_parent.sequence_distance < parent_bridge.sequence_distance:
                        parent_bridge = potential_parent

            # Create node for this bridge
            node = NestingNode(bridge=bridge)
            nodes[bridge] = node

            # Add to parent
            parent_node = nodes.get(parent_bridge, root)
            parent_node.add_child(node)

        return root

    def identify_domain_regions(
        self,
        residues: List['Residue'],
        bridges: List[Bridge],
        hierarchy: NestingNode,
    ) -> List[DomainRegion]:
        """
        Identify all domain regions based on bridge analysis.

        This assigns each residue to a domain region, respecting the
        nesting hierarchy.

        Parameters:
        -----------
        residues : List[Residue]
            All residues
        bridges : List[Bridge]
            Detected bridges
        hierarchy : NestingNode
            Nesting hierarchy from build_nesting_hierarchy

        Returns:
        --------
        List of DomainRegion objects
        """
        n = len(residues)
        regions = []

        # Track which residues are assigned
        assigned = set()

        # Process hierarchy level by level
        def process_node(node: NestingNode, parent_span: Optional[Tuple[int, int]] = None):
            if node.bridge is not None:
                bridge = node.bridge

                # The bridge defines a discontinuous outer domain
                # Its "flanks" are the regions at the start and end

                # Determine flank boundaries
                # Start flank: from bridge.start_idx to the start of any child bridge (or inner region)
                # End flank: from end of any child bridge to bridge.end_idx

                child_spans = []
                for child in node.children:
                    if child.bridge:
                        child_spans.append((child.bridge.start_idx, child.bridge.end_idx))

                # Sort child spans by start position
                child_spans.sort()

                # Build flank regions for this bridge
                flank_indices = []
                flank_resnums = []

                # Start with the region before first child
                current_pos = bridge.start_idx

                for child_start, child_end in child_spans:
                    # Add region from current_pos to just before child_start
                    for idx in range(current_pos, child_start):
                        if idx not in assigned:
                            flank_indices.append(idx)
                            flank_resnums.append(residues[idx].resnum)
                            assigned.add(idx)
                    current_pos = child_end + 1

                # Add remaining region after last child (or all if no children)
                for idx in range(current_pos, bridge.end_idx + 1):
                    if idx not in assigned:
                        flank_indices.append(idx)
                        flank_resnums.append(residues[idx].resnum)
                        assigned.add(idx)

                if flank_indices:
                    regions.append(DomainRegion(
                        indices=flank_indices,
                        resnums=flank_resnums,
                        region_type='outer_flank',
                        parent_bridge=bridge,
                        nesting_level=node.level,
                    ))

            # Process children
            for child in node.children:
                process_node(child, node.bridge.span if node.bridge else None)

        # Start processing from root
        process_node(hierarchy)

        # Any unassigned residues are continuous regions (not involved in bridges)
        unassigned = [i for i in range(n) if i not in assigned]

        if unassigned:
            # Group into contiguous regions
            current_region = [unassigned[0]]

            for idx in unassigned[1:]:
                if idx - current_region[-1] <= 3:  # Allow small gaps
                    current_region.append(idx)
                else:
                    # Save current region and start new one
                    if len(current_region) >= 5:  # Minimum region size
                        regions.append(DomainRegion(
                            indices=current_region,
                            resnums=[residues[i].resnum for i in current_region],
                            region_type='continuous',
                            nesting_level=0,
                        ))
                    current_region = [idx]

            # Handle last region
            if len(current_region) >= 5:
                regions.append(DomainRegion(
                    indices=current_region,
                    resnums=[residues[i].resnum for i in current_region],
                    region_type='continuous',
                    nesting_level=0,
                ))

        return regions

    def get_nesting_constraints(
        self,
        bridges: List[Bridge],
    ) -> List[Tuple[int, int, int, int]]:
        """
        Get constraints for domain assignment based on nesting.

        Returns a list of (outer_start, inner_start, inner_end, outer_end) tuples.
        Each tuple indicates that the inner region should not be merged with
        the outer region.

        Parameters:
        -----------
        bridges : List[Bridge]
            Detected bridges

        Returns:
        --------
        List of nesting constraints
        """
        constraints = []

        for i, outer in enumerate(bridges):
            for inner in bridges[i+1:]:
                if outer.contains(inner):
                    constraints.append((
                        outer.start_idx,
                        inner.start_idx,
                        inner.end_idx,
                        outer.end_idx,
                    ))

        return constraints


def find_bridges_and_hierarchy(
    residues: List['Residue'],
    graph: nx.Graph,
    min_sequence_separation: int = 50,
    min_bridge_contacts: int = 3,
) -> Tuple[List[Bridge], NestingNode, List[DomainRegion]]:
    """
    Convenience function to run the full bridge detection pipeline.

    Parameters:
    -----------
    residues : List[Residue]
        All residues in the protein
    graph : nx.Graph
        Contact graph
    min_sequence_separation : int
        Minimum sequence distance for bridges (default: 50)
    min_bridge_contacts : int
        Minimum contacts to support a bridge (default: 3)

    Returns:
    --------
    Tuple of (bridges, hierarchy, domain_regions)
    """
    detector = BridgeDetector(
        min_sequence_separation=min_sequence_separation,
        min_bridge_contacts=min_bridge_contacts,
    )

    bridges = detector.detect_bridges(residues, graph)
    hierarchy = detector.build_nesting_hierarchy(bridges, len(residues))
    regions = detector.identify_domain_regions(residues, bridges, hierarchy)

    return bridges, hierarchy, regions
