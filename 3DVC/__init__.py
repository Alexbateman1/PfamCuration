"""
3DVC - 3D Voronoi Chop Domain Predictor

A spatial partitioning approach to protein domain prediction.
Uses Voronoi tessellation in 3D to partition residues into domains,
optimizing seed positions to maximize domain quality.

Key differences from clustering approaches:
- Partitions SPACE, not just residue connections
- Boundaries emerge from spatial optimization
- Handles discontinuous domains naturally (A and A' near same seed)
- No arbitrary distance thresholds for domain membership
"""

from .voronoi_predictor import VoronoiPredictor

__all__ = ["VoronoiPredictor"]
__version__ = "0.1.0"
