"""
ABC (Alex Bateman Chop) - Geometric Domain Prediction from AlphaFold Structures

A contact graph-based approach to protein domain prediction that uses:
- Spatial proximity (Cα-Cα distances)
- AlphaFold confidence scores (pLDDT)
- Predicted Aligned Error (PAE) when available
- Community detection algorithms (Leiden/Louvain)
"""

from .abc_predictor import ABCPredictor
from .contact_graph import ContactGraphBuilder
from .domain_quality import DomainQualityAssessor
from .ndr_detector import NDRDetector
from .visualize import DomainVisualizer

__version__ = "0.1.0"
__author__ = "Alex Bateman"

__all__ = [
    "ABCPredictor",
    "ContactGraphBuilder",
    "DomainQualityAssessor",
    "NDRDetector",
    "DomainVisualizer",
]
