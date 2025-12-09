"""
ABC (Alex Bateman Chop) - Geometric Domain Prediction from AlphaFold Structures

A contact graph-based approach to protein domain prediction that uses:
- Spatial proximity (Cα-Cα distances)
- AlphaFold confidence scores (pLDDT)
- Predicted Aligned Error (PAE) when available
- Community detection algorithms (Leiden/Louvain)

Two prediction approaches:
- ABCPredictor (v1): Clustering-based with repeat domain merging
- ABCPredictorV2: Bridge-based with proper nested domain handling
"""

from .abc_predictor import ABCPredictor
from .abc_predictor_v2 import ABCPredictorV2
from .contact_graph import ContactGraphBuilder
from .domain_quality import DomainQualityAssessor
from .ndr_detector import NDRDetector
from .visualize import DomainVisualizer
from .bridge_detector import BridgeDetector, Bridge, NestingNode

__version__ = "0.2.0"
__author__ = "Alex Bateman"

__all__ = [
    "ABCPredictor",
    "ABCPredictorV2",
    "ContactGraphBuilder",
    "DomainQualityAssessor",
    "NDRDetector",
    "DomainVisualizer",
    "BridgeDetector",
    "Bridge",
    "NestingNode",
]
