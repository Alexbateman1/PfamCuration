"""
CCC - Claude Code Chop Domain Predictor

A spectral graph + dynamic programming approach to protein domain prediction.
Uses the Fiedler vector to identify natural partition points and DP to find
the optimal domain assignment along the sequence.
"""

from .ccc_predictor import CCCPredictor

__all__ = ["CCCPredictor"]
__version__ = "0.1.0"
