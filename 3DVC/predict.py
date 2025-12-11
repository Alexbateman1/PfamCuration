#!/usr/bin/env python3
"""Predict domains for a UniProt accession using 3DVC."""
import sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from voronoi_predictor import VoronoiPredictor

if len(sys.argv) < 2:
    print("Usage: python predict.py <UniProt_accession>")
    sys.exit(1)

acc = sys.argv[1]
print(f"Predicting domains for {acc}...")
result = VoronoiPredictor().predict_from_uniprot(acc)
print(result.summary())
