#!/usr/bin/env python3
"""Predict domains for a UniProt accession using 3DVC."""
import argparse
import sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from voronoi_predictor import VoronoiPredictor, setup_file_logging

parser = argparse.ArgumentParser(description="3DVC domain prediction")
parser.add_argument("accession", help="UniProt accession")
parser.add_argument("--log", "-l", help="Log file path (default: no file logging)")
args = parser.parse_args()

if args.log:
    setup_file_logging(args.log)

result = VoronoiPredictor().predict_from_uniprot(args.accession)
print(result.summary())
