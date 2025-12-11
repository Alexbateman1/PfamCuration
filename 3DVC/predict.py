#!/usr/bin/env python3
"""Predict domains for a UniProt accession using 3DVC."""
import argparse
import sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from voronoi_predictor import VoronoiPredictor, setup_file_logging

parser = argparse.ArgumentParser(description="3DVC domain prediction")
parser.add_argument("accession", help="UniProt accession")
parser.add_argument("--log", "-l", help="Log file path (default: no file logging)")
parser.add_argument("--chimerax", "-c", help="Output ChimeraX script for visualization")
parser.add_argument("--nu-domains", type=float, default=100.0,
                    help="Domain count penalty (default: 100)")
args = parser.parse_args()

if args.log:
    setup_file_logging(args.log)

result = VoronoiPredictor(nu_domains=args.nu_domains).predict_from_uniprot(args.accession)
print(result.summary())

if args.chimerax:
    result.save_chimerax_script(args.chimerax)
    print(f"\nChimeraX script saved to: {args.chimerax}")
    print(f"Usage: open the structure in ChimeraX, then run: open {args.chimerax}")
