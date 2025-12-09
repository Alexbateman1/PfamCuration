#!/usr/bin/env python3
"""
ABC Domain Predictor - Standalone script

Run from anywhere:
    python /path/to/abc_predict.py P12931
    python /path/to/abc_predict.py P12931 --resolution 0.5
"""

import argparse
import sys
from pathlib import Path

# Add the parent directory to path so we can import ABC
script_dir = Path(__file__).parent.parent
sys.path.insert(0, str(script_dir))

from ABC import ABCPredictor


def main():
    parser = argparse.ArgumentParser(
        description="ABC Domain Predictor - Predict domains from AlphaFold structures"
    )
    parser.add_argument("uniprot", help="UniProt accession (e.g., P12931)")
    parser.add_argument("--resolution", "-r", type=float, default=0.5,
                        help="Clustering resolution (default: 0.5)")
    parser.add_argument("--distance", "-d", type=float, default=10.0,
                        help="Contact distance threshold in Angstroms (default: 10)")
    parser.add_argument("--min-size", type=int, default=30,
                        help="Minimum domain size (default: 30)")
    parser.add_argument("--output", "-o", type=str, default=None,
                        help="Output filename prefix (default: UniProt accession)")
    parser.add_argument("--debug", action="store_true",
                        help="Show debug information about clustering and NDR detection")
    parser.add_argument("--debug-region", type=str, default=None,
                        help="Focus debug output on specific region (e.g., '1014-1085')")

    args = parser.parse_args()

    # Initialize predictor
    predictor = ABCPredictor(
        distance_threshold=args.distance,
        min_domain_size=args.min_size,
        resolution=args.resolution,
    )

    # Run prediction (with debug if requested)
    if args.debug:
        prediction = predictor.predict_from_uniprot_debug(args.uniprot, args.debug_region)
    else:
        prediction = predictor.predict_from_uniprot(args.uniprot)

    # Print summary
    print(prediction.summary())

    # Output files
    output_prefix = args.output or args.uniprot
    predictor.visualize(prediction, output_prefix)

    print(f"\nOutput files:")
    print(f"  {output_prefix}.cxc  - ChimeraX commands")
    print(f"  {output_prefix}.pml  - PyMOL script")
    print(f"  {output_prefix}.png  - Architecture diagram")
    print(f"  {output_prefix}.html - Detailed report")


if __name__ == "__main__":
    main()
