#!/usr/bin/env python3
"""
ABC Domain Predictor V2 - Bridge-Based Standalone Script

Uses bridge detection for proper handling of nested/inserted domains.

Run from anywhere:
    python /path/to/abc_predict_v2.py P12931
    python /path/to/abc_predict_v2.py A9WIN4 --min-bridge-sep 40
"""

import argparse
import sys
from pathlib import Path

# Add the parent directory to path so we can import ABC
script_dir = Path(__file__).parent.parent
sys.path.insert(0, str(script_dir))

from ABC.abc_predictor_v2 import ABCPredictorV2


def main():
    parser = argparse.ArgumentParser(
        description="ABC Domain Predictor V2 - Bridge-based nested domain detection"
    )
    parser.add_argument("uniprot", help="UniProt accession (e.g., P12931)")
    parser.add_argument("--distance", "-d", type=float, default=10.0,
                        help="Contact distance threshold in Angstroms (default: 10)")
    parser.add_argument("--min-size", type=int, default=30,
                        help="Minimum domain size (default: 30)")
    parser.add_argument("--min-contact-ratio", type=float, default=1.0,
                        help="Minimum internal/external contact ratio (default: 1.0)")
    parser.add_argument("--min-bridge-sep", type=int, default=50,
                        help="Minimum sequence separation for bridges (default: 50)")
    parser.add_argument("--min-bridge-contacts", type=int, default=3,
                        help="Minimum contacts to support a bridge (default: 3)")
    parser.add_argument("--ndr-cutoff", type=float, default=70.0,
                        help="pLDDT cutoff for NDR (default: 70)")
    parser.add_argument("--output", "-o", type=str, default=None,
                        help="Output filename prefix (default: UniProt accession)")

    args = parser.parse_args()

    # Initialize predictor
    predictor = ABCPredictorV2(
        distance_threshold=args.distance,
        min_domain_size=args.min_size,
        min_contact_ratio=args.min_contact_ratio,
        ndr_plddt_cutoff=args.ndr_cutoff,
        min_bridge_separation=args.min_bridge_sep,
        min_bridge_contacts=args.min_bridge_contacts,
    )

    # Run prediction
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
