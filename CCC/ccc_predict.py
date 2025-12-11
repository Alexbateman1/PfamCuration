#!/usr/bin/env python3
"""
CCC Domain Predictor - Standalone script

Claude Code Chop: Spectral + DP approach to domain prediction.

Usage:
    python ccc_predict.py Q9NQT6
    python ccc_predict.py Q9NQT6 --penalty 20
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from CCC import CCCPredictor


def main():
    parser = argparse.ArgumentParser(
        description="CCC Domain Predictor - Spectral + DP approach"
    )
    parser.add_argument("uniprot", help="UniProt accession (e.g., Q9NQT6)")
    parser.add_argument("--distance", "-d", type=float, default=10.0,
                        help="Contact distance threshold (default: 10)")
    parser.add_argument("--min-size", type=int, default=30,
                        help="Minimum domain size (default: 30)")
    parser.add_argument("--penalty", "-p", type=float, default=2.0,
                        help="Domain penalty for MDL (default: 2). "
                             "Higher = fewer domains")
    parser.add_argument("--no-pae", action="store_true",
                        help="Don't use PAE for scoring")
    parser.add_argument("--cache-dir", type=str, default=None,
                        help="Cache directory (default: ./.ccc_cache)")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Verbose output")

    args = parser.parse_args()

    # Set logging level
    import logging
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Create predictor
    predictor = CCCPredictor(
        distance_threshold=args.distance,
        min_domain_size=args.min_size,
        domain_penalty=args.penalty,
        use_pae=not args.no_pae,
        cache_dir=args.cache_dir,
    )

    # Run prediction
    try:
        result = predictor.predict_from_uniprot(args.uniprot)
        print("\n" + result.summary())

        # Print additional details
        print("\nSpectral Analysis:")
        print(f"  Algebraic connectivity: {result.spectral_info.get('algebraic_connectivity', 'N/A'):.4f}")
        if result.spectral_info.get('candidate_boundaries'):
            print(f"  Candidate boundaries: {result.spectral_info['candidate_boundaries']}")

        print("\nDomain Details:")
        for d in result.domains:
            print(f"  {d.domain_id}: {d.start}-{d.end}")
            for k, v in d.metrics.items():
                if isinstance(v, float):
                    print(f"      {k}: {v:.2f}")
                else:
                    print(f"      {k}: {v}")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
