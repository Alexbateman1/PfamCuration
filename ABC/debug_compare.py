#!/usr/bin/env python3
"""
Debug comparison tool - compares CLI and benchmark prediction paths.

This script runs the same prediction through both code paths and
compares the results to identify any discrepancies.

Usage:
    python debug_compare.py P12931 [--use-dssp]
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


def main():
    parser = argparse.ArgumentParser(description="Compare CLI and benchmark prediction paths")
    parser.add_argument("uniprot", help="UniProt accession to test")
    parser.add_argument("--use-dssp", action="store_true", help="Enable DSSP")
    args = parser.parse_args()

    print("=" * 70)
    print(f"DEBUG COMPARISON: {args.uniprot}")
    print(f"DSSP enabled: {args.use_dssp}")
    print("=" * 70)

    # Import the predictor
    from ABC.abc_predictor import ABCPredictor

    # Create predictor with same settings as benchmark
    print("\n--- Creating ABCPredictor (benchmark settings) ---")
    predictor = ABCPredictor(
        resolution=0.5,
        use_dssp=args.use_dssp,
    )

    print(f"Parameters:")
    for key, value in predictor.parameters.items():
        print(f"  {key}: {value}")
    print(f"  cache_dir: {predictor.cache_dir}")

    # Run prediction
    print(f"\n--- Running prediction for {args.uniprot} ---")
    result = predictor.predict_from_uniprot(args.uniprot)

    # Print full details
    print(f"\n--- Results ---")
    print(f"Domains: {len(result.domains)}")

    for i, d in enumerate(result.domains):
        print(f"\nDomain {i+1}:")
        print(f"  Chopping: {d.to_chopping_string()}")
        print(f"  Segments: {d.segments}")
        print(f"  Size: {d.size}")
        print(f"  Discontinuous: {d.is_discontinuous}")
        print(f"  Residue indices (first 10): {d.residue_indices[:10]}...")
        print(f"  Residue indices (last 10): ...{d.residue_indices[-10:]}")
        if d.quality_metrics:
            print(f"  Quality metrics:")
            for key, value in d.quality_metrics.items():
                if isinstance(value, float):
                    print(f"    {key}: {value:.3f}")
                else:
                    print(f"    {key}: {value}")

    print(f"\nNDR Regions: {len(result.ndr_regions)}")
    for ndr in result.ndr_regions:
        print(f"  {ndr.start}-{ndr.end} ({ndr.reason}, pLDDT={ndr.avg_plddt:.1f})")

    # Now convert to benchmark Domain format and verify
    print("\n--- Conversion to benchmark Domain format ---")
    from ABC.benchmark import Domain as BenchmarkDomain

    benchmark_domains = []
    for d in result.domains:
        bd = BenchmarkDomain(segments=d.segments)
        benchmark_domains.append(bd)
        print(f"ABC Domain: {d.to_chopping_string()} -> Benchmark Domain: {bd}")

    # Verify segments are identical
    print("\n--- Segment verification ---")
    for i, (orig, bench) in enumerate(zip(result.domains, benchmark_domains)):
        if orig.segments == bench.segments:
            print(f"Domain {i+1}: OK (segments match)")
        else:
            print(f"Domain {i+1}: MISMATCH!")
            print(f"  Original: {orig.segments}")
            print(f"  Benchmark: {bench.segments}")


if __name__ == "__main__":
    main()
