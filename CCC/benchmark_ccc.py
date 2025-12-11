#!/usr/bin/env python3
"""
CCC Domain Predictor Benchmarking - uses ABC's benchmark framework.

Usage:
    python benchmark_ccc.py ../ABC/TEST/dev_set.txt --output ccc_results/
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from CCC import CCCPredictor

# Import ABC benchmark infrastructure
from ABC.benchmark import (
    parse_dev_set,
    Domain,
    run_benchmark,
    generate_detailed_report,
    generate_plots,
)


def run_ccc_predictions(annotations, verbose=False):
    """
    Run CCC predictions for all proteins.

    Returns dict mapping accession to list of benchmark Domain objects.
    """
    predictor = CCCPredictor()

    predictions = {}

    for ann in annotations:
        if ann.exclude:
            continue

        print(f"Predicting {ann.uniprot_acc}...")
        try:
            result = predictor.predict_from_uniprot(ann.uniprot_acc)

            # Convert CCC domains to benchmark Domain objects
            domains = []
            for d in result.domains:
                # CCC domains are simple start-end, convert to segment list
                segments = [(d.start, d.end)]
                domains.append(Domain(segments=segments))

            print(f"  {ann.uniprot_acc}: {len(result.domains)} domains predicted")
            for i, d in enumerate(result.domains):
                print(f"    Domain {i+1}: {d.start}-{d.end}")

            predictions[ann.uniprot_acc] = domains

        except Exception as e:
            print(f"  Error: {e}")
            if verbose:
                import traceback
                traceback.print_exc()
            predictions[ann.uniprot_acc] = []

    return predictions


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark CCC domain predictions using ABC benchmark framework"
    )
    parser.add_argument(
        "dev_set",
        help="Path to development set file with ground truth annotations"
    )
    parser.add_argument(
        "--output", "-o",
        default="ccc_benchmark_results",
        help="Output directory for results (default: ccc_benchmark_results)"
    )
    parser.add_argument(
        "--iou-threshold",
        type=float,
        default=0.8,
        help="IoU threshold for correct domain match (default: 0.8)"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Print detailed info for debugging"
    )

    args = parser.parse_args()

    # Parse ground truth
    print(f"Loading ground truth from {args.dev_set}...")
    annotations = parse_dev_set(args.dev_set)

    # Count valid proteins
    valid = [a for a in annotations if not a.exclude]
    print(f"Found {len(valid)} proteins for evaluation ({len(annotations) - len(valid)} excluded)")

    # Run CCC predictions
    print("\nRunning CCC predictions...")
    predictions = run_ccc_predictions(annotations, verbose=args.verbose)

    # Run benchmark
    print("\nRunning benchmark evaluation...")
    results = run_benchmark(annotations, predictions, args.iou_threshold)

    # Print summary
    print("\n" + results.summary_brief())

    # Generate reports and plots
    output_dir = Path(args.output)
    generate_detailed_report(results, output_dir)
    generate_plots(results, output_dir)

    print(f"\nResults saved to {output_dir}/")


if __name__ == "__main__":
    main()
