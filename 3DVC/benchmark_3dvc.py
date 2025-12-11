#!/usr/bin/env python3
"""
3DVC Domain Predictor Benchmarking - uses ABC's benchmark framework.

Usage:
    python benchmark_3dvc.py ../ABC/TEST/dev_set.txt --output 3dvc_results/
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# Import VoronoiPredictor - add current directory to path for direct import
sys.path.insert(0, str(Path(__file__).resolve().parent))
from voronoi_predictor import VoronoiPredictor, setup_file_logging

# Import ABC benchmark infrastructure
from ABC.benchmark import (
    parse_dev_set,
    Domain,
    run_benchmark,
    generate_detailed_report,
    generate_plots,
)


def run_3dvc_predictions(annotations, nu_domains=100.0, verbose=False):
    """
    Run 3DVC predictions for all proteins.

    Returns dict mapping accession to list of benchmark Domain objects.
    """
    predictor = VoronoiPredictor(nu_domains=nu_domains)

    predictions = {}

    for ann in annotations:
        if ann.exclude:
            continue

        print(f"Predicting {ann.uniprot_acc}...")
        try:
            result = predictor.predict_from_uniprot(ann.uniprot_acc)

            # Convert 3DVC domains to benchmark Domain objects
            domains = []
            for d in result.domains:
                domains.append(Domain(segments=d.segments))

            print(f"  {ann.uniprot_acc}: {len(result.domains)} domains predicted")
            for i, d in enumerate(result.domains):
                print(f"    Domain {i+1}: {d.to_chopping_string()}")

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
        description="Benchmark 3DVC domain predictions using ABC benchmark framework"
    )
    parser.add_argument(
        "dev_set",
        help="Path to development set file with ground truth annotations"
    )
    parser.add_argument(
        "--output", "-o",
        default="3dvc_benchmark_results",
        help="Output directory for results (default: 3dvc_benchmark_results)"
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
    parser.add_argument(
        "--log", "-l",
        default="3dvc_benchmark.log",
        help="Log file for detailed scores (default: 3dvc_benchmark.log)"
    )
    parser.add_argument(
        "--nu-domains",
        type=float,
        default=100.0,
        help="Domain count penalty - higher=fewer domains, lower=more splits (default: 100)"
    )

    args = parser.parse_args()

    # Set up file logging
    setup_file_logging(args.log)
    print(f"Logging to {args.log}")

    # Parse ground truth
    print(f"Loading ground truth from {args.dev_set}...")
    annotations = parse_dev_set(args.dev_set)

    # Count valid proteins
    valid = [a for a in annotations if not a.exclude]
    print(f"Found {len(valid)} proteins for evaluation ({len(annotations) - len(valid)} excluded)")

    # Run 3DVC predictions
    print(f"\nRunning 3DVC predictions (nu_domains={args.nu_domains})...")
    predictions = run_3dvc_predictions(annotations, nu_domains=args.nu_domains, verbose=args.verbose)

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
