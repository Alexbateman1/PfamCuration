#!/usr/bin/env python3
"""
Simple Voronoi Domain Predictor Benchmarking - uses ABC's benchmark framework.

This benchmarks the hill-climbing Voronoi approach from test_simple_voronoi.py.

Usage:
    python benchmark_simple_voronoi.py ../ABC/TEST/dev_set.txt --output simple_voronoi_results/
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# Import the simple Voronoi predictor functions
sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_simple_voronoi import (
    load_structure,
    filter_structured,
    predict_domains,
    assignments_to_domains,
)

# Import ABC benchmark infrastructure
from ABC.benchmark import (
    parse_dev_set,
    Domain,
    run_benchmark,
    generate_detailed_report,
    generate_plots,
    BenchmarkResults,
)


def run_simple_voronoi_predictions(
    annotations,
    domain_bonus=105,
    crossing_penalty=100,
    intrusion_penalty=10,
    min_size=30,
    max_k=10,
    restarts=10,
    plddt_cutoff=70.0,
    verbose=False,
):
    """
    Run simple Voronoi predictions for all proteins.

    Returns dict mapping accession to list of benchmark Domain objects.
    """
    print(f"  Parameters: domain_bonus={domain_bonus}, crossing_penalty={crossing_penalty}, "
          f"intrusion_penalty={intrusion_penalty}, min_size={min_size}, max_k={max_k}, "
          f"restarts={restarts}, plddt_cutoff={plddt_cutoff}")

    predictions = {}

    for ann in annotations:
        if ann.exclude:
            continue

        print(f"Predicting {ann.uniprot_acc}...")
        try:
            # Load structure
            coords, plddts, resnums = load_structure(ann.uniprot_acc)

            # Filter by pLDDT
            coords, plddts, resnums = filter_structured(coords, plddts, resnums, plddt_cutoff)

            if len(coords) < min_size:
                print(f"  {ann.uniprot_acc}: Too few residues after filtering ({len(coords)})")
                predictions[ann.uniprot_acc] = []
                continue

            # Run prediction
            best_k, seeds, assignments = predict_domains(
                coords, resnums,
                max_k=max_k,
                domain_bonus=domain_bonus,
                crossing_penalty=crossing_penalty,
                intrusion_penalty=intrusion_penalty,
                min_size=min_size,
                restarts=restarts,
            )

            # Convert to domains
            domain_list = assignments_to_domains(assignments, resnums, min_size)

            # Convert to benchmark Domain objects
            domains = []
            for d in domain_list:
                domains.append(Domain(segments=d['segments']))

            print(f"  {ann.uniprot_acc}: {len(domains)} domains predicted (K={best_k})")
            for i, d in enumerate(domains):
                print(f"    Domain {i+1}: {d}")

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
        description="Benchmark simple Voronoi domain predictions using ABC benchmark framework"
    )
    parser.add_argument(
        "dev_set",
        help="Path to development set file with ground truth annotations"
    )
    parser.add_argument(
        "--output", "-o",
        default="simple_voronoi_benchmark_results",
        help="Output directory for results (default: simple_voronoi_benchmark_results)"
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
    # Simple Voronoi specific parameters
    parser.add_argument(
        "--domain-bonus",
        type=float,
        default=105,
        help="Reward for each valid domain (default: 105)"
    )
    parser.add_argument(
        "--crossing-penalty",
        type=float,
        default=100,
        help="Penalty per real chain crossing (default: 100)"
    )
    parser.add_argument(
        "--intrusion-penalty",
        type=float,
        default=10,
        help="Penalty per brief intrusion (default: 10)"
    )
    parser.add_argument(
        "--min-size",
        type=int,
        default=30,
        help="Minimum domain size (default: 30)"
    )
    parser.add_argument(
        "--max-k",
        type=int,
        default=10,
        help="Maximum number of domains to try (default: 10)"
    )
    parser.add_argument(
        "--restarts",
        type=int,
        default=10,
        help="Number of random restarts for hill-climbing (default: 10)"
    )
    parser.add_argument(
        "--plddt-cutoff",
        type=float,
        default=70.0,
        help="pLDDT cutoff for filtering residues (default: 70.0)"
    )

    args = parser.parse_args()

    # Parse ground truth
    print(f"Loading ground truth from {args.dev_set}...")
    annotations = parse_dev_set(args.dev_set)

    # Count valid proteins
    valid = [a for a in annotations if not a.exclude]
    print(f"Found {len(valid)} proteins for evaluation ({len(annotations) - len(valid)} excluded)")

    # Run predictions
    print(f"\nRunning simple Voronoi predictions...")
    predictions = run_simple_voronoi_predictions(
        annotations,
        domain_bonus=args.domain_bonus,
        crossing_penalty=args.crossing_penalty,
        intrusion_penalty=args.intrusion_penalty,
        min_size=args.min_size,
        max_k=args.max_k,
        restarts=args.restarts,
        plddt_cutoff=args.plddt_cutoff,
        verbose=args.verbose,
    )

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
