#!/usr/bin/env python3
"""
ABC Domain Predictor - Command Line Interface

Usage:
    python -m ABC.cli predict --uniprot P12345
    python -m ABC.cli predict --pdb structure.cif --pae pae.json
    python -m ABC.cli batch --input proteins.txt --output results/
"""

import argparse
import json
import logging
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from ABC import ABCPredictor

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def cmd_predict(args):
    """Run domain prediction on a single protein."""
    predictor = ABCPredictor(
        distance_threshold=args.distance_threshold,
        min_domain_size=args.min_domain_size,
        min_contact_ratio=args.min_contact_ratio,
        ndr_plddt_cutoff=args.ndr_cutoff,
        clustering_method=args.clustering,
        resolution=args.resolution,
        sigma=args.sigma,
        use_pae=not args.no_pae,
        min_bridge_seq_sep=args.min_bridge_sep,
        min_bridge_count=args.min_bridge_count,
    )

    if args.uniprot:
        prediction = predictor.predict_from_uniprot(args.uniprot)
        base_name = args.uniprot
    elif args.pdb:
        prediction = predictor.predict_from_file(
            args.pdb,
            args.pae,
            args.uniprot or Path(args.pdb).stem,
        )
        base_name = Path(args.pdb).stem
    else:
        logger.error("Must provide either --uniprot or --pdb")
        sys.exit(1)

    # Print summary
    print(prediction.summary())

    # Generate outputs
    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        output_path = Path(base_name)

    # Visualize and save outputs
    chimerax = predictor.visualize(prediction, str(output_path))

    print(f"\nOutput files:")
    print(f"  ChimeraX: {output_path}.cxc")
    print(f"  PyMOL: {output_path}.pml")
    print(f"  PNG: {output_path}.png")
    print(f"  HTML: {output_path}.html")

    # Also output JSON if requested
    if args.json:
        json_path = output_path.with_suffix(".json")
        with open(json_path, "w") as f:
            json.dump(prediction.to_dict(), f, indent=2)
        print(f"  JSON: {json_path}")

    # Print ChimeraX commands if verbose
    if args.verbose:
        print(f"\nChimeraX commands:")
        print("-" * 40)
        print(chimerax)


def cmd_batch(args):
    """Run domain prediction on multiple proteins."""
    predictor = ABCPredictor(
        distance_threshold=args.distance_threshold,
        min_domain_size=args.min_domain_size,
        min_contact_ratio=args.min_contact_ratio,
        ndr_plddt_cutoff=args.ndr_cutoff,
        clustering_method=args.clustering,
        resolution=args.resolution,
        min_bridge_seq_sep=args.min_bridge_sep,
        min_bridge_count=args.min_bridge_count,
    )

    # Read protein list
    input_path = Path(args.input)
    if not input_path.exists():
        logger.error(f"Input file not found: {args.input}")
        sys.exit(1)

    with open(input_path) as f:
        proteins = [line.strip() for line in f if line.strip() and not line.startswith("#")]

    logger.info(f"Processing {len(proteins)} proteins")

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for i, protein in enumerate(proteins):
        logger.info(f"[{i+1}/{len(proteins)}] Processing {protein}")

        try:
            prediction = predictor.predict_from_uniprot(protein)

            # Save individual outputs
            output_base = output_dir / protein
            predictor.visualize(prediction, str(output_base))

            results.append({
                "uniprot_acc": protein,
                "success": True,
                "n_domains": len(prediction.domains),
                "domains": [d.to_chopping_string() for d in prediction.domains],
            })

        except Exception as e:
            logger.error(f"Error processing {protein}: {e}")
            results.append({
                "uniprot_acc": protein,
                "success": False,
                "error": str(e),
            })

    # Save summary
    summary_path = output_dir / "batch_results.json"
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2)

    # Print summary
    successful = sum(1 for r in results if r.get("success", False))
    print(f"\nBatch complete: {successful}/{len(proteins)} successful")
    print(f"Results saved to: {summary_path}")


def cmd_compare(args):
    """Compare ABC prediction with other methods (placeholder for future)."""
    print("Comparison feature not yet implemented")
    print("Would compare ABC with: Chainsaw, Merizo, UniDoc")


def main():
    parser = argparse.ArgumentParser(
        description="ABC (Alex Bateman Chop) Domain Predictor",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Predict domains from UniProt accession
  python -m ABC.cli predict --uniprot P12345

  # Predict from local structure file
  python -m ABC.cli predict --pdb structure.cif --pae pae.json

  # Batch processing
  python -m ABC.cli batch --input proteins.txt --output results/

  # Custom parameters
  python -m ABC.cli predict --uniprot P12345 --distance-threshold 12 --resolution 0.8
        """
    )

    # Global options
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    subparsers = parser.add_subparsers(dest="command", help="Command to run")

    # Predict command
    predict_parser = subparsers.add_parser("predict", help="Predict domains for a single protein")
    predict_parser.add_argument("--uniprot", "-u", type=str, help="UniProt accession")
    predict_parser.add_argument("--pdb", "-p", type=str, help="Path to PDB/CIF file")
    predict_parser.add_argument("--pae", type=str, help="Path to PAE JSON file")
    predict_parser.add_argument("--output", "-o", type=str, help="Output file base name")
    predict_parser.add_argument("--json", action="store_true", help="Also output JSON")

    # Algorithm parameters
    predict_parser.add_argument("--distance-threshold", type=float, default=10.0,
                               help="Contact distance threshold (default: 10Å)")
    predict_parser.add_argument("--min-domain-size", type=int, default=30,
                               help="Minimum domain size (default: 30)")
    predict_parser.add_argument("--min-contact-ratio", type=float, default=1.5,
                               help="Minimum internal/external contact ratio for domains (default: 1.5)")
    predict_parser.add_argument("--ndr-cutoff", type=float, default=70.0,
                               help="pLDDT cutoff for NDR (default: 70)")
    predict_parser.add_argument("--clustering", type=str, default="leiden",
                               choices=["leiden", "louvain"],
                               help="Clustering algorithm (default: leiden)")
    predict_parser.add_argument("--resolution", type=float, default=0.2,
                               help="Clustering resolution (default: 0.2)")
    predict_parser.add_argument("--sigma", type=float, default=8.0,
                               help="Gaussian decay sigma (default: 8.0)")
    predict_parser.add_argument("--no-pae", action="store_true",
                               help="Don't use PAE for edge weighting")
    predict_parser.add_argument("--min-bridge-sep", type=int, default=50,
                               help="Minimum sequence separation for bridge detection (default: 50)")
    predict_parser.add_argument("--min-bridge-count", type=int, default=10,
                               help="Minimum bridges to merge discontinuous domains (default: 10)")
    predict_parser.set_defaults(func=cmd_predict)

    # Batch command
    batch_parser = subparsers.add_parser("batch", help="Batch predict for multiple proteins")
    batch_parser.add_argument("--input", "-i", type=str, required=True,
                             help="Input file with UniProt accessions (one per line)")
    batch_parser.add_argument("--output", "-o", type=str, required=True,
                             help="Output directory")
    batch_parser.add_argument("--distance-threshold", type=float, default=10.0,
                             help="Contact distance threshold")
    batch_parser.add_argument("--min-domain-size", type=int, default=30,
                             help="Minimum domain size")
    batch_parser.add_argument("--min-contact-ratio", type=float, default=1.5,
                             help="Minimum internal/external contact ratio for domains (default: 1.5)")
    batch_parser.add_argument("--ndr-cutoff", type=float, default=70.0,
                             help="pLDDT cutoff for NDR")
    batch_parser.add_argument("--clustering", type=str, default="leiden",
                             choices=["leiden", "louvain"],
                             help="Clustering algorithm")
    batch_parser.add_argument("--resolution", type=float, default=0.2,
                             help="Clustering resolution (default: 0.2)")
    batch_parser.add_argument("--min-bridge-sep", type=int, default=50,
                             help="Minimum sequence separation for bridge detection (default: 50)")
    batch_parser.add_argument("--min-bridge-count", type=int, default=10,
                             help="Minimum bridges to merge discontinuous domains (default: 10)")
    batch_parser.set_defaults(func=cmd_batch)

    # Compare command (placeholder)
    compare_parser = subparsers.add_parser("compare", help="Compare with other methods")
    compare_parser.add_argument("--uniprot", "-u", type=str, required=True,
                               help="UniProt accession")
    compare_parser.set_defaults(func=cmd_compare)

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    args.func(args)


if __name__ == "__main__":
    main()
