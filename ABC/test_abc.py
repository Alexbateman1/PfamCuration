#!/usr/bin/env python3
"""
Test script for ABC Domain Predictor

Tests the predictor on a set of multi-domain proteins with known structures.
Compares results visually and generates output for manual inspection.

Test proteins selected to cover:
- 2-domain proteins
- 3+ domain proteins
- Discontinuous domains
- NDR/linker regions
- Cases where other methods may fail
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from ABC import ABCPredictor

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


# Test proteins with expected domain architectures
TEST_PROTEINS = [
    {
        "uniprot_acc": "P00520",  # ABL1 kinase - 2 domains
        "name": "ABL1 (Abelson tyrosine kinase 1)",
        "expected_domains": 2,
        "notes": "Classic 2-domain kinase with SH3 and kinase domains",
    },
    {
        "uniprot_acc": "P12931",  # SRC kinase - 3 domains
        "name": "SRC (Proto-oncogene tyrosine-protein kinase Src)",
        "expected_domains": 3,
        "notes": "SH3 + SH2 + kinase domain architecture",
    },
    {
        "uniprot_acc": "P04637",  # p53 - multiple domains with IDR
        "name": "TP53 (Tumor protein p53)",
        "expected_domains": 3,
        "notes": "Has intrinsically disordered regions - good NDR test",
    },
    {
        "uniprot_acc": "P0A9Q1",  # AraC - discontinuous domain
        "name": "AraC (Arabinose operon regulatory protein)",
        "expected_domains": 2,
        "notes": "Has discontinuous DNA-binding domain",
    },
    {
        "uniprot_acc": "P08238",  # HSP90 - multi-domain with linkers
        "name": "HSP90 (Heat shock protein 90)",
        "expected_domains": 3,
        "notes": "N-terminal, middle, and C-terminal domains with linkers",
    },
]


def test_single_protein(
    predictor: ABCPredictor,
    uniprot_acc: str,
    output_dir: Path,
    expected_domains: int = None,
) -> Dict:
    """
    Test prediction on a single protein.

    Returns:
    --------
    Dict with prediction results and comparison metrics
    """
    logger.info(f"Testing protein: {uniprot_acc}")

    try:
        # Run prediction
        prediction = predictor.predict_from_uniprot(uniprot_acc)

        # Generate outputs
        output_base = output_dir / uniprot_acc
        chimerax_cmd = predictor.visualize(prediction, str(output_base))

        # Print summary
        print(f"\n{'='*60}")
        print(prediction.summary())
        print(f"\nChimeraX commands saved to: {output_base}.cxc")
        print(f"PyMOL script saved to: {output_base}.pml")
        print(f"Architecture diagram saved to: {output_base}.png")
        print(f"HTML report saved to: {output_base}.html")

        # Compare with expected
        result = {
            "uniprot_acc": uniprot_acc,
            "success": True,
            "n_domains": len(prediction.domains),
            "n_ndr": len(prediction.ndr_regions),
            "domains": [d.to_chopping_string() for d in prediction.domains],
            "domain_sizes": [d.size for d in prediction.domains],
            "quality_scores": [d.quality_metrics.get("quality_score", 0) for d in prediction.domains],
        }

        if expected_domains is not None:
            result["expected_domains"] = expected_domains
            result["domain_match"] = len(prediction.domains) == expected_domains
            if not result["domain_match"]:
                logger.warning(f"Domain count mismatch: expected {expected_domains}, got {len(prediction.domains)}")

        return result

    except Exception as e:
        logger.error(f"Error processing {uniprot_acc}: {e}")
        return {
            "uniprot_acc": uniprot_acc,
            "success": False,
            "error": str(e),
        }


def run_multi_threshold_test(
    predictor: ABCPredictor,
    uniprot_acc: str,
    output_dir: Path,
) -> Dict:
    """
    Run predictions at multiple thresholds to test stability.
    """
    import tempfile
    import urllib.request

    logger.info(f"Multi-threshold test for: {uniprot_acc}")

    # Download structure once
    cif_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_acc}-F1-model_v4.cif"
    pae_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_acc}-F1-predicted_aligned_error_v4.json"

    with tempfile.TemporaryDirectory() as tmpdir:
        cif_path = Path(tmpdir) / f"{uniprot_acc}.cif"
        pae_path = Path(tmpdir) / f"{uniprot_acc}_pae.json"

        urllib.request.urlretrieve(cif_url, cif_path)
        try:
            urllib.request.urlretrieve(pae_url, pae_path)
        except Exception:
            pae_path = None

        results = predictor.multi_threshold_analysis(
            str(cif_path),
            str(pae_path) if pae_path else None,
            thresholds=[8.0, 10.0, 12.0, 15.0],
        )

    # Analyze stability
    print(f"\n{'='*60}")
    print(f"Multi-threshold analysis for {uniprot_acc}")
    print("-" * 60)

    domain_counts = {}
    for threshold, pred in results.items():
        n_domains = len(pred.domains)
        domain_counts[threshold] = n_domains
        print(f"Threshold {threshold:5.1f}Å: {n_domains} domains")
        for d in pred.domains:
            print(f"    D{d.domain_id}: {d.to_chopping_string()} (size={d.size})")

    # Check stability
    counts = list(domain_counts.values())
    is_stable = len(set(counts)) <= 2  # At most 2 different counts

    print("-" * 60)
    print(f"Stable across thresholds: {'Yes' if is_stable else 'No'}")

    return {
        "uniprot_acc": uniprot_acc,
        "domain_counts": domain_counts,
        "stable": is_stable,
    }


def parameter_sensitivity_test(
    uniprot_acc: str,
    output_dir: Path,
) -> Dict:
    """
    Test sensitivity to different parameters.
    """
    import tempfile
    import urllib.request

    logger.info(f"Parameter sensitivity test for: {uniprot_acc}")

    # Download structure
    cif_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_acc}-F1-model_v4.cif"
    pae_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_acc}-F1-predicted_aligned_error_v4.json"

    with tempfile.TemporaryDirectory() as tmpdir:
        cif_path = Path(tmpdir) / f"{uniprot_acc}.cif"
        pae_path = Path(tmpdir) / f"{uniprot_acc}_pae.json"

        urllib.request.urlretrieve(cif_url, cif_path)
        try:
            urllib.request.urlretrieve(pae_url, pae_path)
        except Exception:
            pae_path = None

        results = {}

        # Test different resolutions
        print(f"\n{'='*60}")
        print(f"Parameter sensitivity for {uniprot_acc}")
        print("-" * 60)

        resolutions = [0.5, 1.0, 1.5, 2.0]
        print("\nClustering resolution sensitivity:")
        for res in resolutions:
            predictor = ABCPredictor(resolution=res)
            pred = predictor.predict_from_file(str(cif_path), str(pae_path) if pae_path else None)
            results[f"resolution_{res}"] = len(pred.domains)
            print(f"  Resolution {res}: {len(pred.domains)} domains")

        # Test NDR cutoffs
        ndr_cutoffs = [60, 70, 80]
        print("\nNDR pLDDT cutoff sensitivity:")
        for cutoff in ndr_cutoffs:
            predictor = ABCPredictor(ndr_plddt_cutoff=cutoff)
            pred = predictor.predict_from_file(str(cif_path), str(pae_path) if pae_path else None)
            results[f"ndr_cutoff_{cutoff}"] = len(pred.ndr_regions)
            print(f"  Cutoff {cutoff}: {len(pred.ndr_regions)} NDR regions")

        # Test with/without PAE
        print("\nPAE usage:")
        for use_pae in [True, False]:
            predictor = ABCPredictor(use_pae=use_pae)
            pred = predictor.predict_from_file(str(cif_path), str(pae_path) if pae_path else None)
            results[f"use_pae_{use_pae}"] = len(pred.domains)
            print(f"  PAE={use_pae}: {len(pred.domains)} domains")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Test ABC Domain Predictor on example proteins"
    )
    parser.add_argument(
        "--uniprot",
        type=str,
        help="Test specific UniProt accession instead of test set"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="./abc_test_results",
        help="Output directory for results"
    )
    parser.add_argument(
        "--multi-threshold",
        action="store_true",
        help="Run multi-threshold analysis"
    )
    parser.add_argument(
        "--sensitivity",
        action="store_true",
        help="Run parameter sensitivity analysis"
    )
    parser.add_argument(
        "--all-tests",
        action="store_true",
        help="Run all test proteins"
    )
    parser.add_argument(
        "--distance-threshold",
        type=float,
        default=10.0,
        help="Distance threshold for contacts (default: 10Å)"
    )
    parser.add_argument(
        "--resolution",
        type=float,
        default=1.0,
        help="Clustering resolution parameter (default: 1.0)"
    )
    parser.add_argument(
        "--min-domain-size",
        type=int,
        default=30,
        help="Minimum domain size (default: 30)"
    )
    parser.add_argument(
        "--ndr-cutoff",
        type=float,
        default=70.0,
        help="pLDDT cutoff for NDR detection (default: 70)"
    )

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Initialize predictor
    predictor = ABCPredictor(
        distance_threshold=args.distance_threshold,
        resolution=args.resolution,
        min_domain_size=args.min_domain_size,
        ndr_plddt_cutoff=args.ndr_cutoff,
    )

    print(f"ABC Domain Predictor Test Suite")
    print(f"================================")
    print(f"Parameters:")
    for key, value in predictor.parameters.items():
        print(f"  {key}: {value}")
    print()

    # Run tests
    all_results = []

    if args.uniprot:
        # Test single protein
        if args.multi_threshold:
            result = run_multi_threshold_test(predictor, args.uniprot, output_dir)
        elif args.sensitivity:
            result = parameter_sensitivity_test(args.uniprot, output_dir)
        else:
            result = test_single_protein(predictor, args.uniprot, output_dir)
        all_results.append(result)

    elif args.all_tests:
        # Run all test proteins
        for protein in TEST_PROTEINS:
            print(f"\n{'#'*60}")
            print(f"# {protein['name']}")
            print(f"# {protein['notes']}")
            print(f"{'#'*60}")

            result = test_single_protein(
                predictor,
                protein["uniprot_acc"],
                output_dir,
                protein.get("expected_domains"),
            )
            result["name"] = protein["name"]
            result["notes"] = protein["notes"]
            all_results.append(result)

    else:
        # Run on first test protein as demo
        protein = TEST_PROTEINS[0]
        print(f"Running demo on: {protein['name']}")
        result = test_single_protein(
            predictor,
            protein["uniprot_acc"],
            output_dir,
            protein.get("expected_domains"),
        )
        all_results.append(result)

    # Save results summary
    results_file = output_dir / "test_results.json"
    with open(results_file, "w") as f:
        json.dump(all_results, f, indent=2)

    print(f"\n{'='*60}")
    print(f"Results saved to: {results_file}")

    # Print summary
    if len(all_results) > 1:
        print(f"\nSummary:")
        print(f"--------")
        successful = sum(1 for r in all_results if r.get("success", False))
        print(f"Successful: {successful}/{len(all_results)}")

        if any("domain_match" in r for r in all_results):
            matches = sum(1 for r in all_results if r.get("domain_match", False))
            print(f"Domain count matches: {matches}/{len([r for r in all_results if 'domain_match' in r])}")


if __name__ == "__main__":
    main()
