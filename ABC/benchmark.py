#!/usr/bin/env python3
"""
ABC Domain Predictor Benchmarking Suite

Evaluates domain predictions against ground truth using metrics from
Merizo and Chainsaw papers:
- IoU (Intersection-over-Union) per domain
- Boundary accuracy (MCC) at multiple tolerances
- Domain count analysis
- Error classification

Usage:
    python benchmark.py dev_set.txt [--output results/] [--run-predictions]
"""

import argparse
import json
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

# Add parent directory to path for ABC imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


@dataclass
class Domain:
    """A domain with potentially discontinuous segments."""
    segments: List[Tuple[int, int]]  # List of (start, end) tuples

    @property
    def residues(self) -> Set[int]:
        """Get all residue positions in this domain."""
        res = set()
        for start, end in self.segments:
            res.update(range(start, end + 1))
        return res

    @property
    def size(self) -> int:
        """Total number of residues."""
        return len(self.residues)

    @property
    def boundaries(self) -> List[int]:
        """Get all boundary positions (starts and ends)."""
        bounds = []
        for start, end in self.segments:
            bounds.extend([start, end])
        return bounds

    def __repr__(self):
        return "_".join(f"{s}-{e}" for s, e in self.segments)


@dataclass
class ProteinAnnotation:
    """Ground truth or predicted annotation for a protein."""
    uniprot_acc: str
    domains: List[Domain]
    exclude: bool = False  # True if marked with ? (complex case)

    @property
    def has_domains(self) -> bool:
        return len(self.domains) > 0

    @property
    def all_boundaries(self) -> List[int]:
        """Get all boundaries across all domains."""
        bounds = []
        for domain in self.domains:
            bounds.extend(domain.boundaries)
        return sorted(bounds)


@dataclass
class DomainMatch:
    """Match between predicted and ground truth domain."""
    pred_domain: Domain
    true_domain: Domain
    iou: float

    @property
    def is_correct(self, threshold: float = 0.8) -> bool:
        return self.iou >= threshold


@dataclass
class ProteinResult:
    """Evaluation results for a single protein."""
    uniprot_acc: str
    true_domains: List[Domain]
    pred_domains: List[Domain]
    matches: List[DomainMatch]
    weighted_iou: float
    domain_count_error: int  # pred - true
    error_type: str  # "correct", "under_split", "over_split", "boundary_error", etc.
    missing_domains: List[Domain] = field(default_factory=list)
    false_domains: List[Domain] = field(default_factory=list)

    @property
    def score(self) -> float:
        """
        Combined score (0-1) where 1 is perfect prediction.

        Combines:
        - IoU (70% weight): measures boundary accuracy
        - Domain count accuracy (30% weight): penalizes wrong number of domains
        """
        # IoU component (0-1)
        iou_score = self.weighted_iou

        # Domain count component (0-1)
        # Penalize based on relative error
        n_true = max(len(self.true_domains), 1)  # Avoid division by zero
        count_error_ratio = abs(self.domain_count_error) / n_true
        count_score = max(0, 1 - count_error_ratio)

        # Combined score
        return 0.7 * iou_score + 0.3 * count_score

    @property
    def status(self) -> str:
        """Simple status: GOOD, OK, or BAD based on score."""
        if self.score >= 0.9:
            return "GOOD"
        elif self.score >= 0.7:
            return "OK"
        else:
            return "BAD"


@dataclass
class BenchmarkResults:
    """Complete benchmark results."""
    protein_results: List[ProteinResult]
    boundary_mcc: Dict[int, float]  # tolerance -> MCC
    overall_iou: float
    domain_mae: float
    under_splits: int
    over_splits: int
    correct_domain_count: int  # proteins with right number of domains
    perfect_predictions: int   # proteins with correct boundaries too

    def summary(self) -> str:
        n = len(self.protein_results)
        lines = [
            "=" * 60,
            "ABC Domain Predictor Benchmark Results",
            "=" * 60,
            "",
            f"Proteins evaluated: {n}",
            f"Overall length-weighted IoU: {self.overall_iou:.3f}",
            "",
            "Domain Count Analysis:",
            f"  Correct domain count: {self.correct_domain_count} ({100*self.correct_domain_count/n:.1f}%)",
            f"  Perfect predictions:  {self.perfect_predictions} ({100*self.perfect_predictions/n:.1f}%)",
            f"  Under-splits: {self.under_splits} ({100*self.under_splits/n:.1f}%)",
            f"  Over-splits:  {self.over_splits} ({100*self.over_splits/n:.1f}%)",
            f"  Domain count MAE: {self.domain_mae:.2f}",
            "",
            "Boundary Accuracy (F1):",
        ]
        for tol in sorted(self.boundary_mcc.keys()):
            lines.append(f"  ±{tol:2d} residues: F1 = {self.boundary_mcc[tol]:.3f}")

        return "\n".join(lines)


def parse_dev_set(filepath: str) -> List[ProteinAnnotation]:
    """
    Parse development set file.

    Format:
        ACCESSION domain1 domain2 ...

    Where domains are:
        - "start-end" for continuous domains
        - "start1-end1_start2-end2" for discontinuous domains
        - Empty line after accession = no domains
        - "?" = exclude from evaluation (complex case)
    """
    annotations = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()

            # Skip comments and empty lines
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            acc = parts[0]
            domain_strs = parts[1:] if len(parts) > 1 else []

            # Check for exclusion marker
            exclude = False
            if domain_strs and domain_strs[0] == '?':
                exclude = True
                domain_strs = []

            # Parse domains
            domains = []
            for dom_str in domain_strs:
                segments = []
                for seg_str in dom_str.split('_'):
                    if '-' in seg_str:
                        start, end = seg_str.split('-')
                        segments.append((int(start), int(end)))
                if segments:
                    domains.append(Domain(segments=segments))

            annotations.append(ProteinAnnotation(
                uniprot_acc=acc,
                domains=domains,
                exclude=exclude,
            ))

    return annotations


def calculate_iou(domain1: Domain, domain2: Domain) -> float:
    """Calculate Intersection-over-Union between two domains."""
    res1 = domain1.residues
    res2 = domain2.residues

    intersection = len(res1 & res2)
    union = len(res1 | res2)

    if union == 0:
        return 0.0

    return intersection / union


def match_domains(
    true_domains: List[Domain],
    pred_domains: List[Domain],
    iou_threshold: float = 0.8,
) -> Tuple[List[DomainMatch], List[Domain], List[Domain]]:
    """
    Match predicted domains to ground truth using IoU.

    Uses greedy matching: best IoU pairs first.

    Returns:
        (matches, missing_domains, false_domains)
    """
    if not true_domains and not pred_domains:
        return [], [], []

    if not true_domains:
        return [], [], pred_domains.copy()

    if not pred_domains:
        return [], true_domains.copy(), []

    # Calculate all pairwise IoUs
    iou_matrix = np.zeros((len(true_domains), len(pred_domains)))
    for i, true_dom in enumerate(true_domains):
        for j, pred_dom in enumerate(pred_domains):
            iou_matrix[i, j] = calculate_iou(true_dom, pred_dom)

    # Greedy matching
    matches = []
    matched_true = set()
    matched_pred = set()

    while True:
        # Find best remaining match
        best_iou = 0
        best_i, best_j = -1, -1

        for i in range(len(true_domains)):
            if i in matched_true:
                continue
            for j in range(len(pred_domains)):
                if j in matched_pred:
                    continue
                if iou_matrix[i, j] > best_iou:
                    best_iou = iou_matrix[i, j]
                    best_i, best_j = i, j

        if best_iou == 0 or best_i == -1:
            break

        matches.append(DomainMatch(
            pred_domain=pred_domains[best_j],
            true_domain=true_domains[best_i],
            iou=best_iou,
        ))
        matched_true.add(best_i)
        matched_pred.add(best_j)

    # Identify missing and false domains
    missing = [true_domains[i] for i in range(len(true_domains)) if i not in matched_true]
    false = [pred_domains[j] for j in range(len(pred_domains)) if j not in matched_pred]

    return matches, missing, false


def calculate_weighted_iou(
    matches: List[DomainMatch],
    true_domains: List[Domain],
) -> float:
    """
    Calculate length-weighted average IoU.

    Weights by ground truth domain size.
    Missing domains contribute 0 to their weight.
    """
    if not true_domains:
        return 1.0  # No domains expected, any empty prediction is correct

    total_weight = sum(d.size for d in true_domains)
    if total_weight == 0:
        return 1.0

    weighted_sum = 0.0
    matched_true = {id(m.true_domain) for m in matches}

    for match in matches:
        weighted_sum += match.iou * match.true_domain.size

    # Missing domains contribute 0
    for dom in true_domains:
        if id(dom) not in matched_true:
            weighted_sum += 0.0  # Explicit for clarity

    return weighted_sum / total_weight


def classify_error(
    true_domains: List[Domain],
    pred_domains: List[Domain],
    matches: List[DomainMatch],
    iou_threshold: float = 0.8,
) -> str:
    """
    Classify the type of error for a protein.

    Returns one of:
        - "correct": all domains matched with IoU >= threshold
        - "under_split": fewer predicted domains (merged)
        - "over_split": more predicted domains (split)
        - "boundary_error": same count but poor IoU
        - "mixed": combination of errors
        - "no_domains_correct": correctly predicted no domains
        - "false_positive": predicted domains when none exist
        - "false_negative": missed all domains
    """
    n_true = len(true_domains)
    n_pred = len(pred_domains)

    # Special case: no domains
    if n_true == 0 and n_pred == 0:
        return "no_domains_correct"
    if n_true == 0 and n_pred > 0:
        return "false_positive"
    if n_true > 0 and n_pred == 0:
        return "false_negative"

    # Check if all matches are above threshold
    good_matches = [m for m in matches if m.iou >= iou_threshold]

    if len(good_matches) == n_true == n_pred:
        return "correct"

    # Determine error type
    if n_pred < n_true:
        return "under_split"
    elif n_pred > n_true:
        return "over_split"
    else:
        return "boundary_error"


def calculate_boundary_mcc(
    protein_results: List[ProteinResult],
    tolerance: int,
) -> float:
    """
    Calculate Matthews Correlation Coefficient for boundary detection.

    A boundary is "correct" if within ±tolerance residues of ground truth.
    """
    tp = 0  # True positives: predicted boundary near true boundary
    fp = 0  # False positives: predicted boundary not near any true
    fn = 0  # False negatives: true boundary not near any predicted
    tn = 0  # True negatives: hard to define for boundaries, typically not used

    for result in protein_results:
        true_bounds = set()
        for dom in result.true_domains:
            true_bounds.update(dom.boundaries)

        pred_bounds = set()
        for dom in result.pred_domains:
            pred_bounds.update(dom.boundaries)

        # Check each true boundary
        for tb in true_bounds:
            found = any(abs(tb - pb) <= tolerance for pb in pred_bounds)
            if found:
                tp += 1
            else:
                fn += 1

        # Check each predicted boundary
        for pb in pred_bounds:
            found = any(abs(pb - tb) <= tolerance for tb in true_bounds)
            if not found:
                fp += 1

    # MCC calculation
    # MCC = (TP*TN - FP*FN) / sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN))
    # Since TN is not well-defined, we use a simplified version
    # that focuses on precision and recall

    if tp + fp == 0 or tp + fn == 0:
        return 0.0

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0

    # F1-like MCC approximation when TN is not available
    # Actually, let's compute proper MCC with TN=0 assumption
    numerator = tp  # Simplified since TN=0
    denominator = np.sqrt((tp + fp) * (tp + fn) * max(1, fp) * max(1, fn))

    if denominator == 0:
        return 0.0

    # Use F1 score as a more interpretable alternative
    if precision + recall == 0:
        return 0.0
    f1 = 2 * precision * recall / (precision + recall)

    return f1


def run_predictions(annotations: List[ProteinAnnotation]) -> Dict[str, List[Domain]]:
    """
    Run ABC predictions for all proteins.

    Returns dict mapping accession to list of predicted domains.
    """
    from ABC.abc_predictor import ABCPredictor

    predictor = ABCPredictor()
    predictions = {}

    for ann in annotations:
        if ann.exclude:
            continue

        print(f"Predicting {ann.uniprot_acc}...")
        try:
            result = predictor.predict_from_uniprot(ann.uniprot_acc)

            # Convert to Domain objects
            domains = []
            for d in result.domains:
                domains.append(Domain(segments=d.segments))

            predictions[ann.uniprot_acc] = domains

        except Exception as e:
            print(f"  Error: {e}")
            predictions[ann.uniprot_acc] = []

    return predictions


def evaluate_protein(
    acc: str,
    true_domains: List[Domain],
    pred_domains: List[Domain],
    iou_threshold: float = 0.8,
) -> ProteinResult:
    """Evaluate predictions for a single protein."""
    matches, missing, false = match_domains(true_domains, pred_domains, iou_threshold)
    weighted_iou = calculate_weighted_iou(matches, true_domains)

    domain_count_error = len(pred_domains) - len(true_domains)
    error_type = classify_error(true_domains, pred_domains, matches, iou_threshold)

    return ProteinResult(
        uniprot_acc=acc,
        true_domains=true_domains,
        pred_domains=pred_domains,
        matches=matches,
        weighted_iou=weighted_iou,
        domain_count_error=domain_count_error,
        error_type=error_type,
        missing_domains=missing,
        false_domains=false,
    )


def run_benchmark(
    annotations: List[ProteinAnnotation],
    predictions: Dict[str, List[Domain]],
    iou_threshold: float = 0.8,
) -> BenchmarkResults:
    """Run full benchmark evaluation."""
    protein_results = []

    for ann in annotations:
        if ann.exclude:
            continue

        if ann.uniprot_acc not in predictions:
            print(f"Warning: No prediction for {ann.uniprot_acc}")
            continue

        result = evaluate_protein(
            ann.uniprot_acc,
            ann.domains,
            predictions[ann.uniprot_acc],
            iou_threshold,
        )
        protein_results.append(result)

    # Calculate overall metrics
    total_weight = sum(
        sum(d.size for d in r.true_domains)
        for r in protein_results
        if r.true_domains
    )

    if total_weight > 0:
        overall_iou = sum(
            r.weighted_iou * sum(d.size for d in r.true_domains)
            for r in protein_results
            if r.true_domains
        ) / total_weight
    else:
        overall_iou = 1.0

    # Domain count analysis
    under_splits = sum(1 for r in protein_results if r.error_type == "under_split")
    over_splits = sum(1 for r in protein_results if r.error_type == "over_split")

    # Correct domain count = right number of domains (includes boundary_error)
    correct_domain_count = sum(
        1 for r in protein_results
        if r.error_type in ["correct", "no_domains_correct", "boundary_error"]
    )
    # Perfect predictions = correct count AND good boundaries
    perfect_predictions = sum(
        1 for r in protein_results
        if r.error_type in ["correct", "no_domains_correct"]
    )

    domain_mae = np.mean([abs(r.domain_count_error) for r in protein_results])

    # Boundary MCC at multiple tolerances
    tolerances = [5, 10, 15, 20, 25, 30]
    boundary_mcc = {}
    for tol in tolerances:
        boundary_mcc[tol] = calculate_boundary_mcc(protein_results, tol)

    return BenchmarkResults(
        protein_results=protein_results,
        boundary_mcc=boundary_mcc,
        overall_iou=overall_iou,
        domain_mae=domain_mae,
        under_splits=under_splits,
        over_splits=over_splits,
        correct_domain_count=correct_domain_count,
        perfect_predictions=perfect_predictions,
    )


def generate_detailed_report(results: BenchmarkResults, output_dir: Path):
    """Generate detailed report files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Sort proteins by score (worst first for easy review)
    sorted_results = sorted(results.protein_results, key=lambda r: r.score)

    # Print console summary sorted by score
    print("\n" + "=" * 80)
    print("Per-Protein Results (sorted by score, worst first)")
    print("=" * 80)
    print(f"{'Accession':<15} {'Score':>6} {'Status':<5} {'IoU':>5} {'#True':>5} {'#Pred':>5} {'Error Type':<15}")
    print("-" * 80)

    for r in sorted_results:
        print(f"{r.uniprot_acc:<15} {r.score:>6.3f} {r.status:<5} {r.weighted_iou:>5.2f} "
              f"{len(r.true_domains):>5} {len(r.pred_domains):>5} {r.error_type:<15}")

    print("-" * 80)

    # Show details for BAD predictions
    bad_results = [r for r in sorted_results if r.status == "BAD"]
    if bad_results:
        print(f"\n{'='*80}")
        print(f"Details for {len(bad_results)} BAD predictions:")
        print("=" * 80)

        for r in bad_results:
            print(f"\n{r.uniprot_acc} (score={r.score:.3f}):")
            print(f"  True:    {'; '.join(str(d) for d in r.true_domains) or 'none'}")
            print(f"  Pred:    {'; '.join(str(d) for d in r.pred_domains) or 'none'}")
            if r.missing_domains:
                print(f"  Missing: {'; '.join(str(d) for d in r.missing_domains)}")
            if r.false_domains:
                print(f"  False:   {'; '.join(str(d) for d in r.false_domains)}")
            if r.matches:
                print(f"  Matches:")
                for m in r.matches:
                    status = "✓" if m.iou >= 0.8 else "✗"
                    print(f"    {status} {m.true_domain} <-> {m.pred_domain} (IoU={m.iou:.3f})")

    # Per-protein results file (sorted by score)
    with open(output_dir / "protein_results.tsv", "w") as f:
        f.write("Accession\tScore\tStatus\tTrue_Domains\tPred_Domains\tTrue_Count\tPred_Count\t"
                "Weighted_IoU\tCount_Error\tError_Type\tMissing\tFalse\n")

        for r in sorted_results:
            true_str = "; ".join(str(d) for d in r.true_domains) or "none"
            pred_str = "; ".join(str(d) for d in r.pred_domains) or "none"
            missing_str = "; ".join(str(d) for d in r.missing_domains) or ""
            false_str = "; ".join(str(d) for d in r.false_domains) or ""

            f.write(f"{r.uniprot_acc}\t{r.score:.3f}\t{r.status}\t{true_str}\t{pred_str}\t"
                    f"{len(r.true_domains)}\t{len(r.pred_domains)}\t"
                    f"{r.weighted_iou:.3f}\t{r.domain_count_error}\t"
                    f"{r.error_type}\t{missing_str}\t{false_str}\n")

    # Summary statistics
    with open(output_dir / "summary.txt", "w") as f:
        f.write(results.summary())

    # Domain matches detail
    with open(output_dir / "domain_matches.tsv", "w") as f:
        f.write("Accession\tTrue_Domain\tPred_Domain\tIoU\tCorrect\n")

        for r in results.protein_results:
            for m in r.matches:
                correct = "yes" if m.iou >= 0.8 else "no"
                f.write(f"{r.uniprot_acc}\t{m.true_domain}\t{m.pred_domain}\t"
                        f"{m.iou:.3f}\t{correct}\n")

    # Boundary accuracy table
    with open(output_dir / "boundary_accuracy.tsv", "w") as f:
        f.write("Tolerance\tMCC_F1\n")
        for tol in sorted(results.boundary_mcc.keys()):
            f.write(f"{tol}\t{results.boundary_mcc[tol]:.3f}\n")

    print(f"\nDetailed reports written to {output_dir}/")


def generate_plots(results: BenchmarkResults, output_dir: Path):
    """Generate visualization plots."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, skipping plots")
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Boundary accuracy vs tolerance
    fig, ax = plt.subplots(figsize=(8, 5))
    tolerances = sorted(results.boundary_mcc.keys())
    mccs = [results.boundary_mcc[t] for t in tolerances]

    ax.plot(tolerances, mccs, 'bo-', linewidth=2, markersize=8)
    ax.set_xlabel('Boundary Tolerance (residues)', fontsize=12)
    ax.set_ylabel('F1 Score', fontsize=12)
    ax.set_title('Boundary Detection Accuracy', fontsize=14)
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3)
    ax.set_xticks(tolerances)

    plt.tight_layout()
    plt.savefig(output_dir / "boundary_accuracy.png", dpi=150)
    plt.close()

    # 2. Error type distribution
    error_counts = defaultdict(int)
    for r in results.protein_results:
        error_counts[r.error_type] += 1

    fig, ax = plt.subplots(figsize=(8, 5))

    labels = list(error_counts.keys())
    counts = [error_counts[l] for l in labels]
    colors = []
    for l in labels:
        if l in ["correct", "no_domains_correct"]:
            colors.append("green")
        elif l == "under_split":
            colors.append("orange")
        elif l == "over_split":
            colors.append("red")
        else:
            colors.append("gray")

    bars = ax.bar(labels, counts, color=colors, edgecolor='black')
    ax.set_ylabel('Number of Proteins', fontsize=12)
    ax.set_title('Error Type Distribution', fontsize=14)
    plt.xticks(rotation=45, ha='right')

    # Add count labels on bars
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                str(count), ha='center', va='bottom', fontsize=10)

    plt.tight_layout()
    plt.savefig(output_dir / "error_distribution.png", dpi=150)
    plt.close()

    # 3. IoU distribution
    ious = [r.weighted_iou for r in results.protein_results]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(ious, bins=20, range=(0, 1), edgecolor='black', alpha=0.7)
    ax.axvline(x=0.8, color='red', linestyle='--', linewidth=2, label='Threshold (0.8)')
    ax.set_xlabel('Length-weighted IoU', fontsize=12)
    ax.set_ylabel('Number of Proteins', fontsize=12)
    ax.set_title('IoU Distribution', fontsize=14)
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_dir / "iou_distribution.png", dpi=150)
    plt.close()

    # 4. Domain count comparison
    fig, ax = plt.subplots(figsize=(6, 6))

    true_counts = [len(r.true_domains) for r in results.protein_results]
    pred_counts = [len(r.pred_domains) for r in results.protein_results]

    max_count = max(max(true_counts), max(pred_counts)) + 1

    ax.scatter(true_counts, pred_counts, alpha=0.6, s=100)
    ax.plot([0, max_count], [0, max_count], 'r--', linewidth=2, label='Perfect')
    ax.set_xlabel('True Domain Count', fontsize=12)
    ax.set_ylabel('Predicted Domain Count', fontsize=12)
    ax.set_title('Domain Count: True vs Predicted', fontsize=14)
    ax.set_xlim(-0.5, max_count)
    ax.set_ylim(-0.5, max_count)
    ax.set_aspect('equal')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / "domain_count_comparison.png", dpi=150)
    plt.close()

    print(f"Plots saved to {output_dir}/")


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark ABC domain predictions against ground truth"
    )
    parser.add_argument(
        "dev_set",
        help="Path to development set file with ground truth annotations"
    )
    parser.add_argument(
        "--predictions",
        help="Path to JSON file with pre-computed predictions (skip running ABC)"
    )
    parser.add_argument(
        "--output", "-o",
        default="benchmark_results",
        help="Output directory for results (default: benchmark_results)"
    )
    parser.add_argument(
        "--run-predictions",
        action="store_true",
        help="Run ABC predictions (otherwise expects --predictions file)"
    )
    parser.add_argument(
        "--save-predictions",
        help="Save predictions to JSON file for future use"
    )
    parser.add_argument(
        "--iou-threshold",
        type=float,
        default=0.8,
        help="IoU threshold for correct domain match (default: 0.8)"
    )

    args = parser.parse_args()

    # Parse ground truth
    print(f"Loading ground truth from {args.dev_set}...")
    annotations = parse_dev_set(args.dev_set)

    # Count valid proteins
    valid = [a for a in annotations if not a.exclude]
    print(f"Found {len(valid)} proteins for evaluation ({len(annotations) - len(valid)} excluded)")

    # Get predictions
    if args.predictions:
        print(f"Loading predictions from {args.predictions}...")
        with open(args.predictions) as f:
            pred_data = json.load(f)

        predictions = {}
        for acc, domains in pred_data.items():
            dom_list = []
            for d in domains:
                if isinstance(d, dict):
                    dom_list.append(Domain(segments=[tuple(s) for s in d["segments"]]))
                else:
                    # Simple format: list of segment tuples
                    dom_list.append(Domain(segments=[tuple(s) for s in d]))
            predictions[acc] = dom_list

    elif args.run_predictions:
        print("Running ABC predictions...")
        predictions = run_predictions(annotations)

        if args.save_predictions:
            # Save for future use
            pred_data = {}
            for acc, domains in predictions.items():
                pred_data[acc] = [{"segments": d.segments} for d in domains]

            with open(args.save_predictions, "w") as f:
                json.dump(pred_data, f, indent=2)
            print(f"Predictions saved to {args.save_predictions}")

    else:
        print("Error: Either --predictions or --run-predictions must be specified")
        sys.exit(1)

    # Run benchmark
    print("\nRunning benchmark evaluation...")
    results = run_benchmark(annotations, predictions, args.iou_threshold)

    # Print summary
    print("\n" + results.summary())

    # Generate reports and plots
    output_dir = Path(args.output)
    generate_detailed_report(results, output_dir)
    generate_plots(results, output_dir)

    print(f"\nResults saved to {output_dir}/")


if __name__ == "__main__":
    main()
