#!/usr/bin/env python3
"""
Perform all-vs-all Pfam AlphaFold model comparison using Reseek.

This script:
1. Creates a Reseek .bca database from Pfam AlphaFold CIF models
2. Runs a single all-vs-all search (much more efficient than individual searches)
3. Parses results into HHsearch-compatible format

Usage:
  python3 pfam_reseek_allvsall.py -n 28000 -dir ResultDB -sensitivity sensitive -threads 32
"""

import argparse
import subprocess
import sys
from pathlib import Path
import re
import pickle


# Configuration
PFAM_MODELS_DIR = Path("/nfs/production/agb/interpro/users/typhaine/interpro-pfam-curation-tools/pfam_add_clan_search/results/chopped_cif")


def run_command(cmd, cwd=None, check=True, verbose=True):
    """Run a shell command and return result."""
    if verbose:
        print(f"  Running: {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        if result.stderr:
            print(f"  Error: {result.stderr}", file=sys.stderr)
        if check:
            sys.exit(1)
    return result


def extract_pfam_accession(filename):
    """Extract Pfam accession from a filename (PF##### pattern)."""
    match = re.search(r'(PF\d{5})', filename)
    if match:
        return match.group(1)
    return None


def build_pfam_model_index(models_dir, num_families=None):
    """Build an index of available Pfam AlphaFold models."""
    print("=== Building Pfam model index ===")

    model_index = {}
    models_dir_path = Path(models_dir)

    if not models_dir_path.exists():
        print(f"Error: Models directory not found: {models_dir}", file=sys.stderr)
        sys.exit(1)

    # Find all CIF files
    cif_files = list(models_dir_path.glob("*.cif"))
    print(f"Found {len(cif_files)} CIF files in {models_dir}")

    for cif_file in cif_files:
        pfam_acc = extract_pfam_accession(cif_file.name)
        if pfam_acc:
            if num_families is None:
                # Include all families
                if pfam_acc not in model_index:
                    model_index[pfam_acc] = cif_file
            else:
                # Only include families within the requested range
                pfam_num = int(pfam_acc[2:])
                if 1 <= pfam_num <= num_families:
                    if pfam_acc not in model_index:
                        model_index[pfam_acc] = cif_file

    if num_families is None:
        print(f"Indexed {len(model_index)} Pfam families (all available)")
    else:
        print(f"Indexed {len(model_index)} Pfam families (PF00001 to PF{num_families:05d})")
    return model_index


def create_reseek_database(model_index, work_dir):
    """Create a Reseek .bca database from Pfam AlphaFold models."""
    print("\n=== Creating Reseek database ===")

    db_marker = work_dir / ".db_complete"
    db_file = work_dir / "pfam_models.bca"

    if db_marker.exists() and db_file.exists():
        print("Reseek database already created, skipping")
        return db_file

    # Create a temporary directory with symlinks to all models
    models_subdir = work_dir / "_models"
    models_subdir.mkdir(exist_ok=True, parents=True)

    print(f"Creating symlinks to {len(model_index)} models...")
    for pfam_acc, model_file in sorted(model_index.items()):
        link_name = models_subdir / f"{pfam_acc}.cif"
        if not link_name.exists():
            link_name.symlink_to(model_file)

    print(f"Converting models to Reseek .bca format...")
    print(f"  Input: {models_subdir}")
    print(f"  Output: {db_file}")

    # Run reseek -convert to create .bca database
    cmd = ['reseek', '-convert', str(models_subdir), '-bca', str(db_file)]
    result = run_command(cmd, cwd=work_dir)

    if result.returncode != 0:
        print("Error: Failed to create Reseek database", file=sys.stderr)
        sys.exit(1)

    db_marker.touch()
    print(f"Database created: {db_file}")

    return db_file


def run_allvsall_search(db_file, output_file, work_dir, sensitivity='sensitive', threads=32):
    """Run Reseek all-vs-all search."""
    print("\n=== Running all-vs-all search ===")

    search_marker = work_dir / ".search_complete"
    raw_output = work_dir / "reseek_allvsall_raw.tsv"

    if search_marker.exists() and raw_output.exists():
        print("All-vs-all search already completed, skipping")
        return raw_output

    print(f"Database: {db_file}")
    print(f"Sensitivity: {sensitivity}")
    print(f"Threads: {threads}")
    print("This may take a while for large databases...")

    # Run reseek all-vs-all search
    # When -search uses a database without specifying -db separately, it does all-vs-all
    cmd = [
        'reseek',
        '-search', str(db_file),
        '-db', str(db_file),  # Search database against itself
        f'-{sensitivity}',
        '-output', str(raw_output),
        '-threads', str(threads),
        '-evalue', '10',
        '-columns', 'query+target+qlo+qhi+tlo+thi+ql+tl+pctid+evalue+aq',
    ]

    result = run_command(cmd, cwd=work_dir, check=False)

    if result.returncode != 0:
        print("Error: All-vs-all search failed", file=sys.stderr)
        sys.exit(1)

    search_marker.touch()
    print(f"Search complete! Raw output: {raw_output}")

    return raw_output


def parse_and_format_results(raw_output, output_file):
    """Parse Reseek output and format as HHsearch-like TSV."""
    print("\n=== Parsing and formatting results ===")

    if not Path(raw_output).exists():
        print(f"Error: Raw output file not found: {raw_output}", file=sys.stderr)
        sys.exit(1)

    results = []
    line_count = 0

    with open(raw_output, 'r') as f:
        for line in f:
            line_count += 1
            if line_count % 100000 == 0:
                print(f"  Processed {line_count} lines, {len(results)} valid hits...")

            line = line.strip()
            if not line or line.startswith('#'):
                continue

            fields = line.split('\t')
            if len(fields) < 11:
                continue

            try:
                query_name = fields[0]
                target_name = fields[1]
                qlo = int(fields[2])
                qhi = int(fields[3])
                tlo = int(fields[4])
                thi = int(fields[5])
                ql = int(fields[6])
                tl = int(fields[7])
                pctid = float(fields[8])
                evalue = float(fields[9])
                aq = float(fields[10])

                # Extract Pfam accessions
                query_pfam = extract_pfam_accession(query_name)
                target_pfam = extract_pfam_accession(target_name)

                if not query_pfam or not target_pfam:
                    continue

                # Calculate alignment length
                alnlen = max(qhi - qlo + 1, thi - tlo + 1)

                # Calculate probability from E-value
                if evalue < 1e-100:
                    prob = 100.0
                elif evalue < 1e-10:
                    prob = 99.9
                elif evalue < 0.01:
                    prob = 99.0 + (1 - evalue * 100)
                elif evalue < 1:
                    prob = 90.0 + (1 - evalue) * 9
                else:
                    prob = 100.0 / (1.0 + evalue)

                # Use alignment quality * 100 as score
                score = aq * 100

                # P-value (use evalue as proxy)
                pvalue = evalue

                # SS not available
                ss = 0.0

                # Query and Template ranges
                query_str = f"{qlo}-{qhi}"
                template_str = f"{tlo}-{thi} ({tl})"

                results.append({
                    'query_family': query_pfam,
                    'hit': target_pfam,
                    'prob': prob,
                    'evalue': evalue,
                    'pvalue': pvalue,
                    'score': score,
                    'ss': ss,
                    'cols': alnlen,
                    'query': query_str,
                    'template': template_str
                })

            except (ValueError, IndexError) as e:
                continue

    print(f"Total lines processed: {line_count}")
    print(f"Valid hits found: {len(results)}")

    # Sort by query family, then by E-value
    print("Sorting results...")
    results.sort(key=lambda x: (x['query_family'], x['evalue']))

    # Write output
    print(f"Writing results to {output_file}...")
    with open(output_file, 'w') as out:
        # Header (HHsearch format)
        out.write("QueryFamily\tHit\tProb\tE-value\tP-value\tScore\tSS\tCols\tQuery\tTemplate\n")

        # Write results
        for r in results:
            out.write('\t'.join([
                r['query_family'],
                r['hit'],
                f"{r['prob']:.1f}",
                f"{r['evalue']:.2e}",
                f"{r['pvalue']:.2e}",
                f"{r['score']:.1f}",
                f"{r['ss']:.1f}",
                str(r['cols']),
                r['query'],
                r['template']
            ]) + '\n')

    size_kb = Path(output_file).stat().st_size / 1024
    print(f"Results written: {output_file} ({size_kb:.1f} KB)")


def main():
    parser = argparse.ArgumentParser(
        description="Perform all-vs-all Pfam comparison using Reseek (efficient version)"
    )
    parser.add_argument(
        '-n', type=int, default=None,
        help="Number of Pfam families to process (PF00001 to PFxxxxx). If not specified, uses all available families."
    )
    parser.add_argument(
        '-dir', type=str, required=True,
        help="Results directory name to create"
    )
    parser.add_argument(
        '-sensitivity', type=str, default='sensitive',
        choices=['fast', 'sensitive', 'verysensitive'],
        help="Reseek sensitivity level (default: sensitive)"
    )
    parser.add_argument(
        '-threads', type=int, default=32,
        help="Number of threads (default: 32)"
    )
    parser.add_argument(
        '--models-dir', type=str, default=str(PFAM_MODELS_DIR),
        help=f"Directory containing Pfam AlphaFold models (default: {PFAM_MODELS_DIR})"
    )
    args = parser.parse_args()

    work_dir = Path(args.dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    if args.n is None:
        print(f"Reseek all-vs-all pipeline for all available Pfam families")
    else:
        print(f"Reseek all-vs-all pipeline for {args.n} Pfam families")
    print(f"Working directory: {work_dir}")
    print(f"Models directory: {args.models_dir}")
    print()

    # Build model index
    model_index = build_pfam_model_index(args.models_dir, args.n)

    if not model_index:
        print("Error: No Pfam models found", file=sys.stderr)
        sys.exit(1)

    # Save model index
    index_file = work_dir / "model_index.pkl"
    with open(index_file, 'wb') as f:
        pickle.dump(model_index, f)
    print(f"Saved model index to {index_file}\n")

    # Create Reseek database
    db_file = create_reseek_database(model_index, work_dir)

    # Run all-vs-all search
    raw_output = run_allvsall_search(
        db_file,
        work_dir / "all_results.tsv",
        work_dir,
        args.sensitivity,
        args.threads
    )

    # Parse and format results
    parse_and_format_results(raw_output, work_dir / "all_results.tsv")

    print()
    print("=" * 50)
    print("=== PIPELINE COMPLETE ===")
    print("=" * 50)
    print(f"Results: {work_dir}/all_results.tsv")


if __name__ == "__main__":
    main()
