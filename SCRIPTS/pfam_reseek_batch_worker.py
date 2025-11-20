#!/usr/bin/env python3
"""
Worker script for Reseek batch processing.

This script:
1. Takes a batch of Pfam families
2. Runs Reseek search for each family against the database
3. Parses and formats results in HHsearch-like format
4. Outputs results as TSV

Usage:
  python3 pfam_reseek_batch_worker.py -batch-id 0 -families families.txt -db pfam.bca -dir results/ -sensitivity sensitive -threads 4
"""

import argparse
import subprocess
import sys
import re
import tempfile
from pathlib import Path


def extract_pfam_accession(name):
    """
    Extract Pfam accession from a name string.

    Examples:
        PF00001 -> PF00001
        PF00001.cif -> PF00001
        AF-P12345-F1-model_v6_PF00001 -> PF00001
    """
    match = re.search(r'(PF\d{5})', name)
    if match:
        return match.group(1)
    return name


def parse_alignment_range(range_str):
    """
    Parse alignment range from Reseek output.

    Examples:
        "1-100" -> (1, 100)
        "50-200" -> (50, 200)
    """
    if '-' in range_str:
        parts = range_str.split('-')
        if len(parts) == 2:
            try:
                return int(parts[0]), int(parts[1])
            except ValueError:
                pass
    return None, None


def run_reseek_search(query_file, db_file, output_file, sensitivity='sensitive', threads=4):
    """
    Run Reseek search for a single query against the database.

    Args:
        query_file: path to query structure file
        db_file: path to Reseek .bca database
        output_file: where to write results
        sensitivity: Reseek sensitivity level (fast/sensitive/verysensitive)
        threads: number of threads

    Returns:
        bool: True if successful, False otherwise
    """
    cmd = [
        'reseek',
        '-search', str(query_file),
        '-db', str(db_file),
        f'-{sensitivity}',
        '-output', str(output_file),
        '-threads', str(threads),
        '-evalue', '10',  # Match HHsearch default
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            print(f"  Error: Reseek search failed: {result.stderr}", file=sys.stderr)
            return False
        return True
    except Exception as e:
        print(f"  Error: Failed to run Reseek: {e}", file=sys.stderr)
        return False


def parse_reseek_output(reseek_output_file, query_pfam_acc):
    """
    Parse Reseek output and convert to HHsearch-like format.

    Reseek output columns (default):
        query, target, identity, alnlen, qstart, qend, tstart, tend, evalue, bitscore

    Returns:
        list of dicts with parsed results
    """
    results = []

    if not Path(reseek_output_file).exists():
        return results

    with open(reseek_output_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            fields = line.split('\t')
            if len(fields) < 10:
                continue

            try:
                query_name = fields[0]
                target_name = fields[1]
                identity = float(fields[2])
                alnlen = int(fields[3])
                qstart = int(fields[4])
                qend = int(fields[5])
                tstart = int(fields[6])
                tend = int(fields[7])
                evalue = float(fields[8])
                bitscore = float(fields[9])

                # Extract Pfam accessions
                query_pfam = extract_pfam_accession(query_name)
                target_pfam = extract_pfam_accession(target_name)

                # Calculate probability (rough approximation from E-value)
                # P = 1 - e^(-E)  but for very small E, P ≈ E
                # For display, convert to percentage-like value
                if evalue < 1e-100:
                    prob = 100.0
                elif evalue < 1:
                    prob = 100.0 * (1.0 - pow(10, min(evalue, 0)))
                else:
                    prob = 100.0 / (1.0 + evalue)

                # P-value (not provided by Reseek, use E-value as proxy)
                pvalue = evalue

                # SS (secondary structure score, not provided by Reseek)
                ss = 0.0

                # Query and Template ranges
                query_range = f"{qstart}-{qend}"
                template_range = f"{tstart}-{tend}"

                # Query length (approximate from alignment range)
                qlen = qend - qstart + 1
                tlen = tend - tstart + 1

                query_str = f"{query_range}"
                template_str = f"{template_range} ({tlen})"

                results.append({
                    'query_family': query_pfam,
                    'hit': target_pfam,
                    'prob': prob,
                    'evalue': evalue,
                    'pvalue': pvalue,
                    'score': bitscore,
                    'ss': ss,
                    'cols': alnlen,
                    'query': query_str,
                    'template': template_str
                })

            except (ValueError, IndexError) as e:
                print(f"  Warning: Failed to parse line: {line[:100]}", file=sys.stderr)
                continue

    return results


def process_batch(families_file, db_file, models_dir, results_dir, batch_id,
                  sensitivity='sensitive', threads=4):
    """
    Process a batch of Pfam families.

    Args:
        families_file: file containing list of Pfam accessions (one per line)
        db_file: Reseek database file
        models_dir: directory containing model CIF files
        results_dir: output directory
        batch_id: batch identifier
        sensitivity: Reseek sensitivity level
        threads: number of threads
    """
    print(f"Processing batch {batch_id}")

    # Read families list
    with open(families_file, 'r') as f:
        families = [line.strip() for line in f if line.strip()]

    print(f"  {len(families)} families to process")

    # Collect all results
    all_results = []

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        for idx, pfam_acc in enumerate(families, 1):
            print(f"  [{idx}/{len(families)}] Searching {pfam_acc}...", end=' ')

            # Find the model file
            query_file = models_dir / f"{pfam_acc}.cif"

            if not query_file.exists():
                print(f"Model not found, skipping")
                continue

            # Run Reseek search
            output_file = tmp_path / f"{pfam_acc}_reseek.tsv"

            success = run_reseek_search(
                query_file, db_file, output_file,
                sensitivity=sensitivity, threads=threads
            )

            if not success:
                print(f"Search failed")
                continue

            # Parse results
            results = parse_reseek_output(output_file, pfam_acc)
            all_results.extend(results)

            print(f"{len(results)} hits")

    # Write batch results
    results_file = Path(results_dir) / f"batch_{batch_id:03d}_results.tsv"

    print(f"\nWriting {len(all_results)} total hits to {results_file}")

    with open(results_file, 'w') as out:
        # Header (HHsearch format)
        out.write("QueryFamily\tHit\tProb\tE-value\tP-value\tScore\tSS\tCols\tQuery\tTemplate\n")

        # Sort by query family, then by E-value
        all_results.sort(key=lambda x: (x['query_family'], x['evalue']))

        # Write results
        for r in all_results:
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

    print(f"Batch {batch_id} complete!")


def main():
    parser = argparse.ArgumentParser(
        description="Worker script for Reseek batch processing"
    )
    parser.add_argument(
        '-batch-id', type=int, required=True,
        help="Batch identifier"
    )
    parser.add_argument(
        '-families', type=str, required=True,
        help="File containing list of Pfam families for this batch"
    )
    parser.add_argument(
        '-db', type=str, required=True,
        help="Reseek .bca database file"
    )
    parser.add_argument(
        '-dir', type=str, required=True,
        help="Results directory"
    )
    parser.add_argument(
        '-sensitivity', type=str, default='sensitive',
        choices=['fast', 'sensitive', 'verysensitive'],
        help="Reseek sensitivity level (default: sensitive)"
    )
    parser.add_argument(
        '-threads', type=int, default=4,
        help="Number of threads (default: 4)"
    )
    args = parser.parse_args()

    # Get the models directory from the database path
    db_path = Path(args.db)
    work_dir = db_path.parent
    models_dir = work_dir / "_models"

    if not models_dir.exists():
        print(f"Error: Models directory not found: {models_dir}", file=sys.stderr)
        sys.exit(1)

    if not Path(args.db).exists():
        print(f"Error: Database file not found: {args.db}", file=sys.stderr)
        sys.exit(1)

    if not Path(args.families).exists():
        print(f"Error: Families file not found: {args.families}", file=sys.stderr)
        sys.exit(1)

    process_batch(
        args.families,
        args.db,
        models_dir,
        args.dir,
        args.batch_id,
        args.sensitivity,
        args.threads
    )


if __name__ == "__main__":
    main()
