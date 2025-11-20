#!/usr/bin/env python3
"""
Build Reseek database from Pfam AlphaFold models and perform all-vs-all comparison.

This script:
1. Creates a Reseek .bca database from Pfam AlphaFold CIF models
2. Submits SLURM jobs to query batches against the database
3. Aggregates all results in HHsearch-like format

Usage:
  python3 pfam_reseek_pipeline.py -n 28000 -dir ResultDB -batch-size 100 -time 60 -mem 16G -cpus 4
"""

import argparse
import subprocess
import sys
from pathlib import Path
import time
import re
import pickle

# Configuration
PFAM_MODELS_DIR = Path("/nfs/production/agb/interpro/users/typhaine/interpro-pfam-curation-tools/pfam_add_clan_search/results/chopped_cif")
WORKER_SCRIPT = Path(__file__).parent / "pfam_reseek_batch_worker.py"


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
    """
    Extract Pfam accession from a filename.

    Examples:
        PF00001_model.cif -> PF00001
        AF-P12345-F1-model_v6_PF00001.cif -> PF00001
    """
    match = re.search(r'(PF\d{5})', filename)
    if match:
        return match.group(1)
    return None


def build_pfam_model_index(models_dir, num_families):
    """
    Build an index of available Pfam AlphaFold models.

    Returns:
        dict: {pfam_acc: model_file_path}
    """
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
            # Only include families within the requested range
            pfam_num = int(pfam_acc[2:])
            if 1 <= pfam_num <= num_families:
                if pfam_acc not in model_index:
                    model_index[pfam_acc] = cif_file
                else:
                    # If duplicate, keep the first one (could add logic to pick best)
                    pass

    print(f"Indexed {len(model_index)} Pfam families (PF00001 to PF{num_families:05d})")
    return model_index


def create_reseek_database(model_index, work_dir):
    """
    Create a Reseek .bca database from Pfam AlphaFold models.

    Args:
        model_index: dict of {pfam_acc: model_file_path}
        work_dir: working directory for database

    Returns:
        Path to the .bca database file
    """
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


def submit_slurm_job(batch_id, pfam_families, results_dir, work_dir, db_file,
                     time_limit, mem, cpus, sensitivity):
    """Submit a SLURM job for a batch and return job ID."""

    # Create a file with the list of families for this batch
    batch_families_file = results_dir / f"batch_{batch_id:03d}_families.txt"
    with open(batch_families_file, 'w') as f:
        for pfam_acc in pfam_families:
            f.write(f"{pfam_acc}\n")

    job_script = f"""#!/bin/bash
#SBATCH --job-name=reseek_query_{batch_id:03d}
#SBATCH --time={time_limit}
#SBATCH --mem={mem}
#SBATCH --cpus-per-task={cpus}
#SBATCH --output={results_dir}/batch_{batch_id:03d}.log
#SBATCH --error={results_dir}/batch_{batch_id:03d}.err

python3 {WORKER_SCRIPT} -batch-id {batch_id} -families {batch_families_file} -db {db_file} -dir {results_dir} -sensitivity {sensitivity} -threads {cpus}
"""

    # Write job script to temp file
    script_file = results_dir / f"batch_{batch_id:03d}.sh"
    with open(script_file, 'w') as f:
        f.write(job_script)
    script_file.chmod(0o755)

    # Submit job
    result = subprocess.run(
        ['sbatch', str(script_file)],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print(f"Error submitting batch {batch_id}: {result.stderr}", file=sys.stderr)
        return None

    # Extract job ID from output
    match = re.search(r'Submitted batch job (\d+)', result.stdout)
    if match:
        job_id = match.group(1)
        print(f"Batch {batch_id:03d} ({len(pfam_families)} families): Submitted as job {job_id}")
        return job_id

    print(f"Error: Could not parse job ID from sbatch output: {result.stdout}", file=sys.stderr)
    return None


def wait_for_jobs(job_ids, check_interval=30):
    """Wait for all SLURM jobs to complete."""
    print()
    print(f"Waiting for {len(job_ids)} jobs to complete...")

    remaining = set(job_ids)
    start_time = time.time()

    while remaining:
        # Check job status
        result = subprocess.run(
            ['squeue', '-j', ','.join(remaining), '-h', '-o', '%i,%T'],
            capture_output=True,
            text=True
        )

        if result.returncode == 0 and result.stdout.strip():
            # Parse job statuses
            active = set()
            for line in result.stdout.strip().split('\n'):
                if line.strip():
                    job_id, status = line.split(',')
                    if status in ['PENDING', 'RUNNING']:
                        active.add(job_id)
            remaining = active
        else:
            # No output means all jobs finished
            remaining = set()

        if remaining:
            elapsed = int(time.time() - start_time)
            print(f"  {len(remaining)} jobs still running ({elapsed}s elapsed)...", end='\r')
            time.sleep(check_interval)

    elapsed = int(time.time() - start_time)
    print(f"\n  All jobs completed ({elapsed}s elapsed)")
    print()


def concatenate_batch_results(results_dir, num_batches):
    """Concatenate all batch TSV files into a master file."""
    print("=== Concatenating batch results ===")

    batch_files = []
    for batch_id in range(num_batches):
        batch_file = results_dir / f"batch_{batch_id:03d}_results.tsv"
        if batch_file.exists():
            batch_files.append(batch_file)
        else:
            print(f"Warning: {batch_file} not found")

    if not batch_files:
        print("Error: No batch result files found", file=sys.stderr)
        return False

    master_file = results_dir / "all_results.tsv"
    try:
        with open(master_file, 'w') as out_f:
            # Write header (HHsearch-like format)
            out_f.write("QueryFamily\tHit\tProb\tE-value\tP-value\tScore\tSS\tCols\tQuery\tTemplate\n")

            # Concatenate all batch files (skip individual headers)
            for batch_file in sorted(batch_files):
                with open(batch_file, 'r') as in_f:
                    next(in_f)  # Skip header
                    for line in in_f:
                        out_f.write(line)

        size_kb = master_file.stat().st_size / 1024
        print(f"Created all_results.tsv ({size_kb:.1f}K)")
        return True
    except Exception as e:
        print(f"Error concatenating results: {e}", file=sys.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Build Reseek database and perform all-vs-all Pfam comparison using SLURM"
    )
    parser.add_argument(
        '-n', type=int, required=True,
        help="Number of Pfam families to process (PF00001 to PFxxxxx)"
    )
    parser.add_argument(
        '-dir', type=str, required=True,
        help="Results directory name to create"
    )
    parser.add_argument(
        '-batch-size', type=int, default=100,
        help="Number of families per query batch (default: 100)"
    )
    parser.add_argument(
        '-time', type=int, default=60,
        help="Time limit per query job in minutes (default: 60)"
    )
    parser.add_argument(
        '-mem', type=str, default='16G',
        help="Memory per job (default: 16G)"
    )
    parser.add_argument(
        '-cpus', type=int, default=4,
        help="CPUs per job (default: 4)"
    )
    parser.add_argument(
        '-sensitivity', type=str, default='sensitive',
        choices=['fast', 'sensitive', 'verysensitive'],
        help="Reseek sensitivity level (default: sensitive)"
    )
    parser.add_argument(
        '--models-dir', type=str, default=str(PFAM_MODELS_DIR),
        help=f"Directory containing Pfam AlphaFold models (default: {PFAM_MODELS_DIR})"
    )
    args = parser.parse_args()

    work_dir = Path(args.dir)
    batch_size = args.batch_size

    # Validate worker script exists
    if not WORKER_SCRIPT.exists():
        print(f"Error: Worker script not found at {WORKER_SCRIPT}", file=sys.stderr)
        sys.exit(1)

    # Create directories
    work_dir.mkdir(parents=True, exist_ok=True)
    results_dir = work_dir / "query_results"
    results_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reseek pipeline for {args.n} Pfam families")
    print(f"Working directory: {work_dir}")
    print(f"Models directory: {args.models_dir}")
    print()

    # Build model index
    model_index = build_pfam_model_index(args.models_dir, args.n)

    if not model_index:
        print("Error: No Pfam models found", file=sys.stderr)
        sys.exit(1)

    # Save model index for worker scripts
    index_file = work_dir / "model_index.pkl"
    with open(index_file, 'wb') as f:
        pickle.dump(model_index, f)
    print(f"Saved model index to {index_file}")
    print()

    # Create Reseek database
    db_file = create_reseek_database(model_index, work_dir)

    # Organize batches
    pfam_families = sorted(model_index.keys())
    num_batches = (len(pfam_families) + batch_size - 1) // batch_size

    print()
    print("=" * 50)
    print(f"Submitting {num_batches} query batches")
    print()

    # Check which batches already have results
    completed_batches = set()
    if results_dir.exists():
        for batch_file in results_dir.glob("batch_*_results.tsv"):
            match = re.search(r'batch_(\d+)_results', batch_file.name)
            if match:
                completed_batches.add(int(match.group(1)))

    if completed_batches:
        print(f"Found {len(completed_batches)} already completed batches")

    # Submit jobs
    print("=== Submitting SLURM query jobs ===")

    job_ids = []
    for batch_id in range(num_batches):
        if batch_id in completed_batches:
            continue

        start_idx = batch_id * batch_size
        end_idx = min((batch_id + 1) * batch_size, len(pfam_families))
        batch_families = pfam_families[start_idx:end_idx]

        job_id = submit_slurm_job(
            batch_id, batch_families, results_dir, work_dir, db_file,
            args.time, args.mem, args.cpus, args.sensitivity
        )
        if job_id:
            job_ids.append(job_id)

    if not job_ids:
        if completed_batches:
            print("All batches already completed!")
        else:
            print("Error: Failed to submit any jobs", file=sys.stderr)
            sys.exit(1)
    else:
        print(f"Submitted {len(job_ids)} new jobs ({len(completed_batches)} already done)")

        # Wait for all jobs
        wait_for_jobs(job_ids)

    # Check for failures
    print("=== Checking job status ===")
    failed_jobs = []
    for batch_id in range(num_batches):
        result_file = results_dir / f"batch_{batch_id:03d}_results.tsv"
        if not result_file.exists():
            failed_jobs.append(batch_id)
            print(f"Batch {batch_id:03d}: FAILED (no results file)")
        else:
            num_hits = sum(1 for line in open(result_file)) - 1  # Subtract header
            print(f"Batch {batch_id:03d}: OK ({num_hits} hits)")

    print()
    if failed_jobs:
        print(f"WARNING: {len(failed_jobs)} batches failed: {failed_jobs}", file=sys.stderr)
        print("Check log files in results directory for details")
        print()

    # Concatenate results
    concatenate_batch_results(results_dir, num_batches)

    print()
    print("=" * 50)
    print("=== PIPELINE COMPLETE ===")
    print("=" * 50)
    print(f"Results: {results_dir}/all_results.tsv")


if __name__ == "__main__":
    main()
