#!/usr/bin/env python3
"""
TED Cluster SEED Curation Helper Script

Loops over TED cluster directories and guides curator through SEED alignment
refinement, then runs pfbuild to create HMM.

Usage:
    python ted_seed_curation.py [--input-dir DIR]
"""

import sys
import os
import subprocess
import shutil
from pathlib import Path


def run_command(cmd, cwd=None, capture_output=False, wait=True):
    """Execute shell command and optionally capture output"""
    try:
        if capture_output:
            result = subprocess.run(cmd, shell=True, cwd=cwd,
                                    capture_output=True, text=True)
            return result.stdout.strip()
        else:
            if wait:
                result = subprocess.run(cmd, shell=True, cwd=cwd)
                return result.returncode == 0
            else:
                subprocess.Popen(cmd, shell=True, cwd=cwd)
                return True
    except Exception as e:
        print(f"Error running command '{cmd}': {e}", file=sys.stderr)
        return None if capture_output else False


def count_sequences(file_path):
    """Count sequences in a Stockholm alignment file"""
    if not Path(file_path).exists():
        return 0

    count = 0
    with open(file_path, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#') and not line.startswith('//'):
                count += 1
    return count


def is_already_built(dir_path):
    """Check if pfbuild has already been run (HMM file exists)"""
    hmm_path = os.path.join(dir_path, 'HMM')
    return os.path.exists(hmm_path)


def get_cluster_directories(input_dir):
    """Get list of cluster directories sorted by name, skipping already built ones"""
    dirs = []
    skipped = 0
    for d in os.listdir(input_dir):
        dir_path = os.path.join(input_dir, d)
        if os.path.isdir(dir_path):
            # Skip DONE, IGNORE, and other special directories
            if d in ['DONE', 'IGNORE', 'REMOVED']:
                continue
            # Check if it has SEED and FA files
            seed_path = os.path.join(dir_path, 'SEED')
            fa_path = os.path.join(dir_path, 'FA')
            if os.path.exists(seed_path) and os.path.exists(fa_path):
                # Skip if already built
                if is_already_built(dir_path):
                    skipped += 1
                    continue
                dirs.append(d)

    if skipped > 0:
        print(f"Skipped {skipped} directories that already have HMM (pfbuild already run)")

    return sorted(dirs)


def curate_seed(cluster_name, start_dir):
    """Guide curator through SEED alignment curation"""

    cluster_dir = os.path.join(start_dir, cluster_name)

    print(f"\n{'='*60}")
    print(f"Working on: {cluster_name}")
    print(f"{'='*60}")

    # Change to cluster directory
    os.chdir(cluster_dir)

    # Count sequences
    seed_count = count_sequences('SEED')
    fa_count = run_command("grep -c '^>' FA", capture_output=True)

    print(f"Sequences in FA: {fa_count}")
    print(f"Sequences in SEED: {seed_count}")

    # Launch belvu for SEED review
    print("\nLaunching belvu for SEED alignment review...")
    print("Edit the alignment in belvu, save, and close when done.")
    belvu_success = run_command("belvu SEED", wait=True)

    if not belvu_success:
        print("Warning: belvu may have encountered an issue", file=sys.stderr)

    # Present editing options
    while True:
        print("\n" + "="*60)
        print("SEED editing options:")
        print("  y/yes     - Happy with SEED, run pfbuild")
        print("  n/no      - Skip to next family")
        print("  i/ignore  - Move to IGNORE directory")
        print("  b/belvu   - Open belvu again")
        print("  w/wholeseq - Run wholeseq on SEED")
        print("  e/extend  - Extend N/C termini")
        print("  p/pad     - Pad alignment ends")
        print("="*60)

        response = input("\nWhat would you like to do? ").strip().lower()

        if response in ['y', 'yes']:
            # Run pfbuild
            print("\nRunning pfbuild -withpfmake...")
            pfbuild_success = run_command("pfbuild -withpfmake", wait=True)

            if pfbuild_success:
                print("pfbuild completed successfully")

                # Move to DONE directory
                os.chdir(start_dir)
                done_dir = Path("DONE")
                if not done_dir.exists():
                    done_dir.mkdir()
                    print("Created DONE directory")

                try:
                    shutil.move(cluster_name, str(done_dir / cluster_name))
                    print(f"Moved {cluster_name} to DONE/")
                except Exception as e:
                    print(f"Warning: Could not move to DONE: {e}", file=sys.stderr)

                return True
            else:
                print("Warning: pfbuild may have failed", file=sys.stderr)
                # Stay in menu to let user decide what to do

        elif response in ['n', 'no']:
            os.chdir(start_dir)
            print("Skipping this family")
            return False

        elif response in ['i', 'ignore']:
            os.chdir(start_dir)
            ignore_dir = Path("IGNORE")
            if not ignore_dir.exists():
                ignore_dir.mkdir()
                print("Created IGNORE directory")

            try:
                shutil.move(cluster_name, str(ignore_dir / cluster_name))
                print(f"Moved {cluster_name} to IGNORE/")
            except Exception as e:
                print(f"Warning: Could not move to IGNORE: {e}", file=sys.stderr)

            return False

        elif response in ['b', 'belvu']:
            print("\nRelaunching belvu...")
            run_command("belvu SEED", wait=True)

        elif response in ['w', 'wholeseq']:
            print("\nRunning wholeseq on SEED...")
            shutil.copy("SEED", "SEED.beforewholeseq")
            wholeseq_cmd = "wholeseq.pl -align SEED.beforewholeseq -m > SEED.tmp && mv SEED.tmp SEED"
            if run_command(wholeseq_cmd, wait=True):
                print("Wholeseq completed")
                print("Relaunching belvu with updated SEED...")
                run_command("belvu SEED", wait=True)
            else:
                print("Warning: wholeseq may have failed", file=sys.stderr)

        elif response in ['e', 'extend']:
            n = input("How many residues to extend N-terminus? [default 0]: ").strip() or "0"
            c = input("How many residues to extend C-terminus? [default 0]: ").strip() or "0"

            print(f"\nExtending SEED N:{n} C:{c}")
            shutil.copy("SEED", "SEED.beforeextend")
            extend_cmd = f"extend.pl -align SEED.beforeextend -n {n} -c {c} -m > SEED.tmp && mv SEED.tmp SEED"
            if run_command(extend_cmd, wait=True):
                print("Extension completed")
                print("Relaunching belvu with updated SEED...")
                run_command("belvu SEED", wait=True)
            else:
                print("Warning: extend may have failed", file=sys.stderr)

        elif response in ['p', 'pad']:
            print("\nPadding SEED alignment...")
            shutil.copy("SEED", "SEED.beforepad")
            pad_cmd = "pad_ends.pl -align SEED.beforepad > SEED.tmp && mv SEED.tmp SEED"
            if run_command(pad_cmd, wait=True):
                print("Padding completed")
                print("Relaunching belvu with updated SEED...")
                run_command("belvu SEED", wait=True)
            else:
                print("Warning: pad_ends may have failed", file=sys.stderr)

        else:
            print("Unrecognized option, please try again.")

    return False


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="TED Cluster SEED Curation Helper"
    )
    parser.add_argument(
        '--input-dir',
        default='.',
        help='Directory containing cluster subdirectories (default: current directory)'
    )

    args = parser.parse_args()

    input_dir = os.path.abspath(args.input_dir)

    if not os.path.exists(input_dir):
        print(f"Error: Directory not found: {input_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"\n{'='*60}")
    print("TED Cluster SEED Curation Helper")
    print(f"{'='*60}")

    # Get cluster directories
    cluster_dirs = get_cluster_directories(input_dir)

    if not cluster_dirs:
        print("No cluster directories found with SEED and FA files.")
        sys.exit(0)

    print(f"Found {len(cluster_dirs)} cluster directories to process")

    # Process each cluster
    completed = 0
    for i, cluster_name in enumerate(cluster_dirs, 1):
        print(f"\n[{i}/{len(cluster_dirs)}]", file=sys.stderr)
        if curate_seed(cluster_name, input_dir):
            completed += 1

    print(f"\n{'='*60}")
    print(f"Complete! Processed {completed} families")
    if completed > 0:
        print("Completed families have been moved to DONE/")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
