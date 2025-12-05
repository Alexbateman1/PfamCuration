#!/usr/bin/env python3
"""
Automated SEED Curation for TED Cluster Alignments

This script automatically processes TED cluster directories to create
good quality SEED alignments by:
1. Trimming ragged alignment ends (trim_ali.py)
2. Padding short truncations (pad_ends.pl)
3. Removing partial sequences (belvu -P)
4. Making non-redundant at 80% or higher (belvu -n)
5. Running pfbuild to create HMM

Usage:
    auto_seed_curation.py [options]

    # Process all directories in current location
    auto_seed_curation.py

    # Process specific directory
    auto_seed_curation.py --input-dir /path/to/clusters

    # Dry run (show what would be done)
    auto_seed_curation.py --dry-run

    # Process only first N directories
    auto_seed_curation.py --limit 10
"""

import argparse
import os
import sys
import subprocess
import shutil
from pathlib import Path


def run_command(cmd, cwd=None, capture_output=False):
    """Execute shell command and return success status or output."""
    try:
        if capture_output:
            result = subprocess.run(
                cmd, shell=True, cwd=cwd,
                capture_output=True, text=True
            )
            return result.stdout.strip(), result.returncode == 0
        else:
            result = subprocess.run(cmd, shell=True, cwd=cwd)
            return result.returncode == 0
    except Exception as e:
        print(f"  Error running command '{cmd}': {e}", file=sys.stderr)
        return (None, False) if capture_output else False


def count_sequences(filepath):
    """Count number of sequences in an alignment file."""
    if not os.path.exists(filepath):
        return 0

    count = 0
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            # Count non-empty, non-comment lines that look like sequence entries
            if line and not line.startswith('#') and not line.startswith('//'):
                # Check if it looks like a sequence line (has accession/coords format)
                if '/' in line and '-' in line.split()[0]:
                    count += 1
    return count


def get_cluster_directories(input_dir):
    """
    Get list of cluster directories that need processing.

    Skips directories that:
    - Already have HMM file (curator has fixed)
    - Are named DONE, IGNORE, REMOVED
    - Don't have a SEED file
    """
    dirs = []
    skipped_hmm = 0
    skipped_no_seed = 0

    for d in sorted(os.listdir(input_dir)):
        dir_path = os.path.join(input_dir, d)

        if not os.path.isdir(dir_path):
            continue

        # Skip special directories
        if d in ['DONE', 'IGNORE', 'REMOVED', '.git']:
            continue

        seed_path = os.path.join(dir_path, 'SEED')
        hmm_path = os.path.join(dir_path, 'HMM')

        # Skip if no SEED file
        if not os.path.exists(seed_path):
            skipped_no_seed += 1
            continue

        # Skip if HMM already exists (curator has already processed)
        if os.path.exists(hmm_path):
            skipped_hmm += 1
            continue

        dirs.append(d)

    if skipped_hmm > 0:
        print(f"Skipped {skipped_hmm} directories with existing HMM")
    if skipped_no_seed > 0:
        print(f"Skipped {skipped_no_seed} directories without SEED file")

    return dirs


def process_cluster(cluster_name, start_dir, dry_run=False, verbose=False):
    """
    Process a single cluster directory to create a good quality SEED.

    Returns:
        True if processing succeeded, False otherwise
    """
    cluster_dir = os.path.join(start_dir, cluster_name)

    print(f"\n{'='*60}")
    print(f"Processing: {cluster_name}")
    print(f"{'='*60}")

    # Work in the cluster directory
    os.chdir(cluster_dir)

    try:
        # Step 1: Move SEED to SEED.orig
        if os.path.exists('SEED.orig'):
            print("  SEED.orig already exists, using existing backup")
        else:
            print("  Moving SEED -> SEED.orig")
            if not dry_run:
                shutil.copy('SEED', 'SEED.orig')

        seed_input = 'SEED.orig'
        initial_count = count_sequences(seed_input)
        print(f"  Initial sequences: {initial_count}")

        if initial_count == 0:
            print("  ERROR: No sequences in SEED.orig, skipping")
            os.chdir(start_dir)
            return False

        # Step 2: Trim alignment ends
        print("  Running trim_ali.py...")
        if not dry_run:
            # Use gradient method for cleaner trimming
            success = run_command(f"trim_ali.py {seed_input} --method gradient > 1")
            if not success or not os.path.exists('1') or os.path.getsize('1') == 0:
                print("  ERROR: trim_ali.py failed")
                os.chdir(start_dir)
                return False

        if verbose and not dry_run:
            print(f"    After trim: {count_sequences('1')} sequences")

        # Step 3: Pad short truncations
        print("  Running pad_ends.pl...")
        if not dry_run:
            success = run_command("pad_ends.pl -align 1 > 2")
            if not success or not os.path.exists('2') or os.path.getsize('2') == 0:
                print("  WARNING: pad_ends.pl may have failed, using trimmed file")
                shutil.copy('1', '2')

        if verbose and not dry_run:
            print(f"    After pad: {count_sequences('2')} sequences")

        # Step 4: Remove partial sequences
        print("  Removing partial sequences (belvu -P)...")
        if not dry_run:
            success = run_command("belvu -P 2 -o mul > 3")
            if not success or not os.path.exists('3') or os.path.getsize('3') == 0:
                print("  WARNING: belvu -P may have failed, using padded file")
                shutil.copy('2', '3')

        if verbose and not dry_run:
            print(f"    After remove partials: {count_sequences('3')} sequences")

        # Step 5: Make non-redundant at 80% (or higher if needed)
        print("  Making non-redundant (belvu -n)...")

        if not dry_run:
            identity_threshold = 80
            max_threshold = 100
            final_count = 0

            while identity_threshold <= max_threshold:
                success = run_command(f"belvu 3 -n {identity_threshold} -o mul > 4")

                if success and os.path.exists('4'):
                    final_count = count_sequences('4')
                    print(f"    At {identity_threshold}% identity: {final_count} sequences")

                    if final_count > 1:
                        break
                    elif final_count == 1:
                        # Only 1 sequence, try higher threshold
                        identity_threshold += 5
                    else:
                        # 0 sequences, something went wrong
                        print("  WARNING: belvu -n produced no sequences")
                        break
                else:
                    print(f"  WARNING: belvu -n {identity_threshold} failed")
                    break

            if final_count < 1:
                print("  ERROR: Could not produce valid alignment")
                os.chdir(start_dir)
                return False

            if final_count == 1:
                print("  WARNING: Only 1 sequence remains even at 100% identity")
                # Still proceed - might be a single-sequence family

        # Step 6: Move result to SEED
        print("  Creating new SEED...")
        if not dry_run:
            if os.path.exists('4') and os.path.getsize('4') > 0:
                shutil.move('4', 'SEED')
            else:
                print("  ERROR: No valid output file to use as SEED")
                os.chdir(start_dir)
                return False

        final_seed_count = count_sequences('SEED') if not dry_run else "N/A"
        print(f"  Final SEED: {final_seed_count} sequences")

        # Step 7: Run pfbuild
        print("  Running pfbuild -withpfmake...")
        if not dry_run:
            success = run_command("pfbuild -withpfmake")
            if not success:
                print("  WARNING: pfbuild may have failed")
            else:
                print("  pfbuild completed successfully")

        # Keep temporary files (1, 2, 3) for debugging

        print(f"  DONE: {cluster_name}")
        os.chdir(start_dir)
        return True

    except Exception as e:
        print(f"  ERROR: {e}", file=sys.stderr)
        os.chdir(start_dir)
        return False


def main():
    parser = argparse.ArgumentParser(
        description='Automated SEED Curation for TED Cluster Alignments',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    %(prog)s                      # Process all directories
    %(prog)s --dry-run            # Show what would be done
    %(prog)s --limit 5            # Process first 5 directories
    %(prog)s --input-dir /path    # Process directories in specific path
        """
    )
    parser.add_argument(
        '--input-dir', '-i',
        default='.',
        help='Directory containing cluster subdirectories (default: current directory)'
    )
    parser.add_argument(
        '--dry-run', '-n',
        action='store_true',
        help='Show what would be done without making changes'
    )
    parser.add_argument(
        '--limit', '-l',
        type=int,
        default=None,
        help='Process only first N directories'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Show detailed progress information'
    )
    parser.add_argument(
        '--continue-on-error',
        action='store_true',
        help='Continue processing other directories if one fails'
    )

    args = parser.parse_args()

    input_dir = os.path.abspath(args.input_dir)

    if not os.path.exists(input_dir):
        print(f"Error: Directory not found: {input_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"\n{'='*60}")
    print("Automated SEED Curation for TED Clusters")
    print(f"{'='*60}")
    print(f"Working directory: {input_dir}")

    if args.dry_run:
        print("DRY RUN - no changes will be made")

    # Get directories to process
    cluster_dirs = get_cluster_directories(input_dir)

    if not cluster_dirs:
        print("\nNo cluster directories found to process.")
        sys.exit(0)

    # Apply limit if specified
    if args.limit:
        cluster_dirs = cluster_dirs[:args.limit]

    print(f"\nFound {len(cluster_dirs)} directories to process")

    # Process each cluster
    processed = 0
    succeeded = 0
    failed = 0

    for i, cluster_name in enumerate(cluster_dirs, 1):
        print(f"\n[{i}/{len(cluster_dirs)}]", file=sys.stderr)

        try:
            success = process_cluster(
                cluster_name, input_dir,
                dry_run=args.dry_run,
                verbose=args.verbose
            )
            processed += 1

            if success:
                succeeded += 1
            else:
                failed += 1
                if not args.continue_on_error and not args.dry_run:
                    print("\nStopping due to error. Use --continue-on-error to continue.")
                    break

        except KeyboardInterrupt:
            print("\n\nInterrupted by user")
            break
        except Exception as e:
            print(f"\nUnexpected error: {e}", file=sys.stderr)
            failed += 1
            if not args.continue_on_error:
                break

    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"Directories processed: {processed}")
    print(f"Succeeded: {succeeded}")
    print(f"Failed: {failed}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
