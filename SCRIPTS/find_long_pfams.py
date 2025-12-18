#!/usr/bin/env python3
"""
Query the pfam_live database to find Pfam families with model_length > 600,
then for each family:
  1. Check out the family with pfco
  2. Run add_image.py
  3. Run swissprot.pl

Usage:
  find_long_pfams.py [--dry-run] [--min-length N]
"""

import argparse
import subprocess
import sys
import os


def query_long_families(min_length=600):
    """
    Query pfam_live database for families with model_length > min_length.

    Returns:
        List of tuples: (pfamA_acc, pfamA_id, type)
    """
    query = f"""SELECT pfamA_acc, pfamA_id, type
                FROM pfamA
                WHERE model_length > {min_length}"""

    cmd = [
        'mysql',
        '--defaults-file=' + os.path.expanduser('~/.my.cnf'),
        'pfam_live',
        '--quick',
        '-N',  # Skip column names
        '-e', query
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: MySQL query failed: {e.stderr}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print("ERROR: mysql command not found. Is MySQL client installed?", file=sys.stderr)
        sys.exit(1)

    families = []
    for line in result.stdout.strip().split('\n'):
        if line:
            parts = line.split('\t')
            if len(parts) >= 3:
                pfamA_acc, pfamA_id, family_type = parts[0], parts[1], parts[2]
                families.append((pfamA_acc, pfamA_id, family_type))

    return families


def run_command(cmd, description, dry_run=False):
    """
    Run a shell command, printing status.

    Args:
        cmd: Command string or list
        description: Human-readable description
        dry_run: If True, only print what would be done

    Returns:
        True if successful (or dry run), False otherwise
    """
    if dry_run:
        print(f"  [DRY RUN] Would run: {cmd}")
        return True

    print(f"  Running: {cmd}")
    try:
        result = subprocess.run(
            cmd,
            shell=isinstance(cmd, str),
            capture_output=True,
            text=True
        )
        if result.returncode != 0:
            print(f"  WARNING: {description} failed with return code {result.returncode}")
            if result.stderr:
                print(f"  stderr: {result.stderr[:500]}")
            return False
        return True
    except Exception as e:
        print(f"  ERROR: {description} failed: {e}")
        return False


def process_family(pfamA_acc, pfamA_id, family_type, dry_run=False):
    """
    Process a single Pfam family:
      1. pfco PFXXXXX
      2. add_image.py PFXXXXX
      3. cd PFXXXXX && swissprot.pl

    Args:
        pfamA_acc: Pfam accession (e.g., PF00001)
        pfamA_id: Pfam identifier (e.g., 7tm_1)
        family_type: Type of the family
        dry_run: If True, only print what would be done
    """
    print(f"\nProcessing {pfamA_acc} ({pfamA_id}) - Type: {family_type}")

    # Step 1: Check out the family
    success = run_command(
        f"pfco {pfamA_acc}",
        f"checkout of {pfamA_acc}",
        dry_run
    )
    if not success and not dry_run:
        print(f"  Skipping remaining steps for {pfamA_acc}")
        return False

    # Step 2: Run add_image.py
    success = run_command(
        f"add_image.py {pfamA_acc}",
        f"add_image.py for {pfamA_acc}",
        dry_run
    )

    # Step 3: Run swissprot.pl inside the family directory
    family_dir = pfamA_acc
    if dry_run:
        print(f"  [DRY RUN] Would run: cd {family_dir} && swissprot.pl && cd ..")
    else:
        if os.path.isdir(family_dir):
            original_dir = os.getcwd()
            try:
                os.chdir(family_dir)
                run_command("swissprot.pl", f"swissprot.pl in {family_dir}", dry_run)
            finally:
                os.chdir(original_dir)
        else:
            print(f"  WARNING: Directory {family_dir} not found, skipping swissprot.pl")
            return False

    return True


def main():
    parser = argparse.ArgumentParser(
        description='Find Pfam families with long models and process them',
        epilog="""
This script queries pfam_live for families with model_length > 600 (or custom value),
then for each family:
  1. Runs: pfco PFXXXXX
  2. Runs: add_image.py PFXXXXX
  3. Runs: swissprot.pl (inside the family directory)

Examples:
  %(prog)s --dry-run           # Show what would be done
  %(prog)s                      # Process all families with model_length > 600
  %(prog)s --min-length 800     # Process families with model_length > 800
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print what would be done without executing commands'
    )
    parser.add_argument(
        '--min-length',
        type=int,
        default=600,
        help='Minimum model length threshold (default: 600)'
    )

    args = parser.parse_args()

    print("=" * 60)
    print("Finding Pfam families with model_length > {}".format(args.min_length))
    print("=" * 60)

    if args.dry_run:
        print("\n*** DRY RUN MODE - No commands will be executed ***\n")

    # Query the database
    print("\nQuerying pfam_live database...")
    families = query_long_families(args.min_length)

    print(f"Found {len(families)} families with model_length > {args.min_length}")

    if not families:
        print("No families to process.")
        return 0

    # Print summary
    print("\nFamilies to process:")
    print("-" * 50)
    for pfamA_acc, pfamA_id, family_type in families:
        print(f"  {pfamA_acc}\t{pfamA_id}\t{family_type}")
    print("-" * 50)

    # Process each family
    success_count = 0
    fail_count = 0

    for pfamA_acc, pfamA_id, family_type in families:
        if process_family(pfamA_acc, pfamA_id, family_type, args.dry_run):
            success_count += 1
        else:
            fail_count += 1

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total families: {len(families)}")
    print(f"Successful: {success_count}")
    print(f"Failed: {fail_count}")

    if args.dry_run:
        print("\n*** This was a dry run - no commands were actually executed ***")

    return 0 if fail_count == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
