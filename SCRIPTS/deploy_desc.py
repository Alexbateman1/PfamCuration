#!/usr/bin/env python3
"""
Move DESC files from the DESC directory to their correct locations 
based on the manifest.json file.

This script:
1. Reads the manifest.json to determine target directories
2. Moves DESC files like P09154_TED01.DESC to the correct location as DESC
3. Backs up existing DESC files to DESC.old
"""

import json
import os
import shutil
from pathlib import Path


def extract_family_prefix(family_id):
    """
    Extract the base family identifier from a family_id.
    
    Examples:
        P39401_TED03/Iterate/Iterate/Iterate -> P39401_TED03
        P37653_TED01/Iterate/Iterate/Iterate -> P37653_TED01
        P11349/Iterate -> P11349
    """
    return family_id.split('/')[0]


def find_desc_file(desc_dir, family_prefix):
    """
    Find the DESC file matching the family prefix.
    Handles both P09154_TED01.DESC and P09154.DESC formats.
    """
    desc_path = desc_dir / f"{family_prefix}.DESC"
    if desc_path.exists():
        return desc_path
    
    # If not found with exact match, try without TED suffix
    if '_TED' in family_prefix:
        base_prefix = family_prefix.split('_TED')[0]
        desc_path = desc_dir / f"{base_prefix}.DESC"
        if desc_path.exists():
            return desc_path
    
    return None


def move_desc_file(desc_file, target_dir, base_dir, dry_run=False):
    """
    Move a DESC file to its target directory, backing up any existing DESC.
    
    Args:
        desc_file: Path to the source DESC file
        target_dir: Target directory path (relative to base_dir)
        base_dir: Base directory for curation
        dry_run: If True, only print what would be done
    """
    target_path = base_dir / target_dir
    target_desc = target_path / "DESC"
    
    if not target_path.exists():
        print(f"  WARNING: Target directory does not exist: {target_path}")
        print(f"  Skipping...")
        return False
    
    # Check if DESC exists in target - we require this
    if not target_desc.exists():
        print(f"  WARNING: No existing DESC file in: {target_path}")
        print(f"  Skipping...")
        return False
    
    # Backup existing DESC
    backup_path = target_path / "DESC.old"
    
    if dry_run:
        print(f"  Would backup: {target_desc} -> {backup_path}")
        print(f"  Would move: {desc_file} -> {target_desc}")
    else:
        shutil.copy2(target_desc, backup_path)
        print(f"  Backed up existing DESC to: DESC.old")
        shutil.copy2(desc_file, target_desc)
        print(f"  Moved: {desc_file.name} -> {target_dir}/DESC")
    
    return True


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Move DESC files to their correct locations based on manifest.json',
        epilog="""
Examples:
  # Navigate to your curation directory, then run:
  cd /homes/agb/Curation/PfamCuration/DATA/EcoliK12
  %(prog)s --manifest manifest.json --dry-run
  %(prog)s --manifest manifest.json

  # Or specify full path to manifest:
  %(prog)s --manifest /path/to/manifest.json --dry-run
  
The script expects:
  - A manifest.json file (required argument)
  - A DESC/ subdirectory in the same location as manifest.json
  - Target family directories relative to current working directory
  - Existing DESC files in target directories (will backup to DESC.old)
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--manifest',
        type=Path,
        required=True,
        help='Path to manifest.json file (REQUIRED)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print what would be done without actually moving files'
    )
    
    args = parser.parse_args()
    
    # Set paths
    base_dir = Path.cwd()
    manifest_path = args.manifest
    
    # DESC directory is in the same location as the manifest
    desc_dir = manifest_path.parent / 'DESC'
    
    # Validate paths
    if not manifest_path.exists():
        print(f"ERROR: Manifest file does not exist: {manifest_path}")
        return 1
    
    if not desc_dir.exists():
        print(f"ERROR: DESC directory does not exist: {desc_dir}")
        print(f"       Expected DESC/ subdirectory next to manifest.json")
        return 1
    
    # Load manifest
    print(f"Base directory (cwd): {base_dir}")
    print(f"Reading manifest from: {manifest_path}")
    print(f"Reading DESC files from: {desc_dir}")
    with open(manifest_path, 'r') as f:
        manifest = json.load(f)
    
    families = manifest.get('families', [])
    print(f"Found {len(families)} families in manifest")
    
    if args.dry_run:
        print("\n*** DRY RUN MODE - No files will be moved ***\n")
    
    # Process each family
    moved_count = 0
    not_found_count = 0
    skipped_count = 0
    
    for family_entry in families:
        family_id = family_entry['family_id']
        family_prefix = extract_family_prefix(family_id)
        
        print(f"\nProcessing: {family_id}")
        print(f"  Looking for: {family_prefix}.DESC")
        
        # Find the DESC file
        desc_file = find_desc_file(desc_dir, family_prefix)
        
        if desc_file is None:
            print(f"  WARNING: DESC file not found for {family_prefix}")
            not_found_count += 1
            continue
        
        print(f"  Found: {desc_file.name}")
        
        # Move the file
        success = move_desc_file(desc_file, family_id, base_dir, args.dry_run)
        if success:
            moved_count += 1
        else:
            skipped_count += 1
    
    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Total families in manifest: {len(families)}")
    print(f"DESC files moved: {moved_count}")
    print(f"DESC files not found in DESC dir: {not_found_count}")
    print(f"Families skipped (no target dir or DESC): {skipped_count}")
    
    if args.dry_run:
        print("\n*** This was a dry run - no files were actually moved ***")
    
    return 0


if __name__ == '__main__':
    exit(main())
