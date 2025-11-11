#!/usr/bin/env python3

"""
A productivity script to help merge curation directories with existing Pfam families.
Identifies potential merges based on gain ratios and guides curators through
the resolution process interactively.

Author: Based on requirements for Pfam curation
Version: v1.0
"""

__version__ = "v1.0"

import os
import sys
import subprocess
import logging
import argparse
import shutil
import re
from datetime import datetime
from pathlib import Path
import time

# Set up logging
def setup_logging():
    """Set up logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(message)s',
        handlers=[
            logging.StreamHandler(sys.stderr),
            logging.FileHandler('triage_merge.log', mode='a')
        ]
    )
    return logging.getLogger(__name__)

def run_command(cmd, check=True, capture_output=False, cwd=None):
    """Run a shell command and handle errors"""
    try:
        if capture_output:
            result = subprocess.run(cmd, shell=True, check=check, capture_output=True, 
                                  text=True, cwd=cwd)
            return result.stdout
        else:
            result = subprocess.run(cmd, shell=True, check=check, cwd=cwd)
            return result.returncode == 0
    except subprocess.CalledProcessError as e:
        if check:
            raise
        return False

def parse_triage_line(line):
    """Parse a line from the triage file"""
    parts = line.strip().split('\t')
    if len(parts) < 15:
        return None
    
    try:
        result = {
            'dir': parts[0],
            'total': int(parts[1]),
            'overlaps': int(parts[2]),
            'non_overlaps': int(parts[3]),
            'ratio': float(parts[4]),
            'overlap_info': []
        }
        
        # Parse overlap information from fields 15 onwards
        for i in range(14, len(parts)):
            if parts[i] and parts[i] != '----' and parts[i] != 'NA':
                result['overlap_info'].append(parts[i])
        
        return result
    except (ValueError, IndexError):
        return None

def get_pfam_accession(family_str, curation_dir):
    """Extract Pfam accession from overlap file or family string"""
    # If it already looks like an accession, return it
    if re.match(r'^PF\d{5}$', family_str):
        return family_str
    
    # Otherwise, parse the overlap file to find the accession
    overlap_file = f"{curation_dir}/overlap"
    if os.path.exists(overlap_file):
        with open(overlap_file, 'r') as f:
            for line in f:
                # Look for pattern like "Adeno_IVa2 PF02456"
                match = re.search(f'{family_str}\\s+(PF\\d{{5}})', line)
                if match:
                    return match.group(1)
    
    # If we can't find it, log a warning
    return None

def find_potential_merges(triage_file):
    """Find curation directories suitable for merging"""
    merges = []
    
    with open(triage_file, 'r') as f:
        for line in f:
            data = parse_triage_line(line)
            if not data:
                continue
            
            # Skip if no overlaps
            if data['overlaps'] == 0:
                continue
            
            # Check for single overlap pattern
            overlap_info = ' '.join(data['overlap_info'])
            
            # Parse all family entries in the overlap info
            # Pattern: [CLAN:]FAMILY COUNT(SIZE)
            # Example: "CL0342:FlgO 13(727)" or "TolB_N 25(3565)"
            family_pattern = r'(?:CL\d+:)?([A-Za-z0-9_]+)\s+\d+\((\d+)\)'
            families_found = re.findall(family_pattern, overlap_info)
            
            # Only process if exactly ONE unique family is found
            unique_families = set(f[0] for f in families_found)
            
            if len(unique_families) == 1:
                # Get the single family and its size
                family_name = families_found[0][0]
                pfam_size = int(families_found[0][1])
                gained = data['non_overlaps']
                
                # Only consider if gain is at least 10% of family size
                if gained >= pfam_size / 10:
                    ratio = gained / pfam_size
                    merges.append({
                        'dir': data['dir'],
                        'family': family_name,
                        'gained': gained,
                        'pfam_size': pfam_size,
                        'ratio': ratio
                    })
            elif len(unique_families) > 1:
                # Log that we're skipping due to multiple families
                pass  # Could add logging here if desired
    
    # Sort by ratio (highest first)
    merges.sort(key=lambda x: x['ratio'], reverse=True)
    return merges

def check_overlap(curation_dir, logger):
    """Check overlaps for a curation directory"""
    logger.info(f"\nChecking overlaps for {curation_dir}")
    
    # Run overlap check - use check=False as it returns 1 when overlaps are found
    run_command(f"pqc-overlap-rdb -no_sigP {curation_dir}", check=False)
    
    # Check if overlap file exists and has content
    overlap_file = f"{curation_dir}/overlap"
    if not os.path.exists(overlap_file):
        return False
    
    if os.path.getsize(overlap_file) == 0:
        logger.info("No overlaps found - already resolved")
        return False
    
    # Display overlap content
    with open(overlap_file, 'r') as f:
        content = f.read()
        logger.info("Overlap content:")
        logger.info("-" * 60)
        logger.info(content)
        logger.info("-" * 60)
    
    return True

def edit_seed_alignment(seed_file, logger):
    """Interactive editing of SEED alignment"""
    while True:
        # Display the alignment
        run_command(f"belvu {seed_file}", check=False)
        
        print("\nSEED editing options:", file=sys.stderr)
        print("  y/yes     - Accept current alignment", file=sys.stderr)
        print("  w/wholeseq - Run wholeseq alignment", file=sys.stderr)
        print("  e/extend  - Extend alignment boundaries", file=sys.stderr)
        print("  p/pad     - Pad alignment ends", file=sys.stderr)
        print("  n/no      - Cancel editing", file=sys.stderr)
        
        reply = input("Choice: ").strip().lower()
        
        if reply in ["y", "yes"]:
            return True
        
        elif reply in ["n", "no"]:
            return False
        
        elif reply in ["w", "wholeseq"]:
            backup = f"{seed_file}.before_wholeseq"
            shutil.copy(seed_file, backup)
            output = run_command(f"wholeseq.pl -align {backup} -m", capture_output=True)
            with open(seed_file, 'w') as f:
                f.write(output)
            logger.info("Applied wholeseq alignment")
        
        elif reply in ["e", "extend"]:
            n = input("N-terminal extension (default 0): ").strip() or "0"
            c = input("C-terminal extension (default 0): ").strip() or "0"
            
            backup = f"{seed_file}.before_extend"
            shutil.copy(seed_file, backup)
            output = run_command(f"extend.pl -align {backup} -n {n} -c {c} -m", 
                               capture_output=True)
            with open(seed_file, 'w') as f:
                f.write(output)
            logger.info(f"Extended alignment N:{n} C:{c}")
        
        elif reply in ["p", "pad"]:
            backup = f"{seed_file}.before_pad"
            shutil.copy(seed_file, backup)
            output = run_command(f"pad_ends.pl -align {backup}", capture_output=True)
            with open(seed_file, 'w') as f:
                f.write(output)
            logger.info("Padded alignment ends")

def merge_families(curation_dir, pfam_family, logger):
    """Merge curation directory with Pfam family"""
    # Get the actual Pfam accession if we have an ID
    pfam_acc = get_pfam_accession(pfam_family, curation_dir)
    if not pfam_acc:
        logger.error(f"Could not determine Pfam accession for {pfam_family}")
        return False
    
    logger.info(f"\nMerging {curation_dir} with {pfam_family} ({pfam_acc})")
    
    # Check out Pfam family if not already present
    pfam_dir = f"{curation_dir}/{pfam_acc}"
    if not os.path.exists(pfam_dir):
        logger.info(f"Checking out {pfam_acc}")
        run_command(f"pfco {pfam_acc}", cwd=curation_dir)
    
    # Check if merge is already in progress
    if os.path.exists(f"{pfam_dir}/.merge_in_progress"):
        logger.info(f"*** Merge already in progress for {pfam_acc} ***")
        logger.info(f"Found existing merge in {pfam_dir}")
        
        # Check if PFAMOUT exists and is newer than SEED
        if check_pfbuild_completion(pfam_dir, logger):
            logger.info("pfbuild has completed")
            reply = input("Complete the merge now? [y/n]: ").strip().lower()
            if reply in ["y", "yes"]:
                return complete_merge(pfam_dir, curation_dir, logger)
            else:
                logger.info("Skipping completion for now")
                return True
        else:
            logger.info("pfbuild still running or not started")
            reply = input("Wait for pfbuild to complete or restart merge? [wait/restart]: ").strip().lower()
            if reply == "restart":
                # Remove the marker and continue with new merge
                os.unlink(f"{pfam_dir}/.merge_in_progress")
                os.unlink(f"{pfam_dir}/.merge_info")
                logger.info("Restarting merge process")
            else:
                logger.info("Waiting for pfbuild to complete")
                return True
    
    # Count sequences in each SEED file
    curation_seed_count = 0
    pfam_seed_count = 0
    
    if os.path.exists(f"{curation_dir}/SEED"):
        with open(f"{curation_dir}/SEED", 'r') as f:
            curation_seed_count = sum(1 for line in f if line.strip() and not line.startswith('#'))
    
    if os.path.exists(f"{pfam_dir}/SEED"):
        with open(f"{pfam_dir}/SEED", 'r') as f:
            pfam_seed_count = sum(1 for line in f if line.strip() and not line.startswith('#'))
    
    logger.info(f"Curation SEED: {curation_seed_count} sequences")
    logger.info(f"Pfam {pfam_acc} SEED: {pfam_seed_count} sequences")
    logger.info(f"Total sequences to merge: {curation_seed_count + pfam_seed_count}")
    
    # Create merged SEED
    logger.info("Creating merged SEED alignment")
    merged_seed = f"{pfam_dir}/SEED_merged"
    run_command(f"merge_alignment.pl -fasta {pfam_dir}/SEED -fasta2 {curation_dir}/SEED -m > {merged_seed}", 
                check=False)
    
    # Check if merge was successful
    if not os.path.exists(merged_seed) or os.path.getsize(merged_seed) == 0:
        logger.error("Failed to create merged SEED")
        return False
    
    # Count sequences in merged SEED
    merged_seed_count = 0
    with open(merged_seed, 'r') as f:
        merged_seed_count = sum(1 for line in f if line.strip() and not line.startswith('#'))
    logger.info(f"Merged SEED contains {merged_seed_count} sequences")
    
    # Allow user to edit the merged SEED
    if edit_seed_alignment(merged_seed, logger):
        # Replace original SEED with edited version
        shutil.move(f"{pfam_dir}/SEED", f"{pfam_dir}/SEED.backup")
        shutil.move(merged_seed, f"{pfam_dir}/SEED")
        
        # Run pfbuild with pfmake
        logger.info("Running pfbuild -withpfmake (this will run in background)")
        run_command(f"pfbuild -withpfmake", cwd=pfam_dir, check=False)
        
        # Mark as in progress
        Path(f"{pfam_dir}/.merge_in_progress").touch()
        with open(f"{pfam_dir}/.merge_info", 'w') as f:
            f.write(f"curation_dir:{curation_dir}\n")
            f.write(f"pfam_acc:{pfam_acc}\n")
            f.write(f"timestamp:{datetime.now()}\n")
        
        return True
    else:
        logger.info("Merge cancelled")
        if os.path.exists(merged_seed):
            os.unlink(merged_seed)
        return False

def check_pfbuild_completion(pfam_dir, logger):
    """Check if pfbuild has completed"""
    # Check for PFAMOUT file that's newer than SEED
    pfamout = f"{pfam_dir}/PFAMOUT"
    seed = f"{pfam_dir}/SEED"
    
    if os.path.exists(pfamout) and os.path.exists(seed):
        if os.path.getmtime(pfamout) > os.path.getmtime(seed):
            return True
    
    # Check if pfbuild job is still running
    result = run_command("squeue -u $USER | grep pfbuild", capture_output=True, check=False)
    if pfam_dir.split('/')[-1] in result:
        return False  # Still running
    
    return False

def complete_merge(pfam_dir, curation_dir, logger):
    """Complete the merge process after pfbuild finishes"""
    logger.info(f"\nCompleting merge for {pfam_dir}")
    
    # Check overlaps - use check=False as it returns 1 when overlaps are found
    run_command(f"pqc-overlap-rdb {pfam_dir}", check=False)
    
    overlap_file = f"{pfam_dir}/overlap"
    has_overlaps = os.path.exists(overlap_file) and os.path.getsize(overlap_file) > 0
    if has_overlaps:
        logger.info("Warning: Pfam family still has overlaps after merge")
        with open(overlap_file, 'r') as f:
            logger.info(f.read())
    
    # Run pqc-missing to check what sequences are missing
    # Must be run from parent directory with just the accession
    pfam_acc = os.path.basename(pfam_dir)
    parent_dir = os.path.dirname(pfam_dir)
    logger.info(f"\nRunning pqc-missing on {pfam_acc} from {parent_dir}")
    run_command(f"pqc-missing {pfam_acc}", cwd=parent_dir, check=False)
    
    # Display missing/found summary if files exist
    missing_count = 0
    found_count = 0
    
    if os.path.exists(f"{pfam_dir}/missing"):
        with open(f"{pfam_dir}/missing", 'r') as f:
            missing_count = sum(1 for line in f if line.strip())
    
    if os.path.exists(f"{pfam_dir}/found"):
        with open(f"{pfam_dir}/found", 'r') as f:
            found_count = sum(1 for line in f if line.strip())
    
    logger.info(f"Found: {found_count} sequences")
    logger.info(f"Missing: {missing_count} sequences")
    logger.info(f"Net gain: {found_count - missing_count} sequences")
    
    # Ask if curator wants to proceed with check-in
    if has_overlaps:
        reply = input("\nProceed with check-in despite overlaps? [y/n/o(verlap)/i(gnore)/s(kip)]: ").strip().lower()
    else:
        reply = input("\nProceed with check-in? [y/n/o(verlap)/i(gnore)/s(kip)]: ").strip().lower()
    
    pfci_success = False
    
    if reply in ["y", "yes"]:
        # Check in the Pfam family
        logger.info("Checking in merged Pfam family")
        # Run pfci from within the parent directory with just the accession
        pfci_success = run_command(f"pfci -m 'Merged in additional sequences from curation' {pfam_acc}", 
                                   cwd=parent_dir, check=False)
        
        if pfci_success:
            destination = "done"
        else:
            logger.error("pfci failed - check error messages above")
            destination = input("Move to: done / overlap / ignore / skip? ").strip().lower()
    elif reply in ["o", "overlap"]:
        destination = "overlap"
    elif reply in ["i", "ignore"]:
        destination = "ignore"
    elif reply in ["s", "skip", "n", "no"]:
        destination = "skip"
    else:
        logger.info(f"Unknown response '{reply}', treating as skip")
        destination = "skip"
    
    # Move curation directory based on destination
    if destination in ["done", "overlap", "ignore"]:
        # Create destination directory if needed
        dest_dir = destination.upper()
        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)
        
        # Move curation directory
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        dest = f"{dest_dir}/{os.path.basename(curation_dir)}_{timestamp}"
        
        if os.path.exists(curation_dir):
            shutil.move(curation_dir, dest)
            logger.info(f"Moved {curation_dir} to {dest}")
        else:
            logger.warning(f"Curation directory {curation_dir} not found")
    elif destination == "skip":
        logger.info("Not moving curation directory")
    
    # Clean up progress markers
    for marker in ['.merge_in_progress', '.merge_info']:
        marker_file = f"{pfam_dir}/{marker}"
        if os.path.exists(marker_file):
            os.unlink(marker_file)
    
    return pfci_success

def abandon_merge(pfam_dir, curation_dir, logger):
    """Abandon a merge and move curation directory to DONE"""
    logger.info(f"Abandoning merge for {pfam_dir}")
    
    # Clean up progress markers
    for marker in ['.merge_in_progress', '.merge_info']:
        marker_file = f"{pfam_dir}/{marker}"
        if os.path.exists(marker_file):
            os.unlink(marker_file)
            logger.info(f"Removed {marker}")
    
    # Move curation directory to DONE
    done_dir = "DONE"
    if not os.path.exists(done_dir):
        os.makedirs(done_dir)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    dest = f"{done_dir}/{os.path.basename(curation_dir)}_abandoned_{timestamp}"
    
    if os.path.exists(curation_dir):
        shutil.move(curation_dir, dest)
        logger.info(f"Moved {curation_dir} to {dest}")
    else:
        logger.warning(f"Curation directory {curation_dir} not found - may have been moved already")
    
    return True

def handle_threshold_change(target_dir, is_pfam, threshold, logger):
    """Handle raising gathering threshold"""
    logger.info(f"Raising threshold to {threshold} in {target_dir}")
    
    # Run pfmake with new threshold
    run_command(f"pfmake -t {threshold}", cwd=target_dir)
    
    if is_pfam:
        # Check in the Pfam family with new threshold
        logger.info("Checking in Pfam family with new threshold")
        # Extract the PF accession from the path (last part)
        pfam_acc = os.path.basename(target_dir)
        # Get the parent directory (curation directory)
        curation_dir = os.path.dirname(target_dir)
        # Run pfci from within the curation directory
        run_command(f"pfci -m 'Raised threshold to resolve overlaps' {pfam_acc}", cwd=curation_dir)
    
    return True

def process_curation_directory(merge_info, logger):
    """Process a single curation directory for merging"""
    curation_dir = merge_info['dir']
    pfam_family = merge_info['family']
    
    # Get the actual Pfam accession
    pfam_acc = get_pfam_accession(pfam_family, curation_dir)
    if not pfam_acc:
        logger.error(f"Could not determine Pfam accession for {pfam_family}")
        return True  # Skip this directory
    
    logger.info("=" * 70)
    logger.info(f"Processing: {curation_dir}")
    logger.info(f"Target Pfam: {pfam_family} ({pfam_acc})")
    logger.info(f"Potential gain: {merge_info['gained']} sequences")
    logger.info(f"Pfam size: {merge_info['pfam_size']} sequences")
    logger.info(f"Gain ratio: {merge_info['ratio']:.2%}")
    logger.info("=" * 70)
    
    # First check if overlaps still exist
    if not check_overlap(curation_dir, logger):
        logger.info("Overlaps already resolved, skipping")
        return True
    
    while True:
        print("\nOptions:", file=sys.stderr)
        print("1 - Merge families", file=sys.stderr)
        print("2 - Resolve overlap in Pfam family by raising threshold", file=sys.stderr)
        print("3 - Resolve overlap in curation directory by raising threshold", file=sys.stderr)
        print("4 - Resolve overlap by editing Pfam SEED", file=sys.stderr)
        print("5 - Resolve overlap by editing curation SEED", file=sys.stderr)
        print("6 - Skip this directory", file=sys.stderr)
        print("7 - Move to IGNORE directory", file=sys.stderr)
        
        choice = input("Choice [1-7]: ").strip()
        
        if choice == "1":
            if merge_families(curation_dir, pfam_family, logger):
                return True
            # If merge was cancelled, continue showing options
        
        elif choice == "2":
            pfam_dir = f"{curation_dir}/{pfam_acc}"
            if not os.path.exists(pfam_dir):
                logger.info(f"Checking out {pfam_acc}")
                run_command(f"pfco {pfam_acc}", cwd=curation_dir)
            
            threshold = input("New threshold value: ").strip()
            if threshold:
                handle_threshold_change(pfam_dir, True, threshold, logger)
                # Re-check overlaps
                check_overlap(curation_dir, logger)
        
        elif choice == "3":
            threshold = input("New threshold value: ").strip()
            if threshold:
                handle_threshold_change(curation_dir, False, threshold, logger)
                # Re-check overlaps
                check_overlap(curation_dir, logger)
        
        elif choice == "4":
            pfam_dir = f"{curation_dir}/{pfam_acc}"
            if not os.path.exists(pfam_dir):
                logger.info(f"Checking out {pfam_acc}")
                run_command(f"pfco {pfam_acc}", cwd=curation_dir)
            
            seed_file = f"{pfam_dir}/SEED"
            if edit_seed_alignment(seed_file, logger):
                logger.info("Running pfbuild -withpfmake")
                run_command(f"pfbuild -withpfmake", cwd=pfam_dir, check=False)
                Path(f"{pfam_dir}/.edit_in_progress").touch()
                return True
        
        elif choice == "5":
            seed_file = f"{curation_dir}/SEED"
            if edit_seed_alignment(seed_file, logger):
                logger.info("Running pfbuild -withpfmake")
                run_command(f"pfbuild -withpfmake", cwd=curation_dir, check=False)
                # Re-check overlaps after build completes
                # Note: We should wait or check later
                return True
        
        elif choice == "6":
            logger.info("Skipping directory")
            return True
        
        elif choice == "7":
            # Move to IGNORE directory
            ignore_dir = "IGNORE"
            if not os.path.exists(ignore_dir):
                os.makedirs(ignore_dir)
            
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            dest = f"{ignore_dir}/{os.path.basename(curation_dir)}_{timestamp}"
            shutil.move(curation_dir, dest)
            logger.info(f"Moved {curation_dir} to {dest}")
            return True
        
        else:
            logger.info("Invalid choice, please try again")

def check_pending_merges(logger):
    """Check for any pending merges from previous runs"""
    pending = []
    
    logger.info("Scanning for pending merges...")
    
    # Directories to skip when scanning
    skip_dirs = {'DONE', 'IGNORE', 'OVERLAP', '.git', '.svn'}
    
    # Look for .merge_in_progress markers - search more thoroughly
    for root, dirs, files in os.walk(".", followlinks=True):
        # Remove directories we should skip from the search
        dirs[:] = [d for d in dirs if d not in skip_dirs]
        
        if '.merge_in_progress' in files:
            pfam_dir = root
            info_file = f"{pfam_dir}/.merge_info"
            
            logger.info(f"Found pending merge marker in: {pfam_dir}")
            
            if os.path.exists(info_file):
                with open(info_file, 'r') as f:
                    info = {}
                    for line in f:
                        if ':' in line:
                            key, value = line.strip().split(':', 1)
                            info[key] = value
                
                if 'curation_dir' in info:
                    pending.append({
                        'pfam_dir': pfam_dir,
                        'curation_dir': info['curation_dir'],
                        'pfam_acc': info.get('pfam_acc', os.path.basename(pfam_dir))
                    })
                    logger.info(f"  Curation dir: {info['curation_dir']}")
                    logger.info(f"  Pfam acc: {info.get('pfam_acc', 'unknown')}")
    
    return pending

def main():
    parser = argparse.ArgumentParser(
        description=f"Productivity script for merging curation directories with Pfam families (v{__version__})",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script identifies curation directories with single Pfam overlaps that could
benefit from merging. It guides curators through various resolution strategies
including merging, threshold adjustments, and SEED editing.

The script processes directories in order of highest gain ratio and provides
interactive options for handling each case.
        """
    )
    
    parser.add_argument('-t', '--triage', type=str, default='triage',
                       help='Triage file to process (default: triage)')
    parser.add_argument('-c', '--continue-pending', action='store_true',
                       help='Continue processing pending merges from previous runs')
    parser.add_argument('-m', '--min-ratio', type=float, default=0.1,
                       help='Minimum gain ratio to consider (default: 0.1)')
    
    args = parser.parse_args()
    
    # Set up logging
    logger = setup_logging()
    
    logger.info(f"Starting triage_merge.py {__version__}")
    logger.info(f"Processing triage file: {args.triage}")
    
    # ALWAYS check for pending merges (not just with -c flag)
    pending = check_pending_merges(logger)
    if pending:
        logger.info(f"Found {len(pending)} pending merges from previous runs")
        for p in pending:
            pfam_dir = p['pfam_dir']
            curation_dir = p['curation_dir']
            pfam_acc = p.get('pfam_acc', os.path.basename(pfam_dir))
            
            logger.info(f"Checking pending merge: {curation_dir} -> {pfam_acc}")
            
            if check_pfbuild_completion(pfam_dir, logger):
                logger.info(f"pfbuild completed for {pfam_dir}")
                
                # Run pqc-missing BEFORE asking about completion
                logger.info(f"\nRunning pqc-missing on {pfam_acc}")
                # Save current directory
                original_dir = os.getcwd()
                # Change to the curation directory
                os.chdir(curation_dir)
                # Run pqc-missing with just the accession
                run_command(f"pqc-missing {pfam_acc}", check=False)
                # Change back to original directory
                os.chdir(original_dir)
                
                # Display missing/found summary
                missing_count = 0
                found_count = 0
                
                if os.path.exists(f"{pfam_dir}/missing"):
                    with open(f"{pfam_dir}/missing", 'r') as f:
                        missing_count = sum(1 for line in f if line.strip())
                
                if os.path.exists(f"{pfam_dir}/found"):
                    with open(f"{pfam_dir}/found", 'r') as f:
                        found_count = sum(1 for line in f if line.strip())
                
                logger.info(f"Found: {found_count} sequences")
                logger.info(f"Missing: {missing_count} sequences")
                logger.info(f"Net gain: {found_count - missing_count} sequences")
                
                reply = input(f"Complete merge for {pfam_acc}? [y/n/abandon]: ").strip().lower()
                if reply in ["y", "yes"]:
                    complete_merge(pfam_dir, curation_dir, logger)
                elif reply == "abandon":
                    abandon_merge(pfam_dir, curation_dir, logger)
                else:
                    logger.info(f"Skipping completion of {pfam_acc} for now")
            else:
                logger.info(f"pfbuild still running or not completed for {pfam_dir}")
                reply = input(f"Continue waiting for {pfam_acc}? [wait/abandon]: ").strip().lower()
                if reply == "abandon":
                    abandon_merge(pfam_dir, curation_dir, logger)
                else:
                    logger.info(f"Continuing to wait for {pfam_acc}")
    
    # Find potential merges
    if not os.path.exists(args.triage):
        logger.error(f"Triage file {args.triage} not found")
        sys.exit(1)
    
    merges = find_potential_merges(args.triage)
    
    # Filter by minimum ratio
    merges = [m for m in merges if m['ratio'] >= args.min_ratio]
    
    if not merges:
        logger.info("No suitable merge candidates found")
        sys.exit(0)
    
    logger.info(f"Found {len(merges)} potential merge candidates")
    logger.info("\nTop candidates:")
    for i, m in enumerate(merges[:10], 1):
        logger.info(f"{i}. {m['dir']} -> {m['family']} (ratio: {m['ratio']:.2%})")
    
    # Process each merge candidate
    for merge in merges:
        try:
            # Skip if this directory has a pending merge
            skip = False
            for p in pending:
                if merge['dir'] in p['curation_dir']:
                    logger.info(f"Skipping {merge['dir']} - has pending merge")
                    skip = True
                    break
            
            if not skip:
                process_curation_directory(merge, logger)
        except KeyboardInterrupt:
            logger.info("\nInterrupted by user")
            break
        except Exception as e:
            logger.error(f"Error processing {merge['dir']}: {e}")
            continue
    
    logger.info("Processing complete")

if __name__ == '__main__':
    main()
