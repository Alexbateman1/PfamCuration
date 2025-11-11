#!/usr/bin/env python3

"""
A script to help iterating Pfam families. Either provide a single
family with -family or a file with a list of pfam names with -list
option.

Note that as families are checked in a file called record is made
that notes the family and number of sequences gained. Families in
this file are skipped if the program is rerun.

Author: Python port of agb's iterate.pl
Version: v41
"""

__version__ = "v41"

import os
import sys
import subprocess
import logging
import argparse
import shutil
import re
import time
import hashlib
import threading
from datetime import datetime
from pathlib import Path

# Set up logging
def setup_logging(process_id):
    """Set up logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(message)s',
        handlers=[
            logging.StreamHandler(sys.stderr),
            logging.FileHandler(f'{process_id}.iteration.log', mode='w')
        ]
    )
    return logging.getLogger(__name__)

def run_command(cmd, check=True, capture_output=False):
    """Run a shell command and handle errors"""
    try:
        if capture_output:
            result = subprocess.run(cmd, shell=True, check=check, capture_output=True, text=True)
            return result.stdout
        else:
            result = subprocess.run(cmd, shell=True, check=check)
            return result.returncode == 0
    except subprocess.CalledProcessError as e:
        if check:
            raise
        return False

def get_user():
    """Get current username"""
    return os.environ.get('USER', 'unknown')

def num_seq(filename):
    """Get number of sequences in an alignment"""
    if not os.path.exists(filename):
        logging.warning(f"failed in num_seq {filename} not found")
        return 0
    
    with open(filename, 'r') as f:
        return sum(1 for line in f if line.strip())

def is_flush(alignment_file):
    """Returns fraction of sequences in alignment that are flush to the ends"""
    if not os.path.exists(alignment_file):
        return 0
    
    total = 0
    flush = 0
    
    with open(alignment_file, 'r') as f:
        for line in f:
            # Match a sequence NAME/S-E SEQUENCE STRING
            parts = line.strip().split()
            if len(parts) >= 2:
                seq_string = parts[1]
                total += 1
                
                # Check if sequence is flush (doesn't start or end with gaps)
                if not (seq_string.startswith('.') or seq_string.startswith('-') or 
                        seq_string.endswith('.') or seq_string.endswith('-')):
                    flush += 1
    
    return flush / total if total > 0 else 0

def get_alignment_hash(filename):
    """Calculate a hash of an alignment file for comparison"""
    if not os.path.exists(filename):
        return None
    
    with open(filename, 'rb') as f:
        return hashlib.md5(f.read()).hexdigest()

def save_alignment_pair(family, before_file, after_file, benchmark_dir="/nfs/production/agb/pfam/users/agb/ALI_BENCHMARK"):
    """Save before and after alignment files if they differ"""
    # Calculate hashes to check if files differ
    before_hash = get_alignment_hash(before_file)
    after_hash = get_alignment_hash(after_file)
    
    if before_hash and after_hash and before_hash != after_hash:
        # Create benchmark directory if it doesn't exist
        os.makedirs(benchmark_dir, exist_ok=True)
        
        # Generate timestamp for unique naming
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Save the alignment pair
        before_dest = os.path.join(benchmark_dir, f"{family}_SEED_before_{timestamp}")
        after_dest = os.path.join(benchmark_dir, f"{family}_SEED_after_{timestamp}")
        
        shutil.copy(before_file, before_dest)
        shutil.copy(after_file, after_dest)
        
        logging.info(f"Saved alignment pair for {family} to benchmark directory")
        logging.info(f"  Before: {before_dest}")
        logging.info(f"  After: {after_dest}")
        return True
    
    return False

def send_terminal_beep():
    """Send a beep to the terminal"""
    print('\a', end='', flush=True)  # ASCII bell character

def get_desc_info(family):
    """Get AC, ID, DE, and TP information from DESC file"""
    desc_file = f"{family}/DESC"
    desc_info = {}
    
    if os.path.exists(desc_file):
        with open(desc_file, 'r') as f:
            for line in f:
                if line.startswith('AC '):
                    desc_info['AC'] = line.strip()
                elif line.startswith('ID '):
                    desc_info['ID'] = line.strip()
                elif line.startswith('DE '):
                    desc_info['DE'] = line.strip()
                elif line.startswith('TP '):
                    desc_info['TP'] = line.strip()
    
    return desc_info

def get_clan_info(family):
    """Get clan information from DESC file if it exists"""
    desc_file = f"{family}/DESC"
    if os.path.exists(desc_file):
        with open(desc_file, 'r') as f:
            for line in f:
                if line.startswith('CL '):
                    return line.strip().split()[1]
    return None

def touch_file(filename):
    """Create an empty file"""
    Path(filename).touch()

def fix_seed(filename):
    """Check if alignment ends on a // and remove it"""
    if not os.path.exists(filename):
        logging.error(f"In fix_seed {filename} does not exist")
        return
    
    with open(filename, 'r') as f:
        content = f.read()
    
    if '//' in content:
        logging.info(f"{filename} has //, removing")
        shutil.copy(filename, f"{filename}.old")
        
        with open(filename, 'w') as f:
            for line in content.splitlines():
                if '//' not in line:
                    f.write(line + '\n')
        
        os.unlink(f"{filename}.old")

def check_dni_ht_flags(family, logger):
    """Check if family has DNI (Do Not Iterate) or HT (High Threshold) flags"""
    desc_file = f"{family}/DESC"
    if not os.path.exists(desc_file):
        return False
    
    with open(desc_file, 'r') as f:
        for line in f:
            if line.startswith('**   DNI'):
                logger.warning(f"{family} has DNI (Do Not Iterate) flag - skipping iteration")
                return True
            elif line.startswith('**   HT'):
                logger.warning(f"{family} has HT (High Threshold) flag - skipping iteration")
                return True
    
    return False

def remove_ed_lines(family, logger):
    """Remove ED lines from DESC file as they can cause issues"""
    desc_file = f"{family}/DESC"
    if not os.path.exists(desc_file):
        return
    
    # Read all lines
    with open(desc_file, 'r') as f:
        lines = f.readlines()
    
    # Filter out ED lines
    filtered_lines = []
    ed_count = 0
    for line in lines:
        if line.startswith('ED '):
            ed_count += 1
        else:
            filtered_lines.append(line)
    
    # Only rewrite if we removed something
    if ed_count > 0:
        logger.info(f"Removing {ed_count} ED lines from DESC file")
        with open(desc_file, 'w') as f:
            f.writelines(filtered_lines)

def remove_family(family, user, archive):
    """Remove a Pfam family from directory"""
    date = datetime.now().strftime("%a_%b_%d_%H:%M:%S_%Z_%Y")
    new_location = f"{archive}/{family}.{user}.{date}"
    
    shutil.move(family, new_location)
    with open("record", "a") as f:
        f.write(f"{family} 0 REMOVED\n")

def pfbuild_family(family, local=False):
    """Run pfbuild for a family"""
    os.chdir(family)
    
    if local:
        run_command("pfbuild -local")
    else:
        run_command("pfbuild")
    
    if os.path.exists(".made"):
        os.unlink(".made")
    
    touch_file(".pfbuild")
    os.chdir("..")

def seedify_alignment(family, logger):
    """
    Create a high-quality SEED alignment using the SEEDify.pl methodology.
    Uses the best scoring sequence as template and preserves SwissProt sequences.
    """
    logger.info("Creating high-quality SEED alignment using SEEDify methodology")
    
    # Step 1: Identify template sequence and coordinates from scores file
    template_acc = None
    start_coord = None
    end_coord = None
    
    if os.path.exists(f"{family}/scores"):
        with open(f"{family}/scores", 'r') as f:
            for line in f:
                # Match pattern like: 353.7 D4GTU5.1/1-165 1-165 6.7e-103
                match = re.search(r'\s+(\S+\.\d+)\/(\d+)-(\d+)\s', line)
                if match:
                    template_acc = match.group(1)
                    start_coord = match.group(2)
                    end_coord = match.group(3)
                    logger.info(f"Template - {template_acc} {start_coord} {end_coord}")
                    break
    
    if not template_acc:
        logger.warning("Could not find template sequence in scores file, using simple method")
        return False
    
    # Step 2: Backup original ALIGN and extract template
    shutil.copy(f"{family}/ALIGN", f"{family}/OLDALIGN")
    
    # Extract template sequence (use check=False as grep may return 1 if no matches)
    run_command(f"grep '{template_acc}' {family}/ALIGN > {family}/tmp.template", check=False)
    
    # Remove template from ALIGN (use check=False as grep -v returns 1 if all lines match)
    run_command(f"grep -v '{template_acc}' {family}/ALIGN > {family}/tmp.ALIGN", check=False)
    
    # Check if tmp.ALIGN is empty (all sequences were the template)
    if not os.path.exists(f"{family}/tmp.ALIGN") or os.path.getsize(f"{family}/tmp.ALIGN") == 0:
        logger.warning("ALIGN only contains the template sequence, cannot build SEED using SEEDify method")
        return False
    
    # Step 3: Create/check sp file for SwissProt sequences
    if not os.path.exists(f"{family}/sp"):
        logger.info("No sp file exists, creating it now")
        run_command(f"cd {family} && swissprot.pl")
    
    # Remove old tmp.sp file if it exists
    if os.path.exists(f"{family}/tmp.sp"):
        os.unlink(f"{family}/tmp.sp")
    
    # Step 4: Extract SwissProt sequences
    # Always create tmp.sp file first (even if empty)
    open(f"{family}/tmp.sp", 'w').close()
    
    if os.path.exists(f"{family}/sp") and os.path.getsize(f"{family}/sp") > 0:
        logger.info("SwissProt sequences exist in this family")
        with open(f"{family}/sp", 'r') as f:
            for line in f:
                match = re.match(r'^(\S+)\.', line)
                if match:
                    acc = match.group(1)
                    # Skip if this is the template sequence (already extracted)
                    if acc == template_acc.split('.')[0]:
                        logger.info(f"Skipping {acc} as it's the template sequence")
                        continue
                    # Use check=False because grep returns 1 if no matches found
                    result = run_command(f"grep '{acc}' {family}/tmp.ALIGN", check=False, capture_output=True)
                    if result:  # Only append if we found something
                        with open(f"{family}/tmp.sp", 'a') as sp_file:
                            sp_file.write(result)
                        # Remove from tmp.ALIGN
                        run_command(f"grep -v '{acc}' {family}/tmp.ALIGN > {family}/tmp.ALIGN2", check=False)
                        if os.path.exists(f"{family}/tmp.ALIGN2"):
                            run_command(f"mv {family}/tmp.ALIGN2 {family}/tmp.ALIGN")
    else:
        logger.info("No SwissProt sequences in this family")
    
    # Step 5: Check ALIGN size and potentially downsample
    align_size = num_seq(f"{family}/ALIGN")
    if align_size > 500:
        logger.info(f"Down sampling alignment of {align_size} to 500 sequences")
        run_command(f"shuf -n 500 {family}/tmp.ALIGN > {family}/tmp.rest")
    else:
        shutil.copy(f"{family}/tmp.ALIGN", f"{family}/tmp.rest")
    
    # Step 6: Reassemble new SEED with template, swissprot seqs and the rest
    run_command(f"cat {family}/tmp.template {family}/tmp.sp {family}/tmp.rest > {family}/tmp")
    
    # Step 7: Extend alignment
    run_command(f"extend.pl -align {family}/tmp -n 50 -c 50 -m > {family}/SEED2")
    if not os.path.exists(f"{family}/SEED2") or os.path.getsize(f"{family}/SEED2") == 0:
        logger.warning("Failed to build SEED2")
        return False
    
    # Step 8: For repeat families, we may need to handle multiple occurrences differently
    # Check if this is a repeat family
    desc_info = get_desc_info(family)
    is_repeat = desc_info.get('TP', '').endswith('Repeat')
    
    if is_repeat:
        logger.info("This is a repeat family, using alternative approach")
        # For repeat families, skip the split_align step and go directly to removing partials
        shutil.copy(f"{family}/SEED2", f"{family}/SEED3")
    else:
        # For non-repeat families, use split_align to trim to template coordinates
        # First check if we have multiple template sequences
        template_count = 0
        template_lines = []
        with open(f"{family}/SEED2", 'r') as f:
            for line in f:
                if line.startswith(template_acc.split('.')[0]):  # Match by accession without version
                    template_count += 1
                    template_lines.append(line.strip()[:100] + "...")
        
        if template_count > 1:
            logger.info(f"Found {template_count} sequences matching template {template_acc}:")
            for tl in template_lines:
                logger.info(f"  {tl}")
            
            # Try to find the best matching template based on coordinates
            logger.info(f"Looking for template with coordinates close to {start_coord}-{end_coord}")
            
            # For now, just use split_align and let it fail if needed
            run_command(f"split_align.pl -align {family}/SEED2 -acc {template_acc.split('.')[0]} -s {start_coord} -e {end_coord} > {family}/SEED3", check=False)
            
            if not os.path.exists(f"{family}/SEED3") or os.path.getsize(f"{family}/SEED3") == 0:
                logger.warning("split_align.pl failed, using whole extended alignment")
                shutil.copy(f"{family}/SEED2", f"{family}/SEED3")
        else:
            # Normal case - single template
            run_command(f"split_align.pl -align {family}/SEED2 -acc {template_acc} -s {start_coord} -e {end_coord} > {family}/SEED3")
    
    if not os.path.exists(f"{family}/SEED3") or os.path.getsize(f"{family}/SEED3") == 0:
        logger.warning("Failed to build SEED3")
        return False
    
    # Step 9: Remove partial sequences
    # Check if we already determined this is a repeat family
    if 'is_repeat' not in locals():
        desc_info = get_desc_info(family)
        is_repeat = desc_info.get('TP', '').endswith('Repeat')
    
    # For repeat families, skip this step as it often removes everything
    if is_repeat:
        logger.info("Skipping partial sequence removal for repeat family")
        shutil.copy(f"{family}/SEED3", f"{family}/SEED4")
    else:
        run_command(f"belvu -P -o mul {family}/SEED3 | grep -v '//' > {family}/SEED4", check=False)
        
        # Check if removing partials left us with no sequences
        if not os.path.exists(f"{family}/SEED4") or os.path.getsize(f"{family}/SEED4") == 0:
            logger.warning("Removing partial sequences left no sequences, using original")
            shutil.copy(f"{family}/SEED3", f"{family}/SEED4")
    
    if not os.path.exists(f"{family}/SEED4") or os.path.getsize(f"{family}/SEED4") == 0:
        logger.warning("Failed to build SEED4")
        return False
    
    # Step 10: Make non-redundant at 80%
    run_command(f"belvu -n 80 {family}/SEED4 -o mul | grep -v '//' > {family}/NEWSEED")
    
    # Step 11: Replace old SEED with new one
    if os.path.exists(f"{family}/NEWSEED") and os.path.getsize(f"{family}/NEWSEED") > 0:
        shutil.move(f"{family}/SEED", f"{family}/OLDSEED")
        shutil.move(f"{family}/NEWSEED", f"{family}/SEED")
        logger.info("Successfully created high-quality SEED alignment")
        
        # Clean up temporary files
        for tmpfile in ['tmp.template', 'tmp.ALIGN', 'tmp.sp', 'tmp.rest', 'tmp', 'SEED2', 'SEED3', 'SEED4']:
            if os.path.exists(f"{family}/{tmpfile}"):
                os.unlink(f"{family}/{tmpfile}")
        
        return True
    else:
        logger.warning("Failed to create NEWSEED")
        return False

def iterate_family(family, seed_max, logger, ignore):
    """
    Main workhorse subroutine. Does the following:
    1) Make a new SEED from the ALIGN using belvu to remove partial sequences
       and make alignment non-redundant at decreasing sequence identities until
       you have less than seed_max
    2) Set a pfbuild running
    """
    if os.path.isdir(family):
        logger.info(f"{family} exists in current dir, skipping seed building step")
        return
    
    logger.info(f"Working on {family}")
    
    # Try to check family out
    if not run_command(f"pfco {family}", check=False):
        logger.warning(f"failed to pfco {family}")
        with open("record", "a") as f:
            f.write(f"{family} 0 NO_CHECK_OUT\n")
        ignore.add(family)
        return
    
    # Check for DNI or HT flags before proceeding
    if check_dni_ht_flags(family, logger):
        ignore.add(family)
        return
    
    # Remove ED lines from DESC file as they can cause issues
    remove_ed_lines(family, logger)
    
    logger.info(f"Making family {family} - SEED building step")
    
    # Display DESC file information
    desc_info = get_desc_info(family)
    if desc_info:
        logger.info("=" * 60)
        for key in ['AC', 'ID', 'DE', 'TP']:
            if key in desc_info:
                logger.info(desc_info[key])
        logger.info("=" * 60)
    
    # Try the SEEDify methodology first
    if seedify_alignment(family, logger):
        # Successfully created SEED using SEEDify method
        old_seed_size = num_seq(f"{family}/OLDSEED")
        seed_size = num_seq(f"{family}/SEED")
        logger.info(f"{family}: oldseed:{old_seed_size} newseed:{seed_size} (using SEEDify method)")
        return
    
    # Fall back to original method if SEEDify fails
    logger.info("Falling back to original SEED building method")
    
    # Copy ALIGN as backup
    shutil.copy(f"{family}/ALIGN", f"{family}/OLDALIGN")
    
    # Start with 80% non-redundancy threshold
    nr_thresh = 80
    
    # Check if ALIGN file is empty or has issues
    if not os.path.exists(f"{family}/ALIGN") or os.path.getsize(f"{family}/ALIGN") == 0:
        logger.warning(f"ALIGN file is empty or missing for {family}")
        return
    
    # Try to run belvu, catching errors
    if not run_command(f"belvu -o mul -n {nr_thresh} {family}/ALIGN | grep -v '//' > {family}/SEED.2", check=False):
        logger.warning(f"belvu failed on {family}/ALIGN - possibly corrupted alignment")
        # Try to copy ALIGN directly as SEED.2
        shutil.copy(f"{family}/ALIGN", f"{family}/SEED.2")
    
    seed_size = num_seq(f"{family}/SEED.2")
    
    if not seed_size:
        logger.warning("Skipping family as it seems to have no sequences in SEED")
        return
    
    # Progressively lower threshold until we have less than seed_max sequences
    while seed_size > seed_max:
        logger.info(f"New seed has {seed_size} sequences using {nr_thresh}%.")
        nr_thresh -= 10
        if not run_command(f"belvu -o mul -n {nr_thresh} {family}/ALIGN | grep -v '//' > {family}/SEED.2", check=False):
            logger.warning(f"belvu failed at {nr_thresh}% threshold")
            break
        seed_size = num_seq(f"{family}/SEED.2")
    
    # Handle case where we end up with only 1 sequence
    if seed_size == 1:
        logger.info("Only 1 sequence in SEED, trying higher identity thresholds")
        nr_thresh = 85  # Start at 85%
        
        while nr_thresh <= 100 and seed_size == 1:
            if not run_command(f"belvu -o mul -n {nr_thresh} {family}/ALIGN | grep -v '//' > {family}/SEED.2", check=False):
                logger.warning(f"belvu failed at {nr_thresh}% threshold")
                nr_thresh += 5
                continue
                
            seed_size = num_seq(f"{family}/SEED.2")
            
            if seed_size > 1:
                logger.info(f"Found {seed_size} sequences at {nr_thresh}% identity")
                break
            
            nr_thresh += 5
        
        if seed_size == 1:
            logger.warning("Still only 1 sequence after trying all thresholds, removing family")
            remove_family(family, get_user(), "/homes/agb/Curation/ITERATION_ARCHIVE")
            ignore.add(family)
            return
    
    logger.info(f"New seed has {seed_size} sequences using {nr_thresh}%. Great!")
    
    # Realign with MAFFT
    if not os.path.exists(f"{family}/OLDSEED"):
        shutil.move(f"{family}/SEED", f"{family}/OLDSEED")
    if not os.path.exists(f"{family}/OLDALIGN"):
        shutil.move(f"{family}/ALIGN", f"{family}/OLDALIGN")
    logger.info("Realign with MAFFT into SEED")
    run_command(f"create_alignment.pl -fasta {family}/SEED.2 -m > {family}/SEED")
    
    old_seed_size = num_seq(f"{family}/OLDSEED")
    logger.info(f"{family}: oldseed:{old_seed_size} newseed:{seed_size} at {nr_thresh} percent")

def check_seed(family, logger, ignore, auto=False):
    """
    Check SEED quality and allow user to modify if needed
    Returns True if family has been dealt with, False if should be run again
    """
    # Check if family directory still exists
    if not os.path.isdir(family):
        logger.info(f"{family} directory no longer exists, skipping")
        return True
        
    if os.path.exists(f"{family}/.goodseed"):
        logger.info(f"{family} SEED has been checked, skipping seed checking step")
        return True
    
    # Display DESC file information
    desc_info = get_desc_info(family)
    if desc_info:
        logger.info("=" * 60)
        for key in ['AC', 'ID', 'DE', 'TP']:
            if key in desc_info:
                logger.info(desc_info[key])
        logger.info("=" * 60)
    
    # Print useful stats
    seed_size = num_seq(f"{family}/SEED")
    old_seed_size = num_seq(f"{family}/OLDSEED")
    logger.info(f"{family}: oldseed:{old_seed_size} newseed:{seed_size}")
    fraction_flush = is_flush(f"{family}/SEED")
    logger.info(f"Fraction flush:{fraction_flush}")
    
    if auto:
        logger.info(f"{family} Auto mode - skipping manual SEED check")
        pfbuild_family(family)
        touch_file(f"{family}/.goodseed")
        return True
    
    # Save SEED before showing to curator
    seed_before = f"{family}/SEED.before_check"
    shutil.copy(f"{family}/SEED", seed_before)
    
    logger.info(f"{family} Please check SEED quality")
    run_command(f"belvu -t {family} {family}/SEED", check=False)
    
    print("Are you happy to pfbuild? [y/n/wholeseq/extend/pad/remove]", file=sys.stderr)
    reply = input().strip().lower()  # Convert to lowercase for comparison
    
    # Check if SEED was modified and save if different
    save_alignment_pair(family, seed_before, f"{family}/SEED")
    
    # Clean up temporary file
    if os.path.exists(seed_before):
        os.unlink(seed_before)
    
    if reply in ["y", "yes"]:
        logger.info(f"SEED is good. Start to pfbuild {family}")
        pfbuild_family(family)
        touch_file(f"{family}/.goodseed")
        return True
    
    elif reply in ["n", "no"]:
        logger.info(f"SEED is not good. Skip {family} to next family")
        return True
    
    elif reply in ["wholeseq", "w"]:
        shutil.copy(f"{family}/SEED", f"{family}/SEED.beforewholeseq")
        run_command(f"wholeseq.pl -align {family}/SEED.beforewholeseq -m > {family}/SEED")
        return False
    
    elif reply in ["extend", "e"]:
        n = input("How many residues should I extend N-terminus [default 0]\n").strip() or "0"
        c = input("How many residues should I extend C-terminus [default 0]\n").strip() or "0"
        
        logger.info(f"extending SEED N:{n} C:{c}")
        shutil.copy(f"{family}/SEED", f"{family}/SEED.beforeextend")
        run_command(f"extend.pl -align {family}/SEED.beforeextend -n {n} -c {c} -m > {family}/SEED")
        return False
    
    elif reply in ["pad", "p"]:
        logger.info("Padding SEED alignment")
        shutil.copy(f"{family}/SEED", f"{family}/SEED.beforepad")
        # Capture output from pad_ends.pl and write to SEED
        output = run_command(f"pad_ends.pl -align {family}/SEED.beforepad", capture_output=True)
        with open(f"{family}/SEED", 'w') as f:
            f.write(output)
        return False
    
    elif reply in ["remove", "r"]:
        remove_family(family, get_user(), "/homes/agb/Curation/ITERATION_ARCHIVE")
        ignore.add(family)  # Add to ignore set immediately
        return True
    
    else:
        logger.info("Did not recognise response. Try again ...")
        return False

def make_family(family, conservative, rebuild, logger, ignore):
    """
    Do a pfmake and QC including all overlap checks
    """
    # Check if family directory still exists
    if not os.path.isdir(family):
        logger.info(f"{family} directory no longer exists, skipping")
        return True
        
    if os.path.exists(f"{family}/.made"):
        logger.info(f"{family} has already been made, skipping make family step")
        return True
    
    logger.info(f"Making family {family} - pfmake and QC steps")
    
    # Check if pfbuild has run
    if (os.path.exists(f"{family}/PFAMOUT") and 
        os.path.getmtime(f"{family}/SEED") < os.path.getmtime(f"{family}/PFAMOUT")):
        
        if os.path.exists(f"{family}/.pfbuild"):
            os.unlink(f"{family}/.pfbuild")
        
        os.chdir(family)
        
        if not conservative:
            logger.info("pfbuild has run: Run pfmake -e 0.01")
            run_command("pfmake -e 0.01")
        else:
            logger.info("pfbuild has run: Run pfmake")
            run_command("pfmake")
        
        os.chdir("..")
        
        # Check if family is in a clan and run clan overlap check
        clan = get_clan_info(family)
        if clan:
            logger.info(f"{family} is a member of clan {clan} - checking clan overlaps")
            # First run with -ignore_cl to get all overlaps including clan members
            run_command(f"pqc-overlap-rdb -filter -ignore_cl {family}", check=False)
            
            # Move the clan overlap file
            if os.path.exists(f"{family}/overlap"):
                shutil.move(f"{family}/overlap", f"{family}/overlap_clan")
            
            # Re-run normal overlap check with -compete flag for clan families
            run_command(f"pqc-overlap-rdb -filter -compete {family}", check=False)
        else:
            # For non-clan families, just run regular overlap check
            run_command(f"pqc-overlap-rdb -filter {family}", check=False)
        
        # Run other QC checks
        for script in ["pqc-missing"]:
            run_command(f"{script} {family}", check=False)
        
        fraction_flush = is_flush(f"{family}/SEED")
        logger.info(f"Fraction flush:{fraction_flush}")
        touch_file(f"{family}/.made")
        return True
    
    elif rebuild:
        logger.info(f"{family} Appears that pfbuild failed to run. Try and pfbuild again")
        pfbuild_family(family)
    else:
        logger.info(f"pfbuild has not finished for {family}")
        return False

def check_in(family, logger, ignore, auto=False):
    """Check in family after validation - now simplified since overlaps are already checked"""
    # Check if family directory still exists
    if not os.path.isdir(family):
        logger.info(f"{family} directory no longer exists, skipping")
        return False
        
    if os.path.exists(f"{family}/.checkin"):
        logger.info(f"{family} is in process of being checked in! skipping")
        return False
    
    if not os.path.exists(f"{family}/.made"):
        logger.info(f"{family} has not had pfmake run yet! skipping")
        return False
    
    logger.info(f"{family} - start checks for check in")
    
    # Get stats
    found = sum(1 for _ in open(f"{family}/found"))
    missing = sum(1 for _ in open(f"{family}/missing"))
    overall_gain = found - missing
    
    oldseedsize = num_seq(f"{family}/OLDSEED")
    newseedsize = num_seq(f"{family}/SEED")
    seed_gain = newseedsize - oldseedsize
    seed_gain_percent = (seed_gain / oldseedsize * 100) if oldseedsize > 0 else 0
    
    oldalignsize = num_seq(f"{family}/OLDALIGN")
    newalignsize = num_seq(f"{family}/ALIGN")
    align_gain = newalignsize - oldalignsize
    
    logger.info(f"Overall sequences gain is {overall_gain} - found:{found} missing:{missing}")
    logger.info(f"Overall seed size gain is {seed_gain} - oldseed:{oldseedsize} newseed:{newseedsize} ({seed_gain_percent:.1f}%)")
    logger.info(f"Overall align size gain is {align_gain} - oldalign:{oldalignsize} newalign:{newalignsize}")
    
    # Check if family is in a clan and report clan overlaps if they exist
    clan = get_clan_info(family)
    if clan:
        logger.info(f"{family} is a member of clan {clan}")
        
        # Always report clan overlap count, even if zero
        clan_overlap_count = 0
        if os.path.exists(f"{family}/overlap_clan"):
            with open(f"{family}/overlap_clan", 'r') as f:
                clan_overlap_count = sum(1 for line in f if line.strip())
        
        if clan_overlap_count > 0:
            logger.info(f"Found {clan_overlap_count} overlaps with clan members (may explain some gains)")
        else:
            logger.info("Found 0 overlaps with clan members")
    
    # Check regular overlaps
    overlap_count = 0
    if os.path.exists(f"{family}/overlap") and os.path.getsize(f"{family}/overlap") > 0:
        with open(f"{family}/overlap", 'r') as f:
            overlap_count = sum(1 for line in f if line.strip())
        if overlap_count > 0:
            logger.info(f"{family} has {overlap_count} overlaps (excluding clan members)")
    else:
        logger.info(f"{family} has 0 overlaps")
    
    # Auto mode decision logic
    if auto:
        logger.info("Auto mode - making automatic decision")
        
        # Auto remove cases
        if align_gain < 0:
            logger.info("Auto-removing family: ALIGN has lost members")
            remove_family(family, get_user(), "/homes/agb/Curation/ITERATION_ARCHIVE")
            ignore.add(family)
            return False
        
        if overlap_count > 0:
            logger.info("Auto-removing family: has overlaps")
            remove_family(family, get_user(), "/homes/agb/Curation/ITERATION_ARCHIVE")
            ignore.add(family)
            return False
        
        # Auto check-in cases
        if overall_gain == 0 and seed_gain_percent >= 10:
            logger.info(f"Auto check-in: No ALIGN change but SEED gained {seed_gain_percent:.1f}%")
            auto_reply = "y"
        elif align_gain > 0 and seed_gain > 0:
            logger.info("Auto check-in: Both ALIGN and SEED have gained members")
            auto_reply = "y"
        else:
            logger.info("Auto mode: No clear decision, removing family")
            remove_family(family, get_user(), "/homes/agb/Curation/ITERATION_ARCHIVE")
            ignore.add(family)
            return False
    
    if is_flush(f"{family}/SEED") != 1:
        logger.info(f"{family}/SEED is not flush")
        if auto:
            logger.info("Auto mode: SEED not flush, skipping family")
            return False
        
        logger.info(f"Please fix {family}/SEED")
        run_command(f"belvu -t {family} {family}/SEED", check=False)
        fix_seed(f"{family}/SEED")
        
        reply = input("Are you happy to pfbuild this family? [y/n/remove]").strip().lower()
        if reply in ["y", "yes"]:
            if not os.path.exists(f"{family}/.pfbuild"):
                pfbuild_family(family)
        elif reply in ["remove", "r"]:
            remove_family(family, get_user(), "/homes/agb/Curation/ITERATION_ARCHIVE")
            ignore.add(family)
            return False
        return False
    
    if overall_gain < 0 and not auto:
        reply = input("Do you wish to remove this family from further analysis? [y/n]\n").strip().lower()
        if reply in ["y", "yes"]:
            remove_family(family, get_user(), "/homes/agb/Curation/ITERATION_ARCHIVE")
            ignore.add(family)
        return False
    
    # Display overlap details if present (only in manual mode)
    if overlap_count > 0 and not auto:
        logger.info("\n=== OVERLAP DETAILS ===")
        with open(f"{family}/overlap", 'r') as f:
            overlap_content = f.read()
            print(overlap_content)
        logger.info("=== END OVERLAP DETAILS ===\n")
    
    # Main check-in loop
    while True:
        if auto:
            reply = auto_reply
        else:
            if overlap_count > 0:
                reply = input("Are you happy to check this family in? [y/n/pfbuild/remove/threshold/ed]\n").strip().lower()
            else:
                reply = input("Are you happy to check this family in? [y/n/pfbuild/remove]\n").strip().lower()
        
        if reply in ["y", "yes"]:
            # Only try to check in if there are no overlaps, or warn the user
            if overlap_count > 0 and not auto:
                logger.warning(f"Cannot check in family with {overlap_count} overlaps")
                logger.info("Please resolve overlaps first using 'threshold' or 'ed' options")
                continue
            
            touch_file(f"{family}/.checkin")
            
            # Try to check in with timeout and beep functionality
            logger.info(f"Attempting to check in {family}...")
            
            # Start a background timer for the beep
            beep_timer = threading.Timer(30.0, send_terminal_beep)
            beep_timer.start()
            
            try:
                # Try to check in, but handle failure
                if run_command(f"pfci -m 'Iterated family and gained {overall_gain} sequences' {family}", check=False):
                    beep_timer.cancel()  # Cancel beep if successful
                    
                    with open("record", "a") as f:
                        f.write(f"{family} {overall_gain} ({align_gain})\n")
                    
                    ignore.add(family)
                    shutil.rmtree(family)
                    break
                else:
                    beep_timer.cancel()  # Cancel beep
                    logger.error("pfci failed - likely due to overlaps or other issues")
                    if auto:
                        logger.info("Auto mode: pfci failed, skipping family")
                        return False
                    logger.info("Please fix the issues and try again")
                    os.unlink(f"{family}/.checkin")  # Remove the checkin flag
                    continue
            except Exception as e:
                beep_timer.cancel()
                raise e
        
        elif reply in ["n", "no"]:
            logger.info(f"Skipping {family} - not checking in")
            break
            
        elif reply in ["pfbuild", "p"]:
            fix_seed(f"{family}/SEED")
            if not os.path.exists(f"{family}/.pfbuild"):
                pfbuild_family(family)
            if os.path.exists(f"{family}/.made"):
                os.unlink(f"{family}/.made")
            break
        
        elif reply in ["remove", "r"]:
            remove_family(family, get_user(), "/homes/agb/Curation/ITERATION_ARCHIVE")
            ignore.add(family)
            break
        
        elif reply in ["threshold", "t"] and overlap_count > 0:
            new_threshold = input("Enter new threshold (e.g., -T 25 -t 20): ").strip()
            if new_threshold:
                logger.info(f"Rethresholding with: pfmake {new_threshold}")
                os.chdir(family)
                run_command(f"pfmake {new_threshold}")
                os.chdir("..")
                
                # Re-run overlap checks after rethresholding
                if clan:
                    # First check with -ignore_cl
                    run_command(f"pqc-overlap-rdb -filter -ignore_cl {family}", check=False)
                    if os.path.exists(f"{family}/overlap"):
                        shutil.move(f"{family}/overlap", f"{family}/overlap_clan")
                    # Then check with -compete
                    run_command(f"pqc-overlap-rdb -filter -compete {family}", check=False)
                else:
                    run_command(f"pqc-overlap-rdb -filter {family}", check=False)
                
                # Update overlap count
                overlap_count = 0
                if os.path.exists(f"{family}/overlap") and os.path.getsize(f"{family}/overlap") > 0:
                    with open(f"{family}/overlap", 'r') as f:
                        overlap_count = sum(1 for line in f if line.strip())
                logger.info(f"After rethresholding: {overlap_count} overlaps")
                
                # Display new overlaps if any
                if overlap_count > 0:
                    logger.info("\n=== NEW OVERLAP DETAILS ===")
                    with open(f"{family}/overlap", 'r') as f:
                        overlap_content = f.read()
                        print(overlap_content)
                    logger.info("=== END OVERLAP DETAILS ===\n")
                
                # Continue loop to ask again
                continue
        
        elif reply in ["ed", "e"] and overlap_count > 0:
            logger.info("Adding ED lines...")
            os.chdir(family)
            run_command("add_ED.pl")
            run_command("pfmake")
            os.chdir("..")
            
            # Re-run overlap checks after adding ED lines
            if clan:
                # First check with -ignore_cl
                run_command(f"pqc-overlap-rdb -filter -ignore_cl {family}", check=False)
                if os.path.exists(f"{family}/overlap"):
                    shutil.move(f"{family}/overlap", f"{family}/overlap_clan")
                # Then check with -compete
                run_command(f"pqc-overlap-rdb -filter -compete {family}", check=False)
            else:
                run_command(f"pqc-overlap-rdb -filter {family}", check=False)
            
            # Update overlap count
            overlap_count = 0
            if os.path.exists(f"{family}/overlap") and os.path.getsize(f"{family}/overlap") > 0:
                with open(f"{family}/overlap", 'r') as f:
                    overlap_count = sum(1 for line in f if line.strip())
            logger.info(f"After adding ED lines: {overlap_count} overlaps")
            
            # Display new overlaps if any
            if overlap_count > 0:
                logger.info("\n=== NEW OVERLAP DETAILS ===")
                with open(f"{family}/overlap", 'r') as f:
                    overlap_content = f.read()
                    print(overlap_content)
                logger.info("=== END OVERLAP DETAILS ===\n")
            
            # Continue loop to ask again
            continue
        
        else:
            if not auto:  # Don't loop in auto mode
                logger.info("Did not recognise response. Try again...")
                continue
        
        # Only break if we're in auto mode or handled a valid response
        if auto:
            break

def main():
    parser = argparse.ArgumentParser(
        description=f"Script to help iterate Pfam families (v{__version__})",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script works on a list of Pfam families and does three steps on each:
1) Check out family and build a new seed from the align file with belvu.
2) Check and correct SEED
3) Run pfbuild on family
4) Runs pfmake and runs QC steps on the family asks if OK to check in

A file called record is made in the current directory which lists
the number of sequences gained per family.
        """
    )
    
    parser.add_argument('-version', '--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument('-list', type=str, help='Run program on a list of families')
    parser.add_argument('-family', type=str, help='Run program on a single family')
    parser.add_argument('-seed_max', type=int, default=200, 
                        help='Max number of sequences to allow in seed (default: 200)')
    parser.add_argument('-delete', action='store_true', 
                        help='Remove family from further analysis')
    parser.add_argument('-conservative', action='store_true',
                        help='Uses GA line thresholds rather than -e 0.01 default for pfmakes')
    parser.add_argument('-rebuild', action='store_true',
                        help='Will rerun pfbuilds where needed')
    parser.add_argument('-local', action='store_true',
                        help='Will run pfbuild locally - good when queues are busy')
    parser.add_argument('-auto', action='store_true',
                        help='Automatic mode - makes decisions based on statistics')
    
    args = parser.parse_args()
    
    # Set up logging
    logger = setup_logging(os.getpid())
    
    # Print version at startup
    logger.info(f"iterate.py version {__version__}")
    
    if args.auto:
        logger.info("Running in AUTO mode - will make automatic decisions")
    
    # Get user and archive location
    user = get_user()
    archive = "/homes/agb/Curation/ITERATION_ARCHIVE"
    
    # Determine families to process
    families = []
    if args.family:
        families.append(args.family)
    
    if args.list:
        with open(args.list, 'r') as f:
            for line in f:
                if line.startswith('PF'):
                    family = line.strip().split()[0]
                    families.append(family)
                else:
                    logger.info(f"{line.strip()} appears not to be a Pfam accession. Skipping.")
    
    if not families:
        parser.print_help()
        sys.exit(0)
    
    # Get list of families to ignore
    ignore = set()
    if os.path.exists("record"):
        with open("record", 'r') as f:
            for line in f:
                parts = line.strip().split()
                if parts:
                    ignore.add(parts[0])
    
    # Handle delete option
    if args.delete and args.family:
        logger.info(f"Deleting {args.family}")
        remove_family(args.family, user, archive)
        sys.exit(0)
    
    # Process families
    # Step 1: Create new SEED files
    for family in families:
        if family not in ignore:
            iterate_family(family, args.seed_max, logger, ignore)
        else:
            logger.info(f"{family} is in ignore list")
    
    # Step 2: Check SEED files
    for family in families:
        if family not in ignore:
            while not check_seed(family, logger, ignore, auto=args.auto):
                pass
        else:
            logger.info(f"{family} is in ignore list")
    
    # Step 3: Make families
    for family in families:
        if family not in ignore:
            make_family(family, args.conservative, args.rebuild, logger, ignore)
        else:
            logger.info(f"{family} is in ignore list")
    
    # Step 4: Check in families
    for family in families:
        if family not in ignore:
            check_in(family, logger, ignore, auto=args.auto)
        else:
            logger.info(f"{family} is in ignore list")

if __name__ == '__main__':
    main()
