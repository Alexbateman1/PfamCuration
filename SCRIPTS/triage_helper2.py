#!/usr/bin/env python3
"""
Pfam Triage Curation Helper Script - Version 2
Extends triage_helper.py to allow resolution of overlaps in deeper Iterate directories.

Key features:
- Focuses on non-overlapping directories (like triage_helper.py)
- Also considers deeper overlapping directories if overlaps < 5% of total
- Offers overlap resolution options including clan joining
- Only allows resolution if single overlapping family OR all overlaps in same clan
"""

import sys
import os
import subprocess
import shutil
import re
from pathlib import Path
from collections import defaultdict
from datetime import datetime


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
                subprocess.run(cmd, shell=True, cwd=cwd)
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


def count_seed_sequences(dir_path):
    """Count sequences in SEED file for a directory"""
    seed_path = Path(dir_path) / 'SEED'
    return count_sequences(str(seed_path))


def parse_overlap_info(overlap_details):
    """
    Parse overlap details from triage file to extract family and clan information.

    Input format examples:
        "CL0351:CHCH 1050(8242)" - family CHCH in clan CL0351, 1050 overlaps, 8242 total size
        "Flavoprotein 2(21413)" - family Flavoprotein not in a clan, 2 overlaps, 21413 total size

    Returns list of dicts with keys: family, clan (or None), overlap_count, pfam_size
    """
    overlaps = []

    for detail in overlap_details:
        if not detail or detail == '----' or detail == 'NA':
            continue

        # Pattern: [CLAN:]FAMILY COUNT(SIZE)
        # Examples: "CL0351:CHCH 1050(8242)" or "Flavoprotein 2(21413)"
        match = re.match(r'^(?:(CL\d+):)?(\S+)\s+(\d+)\((\d+)\)$', detail)
        if match:
            clan = match.group(1)  # May be None
            family = match.group(2)
            overlap_count = int(match.group(3))
            pfam_size = int(match.group(4))

            overlaps.append({
                'family': family,
                'clan': clan,
                'overlap_count': overlap_count,
                'pfam_size': pfam_size,
                'raw': detail
            })

    return overlaps


def check_resolution_eligibility(overlap_info_list):
    """
    Check if overlap resolution should be offered.

    Rules:
    - Only offer if there's exactly ONE overlapping family, OR
    - All overlapping families are in the SAME clan

    Returns: (eligible, reason, clan_to_join)
        - eligible: True/False
        - reason: explanation string
        - clan_to_join: clan accession if all overlaps are in same clan, else None
    """
    if not overlap_info_list:
        return False, "No overlap information", None

    # Get unique families
    families = set(oi['family'] for oi in overlap_info_list)

    if len(families) == 1:
        # Single family - eligible
        oi = overlap_info_list[0]
        if oi['clan']:
            return True, f"Single overlapping family ({oi['family']}) in clan {oi['clan']}", oi['clan']
        else:
            return True, f"Single overlapping family ({oi['family']}) not in a clan", None

    # Multiple families - check if all in same clan
    clans = set(oi['clan'] for oi in overlap_info_list if oi['clan'])
    families_without_clan = [oi['family'] for oi in overlap_info_list if not oi['clan']]

    if families_without_clan:
        # Some families not in any clan
        return False, f"Multiple overlapping families, some not in clans: {', '.join(families_without_clan)}", None

    if len(clans) == 1:
        # All families in the same clan
        clan = list(clans)[0]
        return True, f"All {len(families)} overlapping families are in clan {clan}", clan

    # Multiple different clans
    return False, f"Overlapping families are in different clans: {', '.join(clans)}", None


def parse_triage_file_for_overlaps(triage_file, sp_only=False, min_seed=1, max_overlap_fraction=0.05):
    """
    Parse triage file and identify directories with resolvable overlaps.

    For each root directory, identifies:
    1. The best non-overlapping version (like triage_helper.py)
    2. Deeper versions with overlaps that might be resolvable (< max_overlap_fraction)

    Args:
        triage_file: Path to triage file
        sp_only: If True, only include families with SwissProt proteins
        min_seed: Minimum number of sequences required in SEED
        max_overlap_fraction: Maximum fraction of overlaps to consider for resolution (default 5%)

    Returns:
        dict with keys:
        - 'best_dirs': list of best non-overlapping directories (like triage_helper.py)
        - 'resolvable_overlaps': list of deeper directories with resolvable overlaps
    """

    # Group all entries by root directory
    root_groups = defaultdict(list)

    with open(triage_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) < 5:
                continue

            dir_name = parts[0]
            total_seqs = int(parts[1]) if parts[1].isdigit() else 0
            overlaps = int(parts[2]) if parts[2].isdigit() else 0
            non_overlap_seqs = int(parts[3]) if parts[3].isdigit() else 0

            # Column 8 contains SwissProt count (index 7 in 0-based indexing)
            swissprot_count = 0
            if len(parts) > 7 and parts[7].isdigit():
                swissprot_count = int(parts[7])

            # Skip if -sp option is set and no SwissProt proteins
            if sp_only and swissprot_count == 0:
                continue

            # Get root directory name (before first /Iterate)
            root = dir_name.split('/')[0]

            # Parse overlap details from fields 15+ (index 14+)
            overlap_details = parts[14:] if len(parts) > 14 else []
            overlap_info = parse_overlap_info(overlap_details)

            # Calculate overlap fraction
            overlap_fraction = overlaps / total_seqs if total_seqs > 0 else 0

            # Store all info
            entry = {
                'dir': dir_name,
                'total_seqs': total_seqs,
                'overlaps': overlaps,
                'non_overlap_seqs': non_overlap_seqs,
                'overlap_fraction': overlap_fraction,
                'overlap_details': overlap_details,
                'overlap_info': overlap_info,
                'depth': dir_name.count('/Iterate'),
                'root': root,
                'swissprot_count': swissprot_count
            }

            root_groups[root].append(entry)

    best_dirs = []
    resolvable_overlaps = []

    for root, entries in root_groups.items():
        # Sort by depth
        entries.sort(key=lambda x: x['depth'])

        # Find all non-overlapping directories
        non_overlap = [e for e in entries if e['overlaps'] == 0]

        # Find deeper directories with small overlap fractions
        overlap_entries = [e for e in entries if e['overlaps'] > 0]

        if non_overlap:
            # Count SEED sequences for each non-overlapping entry
            for entry in non_overlap:
                entry['seed_count'] = count_seed_sequences(entry['dir'])

            # Choose the best one using priority:
            # 1. Most sequences in ALIGN (total_seqs)
            # 2. Most sequences in SEED (seed_count)
            # 3. Deeper directories (higher depth)
            best = max(non_overlap, key=lambda x: (x['total_seqs'], x.get('seed_count', 0), x['depth']))
            best['root'] = root

            # Check for deeper overlapping versions that might be worth resolving
            deeper_resolvable = []
            for e in overlap_entries:
                if e['depth'] > best['depth']:
                    # Check if overlap fraction is small enough
                    if e['overlap_fraction'] <= max_overlap_fraction:
                        # Check if resolution is possible
                        eligible, reason, clan = check_resolution_eligibility(e['overlap_info'])
                        if eligible:
                            e['resolution_eligible'] = True
                            e['resolution_reason'] = reason
                            e['clan_to_join'] = clan
                            e['seed_count'] = count_seed_sequences(e['dir'])
                            deeper_resolvable.append(e)

            best['deeper_overlaps'] = []
            best['deeper_resolvable'] = deeper_resolvable

            # Also track non-resolvable deeper overlaps for warning
            for e in overlap_entries:
                if e['depth'] > best['depth'] and e not in deeper_resolvable:
                    overlap_summary = ' '.join(e['overlap_details'])
                    best['deeper_overlaps'].append(f"{e['dir']}: {overlap_summary}")

            # Filter by minimum SEED sequences
            if best.get('seed_count', 0) >= min_seed:
                best_dirs.append(best)

            # Add resolvable overlaps to list
            for entry in deeper_resolvable:
                entry['root'] = root
                entry['best_non_overlap'] = best
                if entry.get('seed_count', 0) >= min_seed:
                    resolvable_overlaps.append(entry)

        else:
            # No non-overlapping version exists - check if any overlapping entry is resolvable
            for e in overlap_entries:
                if e['overlap_fraction'] <= max_overlap_fraction:
                    eligible, reason, clan = check_resolution_eligibility(e['overlap_info'])
                    if eligible:
                        e['resolution_eligible'] = True
                        e['resolution_reason'] = reason
                        e['clan_to_join'] = clan
                        e['root'] = root
                        e['best_non_overlap'] = None
                        e['seed_count'] = count_seed_sequences(e['dir'])
                        if e.get('seed_count', 0) >= min_seed:
                            resolvable_overlaps.append(e)

            if not any(e in resolvable_overlaps for e in overlap_entries):
                print(f"[SKIP] {root}: No non-overlapping versions and no resolvable overlaps", file=sys.stderr)

    # Sort best_dirs by number of sequences (largest first)
    best_dirs.sort(key=lambda x: x['total_seqs'], reverse=True)

    # Sort resolvable_overlaps by total_seqs (largest first)
    resolvable_overlaps.sort(key=lambda x: x['total_seqs'], reverse=True)

    return {
        'best_dirs': best_dirs,
        'resolvable_overlaps': resolvable_overlaps
    }


def check_desc_for_clan(dir_path):
    """Check if DESC file has a CL line, return the clan accession or None"""
    desc_path = Path(dir_path) / 'DESC'
    if desc_path.exists():
        with open(desc_path, 'r') as f:
            for line in f:
                if line.startswith('CL '):
                    # Extract clan accession (format: "CL   CL0351")
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        return parts[1]
    return None


def add_clan_to_desc(dir_path, clan_accession):
    """
    Add a CL line to DESC file to make this family a member of the clan.

    The CL line should be added after the AC line.
    Format: "CL   CL0351"

    Returns True if successful, False otherwise.
    """
    desc_path = Path(dir_path) / 'DESC'
    if not desc_path.exists():
        print(f"Error: DESC file not found in {dir_path}", file=sys.stderr)
        return False

    # Check if already has a CL line
    existing_clan = check_desc_for_clan(dir_path)
    if existing_clan:
        if existing_clan == clan_accession:
            print(f"DESC already has CL line for {clan_accession}")
            return True
        else:
            print(f"Warning: DESC has CL line for {existing_clan}, not {clan_accession}", file=sys.stderr)
            response = input(f"Replace {existing_clan} with {clan_accession}? [y/n]: ").strip().lower()
            if response not in ['y', 'yes']:
                return False

    # Read current DESC
    with open(desc_path, 'r') as f:
        lines = f.readlines()

    # Find where to insert CL line (after AC line) and remove existing CL if present
    new_lines = []
    cl_inserted = False

    for line in lines:
        if line.startswith('CL '):
            # Skip existing CL line (we'll add the new one)
            continue

        new_lines.append(line)

        # Insert CL line after AC line
        if line.startswith('AC ') and not cl_inserted:
            new_lines.append(f"CL   {clan_accession}\n")
            cl_inserted = True

    if not cl_inserted:
        print(f"Warning: Could not find AC line to insert CL after", file=sys.stderr)
        # Insert at the beginning as fallback
        new_lines.insert(0, f"CL   {clan_accession}\n")

    # Write updated DESC
    with open(desc_path, 'w') as f:
        f.writelines(new_lines)

    print(f"Added CL line for {clan_accession} to DESC")
    return True


def update_se_line(dir_path, se_prefix):
    """Update SE line in DESC file with new prefix and strip /Iterate"""
    desc_path = Path(dir_path) / 'DESC'
    if desc_path.exists():
        lines = []
        updated = False
        with open(desc_path, 'r') as f:
            for line in f:
                if line.startswith('SE '):
                    # Extract the part after the colon
                    if ':' in line:
                        parts = line.strip().split(':', 1)
                        if len(parts) == 2:
                            # Get the identifier and remove any /Iterate parts
                            identifier = parts[1]
                            identifier = identifier.replace('/Iterate', '')
                            new_line = f"SE   {se_prefix}:{identifier}\n"
                            lines.append(new_line)
                            updated = True
                        else:
                            lines.append(line)
                    else:
                        lines.append(line)
                else:
                    lines.append(line)

        if updated:
            with open(desc_path, 'w') as f:
                f.writelines(lines)
            return True
    return False


def edit_seed_alignment(seed_file):
    """
    Interactive editing of SEED alignment.
    Ported from triage_merge.py with modifications.

    Returns True if alignment accepted, False if cancelled.
    """
    while True:
        # Display the alignment
        run_command(f"belvu {seed_file}", wait=True)

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
            if output:
                with open(seed_file, 'w') as f:
                    f.write(output)
                print("Applied wholeseq alignment")
            else:
                print("Warning: wholeseq may have failed", file=sys.stderr)

        elif reply in ["e", "extend"]:
            n = input("N-terminal extension (default 0): ").strip() or "0"
            c = input("C-terminal extension (default 0): ").strip() or "0"

            backup = f"{seed_file}.before_extend"
            shutil.copy(seed_file, backup)
            output = run_command(f"extend.pl -align {backup} -n {n} -c {c} -m",
                               capture_output=True)
            if output:
                with open(seed_file, 'w') as f:
                    f.write(output)
                print(f"Extended alignment N:{n} C:{c}")
            else:
                print("Warning: extend may have failed", file=sys.stderr)

        elif reply in ["p", "pad"]:
            backup = f"{seed_file}.before_pad"
            shutil.copy(seed_file, backup)
            output = run_command(f"pad_ends.pl -align {backup}", capture_output=True)
            if output:
                with open(seed_file, 'w') as f:
                    f.write(output)
                print("Padded alignment ends")
            else:
                print("Warning: pad_ends may have failed", file=sys.stderr)


def check_overlap(curation_dir, quiet=True):
    """
    Check overlaps for a curation directory.

    Note: pqc-overlap-rdb should be run from the parent directory
    and the family directory passed as argument.

    When curation_dir is '.', we're already in the family directory,
    so we need to run from parent directory.

    Args:
        curation_dir: Directory to check
        quiet: If True, suppress pqc-overlap-rdb output (default True)

    Returns: (has_overlaps, overlap_content)
    """
    curation_path = Path(curation_dir).resolve()

    # Build command with output suppression if quiet
    redirect = " > /dev/null 2>&1" if quiet else ""

    if curation_dir == '.':
        # We're in the family directory, run from parent
        parent_dir = str(curation_path.parent)
        family_name = curation_path.name
        run_command(f"pqc-overlap-rdb -no_sigP {family_name}{redirect}", cwd=parent_dir, wait=True)
        overlap_file = curation_path / 'overlap'
    else:
        # Run from current directory targeting the specified directory
        run_command(f"pqc-overlap-rdb -no_sigP {curation_dir}{redirect}", wait=True)
        overlap_file = Path(curation_dir) / 'overlap'

    # Check if overlap file exists and has content
    if not overlap_file.exists():
        return False, ""

    if overlap_file.stat().st_size == 0:
        return False, ""

    with open(overlap_file, 'r') as f:
        content = f.read()

    return True, content


def handle_threshold_change(curation_dir, threshold):
    """
    Raise gathering threshold in the curation directory.

    Note: pfmake -t sets the threshold and rebuilds. After this,
    we need to re-check overlaps. The -t flag sets the threshold
    for the pfmake command which updates the model.
    """
    print(f"Raising threshold to {threshold} in {curation_dir}")

    # Run pfmake with new threshold
    # pfmake -t THRESHOLD sets a new gathering threshold
    success = run_command(f"pfmake -t {threshold}", cwd=curation_dir)

    if success:
        print(f"Threshold set to {threshold}")
        return True
    else:
        print("Warning: pfmake may have failed", file=sys.stderr)
        return False


def move_to_directory(root_dir, destination, start_dir):
    """
    Move root directory to specified destination (DONE, IGNORE, OVERLAP).
    Always operates on root directory, not Iterate subdirectories.
    """
    os.chdir(start_dir)

    dest_path = Path(destination.upper())
    if not dest_path.exists():
        dest_path.mkdir()
        print(f"Created {destination.upper()} directory")

    # Add timestamp to avoid conflicts
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    dest_name = f"{root_dir}_{timestamp}" if (dest_path / root_dir).exists() else root_dir

    try:
        shutil.move(root_dir, str(dest_path / dest_name))
        print(f"Moved {root_dir} to {destination.upper()}/")
        return True
    except Exception as e:
        print(f"Warning: Could not move {root_dir} to {destination.upper()}: {e}", file=sys.stderr)
        return False


def display_family_info(dir_path):
    """Display relevant information files for the family"""

    # Display ted.png if it exists
    ted_png = Path(dir_path) / 'ted.png'
    if ted_png.exists():
        print("\nDisplaying ted.png...")
        run_command(f"imgcat {ted_png}", wait=False)

    # Display ted_web.png if it exists
    ted_web_png = Path(dir_path) / 'ted_web.png'
    if ted_web_png.exists():
        print("\nDisplaying ted_web.png...")
        run_command(f"imgcat {ted_web_png}", wait=False)

    # sp.seq_info content
    sp_file = Path(dir_path) / 'sp.seq_info'
    if sp_file.exists():
        print("\n--- sp.seq_info content (copy to Claude) ---", file=sys.stderr)
        with open(sp_file, 'r') as f:
            print(f.read(), file=sys.stderr)
        print("--- End of sp.seq_info ---", file=sys.stderr)
    else:
        print("\nNote: sp.seq_info file not found", file=sys.stderr)

    # species summary
    species_file = Path(dir_path) / 'species'
    print("\n--- species summary (copy to Claude) ---", file=sys.stderr)
    if species_file.exists():
        with open(species_file, 'r') as f:
            print(f.read(), file=sys.stderr)
    else:
        print("Note: species file not found", file=sys.stderr)
    print("--- End of species summary ---", file=sys.stderr)

    # Domain architectures
    arch_file = Path(dir_path) / 'arch'
    print("\n--- Domain architectures (copy to Claude) ---", file=sys.stderr)
    if arch_file.exists():
        with open(arch_file, 'r') as f:
            print(f.read(), file=sys.stderr)
    else:
        print("Note: arch file not found", file=sys.stderr)
    print("--- End of Domain architectures ---", file=sys.stderr)

    # PaperBLAST literature results
    paperblast_file = Path(dir_path) / 'paperblast'
    if paperblast_file.exists() and paperblast_file.stat().st_size > 0:
        print("\n--- PaperBLAST literature results (copy to Claude) ---", file=sys.stderr)
        with open(paperblast_file, 'r') as f:
            lines = f.readlines()
            if len(lines) > 100:
                print(''.join(lines[:100]), file=sys.stderr)
                print(f"... (truncated - showing first 100 of {len(lines)} lines)", file=sys.stderr)
            else:
                print(''.join(lines), file=sys.stderr)
        print("--- End of PaperBLAST ---", file=sys.stderr)

    # TED domain information
    ted_file = Path(dir_path) / 'TED'
    if ted_file.exists():
        print("\n--- TED domain information (copy to Claude) ---", file=sys.stderr)
        with open(ted_file, 'r') as f:
            lines = f.readlines()
            if len(lines) > 100:
                print(''.join(lines[:100]), file=sys.stderr)
                print(f"... (truncated - showing first 100 of {len(lines)} lines)", file=sys.stderr)
            else:
                print(''.join(lines), file=sys.stderr)
        print("--- End of TED ---", file=sys.stderr)

    # STRING protein interaction network
    string_file = Path(dir_path) / 'STRING'
    if string_file.exists() and string_file.stat().st_size > 0:
        print("\n--- STRING protein interactions (copy to Claude) ---", file=sys.stderr)
        with open(string_file, 'r') as f:
            lines = f.readlines()
            if len(lines) > 100:
                print(''.join(lines[:100]), file=sys.stderr)
                print(f"... (truncated)", file=sys.stderr)
            else:
                print(''.join(lines), file=sys.stderr)
        print("--- End of STRING ---", file=sys.stderr)

    # Foldseek structural matches
    foldseek_file = Path(dir_path) / 'foldseek'
    if foldseek_file.exists():
        with open(foldseek_file, 'r') as f:
            lines = f.readlines()
            if len(lines) > 1:  # More than just header
                print("\n--- Foldseek structural matches (copy to Claude) ---", file=sys.stderr)
                if len(lines) > 100:
                    print(''.join(lines[:100]), file=sys.stderr)
                    print(f"... (truncated)", file=sys.stderr)
                else:
                    print(''.join(lines), file=sys.stderr)
                print("--- End of Foldseek ---", file=sys.stderr)

    # Current DESC file
    desc_file = Path(dir_path) / 'DESC'
    print("\n--- Current DESC file (copy to Claude) ---", file=sys.stderr)
    if desc_file.exists():
        with open(desc_file, 'r') as f:
            print(f.read(), file=sys.stderr)
    else:
        print("Note: DESC file not found", file=sys.stderr)
    print("--- End of DESC ---", file=sys.stderr)


def resolve_overlap(entry, start_dir, se_prefix=None, use_nano=False):
    """
    Guide curator through overlap resolution for a single directory.

    Resolution options:
    1. Edit SEED to remove overlapping sequences
    2. Raise threshold
    3. Join clan (if overlapping family is in a clan)
    4. Skip
    5. Move to IGNORE

    Returns True if resolved and ready for pfnew, False otherwise.
    """
    dir_name = entry['dir']
    root_dir = entry['root']

    print(f"\n{'='*70}")
    print(f"OVERLAP RESOLUTION: {dir_name}")
    print(f"{'='*70}")
    print(f"Total sequences: {entry['total_seqs']}")
    print(f"Overlapping sequences: {entry['overlaps']} ({entry['overlap_fraction']*100:.1f}%)")
    print(f"Non-overlapping sequences: {entry['non_overlap_seqs']}")
    if 'seed_count' in entry:
        print(f"SEED sequences: {entry['seed_count']}")
    if entry.get('swissprot_count', 0) > 0:
        print(f"SwissProt proteins: {entry['swissprot_count']}")

    print(f"\nResolution reason: {entry.get('resolution_reason', 'Unknown')}")

    # Show overlap details
    print("\nOverlap details:")
    for oi in entry.get('overlap_info', []):
        clan_str = f" (clan {oi['clan']})" if oi['clan'] else ""
        print(f"  - {oi['family']}{clan_str}: {oi['overlap_count']} overlaps (Pfam size: {oi['pfam_size']})")

    clan_to_join = entry.get('clan_to_join')

    # Check if directory exists
    os.chdir(start_dir)
    if not Path(dir_name).exists():
        print(f"Error: Directory {dir_name} not found, skipping", file=sys.stderr)
        return False

    os.chdir(dir_name)

    # Run initial overlap check
    has_overlaps, overlap_content = check_overlap('.')
    if has_overlaps:
        print("\nCurrent overlap file:")
        print("-" * 40)
        print(overlap_content)
        print("-" * 40)
    else:
        print("\nNo overlaps found - already resolved!")
        os.chdir(start_dir)
        return True  # Ready to continue

    # Main resolution loop
    while True:
        print("\n" + "="*60)
        print("Overlap resolution options:")
        print("  1 - Edit SEED to remove overlapping sequences")
        print("  2 - Raise gathering threshold")
        if clan_to_join:
            print(f"  3 - Join clan {clan_to_join} (resolves overlap with clan members)")
        print("  s - Skip this directory (come back later)")
        print("  i - Move to IGNORE directory")
        print("="*60)

        choice = input("Choice: ").strip().lower()

        if choice == '1':
            # Edit SEED
            seed_file = Path('.') / 'SEED'
            if not seed_file.exists():
                print("Error: SEED file not found", file=sys.stderr)
                continue

            if edit_seed_alignment(str(seed_file)):
                print("Running pfbuild -withpfmake...")
                run_command("pfbuild -withpfmake", wait=False)
                print("pfbuild started in background. Come back to check overlaps later.")
                os.chdir(start_dir)
                return False  # Not ready yet, pfbuild running

        elif choice == '2':
            # Raise threshold
            threshold = input("New threshold value: ").strip()
            if threshold:
                try:
                    float(threshold)  # Validate it's a number
                    handle_threshold_change('.', threshold)

                    # Check overlaps again after threshold change
                    print("\nRe-checking overlaps...")
                    has_overlaps, overlap_content = check_overlap('.')
                    if has_overlaps:
                        print("Overlaps remaining:")
                        print(overlap_content)
                    else:
                        print("All overlaps resolved!")
                        os.chdir(start_dir)
                        return True
                except ValueError:
                    print("Invalid threshold value, please enter a number")

        elif choice == '3' and clan_to_join:
            # Join clan
            print(f"\nJoining clan {clan_to_join}...")
            print("Note: This makes the family a clan member. Overlaps with other")
            print("      clan members will be allowed, but the family will 'compete'")
            print("      with other clan members for sequences.")

            confirm = input(f"Proceed with joining clan {clan_to_join}? [y/n]: ").strip().lower()
            if confirm in ['y', 'yes']:
                if add_clan_to_desc('.', clan_to_join):
                    # Re-run overlap check with -compete flag for clan members
                    # Need to run from parent directory, suppress stdout/stderr
                    print("\nRe-checking overlaps (with clan membership)...")
                    curation_path = Path('.').resolve()
                    parent_dir = str(curation_path.parent)
                    family_name = curation_path.name
                    run_command(f"pqc-overlap-rdb -no_sigP -compete {family_name} > /dev/null 2>&1",
                               cwd=parent_dir, wait=True)

                    overlap_file = Path('overlap')
                    if overlap_file.exists() and overlap_file.stat().st_size > 0:
                        with open(overlap_file, 'r') as f:
                            overlap_content = f.read()
                        print("Overlaps remaining (non-clan):")
                        print("-" * 40)
                        print(overlap_content)
                        print("-" * 40)
                    else:
                        print("All overlaps resolved (clan overlaps are now allowed)!")
                        os.chdir(start_dir)
                        return True

        elif choice == 's':
            print("Skipping this directory")
            os.chdir(start_dir)
            return False

        elif choice == 'i':
            os.chdir(start_dir)
            move_to_directory(root_dir, 'IGNORE', start_dir)
            return False

        else:
            print("Unrecognized option, please try again.")


def continue_to_pfnew(entry, start_dir, se_prefix=None, use_nano=False):
    """
    Continue with annotation and pfnew after overlap resolution.
    Based on curate_family() from triage_helper.py.
    """
    dir_name = entry['dir']
    root_dir = entry['root']

    os.chdir(start_dir)
    if not Path(dir_name).exists():
        print(f"Error: Directory {dir_name} not found", file=sys.stderr)
        return False

    os.chdir(dir_name)

    # Update SE line if prefix provided
    if se_prefix:
        if update_se_line('.', se_prefix):
            print(f"Updated SE line with prefix '{se_prefix}'")

    # Count sequences
    seed_count = count_sequences('SEED')
    align_count = count_sequences('ALIGN')
    print(f"\nSequences in SEED: {seed_count}")
    print(f"Sequences in ALIGN: {align_count}")

    # Launch belvu for SEED review first
    print("\nLaunching belvu for SEED alignment review...")
    run_command("belvu SEED", wait=True)

    # SEED editing options
    while True:
        print("\n" + "="*60)
        print("SEED editing options:")
        print("  y/yes     - Happy with SEED, continue to annotation")
        print("  n/no      - Skip to next family")
        print("  i/ignore  - Move to IGNORE directory")
        print("  w/wholeseq - Run wholeseq on SEED")
        print("  e/extend  - Extend N/C termini")
        print("  p/pad     - Pad alignment ends")
        print("="*60)

        response = input("\nWhat would you like to do? ").strip().lower()

        if response in ['y', 'yes']:
            print("\nSEED looks good, continuing to annotation...")
            break

        elif response in ['n', 'no']:
            os.chdir(start_dir)
            print("Skipping this family")
            return False

        elif response in ['i', 'ignore']:
            os.chdir(start_dir)
            move_to_directory(root_dir, 'IGNORE', start_dir)
            return False

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

    # Display family information
    display_family_info('.')

    # Ask if wish to continue to annotation
    print("\nPlease create annotation using Claude with the information above.")
    while True:
        response = input("Do you wish to continue to edit DESC? (y/n/i for ignore): ").strip().lower()
        if response in ['y', 'n', 'i']:
            break
        print("Please enter 'y', 'n', or 'i'")

    if response == 'i':
        os.chdir(start_dir)
        move_to_directory(root_dir, 'IGNORE', start_dir)
        return False

    if response == 'n':
        os.chdir(start_dir)
        print("Skipping this family")
        return False

    # Open editor for DESC editing
    editor = "nano" if use_nano else "emacs -nw"
    print(f"\nOpening {editor.split()[0]} for DESC file editing...")
    run_command(f"{editor} DESC", wait=True)

    # Ask if wish to pfnew
    while True:
        response = input("\nDo you wish to pfnew this family? (y/n/e to edit DESC again/d for DUF/i for ignore): ").strip().lower()
        if response in ['y', 'n', 'e', 'd', 'i']:
            break
        print("Please enter 'y', 'n', 'e', 'd', or 'i'")

    if response == 'd':
        # Run nextDUF.pl to add DUF number
        print("\nRunning nextDUF.pl to assign DUF number...")
        run_command("nextDUF.pl", wait=True)
        # Recursively ask again
        return continue_to_pfnew(entry, start_dir, se_prefix, use_nano)

    if response == 'e':
        # Edit DESC file again
        print(f"\nOpening {editor.split()[0]} for DESC file editing...")
        run_command(f"{editor} DESC", wait=True)
        # Recursively ask again
        return continue_to_pfnew(entry, start_dir, se_prefix, use_nano)

    if response == 'i':
        os.chdir(start_dir)
        move_to_directory(root_dir, 'IGNORE', start_dir)
        return False

    if response == 'y':
        # Check for clan
        has_clan = check_desc_for_clan('.')

        # Go up one directory to parent
        os.chdir('..')

        # Get just the last directory name for pfnew
        last_dir = dir_name.split('/')[-1]

        # Run pfnew with appropriate flags
        if has_clan:
            cmd = f"pfnew {last_dir} -add_to_clan"
            print(f"Running: {cmd}")
        else:
            cmd = f"pfnew {last_dir}"
            print(f"Running: {cmd}")

        pfnew_success = run_command(cmd, wait=True)

        if pfnew_success:
            print(f"Successfully added {last_dir} to Pfam")
            os.chdir(start_dir)
            move_to_directory(root_dir, 'DONE', start_dir)
            return True
        else:
            print(f"\nERROR: pfnew failed for {last_dir}", file=sys.stderr)
            print("Please review the error above.", file=sys.stderr)

            while True:
                response = input("\nWhat would you like to do? (c=continue/o=move to OVERLAP/i=move to IGNORE): ").strip().lower()
                if response in ['c', 'o', 'i']:
                    break
                print("Please enter 'c', 'o', or 'i'")

            os.chdir(start_dir)

            if response == 'o':
                move_to_directory(root_dir, 'OVERLAP', start_dir)
            elif response == 'i':
                move_to_directory(root_dir, 'IGNORE', start_dir)

            return False

    else:
        os.chdir(start_dir)
        print("Skipped pfnew for this family")
        return False


def process_resolvable_overlaps(resolvable_overlaps, start_dir, se_prefix=None, use_nano=False):
    """
    Process all resolvable overlap directories.

    For each directory:
    1. Attempt to resolve overlaps
    2. If resolved, continue to annotation and pfnew
    """
    if not resolvable_overlaps:
        print("\nNo resolvable overlaps found")
        return 0

    print(f"\nFound {len(resolvable_overlaps)} directories with potentially resolvable overlaps")
    print("\nResolvable overlap directories:")
    for i, entry in enumerate(resolvable_overlaps, 1):
        clan_info = f" [can join {entry['clan_to_join']}]" if entry.get('clan_to_join') else ""
        print(f"  {i}. {entry['dir']} - {entry['overlaps']} overlaps ({entry['overlap_fraction']*100:.1f}%){clan_info}")

    families_added = 0

    for i, entry in enumerate(resolvable_overlaps, 1):
        print(f"\n[{i}/{len(resolvable_overlaps)}]")

        # First resolve overlaps
        if resolve_overlap(entry, start_dir, se_prefix, use_nano):
            # Overlaps resolved, continue to pfnew
            if continue_to_pfnew(entry, start_dir, se_prefix, use_nano):
                families_added += 1

    return families_added


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Pfam Triage Curation Helper v2 - with overlap resolution",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script extends triage_helper.py to handle deeper Iterate directories
that have overlaps. It offers overlap resolution options including:

  - Editing SEED to remove overlapping sequences
  - Raising the gathering threshold
  - Joining a clan (if overlapping families are all in the same clan)

The script only attempts resolution when:
  - Overlaps are less than 5% of total sequences (configurable)
  - There is exactly ONE overlapping family, OR
  - All overlapping families are in the SAME clan

Example usage:
  python triage_helper2.py triage -s TED
  python triage_helper2.py triage --overlap-only
  python triage_helper2.py triage --max-overlap 0.10  # Allow up to 10% overlaps
        """
    )

    parser.add_argument('triage_file', help='Path to triage file')
    parser.add_argument('-s', '--se-prefix', type=str, default=None,
                       help='Update SE line in DESC with this prefix')
    parser.add_argument('-nano', '--use-nano', action='store_true',
                       help='Use nano instead of emacs for editing')
    parser.add_argument('-sp', '--swissprot-only', action='store_true',
                       help='Only process families with SwissProt proteins')
    parser.add_argument('--min-seed', type=int, default=1,
                       help='Minimum number of sequences required in SEED (default: 1)')
    parser.add_argument('--max-overlap', type=float, default=0.05,
                       help='Maximum overlap fraction to consider for resolution (default: 0.05 = 5%%)')
    parser.add_argument('--overlap-only', action='store_true',
                       help='Only process directories with resolvable overlaps (skip non-overlapping)')

    args = parser.parse_args()

    if not os.path.exists(args.triage_file):
        print(f"Error: Triage file '{args.triage_file}' not found", file=sys.stderr)
        sys.exit(1)

    print(f"Processing {args.triage_file}...", file=sys.stderr)
    print(f"Maximum overlap fraction for resolution: {args.max_overlap*100:.0f}%", file=sys.stderr)

    if args.swissprot_only:
        print("Only processing families with SwissProt proteins", file=sys.stderr)

    # Parse triage file
    results = parse_triage_file_for_overlaps(
        args.triage_file,
        sp_only=args.swissprot_only,
        min_seed=args.min_seed,
        max_overlap_fraction=args.max_overlap
    )

    best_dirs = results['best_dirs']
    resolvable_overlaps = results['resolvable_overlaps']

    print(f"\nFound {len(best_dirs)} families with non-overlapping versions", file=sys.stderr)
    print(f"Found {len(resolvable_overlaps)} directories with potentially resolvable overlaps", file=sys.stderr)

    # Count how many best_dirs have deeper resolvable overlaps
    dirs_with_deeper_resolvable = sum(1 for d in best_dirs if d.get('deeper_resolvable'))
    if dirs_with_deeper_resolvable > 0:
        print(f"  ({dirs_with_deeper_resolvable} families have deeper resolvable overlaps)", file=sys.stderr)

    start_dir = os.getcwd()
    families_added = 0

    if args.overlap_only:
        # Only process resolvable overlaps
        print("\n--- Processing resolvable overlaps only ---")
        families_added = process_resolvable_overlaps(
            resolvable_overlaps, start_dir, args.se_prefix, args.use_nano
        )
    else:
        # Process non-overlapping directories first, then offer to resolve deeper overlaps
        print("\n--- Processing non-overlapping directories ---")

        # This would require importing/copying curate_family from triage_helper.py
        # For now, just report what would be processed
        print(f"\nTo process non-overlapping directories, use triage_helper.py")
        print(f"This script focuses on overlap resolution.")
        print(f"\nWould you like to process the resolvable overlaps? [y/n]")

        response = input().strip().lower()
        if response in ['y', 'yes']:
            families_added = process_resolvable_overlaps(
                resolvable_overlaps, start_dir, args.se_prefix, args.use_nano
            )

    print(f"\n{'='*60}")
    print(f"Complete! Added {families_added} families to Pfam")


if __name__ == "__main__":
    main()
