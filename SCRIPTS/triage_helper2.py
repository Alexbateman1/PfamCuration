#!/usr/bin/env python3
"""
Pfam Triage Curation Helper Script - Version 2

This script works exactly like triage_helper.py for building families, with the
added functionality that curators can optionally try to resolve overlaps in
deeper Iterate directories when:
  - The overlap fraction is small (< 5% by default)
  - There is a single overlapping family, OR all overlapping families are in the same clan

Workflow:
1. Parse triage file and find best non-overlapping directory for each family
2. For families with deeper versions that have small/resolvable overlaps, offer a choice:
   - Try to resolve overlaps in the deeper version (more sequences), OR
   - Use the shallower non-overlapping version
3. If overlap resolution fails or is skipped, fall back to the non-overlapping version
4. Continue with normal curation workflow (SEED review, annotation, pfnew)
"""

import sys
import os
import subprocess
import shutil
import re
import json
from pathlib import Path
from collections import defaultdict
from datetime import datetime


def filter_swissprot_content(content):
    """
    Filter SwissProt entry content to reduce token usage by removing:
    - References containing in RP lines: NUCLEOTIDE SEQUENCE, LARGE SCALE ANALYSIS,
      TISSUE SPECIFICITY, or VARIANT
    - DR lines except for STRING, InterPro, Pfam, SMART, Gene3D
    - FT sections for COMPBIAS, MOD_RES, VAR_SEQ, CONFLICT, VARIANT,
      HELIX, STRAND, TURN, MUTAGEN, and disordered REGION

    Args:
        content: Raw SwissProt entry content as a string

    Returns:
        Filtered content string
    """
    lines = content.split('\n')
    filtered_lines = []

    # Databases to keep in DR lines
    keep_databases = {'STRING', 'InterPro', 'Pfam', 'SMART', 'Gene3D'}

    # FT feature types to skip
    skip_ft_types = {'COMPBIAS', 'MOD_RES', 'VAR_SEQ', 'CONFLICT', 'VARIANT',
                     'HELIX', 'STRAND', 'TURN', 'MUTAGEN'}

    # RP line patterns that trigger reference exclusion
    skip_rp_patterns = [
        'NUCLEOTIDE SEQUENCE',
        'LARGE SCALE ANALYSIS',
        'TISSUE SPECIFICITY',
    ]

    i = 0
    while i < len(lines):
        line = lines[i]

        # Handle reference blocks (RN starts a reference)
        if line.startswith('RN   '):
            # Collect the entire reference block
            ref_block = [line]
            i += 1
            while i < len(lines) and lines[i][:2] in ('RP', 'RC', 'RX', 'RG', 'RA', 'RT', 'RL'):
                ref_block.append(lines[i])
                i += 1

            # Check if any RP line contains patterns we want to skip
            skip_reference = False
            for ref_line in ref_block:
                if ref_line.startswith('RP   '):
                    # Check for skip patterns
                    for pattern in skip_rp_patterns:
                        if pattern in ref_line:
                            skip_reference = True
                            break
                    # Also check if RP line is just "VARIANT" (e.g., "RP   VARIANT ARG-763.")
                    if ref_line.startswith('RP   VARIANT'):
                        skip_reference = True
                    if skip_reference:
                        break

            # Only include reference if it doesn't match skip patterns
            if not skip_reference:
                filtered_lines.extend(ref_block)
            continue

        # Handle DR (database reference) lines
        if line.startswith('DR   '):
            # Extract database name (first field after DR)
            # Format: DR   DATABASE; ...
            parts = line[5:].split(';')
            if parts:
                db_name = parts[0].strip()
                if db_name in keep_databases:
                    filtered_lines.append(line)
            i += 1
            continue

        # Handle FT (feature table) lines
        if line.startswith('FT   '):
            # In UniProt flat file format:
            # - New feature lines: "FT   TYPE            position" - feature type at column 5
            # - Continuation lines: "FT                   /qualifier" or "FT                   text" - whitespace at column 5
            # Check if this starts a new feature (character at position 5 is not a space)
            is_continuation = len(line) > 5 and line[5] == ' '

            if not is_continuation:
                # This is a new feature type line
                ft_content = line[5:].strip()
                parts = ft_content.split()
                if parts:
                    feature_type = parts[0]

                    # Check if it's a feature type to skip
                    if feature_type in skip_ft_types:
                        # Skip this feature and all its continuation lines
                        i += 1
                        while i < len(lines) and lines[i].startswith('FT   '):
                            if len(lines[i]) > 5 and lines[i][5] == ' ':
                                # Continuation line, skip it
                                i += 1
                            else:
                                # New feature type, stop skipping
                                break
                        continue

                    elif feature_type == 'REGION':
                        # Need to check if it's a disordered region by looking ahead
                        j = i + 1
                        feature_lines = [line]
                        while j < len(lines) and lines[j].startswith('FT   '):
                            if len(lines[j]) > 5 and lines[j][5] == ' ':
                                # Continuation line
                                feature_lines.append(lines[j])
                                j += 1
                            else:
                                # New feature type
                                break

                        # Check if any line contains "Disordered"
                        feature_text = '\n'.join(feature_lines)
                        if 'Disordered' in feature_text:
                            # Skip the entire disordered REGION block
                            i = j
                            continue
                        # Otherwise, keep REGION - fall through to append

            # Keep this FT line (either a kept feature type or a continuation of a kept feature)
            filtered_lines.append(line)
            i += 1
            continue

        # Keep all other lines
        filtered_lines.append(line)
        i += 1

    return '\n'.join(filtered_lines)


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

    Format: "[CLAN:]FAMILY COUNT(SIZE)" e.g. "CL0351:CHCH 1050(8242)" or "Flavoprotein 2(21413)"
    """
    overlaps = []
    debug = os.environ.get('DEBUG_TRIAGE', '')

    for detail in overlap_details:
        detail = detail.strip() if detail else ""
        if not detail or detail == '----' or detail == 'NA':
            continue

        # More flexible regex to handle whitespace variations
        match = re.match(r'^(?:(CL\d+):)?(\S+)\s+(\d+)\((\d+)\)\s*$', detail)
        if match:
            overlaps.append({
                'family': match.group(2),
                'clan': match.group(1),  # May be None
                'overlap_count': int(match.group(3)),
                'pfam_size': int(match.group(4)),
                'raw': detail
            })
        elif debug:
            print(f"[DEBUG] parse_overlap_info: no match for '{detail}'", file=sys.stderr)

    return overlaps


def check_resolution_eligibility(overlap_info_list):
    """
    Check if overlap resolution should be offered.

    Rules:
    - Only offer if there's exactly ONE overlapping family, OR
    - All overlapping families are in the SAME clan

    Returns: (eligible, reason, clan_to_join)
    """
    if not overlap_info_list:
        return False, "No overlap information", None

    families = set(oi['family'] for oi in overlap_info_list)

    if len(families) == 1:
        oi = overlap_info_list[0]
        if oi['clan']:
            return True, f"Single family ({oi['family']}) in clan {oi['clan']}", oi['clan']
        else:
            return True, f"Single family ({oi['family']}) not in a clan", None

    # Multiple families - check if all in same clan
    clans = set(oi['clan'] for oi in overlap_info_list if oi['clan'])
    families_without_clan = [oi['family'] for oi in overlap_info_list if not oi['clan']]

    if families_without_clan:
        return False, f"Multiple families, some not in clans", None

    if len(clans) == 1:
        clan = list(clans)[0]
        return True, f"All {len(families)} families in clan {clan}", clan

    return False, f"Families in different clans", None


def parse_triage_file(triage_file, sp_only=False, min_seed=1, min_align=0, max_overlap_fraction=0.05):
    """
    Parse triage file and identify best directories to work on.

    For each root directory, identifies:
    - The best non-overlapping version (primary choice)
    - Any deeper versions with small, resolvable overlaps (optional alternative)

    Args:
        triage_file: Path to triage file
        sp_only: If True, only include families with SwissProt proteins
        min_seed: Minimum number of sequences required in SEED (default 1)
        min_align: Minimum number of sequences required in ALIGN (default 0, no filtering)
        max_overlap_fraction: Maximum overlap fraction for resolution (default 0.05)
    """
    root_groups = defaultdict(list)

    with open(triage_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Split and strip each field to handle any whitespace issues
            parts = [p.strip() for p in line.split('\t')]
            if len(parts) < 5:
                continue

            dir_name = parts[0]

            # Parse numeric fields safely
            try:
                total_seqs = int(parts[1]) if parts[1].isdigit() else 0
                overlaps = int(parts[2]) if parts[2].isdigit() else 0
                non_overlap_seqs = int(parts[3]) if parts[3].isdigit() else 0
            except (ValueError, IndexError):
                continue

            swissprot_count = 0
            if len(parts) > 7 and parts[7].isdigit():
                swissprot_count = int(parts[7])

            if sp_only and swissprot_count == 0:
                continue

            root = dir_name.split('/')[0]
            overlap_details = parts[14:] if len(parts) > 14 else []
            overlap_info = parse_overlap_info(overlap_details)
            overlap_fraction = overlaps / total_seqs if total_seqs > 0 else 0

            debug = os.environ.get('DEBUG_TRIAGE', '')
            if debug and overlaps > 0:
                print(f"[DEBUG] Parsing {dir_name}: parts count={len(parts)}, "
                      f"overlap_details={overlap_details}, overlap_info={len(overlap_info)} parsed",
                      file=sys.stderr)

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
    debug = os.environ.get('DEBUG_TRIAGE', '')

    for root, entries in root_groups.items():
        if debug:
            print(f"[DEBUG] Processing root: {root}", file=sys.stderr)
            for e in entries:
                print(f"[DEBUG]   {e['dir']}: depth={e['depth']}, overlaps={e['overlaps']}, "
                      f"total={e['total_seqs']}, overlap_info={e['overlap_info']}", file=sys.stderr)

        # Find non-overlapping directories
        non_overlap = [e for e in entries if e['overlaps'] == 0]

        if not non_overlap:
            print(f"[SKIP] {root}: No non-overlapping versions found", file=sys.stderr)
            continue

        # Count SEED sequences
        for entry in non_overlap:
            entry['seed_count'] = count_seed_sequences(entry['dir'])

        # Choose best non-overlapping version
        best = max(non_overlap, key=lambda x: (x['total_seqs'], x.get('seed_count', 0), x['depth']))
        best['root'] = root

        if debug:
            print(f"[DEBUG] Best non-overlapping: {best['dir']} (depth={best['depth']})", file=sys.stderr)

        # Find deeper versions with small, resolvable overlaps
        deeper_resolvable = []
        deeper_warnings = []

        for e in entries:
            if e['depth'] > best['depth'] and e['overlaps'] > 0:
                if debug:
                    print(f"[DEBUG] Checking {e['dir']}: overlap_fraction={e['overlap_fraction']:.4f}, "
                          f"max={max_overlap_fraction}, overlap_info_count={len(e['overlap_info'])}", file=sys.stderr)

                if e['overlap_fraction'] <= max_overlap_fraction:
                    # Check if we have overlap info to work with
                    if not e['overlap_info']:
                        if debug:
                            print(f"[DEBUG]   No overlap_info parsed - cannot determine eligibility", file=sys.stderr)
                        deeper_warnings.append(f"{e['dir']}: Overlap info not parsed")
                        continue

                    eligible, reason, clan = check_resolution_eligibility(e['overlap_info'])
                    if debug:
                        print(f"[DEBUG]   eligible={eligible}, reason={reason}", file=sys.stderr)

                    if eligible:
                        e['seed_count'] = count_seed_sequences(e['dir'])
                        e['resolution_reason'] = reason
                        e['clan_to_join'] = clan
                        deeper_resolvable.append(e)
                    else:
                        deeper_warnings.append(f"{e['dir']}: {reason}")
                else:
                    overlap_summary = ' '.join(e['overlap_details'])
                    deeper_warnings.append(f"{e['dir']}: {overlap_summary}")

        if debug:
            print(f"[DEBUG] deeper_resolvable: {[e['dir'] for e in deeper_resolvable]}", file=sys.stderr)
            print(f"[DEBUG] deeper_warnings: {deeper_warnings}", file=sys.stderr)

        best['deeper_resolvable'] = deeper_resolvable
        best['deeper_overlaps'] = deeper_warnings

        # Count ALIGN sequences for min_align filtering
        best['align_count'] = count_align_sequences(best['dir'])

        # Filter by minimum SEED and minimum ALIGN
        if best.get('seed_count', 0) >= min_seed:
            if min_align > 0 and best.get('align_count', 0) < min_align:
                print(f"[SKIP] {root}: ALIGN has only {best.get('align_count', 0)} sequences (min: {min_align})", file=sys.stderr)
                continue
            best_dirs.append(best)

    # Sort by total sequences (largest first), then group by protein
    best_dirs.sort(key=lambda x: x['total_seqs'], reverse=True)

    # Group domains from same protein together
    protein_groups = defaultdict(list)
    for entry in best_dirs:
        root = entry['root']
        if '_TED' in root:
            protein_acc = root.split('_TED')[0]
        else:
            protein_acc = root
        protein_groups[protein_acc].append(entry)

    protein_max_size = {p: max(e['total_seqs'] for e in entries)
                        for p, entries in protein_groups.items()}
    sorted_proteins = sorted(protein_groups.keys(),
                           key=lambda p: protein_max_size[p], reverse=True)

    grouped_dirs = []
    for protein_acc in sorted_proteins:
        domains = protein_groups[protein_acc]
        domains.sort(key=lambda x: (-x['total_seqs'], x['root']))
        grouped_dirs.extend(domains)

    return grouped_dirs


def check_desc_for_clan(dir_path):
    """Check if DESC file has a CL line, return the clan accession or None"""
    desc_path = Path(dir_path) / 'DESC'
    if desc_path.exists():
        with open(desc_path, 'r') as f:
            for line in f:
                if line.startswith('CL '):
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        return parts[1]
    return None


def add_clan_to_desc(dir_path, clan_accession):
    """Add a CL line to DESC file after the TP line."""
    desc_path = Path(dir_path) / 'DESC'
    if not desc_path.exists():
        print(f"Error: DESC file not found in {dir_path}", file=sys.stderr)
        return False

    existing_clan = check_desc_for_clan(dir_path)
    if existing_clan:
        if existing_clan == clan_accession:
            print(f"DESC already has CL line for {clan_accession}")
            return True
        else:
            response = input(f"DESC has CL {existing_clan}, replace with {clan_accession}? [y/n]: ").strip().lower()
            if response not in ['y', 'yes']:
                return False

    with open(desc_path, 'r') as f:
        lines = f.readlines()

    new_lines = []
    cl_inserted = False
    for line in lines:
        if line.startswith('CL '):
            continue  # Remove existing CL line
        new_lines.append(line)
        # Insert CL after TP line
        if not cl_inserted and line.startswith('TP '):
            new_lines.append(f"CL   {clan_accession}\n")
            cl_inserted = True

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
                    if ':' in line:
                        parts = line.strip().split(':', 1)
                        if len(parts) == 2:
                            identifier = parts[1].replace('/Iterate', '')
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


def check_overlap_quiet(curation_dir):
    """Run pqc-overlap-rdb and return (has_overlaps, overlap_content), suppressing output."""
    curation_path = Path(curation_dir).resolve()

    if curation_dir == '.':
        parent_dir = str(curation_path.parent)
        family_name = curation_path.name
        run_command(f"pqc-overlap-rdb -no_sigP {family_name} > /dev/null 2>&1",
                   cwd=parent_dir, wait=True)
        overlap_file = curation_path / 'overlap'
    else:
        run_command(f"pqc-overlap-rdb -no_sigP {curation_dir} > /dev/null 2>&1", wait=True)
        overlap_file = Path(curation_dir) / 'overlap'

    if not overlap_file.exists() or overlap_file.stat().st_size == 0:
        return False, ""

    with open(overlap_file, 'r') as f:
        content = f.read()
    return True, content


def get_pfam_accession(family_name, curation_dir):
    """
    Extract Pfam accession (PF#####) from overlap file for a family name.

    Args:
        family_name: Short name of the Pfam family (e.g., 'TF_AP-2')
        curation_dir: Path to the curation directory containing overlap file

    Returns:
        Pfam accession (e.g., 'PF03299') or None if not found
    """
    # If it already looks like an accession, return it
    if re.match(r'^PF\d{5}$', family_name):
        return family_name

    # Parse the overlap file to find the accession
    overlap_file = Path(curation_dir) / 'overlap'
    if overlap_file.exists():
        with open(overlap_file, 'r') as f:
            for line in f:
                # Look for pattern like "TF_AP-2 PF03299" in the overlap line
                match = re.search(rf'{re.escape(family_name)}\s+(PF\d{{5}})', line)
                if match:
                    return match.group(1)

    return None


def raise_pfam_threshold(family_name, curation_dir, threshold):
    """
    Raise the gathering threshold in an existing Pfam family.

    This checks out the Pfam family, runs pfmake with the new threshold,
    and checks it back in.

    Args:
        family_name: Short name of the Pfam family
        curation_dir: Path to the curation directory
        threshold: New threshold value

    Returns:
        True if successful, False otherwise
    """
    # Get Pfam accession
    pfam_acc = get_pfam_accession(family_name, curation_dir)
    if not pfam_acc:
        print(f"Error: Could not determine Pfam accession for {family_name}", file=sys.stderr)
        return False

    print(f"\nRaising threshold for {family_name} ({pfam_acc})")

    # Check out the Pfam family
    pfam_dir = Path(curation_dir) / pfam_acc
    if not pfam_dir.exists():
        print(f"Checking out {pfam_acc}...")
        success = run_command(f"pfco {pfam_acc}", cwd=curation_dir, wait=True)
        if not success:
            print(f"Error: Failed to check out {pfam_acc}", file=sys.stderr)
            return False

    # Run pfmake with new threshold
    print(f"Running pfmake -t {threshold} in {pfam_acc}...")
    success = run_command(f"pfmake -t {threshold}", cwd=str(pfam_dir), wait=True)
    if not success:
        print(f"Error: pfmake failed", file=sys.stderr)
        return False

    # Check in the family
    print(f"Checking in {pfam_acc} with new threshold...")
    success = run_command(
        f"pfci -m 'Raised threshold to {threshold} to resolve overlaps' {pfam_acc}",
        cwd=curation_dir, wait=True
    )
    if not success:
        print(f"Error: pfci failed", file=sys.stderr)
        return False

    print(f"Successfully raised threshold for {pfam_acc}")
    return True


def try_resolve_overlaps(entry, start_dir):
    """
    Attempt to resolve overlaps for a deeper directory.

    Returns:
        True if overlaps resolved successfully
        False if resolution failed/skipped (should fall back to non-overlapping)
    """
    dir_name = entry['dir']
    clan_to_join = entry.get('clan_to_join')

    print(f"\n{'='*60}")
    print(f"ATTEMPTING OVERLAP RESOLUTION: {dir_name}")
    print(f"{'='*60}")
    print(f"Total sequences: {entry['total_seqs']}")
    print(f"Overlapping: {entry['overlaps']} ({entry['overlap_fraction']*100:.1f}%)")
    print(f"Reason eligible: {entry.get('resolution_reason', 'Unknown')}")

    print("\nOverlap details:")
    for oi in entry.get('overlap_info', []):
        clan_str = f" (clan {oi['clan']})" if oi['clan'] else ""
        print(f"  - {oi['family']}{clan_str}: {oi['overlap_count']} overlaps")

    os.chdir(start_dir)
    if not Path(dir_name).exists():
        print(f"Error: Directory {dir_name} not found", file=sys.stderr)
        return False

    os.chdir(dir_name)

    # Check current overlaps
    has_overlaps, overlap_content = check_overlap_quiet('.')
    if not has_overlaps:
        print("No overlaps found - already resolved!")
        os.chdir(start_dir)
        return True

    print(f"\nCurrent overlaps ({entry['overlaps']} sequences):")
    print("-" * 40)
    # Show first 20 lines only
    lines = overlap_content.strip().split('\n')
    for line in lines[:20]:
        print(line)
    if len(lines) > 20:
        print(f"... ({len(lines) - 20} more lines)")
    print("-" * 40)

    # Get the overlapping family name for option 3 (raise Pfam threshold)
    overlap_families = entry.get('overlap_info', [])
    pfam_family_name = overlap_families[0]['family'] if overlap_families else None

    # Resolution menu
    while True:
        print("\nOverlap resolution options:")
        print("  1 - Edit SEED (remove overlapping sequences)")
        print("  2 - Raise threshold in this curation directory")
        if pfam_family_name:
            print(f"  3 - Raise threshold in Pfam family ({pfam_family_name})")
        if clan_to_join:
            print(f"  4 - Join clan {clan_to_join}")
        print("  5 - Add ED lines (exclude overlapping sequences)")
        print("  s - Skip (use non-overlapping version instead)")
        print("  i - Move to IGNORE")

        choice = input("Choice: ").strip().lower()

        if choice == '1':
            # Edit SEED with belvu
            print("\nLaunching belvu to edit SEED...")
            run_command("belvu SEED", wait=True)

            print("\nRunning pfbuild -withpfmake...")
            run_command("pfbuild -withpfmake", wait=False)
            print("pfbuild started. Please come back when complete.")
            print("(Falling back to non-overlapping version for now)")
            os.chdir(start_dir)
            return False

        elif choice == '2':
            threshold = input("New threshold value for curation directory: ").strip()
            if threshold:
                try:
                    float(threshold)
                    print(f"Running pfmake -t {threshold}...")
                    run_command(f"pfmake -t {threshold}", wait=True)

                    print("Re-checking overlaps...")
                    has_overlaps, overlap_content = check_overlap_quiet('.')
                    if not has_overlaps:
                        print("All overlaps resolved!")
                        os.chdir(start_dir)
                        return True
                    else:
                        print(f"Overlaps still present:")
                        lines = overlap_content.strip().split('\n')
                        for line in lines[:10]:
                            print(f"  {line}")
                        if len(lines) > 10:
                            print(f"  ... ({len(lines) - 10} more)")
                except ValueError:
                    print("Invalid threshold")

        elif choice == '3' and pfam_family_name:
            threshold = input(f"New threshold value for {pfam_family_name}: ").strip()
            if threshold:
                try:
                    float(threshold)
                    current_dir = os.getcwd()
                    if raise_pfam_threshold(pfam_family_name, current_dir, threshold):
                        print("Re-checking overlaps...")
                        has_overlaps, overlap_content = check_overlap_quiet('.')
                        if not has_overlaps:
                            print("All overlaps resolved!")
                            os.chdir(start_dir)
                            return True
                        else:
                            print(f"Overlaps still present:")
                            lines = overlap_content.strip().split('\n')
                            for line in lines[:10]:
                                print(f"  {line}")
                            if len(lines) > 10:
                                print(f"  ... ({len(lines) - 10} more)")
                except ValueError:
                    print("Invalid threshold")

        elif choice == '4' and clan_to_join:
            print(f"\nJoining clan {clan_to_join}...")
            print("Note: Overlaps with clan members will be allowed.")

            confirm = input(f"Proceed? [y/n]: ").strip().lower()
            if confirm in ['y', 'yes']:
                if add_clan_to_desc('.', clan_to_join):
                    # Re-check overlaps after adding clan
                    print("Re-checking overlaps after adding clan membership...")
                    curation_path = Path('.').resolve()
                    parent_dir = str(curation_path.parent)
                    family_name = curation_path.name
                    run_command(f"pqc-overlap-rdb -no_sigP {family_name} > /dev/null 2>&1",
                               cwd=parent_dir, wait=True)

                    overlap_file = Path('overlap')
                    if overlap_file.exists() and overlap_file.stat().st_size > 0:
                        with open(overlap_file, 'r') as f:
                            content = f.read()
                        print("Non-clan overlaps remaining:")
                        print(content)
                    else:
                        print("All overlaps resolved (clan overlaps now allowed)!")
                        os.chdir(start_dir)
                        return True

        elif choice == '5':
            # Add ED lines to exclude overlapping sequences
            print("\nAdding ED lines to exclude overlapping sequences...")
            print("Running add_ED.pl...")
            success = run_command("add_ED.pl", wait=True)
            if success:
                print("Running pfmake...")
                run_command("pfmake", wait=True)

                print("Re-checking overlaps...")
                has_overlaps, overlap_content = check_overlap_quiet('.')
                if not has_overlaps:
                    print("All overlaps resolved!")
                    os.chdir(start_dir)
                    return True
                else:
                    print(f"Overlaps still present:")
                    lines = overlap_content.strip().split('\n')
                    for line in lines[:10]:
                        print(f"  {line}")
                    if len(lines) > 10:
                        print(f"  ... ({len(lines) - 10} more)")
            else:
                print("add_ED.pl failed", file=sys.stderr)

        elif choice == 's':
            print("Skipping - will use non-overlapping version")
            os.chdir(start_dir)
            return False

        elif choice == 'i':
            # Move to IGNORE
            os.chdir(start_dir)
            root_dir = entry['root']
            ignore_dir = Path("IGNORE")
            if not ignore_dir.exists():
                ignore_dir.mkdir()
            try:
                shutil.move(root_dir, str(ignore_dir / root_dir))
                print(f"Moved {root_dir} to IGNORE/")
            except Exception as e:
                print(f"Warning: Could not move to IGNORE: {e}", file=sys.stderr)
            return None  # Signal that family was ignored entirely

        else:
            print("Unrecognized option")


def curate_family(entry, start_dir, se_prefix=None, use_nano=False, working_dir=None, skip_to_annotation=False):
    """
    Guide curator through family curation.

    Args:
        entry: The family entry dict
        start_dir: Starting directory
        se_prefix: Optional SE line prefix
        use_nano: Use nano instead of emacs
        working_dir: Override directory to work in (for overlap resolution)
        skip_to_annotation: Skip SEED editing and go directly to annotation (after overlap resolution)
    """
    dir_name = working_dir if working_dir else entry['dir']
    root_dir = entry['root']

    print(f"\n{'='*60}")
    print(f"Working on: {dir_name}")
    print(f"Total sequences: {entry['total_seqs']}")
    if 'seed_count' in entry:
        print(f"SEED sequences: {entry['seed_count']}")
    if entry.get('swissprot_count', 0) > 0:
        print(f"SwissProt proteins: {entry['swissprot_count']}")

    # Warn about deeper overlaps that weren't resolvable (only when not using deeper dir)
    if entry.get('deeper_overlaps') and not working_dir:
        print("\nNote: Deeper iterations have overlaps (not resolvable):")
        for warning in entry['deeper_overlaps'][:3]:
            print(f"  - {warning}")
        if len(entry['deeper_overlaps']) > 3:
            print(f"  ... and {len(entry['deeper_overlaps']) - 3} more")

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

    # Auto-ignore single sequence Iterate directories
    if seed_count <= 1 and '/Iterate' in dir_name:
        print(f"\nIterate directory has only {seed_count} sequence(s), auto-ignoring")
        os.chdir(start_dir)
        ignore_dir = Path("IGNORE")
        if not ignore_dir.exists():
            ignore_dir.mkdir()
        try:
            shutil.move(root_dir, str(ignore_dir / root_dir))
            print(f"Moved {root_dir} to IGNORE/")
        except Exception as e:
            print(f"Warning: Could not move to IGNORE: {e}", file=sys.stderr)
        return False

    # Skip SEED editing if coming from successful overlap resolution
    if not skip_to_annotation:
        # Display images if available
        ted_png = Path('ted.png')
        if ted_png.exists():
            print("\nDisplaying ted.png...")
            run_command(f"imgcat {ted_png}", wait=False)

        pfam_png = Path('pfam.png')
        if pfam_png.exists():
            print("\nDisplaying pfam.png...")
            run_command(f"imgcat {pfam_png}", wait=False)

        ted_web_png = Path('ted_web.png')
        if ted_web_png.exists():
            print("\nDisplaying ted_web.png...")
            run_command(f"imgcat {ted_web_png}", wait=False)

        # Launch belvu for SEED review
        print("\nLaunching belvu for SEED alignment review...")
        run_command("belvu SEED", wait=True)

        # SEED editing loop
        while True:
            print("\n" + "="*60)
            print("SEED editing options:")
            print("  y/yes     - Happy with SEED, continue to annotation")
            print("  n/no      - Skip to next family")
            print("  i/ignore  - Move to IGNORE directory")
            print("  b/build   - Run pfbuild on this family")
            print("  w/wholeseq - Run wholeseq on SEED")
            print("  e/extend  - Extend N/C termini")
            print("  p/pad     - Pad alignment ends")
            print("  r/remove  - Remove family to archive")
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
                ignore_dir = Path("IGNORE")
                if not ignore_dir.exists():
                    ignore_dir.mkdir()
                try:
                    shutil.move(root_dir, str(ignore_dir / root_dir))
                    print(f"Moved {root_dir} to IGNORE/")
                except Exception as e:
                    print(f"Warning: Could not move to IGNORE: {e}", file=sys.stderr)
                return False

            elif response in ['b', 'build']:
                print("\nRunning pfbuild -withpfmake...")
                run_command("pfbuild -withpfmake", wait=False)
                print("pfbuild started, moving to next family...")
                os.chdir(start_dir)
                return False

            elif response in ['w', 'wholeseq']:
                print("\nRunning wholeseq on SEED...")
                shutil.copy("SEED", "SEED.beforewholeseq")
                cmd = "wholeseq.pl -align SEED.beforewholeseq -m > SEED.tmp && mv SEED.tmp SEED"
                if run_command(cmd, wait=True):
                    print("Wholeseq completed")
                    run_command("belvu SEED", wait=True)
                else:
                    print("Warning: wholeseq may have failed", file=sys.stderr)

            elif response in ['e', 'extend']:
                n = input("N-terminus extension [0]: ").strip() or "0"
                c = input("C-terminus extension [0]: ").strip() or "0"
                print(f"\nExtending SEED N:{n} C:{c}")
                shutil.copy("SEED", "SEED.beforeextend")
                cmd = f"extend.pl -align SEED.beforeextend -n {n} -c {c} -m > SEED.tmp && mv SEED.tmp SEED"
                if run_command(cmd, wait=True):
                    print("Extension completed")
                    run_command("belvu SEED", wait=True)
                else:
                    print("Warning: extend may have failed", file=sys.stderr)

            elif response in ['p', 'pad']:
                print("\nPadding SEED alignment...")
                shutil.copy("SEED", "SEED.beforepad")
                cmd = "pad_ends.pl -align SEED.beforepad > SEED.tmp && mv SEED.tmp SEED"
                if run_command(cmd, wait=True):
                    print("Padding completed")
                    run_command("belvu SEED", wait=True)
                else:
                    print("Warning: pad_ends may have failed", file=sys.stderr)

            elif response in ['r', 'remove']:
                os.chdir(start_dir)
                removed_dir = Path("REMOVED")
                if not removed_dir.exists():
                    removed_dir.mkdir()
                try:
                    shutil.move(root_dir, str(removed_dir / root_dir))
                    print(f"Moved {root_dir} to REMOVED/")
                except Exception as e:
                    print(f"Warning: Could not move to REMOVED: {e}", file=sys.stderr)
                return False

            else:
                print("Unrecognized option, please try again.")
    else:
        print("\n(Skipping SEED review - overlaps already resolved)")

    # Display information files for annotation
    display_info_files('.')

    # Ask about annotation
    print("\nPlease create annotation using Claude with the information above.")
    while True:
        response = input("Continue to edit DESC? (y/n/i for ignore): ").strip().lower()
        if response in ['y', 'n', 'i']:
            break
        print("Please enter 'y', 'n', or 'i'")

    if response == 'i':
        os.chdir(start_dir)
        ignore_dir = Path("IGNORE")
        if not ignore_dir.exists():
            ignore_dir.mkdir()
        try:
            shutil.move(root_dir, str(ignore_dir / root_dir))
            print(f"Moved {root_dir} to IGNORE/")
        except Exception as e:
            print(f"Warning: Could not move to IGNORE: {e}", file=sys.stderr)
        return False

    if response == 'n':
        os.chdir(start_dir)
        print("Skipping this family")
        return False

    # Open editor for DESC
    editor = "nano" if use_nano else "emacs -nw"
    print(f"\nOpening {editor.split()[0]} for DESC editing...")
    run_command(f"{editor} DESC", wait=True)

    # pfnew loop
    while True:
        response = input("\npfnew this family? (y/n/e=edit DESC/d=DUF/i=ignore): ").strip().lower()

        if response == 'd':
            print("\nRunning nextDUF.pl...")
            run_command("nextDUF.pl", wait=True)
            continue

        if response == 'e':
            run_command(f"{editor} DESC", wait=True)
            continue

        if response == 'i':
            os.chdir(start_dir)
            ignore_dir = Path("IGNORE")
            if not ignore_dir.exists():
                ignore_dir.mkdir()
            try:
                shutil.move(root_dir, str(ignore_dir / root_dir))
                print(f"Moved {root_dir} to IGNORE/")
            except Exception as e:
                print(f"Warning: Could not move to IGNORE: {e}", file=sys.stderr)
            return False

        if response == 'n':
            os.chdir(start_dir)
            print("Skipped pfnew")
            return False

        if response == 'y':
            has_clan = check_desc_for_clan('.')
            os.chdir('..')
            last_dir = dir_name.split('/')[-1]

            if has_clan:
                cmd = f"pfnew {last_dir} -add_to_clan"
            else:
                cmd = f"pfnew {last_dir}"

            print(f"Forking pfnew: {cmd}")

            # Fork the pfnew process to run in background
            pfnew_log = Path(start_dir) / "PFNEWLOG"
            current_cwd = os.getcwd()

            # Capture all paths needed by child BEFORE fork
            # Use absolute paths to avoid issues after parent changes directories
            curation_dir = str(Path(current_cwd) / last_dir)
            root_dir_abs = str(Path(start_dir) / root_dir)
            done_dir_path = str(Path(start_dir) / "DONE")

            pid = os.fork()
            if pid == 0:
                # Child process - run pfmake then pfnew and log output
                try:
                    # Create DONE directory if needed (child does this now)
                    done_dir = Path(done_dir_path)
                    if not done_dir.exists():
                        done_dir.mkdir()

                    with open(pfnew_log, 'a') as log_file:
                        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                        log_file.write(f"\n{'='*60}\n")
                        log_file.write(f"[{timestamp}] Processing: {last_dir}\n")
                        log_file.write(f"Directory: {curation_dir}\n")
                        log_file.write(f"Command: {cmd}\n")
                        log_file.write(f"{'='*60}\n")
                        log_file.flush()

                        # Run pfmake first to fix any TC/NC line issues
                        log_file.write(f"\nRunning pfmake in {curation_dir}...\n")
                        log_file.flush()
                        pfmake_result = subprocess.run(
                            "pfmake", shell=True, cwd=curation_dir,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
                        )
                        log_file.write(pfmake_result.stdout if pfmake_result.stdout else "(no output)\n")
                        if pfmake_result.returncode != 0:
                            log_file.write(f"\n*** WARNING: pfmake returned {pfmake_result.returncode} ***\n")
                        else:
                            log_file.write(f"pfmake completed successfully\n")
                        log_file.flush()

                        # Now run pfnew from the parent directory
                        log_file.write(f"\nRunning: {cmd} (from {current_cwd})\n")
                        log_file.flush()
                        result = subprocess.run(
                            cmd, shell=True, cwd=current_cwd,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
                        )

                        # Write ALL output from pfnew
                        log_file.write("\n--- pfnew output start ---\n")
                        log_file.write(result.stdout if result.stdout else "(no output)\n")
                        log_file.write("--- pfnew output end ---\n")
                        log_file.write(f"\nReturn code: {result.returncode}\n")

                        if result.returncode != 0:
                            log_file.write(f"\n*** WARNING: pfnew FAILED for {last_dir} ***\n")
                            log_file.write(f"*** Directory NOT moved to DONE - please check and fix ***\n")
                        else:
                            log_file.write(f"\nSUCCESS: {last_dir} added to Pfam\n")
                            # Move to DONE on success
                            try:
                                shutil.move(root_dir_abs, str(done_dir / root_dir))
                                log_file.write(f"Moved {root_dir} to DONE/\n")
                            except Exception as move_err:
                                log_file.write(f"Warning: Could not move to DONE: {move_err}\n")

                        log_file.flush()
                except Exception as e:
                    with open(pfnew_log, 'a') as log_file:
                        log_file.write(f"\n*** ERROR running pfmake/pfnew: {e} ***\n")
                        import traceback
                        log_file.write(traceback.format_exc())
                os._exit(0)
            else:
                # Parent process - continue without waiting
                # DO NOT move directory here - child will do it after pfnew completes
                print(f"pfnew started in background (pid {pid})")
                print(f"Output will be logged to: {pfnew_log}")
                print(f"Directory will be moved to DONE/ after pfnew completes successfully")
                os.chdir(start_dir)
                return True

        print("Please enter y, n, e, d, or i")

    os.chdir(start_dir)
    return False


def display_info_files(dir_path):
    """Display information files for annotation"""
    # 1. DESC file first
    desc_file = Path(dir_path) / 'DESC'
    if desc_file.exists():
        print("\n--- DESC ---", file=sys.stderr)
        with open(desc_file, 'r') as f:
            print(f.read(), file=sys.stderr)
        print("--- End DESC ---", file=sys.stderr)

    # 2. species
    species_file = Path(dir_path) / 'species'
    if species_file.exists():
        print("\n--- species summary ---", file=sys.stderr)
        print("The species distribution showing number of matches at each taxonomic rank.", file=sys.stderr)
        print("", file=sys.stderr)
        with open(species_file, 'r') as f:
            print(f.read(), file=sys.stderr)
        print("--- End species ---", file=sys.stderr)

    # 3. SwissProt entries (sp.seq_info) - filtered to reduce token usage
    sp_file = Path(dir_path) / 'sp.seq_info'
    if sp_file.exists() and sp_file.stat().st_size > 0:
        print("\n--- sp.seq_info (SwissProt entries) ---", file=sys.stderr)
        print("SwissProt/UniProt annotations for sequences in this family.", file=sys.stderr)
        print("(Filtered: removed low-value refs and FT sections to reduce tokens)", file=sys.stderr)
        print("", file=sys.stderr)
        with open(sp_file, 'r') as f:
            raw_content = f.read()
        filtered_content = filter_swissprot_content(raw_content)
        print(filtered_content, file=sys.stderr)
        print("--- End sp.seq_info ---", file=sys.stderr)

    # 4. Foldseek (after sp.seq_info, only if matches beyond header)
    fs_file = Path(dir_path) / 'foldseek'
    if fs_file.exists():
        with open(fs_file, 'r') as f:
            lines = f.readlines()
            if len(lines) > 1:
                print("\n--- Foldseek ---", file=sys.stderr)
                print("Structural similarity matches. Can be suggestive of clan lines to add to DESC if results are consistent,", file=sys.stderr)
                print("or adding a CC line noting structural relationship to relevant Pfam families identified.", file=sys.stderr)
                print("", file=sys.stderr)
                for line in lines[:100]:
                    print(line, end='', file=sys.stderr)
                if len(lines) > 100:
                    print(f"... (truncated)", file=sys.stderr)
                print("--- End Foldseek ---", file=sys.stderr)

    # 5. Domain architectures (arch)
    arch_file = Path(dir_path) / 'arch'
    if arch_file.exists():
        print("\n--- Domain architectures ---", file=sys.stderr)
        print("Example domain architectures. Indented lines show domain regions in the protein.", file=sys.stderr)
        print("Lines starting with 'Iterate' or not starting with PFXXXXX are the domain we are considering.", file=sys.stderr)
        print("Note that the protein descriptions in this section are often unreliable.", file=sys.stderr)
        print("", file=sys.stderr)
        with open(arch_file, 'r') as f:
            print(f.read(), file=sys.stderr)
        print("--- End arch ---", file=sys.stderr)

    # 6. PaperBLAST (truncated)
    pb_file = Path(dir_path) / 'paperblast'
    if pb_file.exists() and pb_file.stat().st_size > 0:
        print("\n--- PaperBLAST ---", file=sys.stderr)
        print("Potentially relevant papers to include. Use Pubmed connector to add highly relevant papers to DESC.", file=sys.stderr)
        print("", file=sys.stderr)
        with open(pb_file, 'r') as f:
            lines = f.readlines()
            for line in lines[:100]:
                print(line, end='', file=sys.stderr)
            if len(lines) > 100:
                print(f"... (truncated, {len(lines)} total lines)", file=sys.stderr)
        print("--- End PaperBLAST ---", file=sys.stderr)

    # 7. UniProt bibliography
    bibl_file = Path(dir_path) / 'uniprot_bibl'
    if bibl_file.exists() and bibl_file.stat().st_size > 0:
        print("\n--- UniProt Bibliography ---", file=sys.stderr)
        print("Potentially relevant papers to include. Use Pubmed connector to add highly relevant papers to DESC.", file=sys.stderr)
        print("", file=sys.stderr)
        with open(bibl_file, 'r') as f:
            lines = f.readlines()
            for line in lines[:100]:
                print(line, end='', file=sys.stderr)
            if len(lines) > 100:
                print(f"... (truncated, {len(lines)} total lines)", file=sys.stderr)
        print("--- End UniProt Bibliography ---", file=sys.stderr)

    # 8. TED (truncated)
    ted_file = Path(dir_path) / 'TED'
    if ted_file.exists():
        print("\n--- TED ---", file=sys.stderr)
        print("TED domain architectures. The percentage value at the end is the percentage overlap", file=sys.stderr)
        print("with the domain/family of interest.", file=sys.stderr)
        print("", file=sys.stderr)
        with open(ted_file, 'r') as f:
            lines = f.readlines()
            for line in lines[:100]:
                print(line, end='', file=sys.stderr)
            if len(lines) > 100:
                print(f"... (truncated)", file=sys.stderr)
        print("--- End TED ---", file=sys.stderr)

    # 9. STRING (truncated)
    string_file = Path(dir_path) / 'STRING'
    if string_file.exists() and string_file.stat().st_size > 0:
        print("\n--- STRING ---", file=sys.stderr)
        print("Pfam families predicted to be functionally related by STRING. Mostly not relevant for inclusion.", file=sys.stderr)
        print("", file=sys.stderr)
        with open(string_file, 'r') as f:
            lines = f.readlines()
            for line in lines[:100]:
                print(line, end='', file=sys.stderr)
            if len(lines) > 100:
                print(f"... (truncated)", file=sys.stderr)
        print("--- End STRING ---", file=sys.stderr)


def count_align_sequences(dir_path):
    """Count sequences in ALIGN file for a directory"""
    align_path = Path(dir_path) / 'ALIGN'
    return count_sequences(str(align_path))


def collect_directory_info(dir_name, start_dir):
    """Collect all information for a directory without interactive curation"""

    # Change to start directory
    os.chdir(start_dir)

    if not Path(dir_name).exists():
        return None, f"Directory {dir_name} not found"

    os.chdir(dir_name)

    # Collect all the information in order
    info_sections = []

    # 1. DESC file first
    desc_file = Path('DESC')
    if desc_file.exists():
        info_sections.append("--- DESC ---")
        with open(desc_file, 'r') as f:
            info_sections.append(f.read())
        info_sections.append("--- End DESC ---")

    # 2. species summary
    species_file = Path('species')
    if species_file.exists():
        info_sections.append("\n--- species summary ---")
        info_sections.append("The species distribution showing number of matches at each taxonomic rank.")
        info_sections.append("")
        with open(species_file, 'r') as f:
            info_sections.append(f.read())
        info_sections.append("--- End species ---")

    # 3. SwissProt entries (sp.seq_info) - filtered to reduce token usage
    sp_file = Path('sp.seq_info')
    if sp_file.exists() and sp_file.stat().st_size > 0:
        info_sections.append("\n--- sp.seq_info (SwissProt entries) ---")
        info_sections.append("SwissProt/UniProt annotations for sequences in this family.")
        info_sections.append("(Filtered: removed low-value refs and FT sections to reduce tokens)")
        info_sections.append("")
        with open(sp_file, 'r') as f:
            raw_content = f.read()
        filtered_content = filter_swissprot_content(raw_content)
        info_sections.append(filtered_content)
        info_sections.append("--- End sp.seq_info ---")

    # 4. Foldseek (after sp.seq_info, only if matches beyond header)
    foldseek_file = Path('foldseek')
    if foldseek_file.exists():
        with open(foldseek_file, 'r') as f:
            lines = f.readlines()
            if len(lines) > 1:
                info_sections.append("\n--- Foldseek ---")
                info_sections.append("Structural similarity matches. Can be suggestive of clan lines to add to DESC if results are consistent,")
                info_sections.append("or adding a CC line noting structural relationship to relevant Pfam families identified.")
                info_sections.append("")
                if len(lines) > 100:
                    info_sections.append(''.join(lines[:100]))
                    info_sections.append(f"... (truncated)")
                else:
                    info_sections.append(''.join(lines))
                info_sections.append("--- End Foldseek ---")

    # 5. Domain architectures
    arch_file = Path('arch')
    if arch_file.exists():
        info_sections.append("\n--- Domain architectures ---")
        info_sections.append("Example domain architectures. Indented lines show domain regions in the protein.")
        info_sections.append("Lines starting with 'Iterate' or not starting with PFXXXXX are the domain we are considering.")
        info_sections.append("Note that the protein descriptions in this section are often unreliable.")
        info_sections.append("")
        with open(arch_file, 'r') as f:
            info_sections.append(f.read())
        info_sections.append("--- End arch ---")

    # 6. PaperBLAST literature results
    paperblast_file = Path('paperblast')
    if paperblast_file.exists() and paperblast_file.stat().st_size > 0:
        info_sections.append("\n--- PaperBLAST ---")
        info_sections.append("Potentially relevant papers to include. Use Pubmed connector to add highly relevant papers to DESC.")
        info_sections.append("")
        with open(paperblast_file, 'r') as f:
            lines = f.readlines()
            if len(lines) > 100:
                info_sections.append(''.join(lines[:100]))
                info_sections.append(f"... (truncated, {len(lines)} total lines)")
            else:
                info_sections.append(''.join(lines))
        info_sections.append("--- End PaperBLAST ---")

    # 7. UniProt bibliography
    bibl_file = Path('uniprot_bibl')
    if bibl_file.exists() and bibl_file.stat().st_size > 0:
        info_sections.append("\n--- UniProt Bibliography ---")
        info_sections.append("Potentially relevant papers to include. Use Pubmed connector to add highly relevant papers to DESC.")
        info_sections.append("")
        with open(bibl_file, 'r') as f:
            lines = f.readlines()
            if len(lines) > 100:
                info_sections.append(''.join(lines[:100]))
                info_sections.append(f"... (truncated, {len(lines)} total lines)")
            else:
                info_sections.append(''.join(lines))
        info_sections.append("--- End UniProt Bibliography ---")

    # 8. TED domain information
    ted_file = Path('TED')
    if ted_file.exists():
        info_sections.append("\n--- TED ---")
        info_sections.append("TED domain architectures. The percentage value at the end is the percentage overlap")
        info_sections.append("with the domain/family of interest.")
        info_sections.append("")
        with open(ted_file, 'r') as f:
            lines = f.readlines()
            if len(lines) > 100:
                info_sections.append(''.join(lines[:100]))
                info_sections.append(f"... (truncated)")
            else:
                info_sections.append(''.join(lines))
        info_sections.append("--- End TED ---")

    # 9. STRING protein interaction network
    string_file = Path('STRING')
    if string_file.exists() and string_file.stat().st_size > 0:
        info_sections.append("\n--- STRING ---")
        info_sections.append("Pfam families predicted to be functionally related by STRING. Mostly not relevant for inclusion.")
        info_sections.append("")
        with open(string_file, 'r') as f:
            lines = f.readlines()
            if len(lines) > 100:
                info_sections.append(''.join(lines[:100]))
                info_sections.append(f"... (truncated)")
            else:
                info_sections.append(''.join(lines))
        info_sections.append("--- End STRING ---")

    # Return to start directory
    os.chdir(start_dir)

    return '\n'.join(info_sections), None


def create_manifest_bundle(directories, manifest_entries=None, output_tarball='triage_bundle.tar.gz'):
    """Create a self-contained bundle for batch DESC generation

    Args:
        directories: List of directory paths to include
        manifest_entries: Optional list of dicts with metadata for each directory
                         (type, root, total_seqs, overlaps, etc.)
        output_tarball: Name of output tarball file
    """

    start_dir = os.getcwd()
    bundle_dir = 'pfam_batch'

    # Create bundle directory
    if Path(bundle_dir).exists():
        shutil.rmtree(bundle_dir)
    Path(bundle_dir).mkdir()

    # Build lookup for manifest entries if provided
    entry_lookup = {}
    if manifest_entries:
        for entry in manifest_entries:
            entry_lookup[entry['dir']] = entry

    manifest = {
        'output_dir': 'desc_output',
        'families': []
    }

    print(f"Creating manifest bundle for {len(directories)} directories...", file=sys.stderr)

    for dir_name in directories:
        entry_meta = entry_lookup.get(dir_name, {})
        entry_type = entry_meta.get('type', 'unknown')
        print(f"Processing {dir_name} ({entry_type})...", file=sys.stderr)

        # Collect information
        info_content, error = collect_directory_info(dir_name, start_dir)

        if error:
            print(f"  Error: {error}", file=sys.stderr)
            continue

        # Create family directory in bundle
        family_dir = Path(bundle_dir) / dir_name
        family_dir.mkdir(parents=True, exist_ok=True)

        # Write triage output
        triage_output_file = family_dir / 'triage_output.txt'
        with open(triage_output_file, 'w') as f:
            f.write(info_content)

        # Add to manifest with relative path and metadata
        family_entry = {
            'family_id': dir_name,
            'triage_file': f"{dir_name}/triage_output.txt"
        }

        # Include metadata if available
        if entry_meta:
            family_entry['type'] = entry_meta.get('type', 'unknown')
            family_entry['root'] = entry_meta.get('root', '')
            family_entry['total_seqs'] = entry_meta.get('total_seqs', 0)
            if entry_meta.get('type') == 'deeper_with_overlaps':
                family_entry['overlaps'] = entry_meta.get('overlaps', 0)
                family_entry['overlap_fraction'] = entry_meta.get('overlap_fraction', 0)
                family_entry['resolution_reason'] = entry_meta.get('resolution_reason', '')

        manifest['families'].append(family_entry)

        print(f"  Collected information", file=sys.stderr)

    # Write manifest
    manifest_file = Path(bundle_dir) / 'manifest.json'
    with open(manifest_file, 'w') as f:
        json.dump(manifest, f, indent=2)

    print(f"\nCreated manifest with {len(manifest['families'])} families", file=sys.stderr)

    # Create tarball
    print(f"Creating tarball {output_tarball}...", file=sys.stderr)

    # Use tar to create archive with relative paths
    tar_cmd = f"tar -czf {output_tarball} -C {bundle_dir} ."
    if run_command(tar_cmd, wait=True):
        print(f"Created {output_tarball}", file=sys.stderr)

        # Show tarball contents for verification
        print(f"\nTarball contents:", file=sys.stderr)
        run_command(f"tar -tzf {output_tarball} | head -20", wait=True)

        # Clean up bundle directory
        shutil.rmtree(bundle_dir)
        print(f"\nBundle creation complete!", file=sys.stderr)
        print(f"Transfer {output_tarball} to your laptop and extract with:", file=sys.stderr)
        print(f"  tar -xzf {output_tarball}", file=sys.stderr)
    else:
        print(f"Error creating tarball", file=sys.stderr)


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Pfam Triage Curation Helper v2 - like triage_helper.py with overlap resolution",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script works like triage_helper.py but additionally offers the option to
resolve overlaps in deeper Iterate directories when:
  - Overlaps are < 5% of total sequences
  - Single overlapping family OR all overlapping families in same clan

For each family, if a deeper version with resolvable overlaps exists, you can:
  - Try to resolve overlaps (edit SEED, raise threshold, join clan)
  - Or use the shallower non-overlapping version

Example:
  python triage_helper2.py triage -s TED
        """
    )

    parser.add_argument('triage_file', help='Path to triage file')
    parser.add_argument('-s', '--se-prefix', type=str,
                       help='Update SE line in DESC with this prefix')
    parser.add_argument('-nano', action='store_true',
                       help='Use nano instead of emacs')
    parser.add_argument('-sp', action='store_true',
                       help='Only process families with SwissProt proteins')
    parser.add_argument('--min-seed', type=int, default=1,
                       help='Minimum SEED sequences (default: 1)')
    parser.add_argument('--min-align', type=int, default=0,
                       help='Minimum ALIGN sequences (default: 0, no filtering)')
    parser.add_argument('--max-overlap', type=float, default=0.05,
                       help='Max overlap fraction for resolution (default: 0.05)')
    parser.add_argument('--manifest', action='store_true',
                       help='Create manifest bundle for batch processing instead of interactive curation')

    args = parser.parse_args()

    if not os.path.exists(args.triage_file):
        print(f"Error: Triage file '{args.triage_file}' not found", file=sys.stderr)
        sys.exit(1)

    print(f"Processing {args.triage_file}...", file=sys.stderr)
    if args.sp:
        print("Only processing families with SwissProt proteins", file=sys.stderr)
    if args.min_align > 0:
        print(f"Minimum ALIGN sequences: {args.min_align}", file=sys.stderr)

    best_dirs = parse_triage_file(
        args.triage_file,
        sp_only=args.sp,
        min_seed=args.min_seed,
        min_align=args.min_align,
        max_overlap_fraction=args.max_overlap
    )

    if not best_dirs:
        print("\nNo suitable families found", file=sys.stderr)
        return

    print(f"\nFound {len(best_dirs)} families to potentially curate", file=sys.stderr)

    # Handle manifest mode
    if args.manifest:
        # Include both non-overlapping and deeper resolvable directories
        manifest_entries = []
        for entry in best_dirs:
            # Add the non-overlapping directory
            manifest_entries.append({
                'dir': entry['dir'],
                'type': 'non_overlapping',
                'root': entry['root'],
                'total_seqs': entry['total_seqs']
            })
            # Add any deeper resolvable directories
            for deeper in entry.get('deeper_resolvable', []):
                manifest_entries.append({
                    'dir': deeper['dir'],
                    'type': 'deeper_with_overlaps',
                    'root': entry['root'],
                    'total_seqs': deeper['total_seqs'],
                    'overlaps': deeper['overlaps'],
                    'overlap_fraction': deeper['overlap_fraction'],
                    'resolution_reason': deeper.get('resolution_reason', '')
                })
        directories = [e['dir'] for e in manifest_entries]
        create_manifest_bundle(directories, manifest_entries)
        return

    # Count families with resolvable deeper overlaps
    with_deeper = sum(1 for d in best_dirs if d.get('deeper_resolvable'))
    if with_deeper > 0:
        print(f"  ({with_deeper} have deeper versions with resolvable overlaps)", file=sys.stderr)

    start_dir = os.getcwd()
    families_added = 0

    for i, entry in enumerate(best_dirs, 1):
        print(f"\n[{i}/{len(best_dirs)}]", file=sys.stderr)

        # Check for deeper resolvable overlaps
        deeper = entry.get('deeper_resolvable', [])
        working_dir = None

        if deeper:
            # Sort by total_seqs to get the best deeper option
            deeper.sort(key=lambda x: x['total_seqs'], reverse=True)
            best_deeper = deeper[0]

            print(f"\n{'*'*60}")
            print(f"Family: {entry['root']}")
            print(f"Non-overlapping version: {entry['dir']} ({entry['total_seqs']} seqs)")
            print(f"Deeper version with overlaps: {best_deeper['dir']} ({best_deeper['total_seqs']} seqs)")
            print(f"  - Overlaps: {best_deeper['overlaps']} ({best_deeper['overlap_fraction']*100:.1f}%)")
            print(f"  - {best_deeper['resolution_reason']}")
            print(f"{'*'*60}")

            while True:
                choice = input("\nTry deeper version (d) or use non-overlapping (n)? ").strip().lower()
                if choice in ['d', 'n']:
                    break
                print("Please enter 'd' or 'n'")

            if choice == 'd':
                result = try_resolve_overlaps(best_deeper, start_dir)
                if result is True:
                    # Overlaps resolved, use the deeper directory
                    working_dir = best_deeper['dir']
                    entry['seed_count'] = best_deeper['seed_count']
                    entry['total_seqs'] = best_deeper['total_seqs']
                    # Skip to annotation since overlaps are already resolved
                    if curate_family(entry, start_dir, args.se_prefix, args.nano, working_dir, skip_to_annotation=True):
                        families_added += 1
                    continue
                elif result is None:
                    # Family was moved to IGNORE
                    continue
                else:
                    # Resolution failed, fall back to non-overlapping
                    print(f"\nFalling back to non-overlapping version: {entry['dir']}")
                    working_dir = None

        # Curate the family (normal flow without overlap resolution)
        if curate_family(entry, start_dir, args.se_prefix, args.nano, working_dir):
            families_added += 1

    print(f"\n{'='*60}")
    print(f"Complete! Added {families_added} families to Pfam")
    if families_added > 0:
        print("Processed families have been moved to DONE/")


if __name__ == "__main__":
    main()
