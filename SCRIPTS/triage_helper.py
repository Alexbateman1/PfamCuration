#!/usr/bin/env python3
"""
Pfam Triage Curation Helper Script
Processes triage file to identify best families to curate and guides through annotation
"""

import sys
import os
import subprocess
import shutil
import json
from pathlib import Path
from collections import defaultdict

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
    if not seed_path.exists():
        return 0
    
    count = 0
    with open(seed_path, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#') and not line.startswith('//'):
                count += 1
    return count

def parse_triage_file(triage_file, sp_only=False, min_seed=1):
    """Parse triage file and identify best directories to work on
    
    Args:
        triage_file: Path to triage file
        sp_only: If True, only include families with SwissProt proteins
        min_seed: Minimum number of sequences required in SEED (default 1)
    """
    
    # Group by root directory
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
            
            # Store all info
            entry = {
                'dir': dir_name,
                'total_seqs': total_seqs,
                'overlaps': overlaps,
                'non_overlap_seqs': non_overlap_seqs,
                'overlap_details': parts[14:] if len(parts) > 14 else [],
                'depth': dir_name.count('/Iterate'),
                'root': root,
                'swissprot_count': swissprot_count
            }
            
            root_groups[root].append(entry)
    
    # Find best directory for each root
    best_dirs = []
    
    for root, entries in root_groups.items():
        # Find all non-overlapping directories
        non_overlap = [e for e in entries if e['overlaps'] == 0]
        
        if not non_overlap:
            print(f"[SKIP] {root}: No non-overlapping versions found", file=sys.stderr)
            continue
        
        # Count SEED sequences for each non-overlapping entry
        for entry in non_overlap:
            entry['seed_count'] = count_seed_sequences(entry['dir'])
        
        # Choose the best one using priority:
        # 1. Most sequences in ALIGN (total_seqs)
        # 2. Most sequences in SEED (seed_count) - prefers Iterate dirs with multiple seqs
        # 3. Deeper directories (higher depth) - Iterate versions are often improved
        best = max(non_overlap, key=lambda x: (x['total_seqs'], x['seed_count'], x['depth']))
        
        # Add root to the best entry
        best['root'] = root
        
        # Check if deeper iterations have overlaps
        deeper_overlaps = []
        for e in entries:
            if e['depth'] > best['depth'] and e['overlaps'] > 0:
                overlap_info = ' '.join(e['overlap_details'])
                deeper_overlaps.append(f"{e['dir']}: {overlap_info}")
        
        best['deeper_overlaps'] = deeper_overlaps
        best_dirs.append(best)
    
    # Filter by minimum SEED sequences if requested
    if min_seed > 1:
        before_count = len(best_dirs)
        best_dirs = [entry for entry in best_dirs if entry.get('seed_count', 0) >= min_seed]
        after_count = len(best_dirs)
        if before_count > after_count:
            print(f"[FILTER] Removed {before_count - after_count} families with fewer than {min_seed} SEED sequences", file=sys.stderr)
    
    # Sort by number of sequences (largest first)
    best_dirs.sort(key=lambda x: x['total_seqs'], reverse=True)
    
    # Group domains from the same protein together
    # Extract protein accession from directory names (e.g., A0A073CEA3 from A0A073CEA3_TED01)
    protein_groups = defaultdict(list)
    for entry in best_dirs:
        root = entry['root']
        # Extract protein accession (everything before _TED or other suffix)
        if '_TED' in root:
            protein_acc = root.split('_TED')[0]
        elif '_' in root:
            # Be careful with underscores - some proteins have them in their names
            # Only split on last underscore if it looks like a domain suffix
            parts = root.rsplit('_', 1)
            if len(parts) == 2 and len(parts[1]) <= 10:  # Likely a domain suffix
                protein_acc = parts[0]
            else:
                protein_acc = root
        else:
            protein_acc = root
        
        protein_groups[protein_acc].append(entry)
    
    # Find the maximum size for each protein group (to determine ordering)
    protein_max_size = {}
    for protein_acc, entries in protein_groups.items():
        max_size = max(e['total_seqs'] for e in entries)
        protein_max_size[protein_acc] = max_size
    
    # Sort proteins by their largest domain
    sorted_proteins = sorted(protein_groups.keys(), 
                           key=lambda p: protein_max_size[p], 
                           reverse=True)
    
    # Rebuild the list with grouped proteins
    grouped_dirs = []
    for protein_acc in sorted_proteins:
        # Get all domains for this protein
        protein_domains = protein_groups[protein_acc]
        # Sort domains within the protein by size (largest first) and then by name
        protein_domains.sort(key=lambda x: (-x['total_seqs'], x['root']))
        grouped_dirs.extend(protein_domains)
    
    return grouped_dirs

def check_desc_for_clan(dir_path):
    """Check if DESC file has a CL line"""
    desc_path = Path(dir_path) / 'DESC'
    if desc_path.exists():
        with open(desc_path, 'r') as f:
            for line in f:
                if line.startswith('CL '):
                    return True
    return False

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
                            # Remove all /Iterate occurrences
                            identifier = identifier.replace('/Iterate', '')
                            # Create new SE line with the provided prefix
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

def curate_family(entry, start_dir, se_prefix=None, use_nano=False, skip_to_pfnew=False):
    """Guide curator through family curation"""
    
    dir_name = entry['dir']
    
    # If skip_to_pfnew is True, we're coming back after editing DESC
    if not skip_to_pfnew:
        print(f"\n{'='*60}")
        print(f"Working on: {dir_name}")
        print(f"Total sequences: {entry['total_seqs']}")
        if 'seed_count' in entry:
            print(f"SEED sequences: {entry['seed_count']}")
        if 'swissprot_count' in entry and entry['swissprot_count'] > 0:
            print(f"SwissProt proteins: {entry['swissprot_count']}")
        
        # Warn about deeper overlaps
        if entry['deeper_overlaps']:
            print("\n⚠️  Warning: Deeper iterations have overlaps:")
            for overlap in entry['deeper_overlaps']:
                print(f"  - {overlap}")
        
        # Always start from the original directory
        os.chdir(start_dir)
        
        # Check if directory exists
        if not Path(dir_name).exists():
            print(f"Error: Directory {dir_name} not found, skipping", file=sys.stderr)
            return False
        
        # Change to directory
        os.chdir(dir_name)
    
    if not skip_to_pfnew:
        # Update SE line if prefix provided
        if se_prefix:
            if update_se_line('.', se_prefix):
                print(f"✓ Updated SE line with prefix '{se_prefix}'")
        
        # Count sequences in SEED and ALIGN
        seed_count = count_sequences('SEED')
        align_count = count_sequences('ALIGN')
        
        print(f"\nSequences in SEED: {seed_count}")
        print(f"Sequences in ALIGN: {align_count}")
        
        # Only auto-ignore single sequence SEEDs if this is an Iterate directory
        # Root directories often start with single sequences and get expanded
        if seed_count <= 1 and '/Iterate' in dir_name:
            print(f"\n⚠️  Iterate directory has only {seed_count} sequence(s), auto-ignoring this family", file=sys.stderr)
            os.chdir(start_dir)
            root_dir = dir_name.split('/')[0]
            
            # Create IGNORE directory if needed
            ignore_dir = Path("IGNORE")
            if not ignore_dir.exists():
                ignore_dir.mkdir()
                print("Created IGNORE directory")
            
            # Move root directory to IGNORE
            try:
                shutil.move(root_dir, str(ignore_dir / root_dir))
                print(f"Moved {root_dir} to IGNORE/ (single sequence SEED in Iterate)")
            except Exception as e:
                print(f"⚠ Warning: Could not move {root_dir} to IGNORE: {e}", file=sys.stderr)
            
            return False
    
    # Skip only the initial setup and SEED editing if we're coming back from DESC edit
    if not skip_to_pfnew:
        # Check for ted.png and display it if it exists
        ted_png = Path('ted.png')
        if ted_png.exists():
            print("\nDisplaying ted.png...")
            run_command("display ted.png", wait=False)

        # Launch belvu for SEED review
        print("\nLaunching belvu for SEED alignment review...")
        belvu_success = run_command("belvu SEED", wait=True)
        
        if not belvu_success:
            print("⚠ Warning: belvu may have encountered an issue", file=sys.stderr)
        
        # Present editing options from iterate.py
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
                # Happy with SEED, continue to annotation workflow
                print("\nSEED looks good, continuing to annotation...")
                break
                
            elif response in ['n', 'no']:
                # Skip to next family
                os.chdir(start_dir)
                print("Skipping this family")
                return False
                
            elif response in ['i', 'ignore']:
                # Move to IGNORE directory
                os.chdir(start_dir)
                root_dir = dir_name.split('/')[0]
                
                # Create IGNORE directory if needed
                ignore_dir = Path("IGNORE")
                if not ignore_dir.exists():
                    ignore_dir.mkdir()
                    print("Created IGNORE directory")
                
                # Move root directory to IGNORE
                try:
                    shutil.move(root_dir, str(ignore_dir / root_dir))
                    print(f"Moved {root_dir} to IGNORE/")
                except Exception as e:
                    print(f"⚠ Warning: Could not move {root_dir} to IGNORE: {e}", file=sys.stderr)
                
                return False
                
            elif response in ['b', 'build']:
                # Run pfbuild with -withpfmake flag
                print("\nRunning pfbuild -withpfmake...")
                pfbuild_success = run_command("pfbuild -withpfmake", wait=False)
                if pfbuild_success:
                    print("✓ pfbuild started successfully")
                    print("Moving to next family while pfbuild runs...")
                else:
                    print("⚠ Warning: pfbuild may have failed", file=sys.stderr)
                
                # Return to start directory and skip to next family
                os.chdir(start_dir)
                return False
                
            elif response in ['w', 'wholeseq']:
                print("\nRunning wholeseq on SEED...")
                shutil.copy("SEED", "SEED.beforewholeseq")
                wholeseq_cmd = "wholeseq.pl -align SEED.beforewholeseq -m > SEED.tmp && mv SEED.tmp SEED"
                if run_command(wholeseq_cmd, wait=True):
                    print("✓ Wholeseq completed")
                    # Relaunch belvu to show updated SEED
                    print("Relaunching belvu with updated SEED...")
                    run_command("belvu SEED", wait=True)
                else:
                    print("⚠ Warning: wholeseq may have failed", file=sys.stderr)
                    
            elif response in ['e', 'extend']:
                n = input("How many residues to extend N-terminus? [default 0]: ").strip() or "0"
                c = input("How many residues to extend C-terminus? [default 0]: ").strip() or "0"
                
                print(f"\nExtending SEED N:{n} C:{c}")
                shutil.copy("SEED", "SEED.beforeextend")
                extend_cmd = f"extend.pl -align SEED.beforeextend -n {n} -c {c} -m > SEED.tmp && mv SEED.tmp SEED"
                if run_command(extend_cmd, wait=True):
                    print("✓ Extension completed")
                    # Relaunch belvu to show updated SEED
                    print("Relaunching belvu with updated SEED...")
                    run_command("belvu SEED", wait=True)
                else:
                    print("⚠ Warning: extend may have failed", file=sys.stderr)
                    
            elif response in ['p', 'pad']:
                print("\nPadding SEED alignment...")
                shutil.copy("SEED", "SEED.beforepad")
                pad_cmd = "pad_ends.pl -align SEED.beforepad > SEED.tmp && mv SEED.tmp SEED"
                if run_command(pad_cmd, wait=True):
                    print("✓ Padding completed")
                    # Relaunch belvu to show updated SEED
                    print("Relaunching belvu with updated SEED...")
                    run_command("belvu SEED", wait=True)
                else:
                    print("⚠ Warning: pad_ends may have failed", file=sys.stderr)
                    
            elif response in ['r', 'remove']:
                print("\nRemoving family to archive...")
                os.chdir(start_dir)
                root_dir = dir_name.split('/')[0]
                
                # Create REMOVED directory if needed (similar to archive)
                removed_dir = Path("REMOVED")
                if not removed_dir.exists():
                    removed_dir.mkdir()
                    print("Created REMOVED directory")
                
                # Move root directory to REMOVED
                try:
                    shutil.move(root_dir, str(removed_dir / root_dir))
                    print(f"Moved {root_dir} to REMOVED/")
                except Exception as e:
                    print(f"⚠ Warning: Could not move {root_dir} to REMOVED: {e}", file=sys.stderr)
                
                return False
                
            else:
                print("Unrecognized option, please try again.")
        
        # Now continue with annotation workflow
        
        # Check for sp.seq_info file and output to STDERR
        sp_file = Path('sp.seq_info')
        if sp_file.exists():
            print("\n--- sp.seq_info content (copy to Claude) ---", file=sys.stderr)
            with open(sp_file, 'r') as f:
                print(f.read(), file=sys.stderr)
            print("--- End of sp.seq_info ---", file=sys.stderr)
        else:
            print(f"\nNote: sp.seq_info file not found", file=sys.stderr)
            # Print current path which often contains useful protein accession
            current_path = os.getcwd()
            print(f"Current directory: {current_path}", file=sys.stderr)
            # Extract potential protein accession from path
            path_parts = current_path.split('/')
            if path_parts:
                dir_name_local = path_parts[-1]
                # Check if it looks like a protein accession (e.g., A0A073CEA3_TED01)
                if '_TED' in dir_name_local:
                    protein_acc = dir_name_local.split('_TED')[0]
                    print(f"Potential protein accession for UniProt search: {protein_acc}", file=sys.stderr)
        
        # Check for species file and output to STDERR
        species_file = Path('species')
        print("\n--- species summary (copy to Claude) ---", file=sys.stderr)
        if species_file.exists():
            with open(species_file, 'r') as f:
                print(f.read(), file=sys.stderr)
        else:
            print("Note: species file not found", file=sys.stderr)
        print("--- End of species summary ---", file=sys.stderr)
        
        # Check for arch file and output to STDERR
        arch_file = Path('arch')
        print("\n--- Domain architectures (copy to Claude) ---", file=sys.stderr)
        if arch_file.exists():
            with open(arch_file, 'r') as f:
                print(f.read(), file=sys.stderr)
        else:
            print("Note: arch file not found", file=sys.stderr)
        print("--- End of Domain architectures ---", file=sys.stderr)
        
        # Check for paperblast file and output to STDERR if it has content
        paperblast_file = Path('paperblast')
        if paperblast_file.exists() and paperblast_file.stat().st_size > 0:
            print("\n--- PaperBLAST literature results (copy to Claude) ---", file=sys.stderr)
            with open(paperblast_file, 'r') as f:
                lines = f.readlines()
                # Truncate if too long (similar to TED file)
                if len(lines) > 100:
                    print(''.join(lines[:100]), file=sys.stderr)
                    print(f"... (truncated - showing first 100 of {len(lines)} lines)", file=sys.stderr)
                else:
                    print(''.join(lines), file=sys.stderr)
            print("--- End of PaperBLAST literature results ---", file=sys.stderr)
        
        # Check for TED file and output to STDERR (truncated to 100 lines)
        ted_file = Path('TED')
        print("\n--- TED domain information (copy to Claude) ---", file=sys.stderr)
        if ted_file.exists():
            with open(ted_file, 'r') as f:
                lines = f.readlines()
                if len(lines) > 100:
                    print(''.join(lines[:100]), file=sys.stderr)
                    print(f"... (truncated - showing first 100 of {len(lines)} lines)", file=sys.stderr)
                else:
                    print(''.join(lines), file=sys.stderr)
        else:
            print("Note: TED file not found", file=sys.stderr)
        print("--- End of TED domain information ---", file=sys.stderr)

        # Check for STRING file and output to STDERR
        string_file = Path('STRING')
        if string_file.exists() and string_file.stat().st_size > 0:
            print("\n--- STRING protein interaction network (copy to Claude) ---", file=sys.stderr)
            with open(string_file, 'r') as f:
                lines = f.readlines()
                # Truncate if too long
                if len(lines) > 100:
                    print(''.join(lines[:100]), file=sys.stderr)
                    print(f"... (truncated - showing first 100 of {len(lines)} lines)", file=sys.stderr)
                else:
                    print(''.join(lines), file=sys.stderr)
            print("--- End of STRING protein interaction network ---", file=sys.stderr)

        # Check for foldseek file and output to STDERR if it has more than header line
        foldseek_file = Path('foldseek')
        if foldseek_file.exists():
            with open(foldseek_file, 'r') as f:
                lines = f.readlines()
                # Only output if more than header line
                if len(lines) > 1:
                    print("\n--- Foldseek structural matches (copy to Claude) ---", file=sys.stderr)
                    if len(lines) > 100:
                        print(''.join(lines[:100]), file=sys.stderr)
                        print(f"... (truncated - showing first 100 of {len(lines)} lines)", file=sys.stderr)
                    else:
                        print(''.join(lines), file=sys.stderr)
                    print("--- End of Foldseek structural matches ---", file=sys.stderr)
        
        # Check for DESC file and output to STDERR
        desc_file = Path('DESC')
        print("\n--- Current DESC file (copy to Claude) ---", file=sys.stderr)
        if desc_file.exists():
            with open(desc_file, 'r') as f:
                print(f.read(), file=sys.stderr)
        else:
            print("Note: DESC file not found", file=sys.stderr)
        print("--- End of Current DESC file ---", file=sys.stderr)
        
        # Ask if wish to continue after annotation
        print("\nPlease create annotation using Claude with the sp.seq_info above.")
        while True:
            response = input("Do you wish to continue after annotation? (y/n/i for ignore): ").strip().lower()
            if response in ['y', 'n', 'i']:
                break
            print("Please enter 'y', 'n', or 'i'")
        
        if response == 'i':
            # Move to IGNORE directory
            os.chdir(start_dir)
            root_dir = dir_name.split('/')[0]
            
            # Create IGNORE directory if needed
            ignore_dir = Path("IGNORE")
            if not ignore_dir.exists():
                ignore_dir.mkdir()
                print("Created IGNORE directory")
            
            # Move root directory to IGNORE
            try:
                shutil.move(root_dir, str(ignore_dir / root_dir))
                print(f"Moved {root_dir} to IGNORE/")
            except Exception as e:
                print(f"⚠ Warning: Could not move {root_dir} to IGNORE: {e}", file=sys.stderr)
            
            return False
        
        if response == 'n':
            # Return to start directory
            os.chdir(start_dir)
            print("Skipping this family after annotation")
            return False
        
        # Open editor for DESC editing
        editor = "nano" if use_nano else "emacs -nw"
        print(f"\nOpening {editor.split()[0]} for DESC file editing...")
        editor_cmd = f"{editor} DESC"
        editor_success = run_command(editor_cmd, wait=True)
        
        if not editor_success:
            print(f"⚠ Warning: {editor.split()[0]} may have encountered an issue", file=sys.stderr)
    
    # Ask if wish to pfnew
    while True:
        response = input("\nDo you wish to pfnew this family? (y/n/e to edit DESC again/d for DUF/i for ignore): ").strip().lower()
        if response in ['y', 'n', 'e', 'd', 'i']:
            break
        print("Please enter 'y', 'n', 'e', 'd', or 'i'")
    
    if response == 'd':
        # Run nextDUF.pl to add DUF number
        print("\nRunning nextDUF.pl to assign DUF number...")
        nextduf_success = run_command("nextDUF.pl", wait=True)
        
        if nextduf_success:
            print("✓ DUF number assigned successfully")
        else:
            print("⚠ Warning: nextDUF.pl may have encountered an issue", file=sys.stderr)
        
        # Return to parent directory (we're already in the curation directory)
        # No need to cd as we're already in the right place
        
        # Loop back to ask the pfnew question again
        return curate_family(entry, start_dir, se_prefix, use_nano, skip_to_pfnew=True)
    
    if response == 'e':
        # Edit DESC file again
        editor = "nano" if use_nano else "emacs -nw"
        print(f"\nOpening {editor.split()[0]} for DESC file editing...")
        editor_cmd = f"{editor} DESC"
        editor_success = run_command(editor_cmd, wait=True)
        
        if not editor_success:
            print(f"⚠ Warning: {editor.split()[0]} may have encountered an issue", file=sys.stderr)
        
        # Loop back to ask the pfnew question again
        return curate_family(entry, start_dir, se_prefix, use_nano, skip_to_pfnew=True)
    
    if response == 'i':
        # Move to IGNORE directory
        os.chdir(start_dir)
        root_dir = dir_name.split('/')[0]
        
        # Create IGNORE directory if needed
        ignore_dir = Path("IGNORE")
        if not ignore_dir.exists():
            ignore_dir.mkdir()
            print("Created IGNORE directory")
        
        # Move root directory to IGNORE
        try:
            shutil.move(root_dir, str(ignore_dir / root_dir))
            print(f"Moved {root_dir} to IGNORE/")
        except Exception as e:
            print(f"⚠ Warning: Could not move {root_dir} to IGNORE: {e}", file=sys.stderr)
        
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
            print(f"✓ Successfully added {last_dir} to Pfam")
            
            # Get root directory name
            root_dir = dir_name.split('/')[0]
            
            # Navigate back to start directory to move the root
            os.chdir(start_dir)
            
            # Create DONE directory if needed
            done_dir = Path("DONE")
            if not done_dir.exists():
                done_dir.mkdir()
                print("Created DONE directory")
            
            # Move root directory to DONE
            try:
                shutil.move(root_dir, str(done_dir / root_dir))
                print(f"Moved {root_dir} to DONE/")
            except Exception as e:
                print(f"⚠ Warning: Could not move {root_dir} to DONE: {e}", file=sys.stderr)
            
            return True
        else:
            print(f"\n⚠ ERROR: pfnew failed for {last_dir}", file=sys.stderr)
            print("Please review the error above.", file=sys.stderr)
            
            # Ask what to do with the failed family
            while True:
                response = input("\nWhat would you like to do? (c=continue/o=move to OVERLAP/i=move to IGNORE): ").strip().lower()
                if response in ['c', 'o', 'i']:
                    break
                print("Please enter 'c', 'o', or 'i'")
            
            # Get root directory name
            root_dir = dir_name.split('/')[0]
            
            # Navigate back to start directory
            os.chdir(start_dir)
            
            if response == 'o':
                # Move to OVERLAP directory
                overlap_dir = Path("OVERLAP")
                if not overlap_dir.exists():
                    overlap_dir.mkdir()
                    print("Created OVERLAP directory")
                
                try:
                    shutil.move(root_dir, str(overlap_dir / root_dir))
                    print(f"Moved {root_dir} to OVERLAP/")
                except Exception as e:
                    print(f"⚠ Warning: Could not move {root_dir} to OVERLAP: {e}", file=sys.stderr)
            
            elif response == 'i':
                # Move to IGNORE directory
                ignore_dir = Path("IGNORE")
                if not ignore_dir.exists():
                    ignore_dir.mkdir()
                    print("Created IGNORE directory")
                
                try:
                    shutil.move(root_dir, str(ignore_dir / root_dir))
                    print(f"Moved {root_dir} to IGNORE/")
                except Exception as e:
                    print(f"⚠ Warning: Could not move {root_dir} to IGNORE: {e}", file=sys.stderr)
            
            # For all options (c, o, i), continue to next family
            return False
    else:
        # Return to start directory
        os.chdir(start_dir)
        print("Skipped pfnew for this family")
        return False

def collect_directory_info(dir_name, start_dir):
    """Collect all information for a directory without interactive curation"""
    
    # Change to start directory
    os.chdir(start_dir)
    
    if not Path(dir_name).exists():
        return None, f"Directory {dir_name} not found"
    
    os.chdir(dir_name)
    
    # Collect all the information
    info_sections = []
    
    # sp.seq_info content
    sp_file = Path('sp.seq_info')
    if sp_file.exists():
        info_sections.append("--- sp.seq_info content (copy to Claude) ---")
        with open(sp_file, 'r') as f:
            info_sections.append(f.read())
        info_sections.append("--- End of sp.seq_info ---")
    else:
        info_sections.append("--- sp.seq_info content (copy to Claude) ---")
        info_sections.append("Note: sp.seq_info file not found")
        # Extract potential protein accession
        current_path = os.getcwd()
        info_sections.append(f"Current directory: {current_path}")
        if '_TED' in dir_name:
            protein_acc = dir_name.split('_TED')[0]
            info_sections.append(f"Potential protein accession for UniProt search: {protein_acc}")
        info_sections.append("--- End of sp.seq_info ---")
    
    # species summary
    info_sections.append("\n--- species summary (copy to Claude) ---")
    species_file = Path('species')
    if species_file.exists():
        with open(species_file, 'r') as f:
            info_sections.append(f.read())
    else:
        info_sections.append("Note: species file not found")
    info_sections.append("--- End of species summary ---")
    
    # Domain architectures
    info_sections.append("\n--- Domain architectures (copy to Claude) ---")
    arch_file = Path('arch')
    if arch_file.exists():
        with open(arch_file, 'r') as f:
            info_sections.append(f.read())
    else:
        info_sections.append("Note: arch file not found")
    info_sections.append("--- End of Domain architectures ---")
    
    # PaperBLAST literature results
    paperblast_file = Path('paperblast')
    if paperblast_file.exists() and paperblast_file.stat().st_size > 0:
        info_sections.append("\n--- PaperBLAST literature results (copy to Claude) ---")
        with open(paperblast_file, 'r') as f:
            lines = f.readlines()
            if len(lines) > 100:
                info_sections.append(''.join(lines[:100]))
                info_sections.append(f"... (truncated - showing first 100 of {len(lines)} lines)")
            else:
                info_sections.append(''.join(lines))
        info_sections.append("--- End of PaperBLAST literature results ---")
    
    # TED domain information
    info_sections.append("\n--- TED domain information (copy to Claude) ---")
    ted_file = Path('TED')
    if ted_file.exists():
        with open(ted_file, 'r') as f:
            lines = f.readlines()
            if len(lines) > 100:
                info_sections.append(''.join(lines[:100]))
                info_sections.append(f"... (truncated - showing first 100 of {len(lines)} lines)")
            else:
                info_sections.append(''.join(lines))
    else:
        info_sections.append("Note: TED file not found")
    info_sections.append("--- End of TED domain information ---")

    # STRING protein interaction network
    string_file = Path('STRING')
    if string_file.exists() and string_file.stat().st_size > 0:
        info_sections.append("\n--- STRING protein interaction network (copy to Claude) ---")
        with open(string_file, 'r') as f:
            lines = f.readlines()
            if len(lines) > 100:
                info_sections.append(''.join(lines[:100]))
                info_sections.append(f"... (truncated - showing first 100 of {len(lines)} lines)")
            else:
                info_sections.append(''.join(lines))
        info_sections.append("--- End of STRING protein interaction network ---")

    # Foldseek results
    foldseek_file = Path('foldseek')
    if foldseek_file.exists():
        with open(foldseek_file, 'r') as f:
            lines = f.readlines()
            # Only include if more than header line
            if len(lines) > 1:
                info_sections.append("\n--- Foldseek structural matches (copy to Claude) ---")
                if len(lines) > 100:
                    info_sections.append(''.join(lines[:100]))
                    info_sections.append(f"... (truncated - showing first 100 of {len(lines)} lines)")
                else:
                    info_sections.append(''.join(lines))
                info_sections.append("--- End of Foldseek structural matches ---")
    
    # Current DESC file
    info_sections.append("\n--- Current DESC file (copy to Claude) ---")
    desc_file = Path('DESC')
    if desc_file.exists():
        with open(desc_file, 'r') as f:
            info_sections.append(f.read())
    else:
        info_sections.append("Note: DESC file not found")
    info_sections.append("--- End of Current DESC file ---")
    
    # Return to start directory
    os.chdir(start_dir)
    
    return '\n'.join(info_sections), None


def create_manifest_bundle(directories, output_tarball='triage_bundle.tar.gz'):
    """Create a self-contained bundle for batch DESC generation"""
    
    start_dir = os.getcwd()
    bundle_dir = 'pfam_batch'
    
    # Create bundle directory
    if Path(bundle_dir).exists():
        shutil.rmtree(bundle_dir)
    Path(bundle_dir).mkdir()
    
    manifest = {
        'output_dir': 'desc_output',
        'families': []
    }
    
    print(f"Creating manifest bundle for {len(directories)} families...", file=sys.stderr)
    
    for dir_name in directories:
        print(f"Processing {dir_name}...", file=sys.stderr)
        
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
        
        # Add to manifest with relative path
        manifest['families'].append({
            'family_id': dir_name,
            'triage_file': f"{dir_name}/triage_output.txt"
        })
        
        print(f"  ✓ Collected information", file=sys.stderr)
    
    # Write manifest
    manifest_file = Path(bundle_dir) / 'manifest.json'
    with open(manifest_file, 'w') as f:
        json.dump(manifest, f, indent=2)
    
    print(f"\n✓ Created manifest with {len(manifest['families'])} families", file=sys.stderr)
    
    # Create tarball
    print(f"Creating tarball {output_tarball}...", file=sys.stderr)
    
    # Use tar to create archive with relative paths
    # -C changes to directory before archiving
    # . archives everything in that directory
    tar_cmd = f"tar -czf {output_tarball} -C {bundle_dir} ."
    if run_command(tar_cmd, wait=True):
        print(f"✓ Created {output_tarball}", file=sys.stderr)
        
        # Show tarball contents for verification
        print(f"\nTarball contents:", file=sys.stderr)
        run_command(f"tar -tzf {output_tarball} | head -20", wait=True)
        
        # Clean up bundle directory
        shutil.rmtree(bundle_dir)
        print(f"\n✓ Bundle creation complete!", file=sys.stderr)
        print(f"Transfer {output_tarball} to your laptop and extract with:", file=sys.stderr)
        print(f"  tar -xzf {output_tarball}", file=sys.stderr)
    else:
        print(f"⚠ Error creating tarball", file=sys.stderr)


def main():
    # Check for --create-manifest mode
    if '--create-manifest' in sys.argv:
        # In manifest mode, first argument should be triage file
        if len(sys.argv) < 3:
            print("Error: --create-manifest requires a triage file", file=sys.stderr)
            print("Usage: python triage_helper.py --create-manifest <triage_file> [-sp] [--min-seed N]", file=sys.stderr)
            sys.exit(1)
        
        triage_file = sys.argv[2]  # First arg after --create-manifest
        
        if not os.path.exists(triage_file):
            print(f"Error: Triage file '{triage_file}' not found", file=sys.stderr)
            sys.exit(1)
        
        # Check for -sp option in manifest mode
        sp_only = '-sp' in sys.argv
        if sp_only:
            print("Only including families with SwissProt proteins", file=sys.stderr)
        
        # Check for --min-seed option
        min_seed = 1
        if '--min-seed' in sys.argv:
            try:
                seed_index = sys.argv.index('--min-seed')
                if seed_index + 1 < len(sys.argv):
                    min_seed = int(sys.argv[seed_index + 1])
                    print(f"Minimum SEED sequences: {min_seed}", file=sys.stderr)
                else:
                    print("Error: --min-seed requires a number", file=sys.stderr)
                    sys.exit(1)
            except ValueError:
                print("Error: --min-seed requires a valid number", file=sys.stderr)
                sys.exit(1)
        
        print(f"Processing {triage_file} to identify best directories...", file=sys.stderr)
        
        # Use same selection logic as normal mode
        best_dirs = parse_triage_file(triage_file, sp_only, min_seed)
        
        if not best_dirs:
            print("\nNo suitable families found for manifest creation", file=sys.stderr)
            return
        
        print(f"Found {len(best_dirs)} families to include in manifest", file=sys.stderr)
        
        # Extract directory names from best_dirs
        directories = [entry['dir'] for entry in best_dirs]
        
        create_manifest_bundle(directories)
        return
    
    if len(sys.argv) < 2:
        print("Usage: python triage_helper.py <triage_file> [-s SE_PREFIX] [-nano] [-sp] [--min-seed N]", file=sys.stderr)
        print("   OR: python triage_helper.py --create-manifest <triage_file> [-sp] [--min-seed N]", file=sys.stderr)
        print("", file=sys.stderr)
        print("Normal mode options:", file=sys.stderr)
        print("  -s SE_PREFIX  Optional: Update SE line in DESC with this prefix", file=sys.stderr)
        print("  -nano         Optional: Use nano instead of emacs for editing", file=sys.stderr)
        print("  -sp           Optional: Only process families with SwissProt proteins", file=sys.stderr)
        print("  --min-seed N  Optional: Only process families with at least N sequences in SEED (default: 1)", file=sys.stderr)
        print("", file=sys.stderr)
        print("Manifest mode:", file=sys.stderr)
        print("  --create-manifest <triage_file> [-sp] [--min-seed N]  Create self-contained bundle for best families", file=sys.stderr)
        print("                                                        Uses same selection logic as normal mode", file=sys.stderr)
        sys.exit(1)
    
    triage_file = sys.argv[1]
    
    # Check for -s option
    se_prefix = None
    if '-s' in sys.argv:
        s_index = sys.argv.index('-s')
        if s_index + 1 < len(sys.argv) and not sys.argv[s_index + 1].startswith('-'):
            se_prefix = sys.argv[s_index + 1]
            print(f"Will update SE lines with prefix: {se_prefix}", file=sys.stderr)
        else:
            print("Error: -s option requires a prefix string", file=sys.stderr)
            sys.exit(1)
    
    # Check for -nano option
    use_nano = '-nano' in sys.argv
    if use_nano:
        print("Using nano for editing DESC files", file=sys.stderr)
    
    # Check for -sp option
    sp_only = '-sp' in sys.argv
    if sp_only:
        print("Only processing families with SwissProt proteins", file=sys.stderr)
    
    # Check for --min-seed option
    min_seed = 1
    if '--min-seed' in sys.argv:
        try:
            seed_index = sys.argv.index('--min-seed')
            if seed_index + 1 < len(sys.argv):
                min_seed = int(sys.argv[seed_index + 1])
                print(f"Minimum SEED sequences: {min_seed}", file=sys.stderr)
            else:
                print("Error: --min-seed requires a number", file=sys.stderr)
                sys.exit(1)
        except ValueError:
            print("Error: --min-seed requires a valid number", file=sys.stderr)
            sys.exit(1)
    
    if not os.path.exists(triage_file):
        print(f"Error: Triage file '{triage_file}' not found", file=sys.stderr)
        sys.exit(1)
    
    print(f"Processing {triage_file}...", file=sys.stderr)
    
    # Get list of families to review
    best_dirs = parse_triage_file(triage_file, sp_only, min_seed)
    
    if not best_dirs:
        print("\nNo suitable families found for curation", file=sys.stderr)
        return
    
    print(f"\nFound {len(best_dirs)} families to potentially curate", file=sys.stderr)
    
    if sp_only:
        total_sp = sum(e['swissprot_count'] for e in best_dirs)
        print(f"Total SwissProt proteins in selected families: {total_sp}", file=sys.stderr)
    
    # Report on protein grouping
    protein_counts = defaultdict(int)
    for entry in best_dirs:
        root = entry['root']
        if '_TED' in root:
            protein_acc = root.split('_TED')[0]
        elif '_' in root:
            protein_acc = root.split('_')[0]
        else:
            protein_acc = root
        protein_counts[protein_acc] += 1
    
    multi_domain_proteins = {p: c for p, c in protein_counts.items() if c > 1}
    if multi_domain_proteins:
        print(f"Found {len(multi_domain_proteins)} proteins with multiple domains:", file=sys.stderr)
        for protein, count in sorted(multi_domain_proteins.items(), key=lambda x: -x[1])[:5]:
            print(f"  {protein}: {count} domains", file=sys.stderr)
        if len(multi_domain_proteins) > 5:
            print(f"  ... and {len(multi_domain_proteins) - 5} more", file=sys.stderr)
    
    # Store starting directory
    start_dir = os.getcwd()
    
    # Process each family
    families_added = 0
    for i, entry in enumerate(best_dirs, 1):
        print(f"\n[{i}/{len(best_dirs)}]", file=sys.stderr)
        if curate_family(entry, start_dir, se_prefix, use_nano):
            families_added += 1
    
    print(f"\n{'='*60}")
    print(f"Complete! Added {families_added} families to Pfam")
    if families_added > 0:
        print(f"Processed families have been moved to DONE/")

if __name__ == "__main__":
    main()
