#!/usr/bin/env python3
"""
foldmason_realign.py - Realign Pfam SEED using FoldMason structure-based alignment

This script takes a Pfam SEED alignment (mul format) and creates a structure-based
realignment using FoldMason with AlphaFold models.

Usage:
    foldmason_realign.py [curation_dir] [options]

Example:
    foldmason_realign.py /path/to/PF00001
    foldmason_realign.py .  # Use current directory
    foldmason_realign.py . --par  # Pad and realign (two-pass)
"""

import argparse
import os
import re
import shutil
import sys
import tempfile
import subprocess
import urllib.request
from datetime import datetime
from pathlib import Path


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Realign Pfam SEED using FoldMason structure-based alignment',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    %(prog)s /path/to/PF00001
    %(prog)s .  --refine-iters 3
    %(prog)s . --keep-temp --plddt-warn 70
    %(prog)s . --par  # Two-pass: align, pad ends, realign
        """
    )
    parser.add_argument(
        'curation_dir',
        nargs='?',
        default='.',
        help='Path to Pfam curation directory containing SEED file (default: current dir)'
    )
    parser.add_argument(
        '--output',
        default=None,
        help='Output alignment file (default: SEED.foldmason in curation_dir)'
    )
    parser.add_argument(
        '--foldmason',
        default='foldmason',
        help='Path to foldmason executable (default: foldmason in PATH)'
    )
    parser.add_argument(
        '--keep-temp',
        action='store_true',
        help="Don't delete temporary files (for debugging)"
    )
    parser.add_argument(
        '--refine-iters',
        type=int,
        default=0,
        help='Number of FoldMason refinement iterations (default: 0)'
    )
    parser.add_argument(
        '--plddt-warn',
        type=float,
        default=50.0,
        help='Warn if mean pLDDT below this threshold (default: 50)'
    )
    parser.add_argument(
        '--par',
        action='store_true',
        help='Pad and realign: run FoldMason, pad ends with padd_ends.pl, then realign'
    )

    return parser.parse_args()


def parse_seed_alignment(seed_file):
    """
    Parse a Pfam SEED alignment file in mul format.

    Args:
        seed_file: Path to SEED file

    Returns:
        list of dicts with keys: accession, version_accession, start, end, sequence
    """
    sequences = []

    with open(seed_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.rstrip()
            if not line:
                continue

            # Match lines like: A0A0D4CIU0.1/281-423    PAILILTNSDQIM...
            # Accession may or may not have version number
            match = re.match(r'^(\S+)/(\d+)-(\d+)\s+(.+)$', line)
            if match:
                full_acc = match.group(1)  # e.g., A0A0D4CIU0.1
                start = int(match.group(2))
                end = int(match.group(3))
                sequence = match.group(4)

                # Strip version number for API calls
                base_acc = full_acc.split('.')[0]

                sequences.append({
                    'full_name': f"{full_acc}/{start}-{end}",
                    'version_accession': full_acc,
                    'base_accession': base_acc,
                    'start': start,
                    'end': end,
                    'sequence': sequence
                })
            else:
                print(f"Warning: Could not parse line {line_num}: {line[:50]}...",
                      file=sys.stderr)

    return sequences


def download_alphafold_model(uniprot_acc, output_dir):
    """
    Download AlphaFold model CIF file for a given UniProt accession.
    Tries v6 first, falls back to v4 if not found.

    Args:
        uniprot_acc: UniProt accession (without version number)
        output_dir: Directory to save the file

    Returns:
        tuple: (path to downloaded CIF file, error message or None)
    """
    # Try versions in order of preference (newest first)
    versions = ['v6', 'v4']
    last_error = None

    for version in versions:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_acc}-F1-model_{version}.cif"
        output_file = os.path.join(output_dir, f"AF-{uniprot_acc}-F1-model_{version}.cif")

        # Check if already downloaded
        if os.path.exists(output_file):
            return output_file, None

        try:
            urllib.request.urlretrieve(url, output_file)
            return output_file, None
        except urllib.error.HTTPError as e:
            last_error = f"HTTP {e.code}: {e.reason}"
            # Continue to try next version
        except urllib.error.URLError as e:
            last_error = f"URL Error: {e.reason}"
        except Exception as e:
            last_error = str(e)

    return None, last_error


def chop_structure_to_domain(cif_file, start, end, output_pdb):
    """
    Extract a domain region from a CIF file and save as PDB.

    Args:
        cif_file: Path to input CIF file
        start: Start residue position (1-indexed, UniProt coordinates)
        end: End residue position (1-indexed, UniProt coordinates)
        output_pdb: Path to output PDB file

    Returns:
        tuple: (success boolean, error message or None)
    """
    try:
        from Bio.PDB import MMCIFParser, PDBIO, Select
    except ImportError:
        return False, "Biopython not installed. Please install with: pip install biopython"

    class DomainSelect(Select):
        """Select only residues within the domain boundaries"""
        def __init__(self, start, end):
            self.start = start
            self.end = end

        def accept_residue(self, residue):
            res_id = residue.get_id()[1]
            return self.start <= res_id <= self.end

    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure('protein', cif_file)

        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb, DomainSelect(start, end))

        # Verify the output file has content
        if os.path.getsize(output_pdb) < 100:
            return False, "Output PDB file is too small - domain may be outside structure"

        return True, None
    except Exception as e:
        return False, str(e)


def calculate_mean_plddt(cif_file, start, end):
    """
    Calculate mean pLDDT score for a specific region in the structure.
    pLDDT scores are stored in the B-factor field of AlphaFold models.

    Args:
        cif_file: Path to CIF file
        start: Start residue position (1-indexed)
        end: End residue position (1-indexed)

    Returns:
        float: Mean pLDDT score, or None if calculation fails
    """
    try:
        from Bio.PDB import MMCIFParser
    except ImportError:
        return None

    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure('protein', cif_file)

        chain = next(structure.get_chains())

        plddt_scores = []
        for residue in chain.get_residues():
            res_id = residue.get_id()[1]
            if start <= res_id <= end:
                # pLDDT is stored in the B-factor field
                # Use CA atom if available, otherwise first atom
                if 'CA' in residue:
                    plddt_scores.append(residue['CA'].get_bfactor())
                else:
                    for atom in residue.get_atoms():
                        plddt_scores.append(atom.get_bfactor())
                        break

        if plddt_scores:
            return sum(plddt_scores) / len(plddt_scores)
        return None
    except Exception:
        return None


def check_foldmason_installed(foldmason_path):
    """
    Check if FoldMason is installed and accessible.

    Args:
        foldmason_path: Path to foldmason executable

    Returns:
        tuple: (success boolean, version string or error message)
    """
    try:
        result = subprocess.run(
            [foldmason_path, '--help'],
            capture_output=True,
            text=True,
            timeout=10
        )
        # FoldMason returns 0 for --help
        return True, "FoldMason found"
    except FileNotFoundError:
        return False, f"FoldMason not found at '{foldmason_path}'. Please install or specify path with --foldmason"
    except subprocess.TimeoutExpired:
        return False, "FoldMason timed out"
    except Exception as e:
        return False, str(e)


def run_foldmason(structure_dir, output_prefix, foldmason_path, refine_iters=0, tmp_dir=None):
    """
    Run FoldMason easy-msa on a directory of PDB files.

    Args:
        structure_dir: Directory containing PDB files
        output_prefix: Prefix for output files
        foldmason_path: Path to foldmason executable
        refine_iters: Number of refinement iterations
        tmp_dir: Temporary directory for FoldMason

    Returns:
        tuple: (success boolean, dict with output paths or error message)
    """
    if tmp_dir is None:
        tmp_dir = tempfile.mkdtemp(prefix='foldmason_')

    foldmason_tmp = os.path.join(tmp_dir, 'foldmason_tmp')
    os.makedirs(foldmason_tmp, exist_ok=True)

    # Build command
    cmd = [
        foldmason_path, 'easy-msa',
        structure_dir,
        output_prefix,
        foldmason_tmp
    ]

    if refine_iters > 0:
        cmd.extend(['--refine-iters', str(refine_iters)])

    log_file = f"{output_prefix}.foldmason.log"

    try:
        print(f"Running FoldMason...")
        print(f"  Command: {' '.join(cmd)}")

        with open(log_file, 'w') as log:
            result = subprocess.run(
                cmd,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=3600  # 1 hour timeout
            )

        if result.returncode != 0:
            return False, f"FoldMason failed with return code {result.returncode}. Check {log_file}"

        # FoldMason outputs:
        # - output_prefix (FASTA alignment)
        # - output_prefix_aa.fa (amino acid alignment)
        # - output_prefix_3di.fa (3Di alphabet alignment)
        # - output_prefix.nw or similar (tree)

        outputs = {
            'alignment': output_prefix,
            'log': log_file
        }

        # Find tree file (might have different extensions)
        for ext in ['.nw', '.newick', '_tree.nw']:
            tree_file = f"{output_prefix}{ext}"
            if os.path.exists(tree_file):
                outputs['tree'] = tree_file
                break

        # Check for AA alignment
        aa_file = f"{output_prefix}_aa.fa"
        if os.path.exists(aa_file):
            outputs['aa_alignment'] = aa_file

        return True, outputs

    except subprocess.TimeoutExpired:
        return False, "FoldMason timed out after 1 hour"
    except Exception as e:
        return False, str(e)


def parse_fasta_alignment(fasta_file):
    """
    Parse a FASTA alignment file.

    Args:
        fasta_file: Path to FASTA file

    Returns:
        dict: {sequence_name: aligned_sequence}
    """
    sequences = {}
    current_name = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if current_name:
                    sequences[current_name] = ''.join(current_seq)
                current_name = line[1:].split()[0]  # Take first word after >
                current_seq = []
            else:
                current_seq.append(line)

        if current_name:
            sequences[current_name] = ''.join(current_seq)

    return sequences


def fasta_to_mul(fasta_file, name_mapping, output_file):
    """
    Convert FoldMason FASTA output to mul format with original sequence names.

    Args:
        fasta_file: Path to FoldMason FASTA output
        name_mapping: dict mapping PDB filename stems to original full names
        output_file: Path to output mul file

    Returns:
        tuple: (success boolean, error message or None)
    """
    try:
        alignments = parse_fasta_alignment(fasta_file)

        if not alignments:
            return False, "No sequences found in FoldMason output"

        # Find the maximum name length for formatting
        max_name_len = max(len(name_mapping.get(k, k)) for k in alignments.keys())

        with open(output_file, 'w') as f:
            for pdb_name, sequence in alignments.items():
                # Map back to original name
                original_name = name_mapping.get(pdb_name, pdb_name)

                # Convert gaps: FoldMason uses '-', mul format uses '.'
                sequence = sequence.replace('-', '.')

                # Write in mul format with proper spacing
                f.write(f"{original_name:<{max_name_len}}    {sequence}\n")

        return True, None
    except Exception as e:
        return False, str(e)


def write_log(log_file, seed_file, total_seqs, removed_seqs, plddt_warnings, successful_seqs, par_mode=False):
    """
    Write a log file with information about the realignment.

    Args:
        log_file: Path to log file
        seed_file: Path to original SEED file
        total_seqs: Total number of sequences in SEED
        removed_seqs: List of tuples (name, reason)
        plddt_warnings: List of tuples (name, plddt_score)
        successful_seqs: Number of successfully aligned sequences
        par_mode: Whether pad-and-realign mode was used
    """
    with open(log_file, 'w') as f:
        f.write(f"# foldmason_realign.log\n")
        f.write(f"# Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# Input: {os.path.abspath(seed_file)} ({total_seqs} sequences)\n")
        if par_mode:
            f.write(f"# Mode: Pad and Realign (two-pass)\n")
        f.write(f"\n")

        f.write(f"## Sequences removed (no AlphaFold model available): {len(removed_seqs)}\n")
        if removed_seqs:
            for name, reason in removed_seqs:
                f.write(f"- {name}: {reason}\n")
        else:
            f.write("- None\n")
        f.write(f"\n")

        f.write(f"## Low pLDDT warnings: {len(plddt_warnings)}\n")
        if plddt_warnings:
            for name, plddt in plddt_warnings:
                f.write(f"- {name}: mean pLDDT = {plddt:.1f}\n")
        else:
            f.write("- None\n")
        f.write(f"\n")

        f.write(f"## Summary\n")
        f.write(f"- Input sequences: {total_seqs}\n")
        f.write(f"- Successfully aligned: {successful_seqs}\n")
        f.write(f"- Removed: {len(removed_seqs)}")
        if total_seqs > 0:
            f.write(f" ({100*len(removed_seqs)/total_seqs:.1f}%)")
        f.write(f"\n")


def run_pad_ends(input_file, output_file):
    """
    Run pad_ends.pl to tidy up alignment ends.

    Args:
        input_file: Path to input alignment file
        output_file: Path to output alignment file

    Returns:
        tuple: (success boolean, error message or None)
    """
    cmd = ['pad_ends.pl', '-align', input_file]

    try:
        print(f"Running pad_ends.pl...")
        print(f"  Command: {' '.join(cmd)} > {output_file}")

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300
        )

        if result.returncode != 0:
            return False, f"pad_ends.pl failed: {result.stderr}"

        # Write stdout to output file
        with open(output_file, 'w') as f:
            f.write(result.stdout)

        return True, None

    except FileNotFoundError:
        return False, "pad_ends.pl not found in PATH"
    except subprocess.TimeoutExpired:
        return False, "pad_ends.pl timed out"
    except Exception as e:
        return False, str(e)


def download_models_to_cache(sequences, afdb_dir):
    """
    Download AlphaFold models to cache directory.

    Args:
        sequences: List of sequence dicts from parse_seed_alignment
        afdb_dir: Directory to store cached CIF files

    Returns:
        dict: {base_accession: cif_file_path} for successful downloads
    """
    os.makedirs(afdb_dir, exist_ok=True)

    # Get unique accessions
    unique_accessions = set(seq['base_accession'] for seq in sequences)
    print(f"  {len(unique_accessions)} unique UniProt accessions")

    downloaded_models = {}
    for i, base_acc in enumerate(sorted(unique_accessions), 1):
        # Check if already cached
        cached_v6 = os.path.join(afdb_dir, f"AF-{base_acc}-F1-model_v6.cif")
        cached_v4 = os.path.join(afdb_dir, f"AF-{base_acc}-F1-model_v4.cif")

        if os.path.exists(cached_v6):
            print(f"  [{i}/{len(unique_accessions)}] {base_acc}... cached (v6)")
            downloaded_models[base_acc] = cached_v6
            continue
        elif os.path.exists(cached_v4):
            print(f"  [{i}/{len(unique_accessions)}] {base_acc}... cached (v4)")
            downloaded_models[base_acc] = cached_v4
            continue

        print(f"  [{i}/{len(unique_accessions)}] Downloading {base_acc}...", end=' ')
        cif_file, error = download_alphafold_model(base_acc, afdb_dir)

        if cif_file:
            print("OK")
            downloaded_models[base_acc] = cif_file
        else:
            print(f"FAILED ({error})")

    return downloaded_models


def run_foldmason_alignment(sequences, downloaded_models, output_file, tree_file,
                            foldmason_path, refine_iters, plddt_warn, tmp_dir, pass_name=""):
    """
    Run the FoldMason alignment workflow.

    Args:
        sequences: List of sequence dicts from parse_seed_alignment
        downloaded_models: Dict of {base_accession: cif_file_path}
        output_file: Path to output mul file
        tree_file: Path to output tree file
        foldmason_path: Path to foldmason executable
        refine_iters: Number of refinement iterations
        plddt_warn: pLDDT warning threshold
        tmp_dir: Temporary directory for intermediate files
        pass_name: Name for this pass (for logging)

    Returns:
        tuple: (success, removed_seqs, plddt_warnings, successful_count)
    """
    pdb_dir = os.path.join(tmp_dir, f'pdb_files_{pass_name}' if pass_name else 'pdb_files')
    os.makedirs(pdb_dir, exist_ok=True)

    removed_seqs = []
    plddt_warnings = []
    name_mapping = {}
    successful_count = 0

    # Process each sequence
    if pass_name:
        print(f"\n[{pass_name}] Extracting domain structures...")
    else:
        print(f"\nExtracting domain structures...")

    for seq in sequences:
        base_acc = seq['base_accession']
        full_name = seq['full_name']

        # Check if we have the model
        if base_acc not in downloaded_models:
            removed_seqs.append((full_name, "AlphaFold model not available"))
            continue

        cif_file = downloaded_models[base_acc]

        # Calculate pLDDT for this domain
        plddt = calculate_mean_plddt(cif_file, seq['start'], seq['end'])
        if plddt is not None and plddt < plddt_warn:
            plddt_warnings.append((full_name, plddt))

        # Create a unique filename for this domain
        pdb_stem = f"{seq['version_accession']}_{seq['start']}_{seq['end']}"
        pdb_file = os.path.join(pdb_dir, f"{pdb_stem}.pdb")

        success, error = chop_structure_to_domain(
            cif_file, seq['start'], seq['end'], pdb_file
        )

        if success:
            name_mapping[pdb_stem] = full_name
            successful_count += 1
        else:
            removed_seqs.append((full_name, f"Structure extraction failed: {error}"))

    print(f"  Successfully processed: {successful_count}")
    print(f"  Removed: {len(removed_seqs)}")

    if successful_count < 2:
        return False, removed_seqs, plddt_warnings, successful_count

    # Report pLDDT warnings
    if plddt_warnings:
        print(f"\nWarning: {len(plddt_warnings)} sequences have low pLDDT (< {plddt_warn}):")
        for name, plddt in plddt_warnings[:5]:
            print(f"  - {name}: {plddt:.1f}")
        if len(plddt_warnings) > 5:
            print(f"  ... and {len(plddt_warnings) - 5} more (see log file)")

    # Run FoldMason
    if pass_name:
        print(f"\n[{pass_name}] Running FoldMason alignment...")
    else:
        print(f"\nRunning FoldMason alignment...")

    foldmason_prefix = os.path.join(tmp_dir, f'foldmason_out_{pass_name}' if pass_name else 'foldmason_out')

    success, result = run_foldmason(
        pdb_dir,
        foldmason_prefix,
        foldmason_path,
        refine_iters=refine_iters,
        tmp_dir=tmp_dir
    )

    if not success:
        print(f"Error: {result}", file=sys.stderr)
        return False, removed_seqs, plddt_warnings, successful_count

    print("  FoldMason completed successfully")

    # Convert output to mul format
    print(f"\nConverting to mul format...")

    fasta_file = None
    for candidate in [result.get('alignment'), result.get('aa_alignment'), foldmason_prefix]:
        if candidate and os.path.exists(candidate):
            fasta_file = candidate
            break

    if not fasta_file:
        print("Error: Could not find FoldMason alignment output", file=sys.stderr)
        return False, removed_seqs, plddt_warnings, successful_count

    success, error = fasta_to_mul(fasta_file, name_mapping, output_file)

    if not success:
        print(f"Error converting to mul format: {error}", file=sys.stderr)
        return False, removed_seqs, plddt_warnings, successful_count

    print(f"  Written: {output_file}")

    # Copy tree file if it exists
    if 'tree' in result and os.path.exists(result['tree']):
        shutil.copy(result['tree'], tree_file)
        print(f"  Tree file: {tree_file}")

    return True, removed_seqs, plddt_warnings, successful_count


def main():
    args = parse_arguments()

    # Resolve paths
    curation_dir = os.path.abspath(args.curation_dir)
    seed_file = os.path.join(curation_dir, 'SEED')
    afdb_dir = os.path.join(curation_dir, 'AFDB')

    if args.output:
        output_file = os.path.abspath(args.output)
    else:
        output_file = os.path.join(curation_dir, 'SEED.foldmason')

    tree_file = f"{output_file}.tree"
    log_file = os.path.join(curation_dir, 'foldmason_realign.log')

    # Check SEED file exists
    if not os.path.exists(seed_file):
        print(f"Error: SEED file not found at {seed_file}", file=sys.stderr)
        sys.exit(1)

    # Check FoldMason is installed
    print("Checking FoldMason installation...")
    success, msg = check_foldmason_installed(args.foldmason)
    if not success:
        print(f"Error: {msg}", file=sys.stderr)
        sys.exit(1)
    print(f"  {msg}")

    if args.par:
        print("\n" + "="*60)
        print("PAD AND REALIGN MODE (two-pass)")
        print("="*60)

    # Parse SEED file
    print(f"\nParsing SEED file: {seed_file}")
    sequences = parse_seed_alignment(seed_file)
    print(f"  Found {len(sequences)} sequences")

    if len(sequences) == 0:
        print("Error: No sequences found in SEED file", file=sys.stderr)
        sys.exit(1)

    # Create temporary directory for intermediate files
    if args.keep_temp:
        tmp_dir = tempfile.mkdtemp(prefix='foldmason_realign_')
        print(f"\nTemporary directory (will be kept): {tmp_dir}")
    else:
        tmp_context = tempfile.TemporaryDirectory(prefix='foldmason_realign_')
        tmp_dir = tmp_context.name
        print(f"\nTemporary directory: {tmp_dir}")

    # Download AlphaFold models to AFDB cache directory
    print(f"\nDownloading AlphaFold models to {afdb_dir}...")
    downloaded_models = download_models_to_cache(sequences, afdb_dir)

    if args.par:
        # ===== PASS 1: Initial alignment =====
        print("\n" + "-"*60)
        print("PASS 1: Initial FoldMason alignment")
        print("-"*60)

        pass1_output = os.path.join(tmp_dir, 'pass1.foldmason')
        pass1_tree = os.path.join(tmp_dir, 'pass1.tree')

        success, removed_seqs, plddt_warnings, successful_count = run_foldmason_alignment(
            sequences, downloaded_models, pass1_output, pass1_tree,
            args.foldmason, args.refine_iters, args.plddt_warn, tmp_dir, "Pass1"
        )

        if not success:
            print("\nError: Pass 1 alignment failed", file=sys.stderr)
            write_log(log_file, seed_file, len(sequences), removed_seqs, plddt_warnings,
                     successful_count, par_mode=True)
            print(f"\nLog written to: {log_file}")
            sys.exit(1)

        # ===== PAD ENDS =====
        print("\n" + "-"*60)
        print("PADDING ENDS")
        print("-"*60)

        # Backup original SEED
        seed_backup = os.path.join(tmp_dir, 'SEED.original')
        shutil.copy(seed_file, seed_backup)
        print(f"  Backed up original SEED to: {seed_backup}")

        # Run pad_ends.pl: SEED.foldmason -> SEED
        success, error = run_pad_ends(pass1_output, seed_file)
        if not success:
            print(f"Error: {error}", file=sys.stderr)
            # Restore original SEED
            shutil.copy(seed_backup, seed_file)
            print("  Restored original SEED")
            sys.exit(1)

        print("  padd_ends.pl completed successfully")
        print(f"  Updated SEED file with padded ends")

        # ===== PASS 2: Realignment with padded SEED =====
        print("\n" + "-"*60)
        print("PASS 2: FoldMason realignment with padded ends")
        print("-"*60)

        # Re-parse the updated SEED file (coordinates may have changed)
        print(f"\nRe-parsing padded SEED file...")
        sequences_pass2 = parse_seed_alignment(seed_file)
        print(f"  Found {len(sequences_pass2)} sequences")

        success, removed_seqs, plddt_warnings, successful_count = run_foldmason_alignment(
            sequences_pass2, downloaded_models, output_file, tree_file,
            args.foldmason, args.refine_iters, args.plddt_warn, tmp_dir, "Pass2"
        )

        if not success:
            print("\nError: Pass 2 alignment failed", file=sys.stderr)
            write_log(log_file, seed_file, len(sequences_pass2), removed_seqs, plddt_warnings,
                     successful_count, par_mode=True)
            print(f"\nLog written to: {log_file}")
            sys.exit(1)

        total_seqs = len(sequences_pass2)

    else:
        # ===== SINGLE PASS MODE =====
        success, removed_seqs, plddt_warnings, successful_count = run_foldmason_alignment(
            sequences, downloaded_models, output_file, tree_file,
            args.foldmason, args.refine_iters, args.plddt_warn, tmp_dir
        )

        if not success:
            print("\nError: Alignment failed", file=sys.stderr)
            write_log(log_file, seed_file, len(sequences), removed_seqs, plddt_warnings,
                     successful_count)
            print(f"\nLog written to: {log_file}")
            sys.exit(1)

        total_seqs = len(sequences)

    # Write log file
    write_log(log_file, seed_file, total_seqs, removed_seqs, plddt_warnings,
              successful_count, par_mode=args.par)
    print(f"  Log file: {log_file}")

    # Summary
    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    if args.par:
        print(f"  Mode:             Pad and Realign (two-pass)")
    print(f"  Input sequences:  {total_seqs}")
    print(f"  Aligned:          {successful_count}")
    print(f"  Removed:          {len(removed_seqs)} ({100*len(removed_seqs)/total_seqs:.1f}%)")
    if plddt_warnings:
        print(f"  Low pLDDT:        {len(plddt_warnings)}")
    print(f"\nOutput files:")
    print(f"  Alignment: {output_file}")
    if os.path.exists(tree_file):
        print(f"  Tree:      {tree_file}")
    print(f"  Log:       {log_file}")
    print(f"  AFDB cache: {afdb_dir}")
    if args.par:
        print(f"\nNote: SEED file has been updated with padded ends")

    # Clean up temp directory if not keeping
    if not args.keep_temp:
        tmp_context.cleanup()


if __name__ == '__main__':
    main()
