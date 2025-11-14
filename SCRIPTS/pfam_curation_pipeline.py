#!/usr/bin/env python3
"""
Automated Pfam Curation Directory Builder

This script automates the process of building Pfam curation directories from
UniProt accessions, either from a list file or by downloading from a proteome.

Usage:
    pfam_curation_pipeline.py --list accessions.txt
    pfam_curation_pipeline.py --proteome UP000000625
    pfam_curation_pipeline.py --list accessions.txt --resume
"""

import argparse
import subprocess
import sys
import time
import os
from pathlib import Path
import requests
from typing import List, Set

class CurationPipeline:
    """Manages the Pfam curation directory building pipeline"""
    
    # Hidden files to track progress
    PROGRESS_FILES = {
        'ted_build': '.ted_build_complete',
        'initial_triage': '.initial_triage_complete',
        'iteration': '.iteration_complete',
        'add_author': '.add_author_complete',
        'add_sp': '.add_sp_complete',
        'add_pdb': '.add_pdb_complete',
        'add_abb': '.add_abb_complete',
        'add_paperblast': '.add_paperblast_complete',
        'add_species': '.add_species_complete',
        'add_swiss': '.add_swiss_complete',
        'add_ted': '.add_ted_complete',
        'add_foldseek': '.add_foldseek_complete'
    }
    
    def __init__(self, working_dir: str = '.'):
        self.working_dir = Path(working_dir)
        self.accessions: List[str] = []
        
    def get_proteome_accessions(self, proteome_id: str) -> List[str]:
        """
        Download list of UniProt accessions for a given proteome ID
        Gets all accessions without any size limit.
        
        Args:
            proteome_id: UniProt proteome identifier (e.g., UP000000625)
            
        Returns:
            List of UniProt accessions
        """
        print(f"Downloading accessions for proteome {proteome_id}...")
        print("(This may take a while for large proteomes)")
        
        url = f"https://rest.uniprot.org/uniprotkb/stream"
        params = {
            'query': f'proteome:{proteome_id}',
            'format': 'list',
            'compressed': 'false'  # Ensure we get plain text
        }
        
        try:
            response = requests.get(url, params=params, timeout=300)  # 5 min timeout
            response.raise_for_status()
            
            accessions = [line.strip() for line in response.text.strip().split('\n') if line.strip()]
            
            if not accessions:
                print(f"⚠ Warning: No accessions found for proteome {proteome_id}")
                print("   Please check the proteome ID is correct")
                sys.exit(1)
            
            print(f"✓ Downloaded {len(accessions)} accessions")
            return accessions
            
        except requests.exceptions.Timeout:
            print(f"⚠ Error: Request timed out after 5 minutes", file=sys.stderr)
            print("   The proteome may be very large. Try breaking it into smaller chunks.", file=sys.stderr)
            sys.exit(1)
        except requests.exceptions.RequestException as e:
            print(f"⚠ Error downloading proteome: {e}", file=sys.stderr)
            sys.exit(1)
    
    def load_accessions_from_file(self, filename: str) -> List[str]:
        """
        Load UniProt accessions from a file (one per line)
        
        Args:
            filename: Path to file containing accessions
            
        Returns:
            List of UniProt accessions
        """
        print(f"Loading accessions from {filename}...")
        try:
            with open(filename, 'r') as f:
                accessions = [line.strip() for line in f if line.strip()]
            print(f"Loaded {len(accessions)} accessions")
            return accessions
        except IOError as e:
            print(f"Error reading file: {e}", file=sys.stderr)
            sys.exit(1)
    
    def is_step_complete(self, step: str) -> bool:
        """Check if a pipeline step has been completed"""
        progress_file = self.working_dir / self.PROGRESS_FILES[step]
        return progress_file.exists()
    
    def mark_step_complete(self, step: str):
        """Mark a pipeline step as completed"""
        progress_file = self.working_dir / self.PROGRESS_FILES[step]
        progress_file.touch()
        print(f"✓ Marked {step} as complete")
    
    def wait_for_jobs(self, timeout: int = 3600):
        """
        Wait for all user jobs to complete by monitoring squeue
        
        Args:
            timeout: Maximum time to wait in seconds (default 1 hour)
        """
        print("Waiting for jobs to complete...")
        start_time = time.time()
        
        while True:
            try:
                # Check for any running jobs from this user
                result = subprocess.run(
                    ['squeue', '--me', '--noheader'],
                    capture_output=True,
                    text=True,
                    check=True
                )
                
                if not result.stdout.strip():
                    print("✓ All jobs completed")
                    return True
                
                # Count running jobs
                job_count = len([line for line in result.stdout.strip().split('\n') if line])
                print(f"  {job_count} job(s) still running...", end='\r')
                
                # Check timeout
                if time.time() - start_time > timeout:
                    print(f"\n⚠ Timeout reached after {timeout}s")
                    return False
                
                time.sleep(30)  # Check every 30 seconds
                
            except subprocess.CalledProcessError as e:
                print(f"\n⚠ Error checking job status: {e}", file=sys.stderr)
                return False
    
    def accession_directories_exist(self, accession: str) -> bool:
        """
        Check if any directories already exist for this accession

        Args:
            accession: UniProt accession to check

        Returns:
            True if directories exist, False otherwise
        """
        # Check for directories that start with the accession
        # TED domains typically create directories like: ACCESSION_1-100
        # build.sh creates directories like: ACCESSION
        pattern = f"{accession}*"
        matching_dirs = list(self.working_dir.glob(pattern))
        # Filter to only directories (not files)
        matching_dirs = [d for d in matching_dirs if d.is_dir()]
        return len(matching_dirs) > 0

    def run_ted_build(self):
        """Run ted_build.py for each accession, fallback to build.sh if no TED domains"""
        if self.is_step_complete('ted_build'):
            print("⊳ TED build already completed, skipping...")
            return

        print(f"\n{'='*60}")
        print("STEP 1: Running TED build for all accessions")
        print(f"{'='*60}")

        failed_accessions = []
        skipped_accessions = []

        for i, accession in enumerate(self.accessions, 1):
            # Check if directories already exist for this accession
            if self.accession_directories_exist(accession):
                print(f"[{i}/{len(self.accessions)}] Skipping {accession} (directories already exist)")
                skipped_accessions.append(accession)
                continue

            print(f"[{i}/{len(self.accessions)}] Building domains for {accession}")

            # Count directories before TED build
            dirs_before = set(d.name for d in self.working_dir.iterdir() if d.is_dir())

            try:
                # Try TED build first
                result = subprocess.run(
                    ['ted_build.py', accession],
                    capture_output=True,
                    text=True,
                    cwd=self.working_dir
                )

                # Count directories after TED build
                dirs_after = set(d.name for d in self.working_dir.iterdir() if d.is_dir())
                new_dirs = dirs_after - dirs_before

                # If no new directories created, TED probably found no domains
                if not new_dirs:
                    print(f"  → No TED domains found, using whole sequence (build.sh)")
                    subprocess.run(
                        ['bash', 'build.sh', accession],
                        check=True,
                        cwd=self.working_dir
                    )
                else:
                    print(f"  → Created {len(new_dirs)} TED domain(s)")

            except subprocess.CalledProcessError as e:
                print(f"  ⚠ Failed to build {accession}: {e}")
                failed_accessions.append(accession)
        
        # Print summary
        if skipped_accessions:
            print(f"\n✓ Skipped {len(skipped_accessions)} accessions (directories already exist)")

        if failed_accessions:
            print(f"\n⚠ Warning: {len(failed_accessions)} accessions failed:")
            for acc in failed_accessions:
                print(f"  - {acc}")

        # Wait for all pfbuild jobs to complete
        if not self.wait_for_jobs():
            print("⚠ Warning: Some jobs may still be running")

        self.mark_step_complete('ted_build')
    
    def run_triage(self) -> bool:
        """
        Run triage.pl to generate triage file
        
        Returns:
            True if successful, False otherwise
        """
        print("\nRunning triage.pl...")
        try:
            with open(self.working_dir / 'triage', 'w') as f:
                subprocess.run(
                    ['triage.pl'],
                    stdout=f,
                    check=True,
                    cwd=self.working_dir
                )
            print("✓ Triage completed")
            return True
        except subprocess.CalledProcessError as e:
            print(f"⚠ Error running triage: {e}", file=sys.stderr)
            return False
    
    def run_iteration_cycle(self) -> bool:
        """
        Run a single iteration cycle (iterate.sh + wait + triage)
        
        Returns:
            True if new Iterate directories were created, False otherwise
        """
        print("\nRunning iteration cycle...")

        # Count existing Iterate directories before (at all levels)
        before_count = len(list(self.working_dir.glob('**/Iterate')))

        # Run iterate.sh equivalent
        try:
            with open(self.working_dir / 'triage', 'r') as f:
                triage_lines = f.readlines()

            # Shuffle lines and filter
            import random
            random.shuffle(triage_lines)

            for line in triage_lines:
                # Skip if at max Iterate depth (6 levels)
                if 'Iterate/Iterate/Iterate/Iterate/Iterate/Iterate/Iterate' in line:
                    continue
                
                fields = line.split()
                if len(fields) < 5:
                    continue
                
                try:
                    # Check conditions: $5>0.99 && $4>1
                    if float(fields[4]) > 0.99 and float(fields[3]) > 1:
                        directory = fields[0]
                        print(f"  Iterating {directory}")
                        subprocess.run(
                            ['/homes/agb/Scripts/iterate_inline.pl', directory],
                            check=True,
                            cwd=self.working_dir
                        )
                except (ValueError, subprocess.CalledProcessError) as e:
                    continue
            
        except IOError as e:
            print(f"⚠ Error reading triage file: {e}", file=sys.stderr)
            return False
        
        # Wait for pfbuilds to complete
        if not self.wait_for_jobs():
            print("⚠ Warning: Some jobs may still be running")

        # Run triage again
        self.run_triage()

        # Count Iterate directories after (at all levels)
        after_count = len(list(self.working_dir.glob('**/Iterate')))
        
        new_directories = after_count - before_count
        if new_directories > 0:
            print(f"✓ Created {new_directories} new Iterate directories")
            return True
        else:
            print("✓ No new Iterate directories created")
            return False
    
    def run_iterations(self):
        """Run iteration cycles until no new directories are created"""
        if self.is_step_complete('iteration'):
            print("⊳ Iterations already completed, skipping...")
            return
        
        print(f"\n{'='*60}")
        print("STEP 2: Running initial triage")
        print(f"{'='*60}")
        
        if not self.is_step_complete('initial_triage'):
            self.run_triage()
            self.mark_step_complete('initial_triage')
        
        print(f"\n{'='*60}")
        print("STEP 3: Running iteration cycles")
        print(f"{'='*60}")
        
        iteration_count = 0
        max_iterations = 20  # Safety limit
        
        while iteration_count < max_iterations:
            iteration_count += 1
            print(f"\n--- Iteration cycle {iteration_count} ---")
            
            if not self.run_iteration_cycle():
                print("\n✓ All iterations complete")
                break
        else:
            print(f"\n⚠ Warning: Reached maximum iteration limit ({max_iterations})")
        
        self.mark_step_complete('iteration')
    
    def get_good_quality_directories(self) -> Set[str]:
        """
        Parse triage file and return directories with > 50% non-overlapping sequences
        
        Returns:
            Set of directory names that pass quality filter
        """
        triage_file = self.working_dir / 'triage'
        if not triage_file.exists():
            print("⚠ Warning: triage file not found, processing all directories")
            return set()
        
        good_dirs = set()
        skipped_dirs = []
        
        try:
            with open(triage_file, 'r') as f:
                for line in f:
                    fields = line.strip().split()
                    if len(fields) < 5:
                        continue
                    
                    directory = fields[0]
                    try:
                        fraction = float(fields[4])
                        if fraction > 0.5:
                            good_dirs.add(directory)
                        else:
                            skipped_dirs.append((directory, fraction))
                    except ValueError:
                        continue
            
            print(f"\n{'='*60}")
            print("Quality Filter Results")
            print(f"{'='*60}")
            print(f"✓ Good quality directories (>50% non-overlapping): {len(good_dirs)}")
            print(f"⊗ Skipped directories (≤50% non-overlapping): {len(skipped_dirs)}")
            
            if skipped_dirs and len(skipped_dirs) <= 20:
                print(f"\nSkipped directories:")
                for dir_name, frac in sorted(skipped_dirs, key=lambda x: x[1]):
                    print(f"  - {dir_name} ({frac:.2%} non-overlapping)")
            elif len(skipped_dirs) > 20:
                print(f"\nExample skipped directories (showing worst 10):")
                for dir_name, frac in sorted(skipped_dirs, key=lambda x: x[1])[:10]:
                    print(f"  - {dir_name} ({frac:.2%} non-overlapping)")
            
            return good_dirs
            
        except IOError as e:
            print(f"⚠ Error reading triage file: {e}", file=sys.stderr)
            return set()
    
    def create_filtered_triage(self) -> Path:
        """
        Create a filtered triage file with only high-quality directories
        
        Returns:
            Path to filtered triage file
        """
        triage_file = self.working_dir / 'triage'
        filtered_file = self.working_dir / 'triage.filtered'
        
        if not triage_file.exists():
            print("⚠ Warning: triage file not found")
            return None
        
        good_lines = []
        skipped_count = 0
        
        with open(triage_file, 'r') as f:
            for line in f:
                fields = line.strip().split()
                if len(fields) < 5:
                    continue
                
                try:
                    fraction = float(fields[4])
                    if fraction > 0.5:
                        good_lines.append(line)
                    else:
                        skipped_count += 1
                except ValueError:
                    continue
        
        with open(filtered_file, 'w') as f:
            f.writelines(good_lines)
        
        print(f"Created filtered triage: {len(good_lines)} good, {skipped_count} skipped")
        return filtered_file
    
    def run_annotation_scripts(self):
        """Run all annotation shell scripts on high-quality directories only"""
        print(f"\n{'='*60}")
        print("STEP 4: Running annotation scripts")
        print(f"{'='*60}")
        
        # Get filtered list and display summary
        good_dirs = self.get_good_quality_directories()
        
        if not good_dirs:
            print("⚠ No directories passed quality filter - skipping annotation")
            return
        
        # Create filtered triage file
        filtered_triage = self.create_filtered_triage()
        if not filtered_triage:
            print("⚠ Could not create filtered triage file")
            return
        
        # Dictionary of scripts and their commands
        # Using filtered triage file instead of regular triage
        scripts = {
            'add_author': {
                'cmd': f'awk \'{{system("add_author.pl "$1)}}\' {filtered_triage}',
                'desc': 'Adding authors'
            },
            'add_sp': {
                'cmd': f'''current_dir=$(pwd)
awk '{{
    system("cd "$1" && \\
    if [ -e sp ]; then \\
        echo \\"Skipping "$1": sp file already exists\\"; \\
    else \\
        echo \\"Processing "$1"\\"; \\
        align_lines=$(grep -c \\"^\\" ALIGN 2>/dev/null || echo 0); \\
        n_value=$(( align_lines < 100 ? align_lines : 100 )); \\
        if [ $n_value -gt 0 ]; then \\
            swissprot.pl -n $n_value; \\
        else \\
            echo \\"Warning: No ALIGN file or empty ALIGN in "$1"\\"; \\
        fi; \\
    fi && \\
    cd '"\\"$current_dir\\""'")
}}' {filtered_triage}''',
                'desc': 'Adding SwissProt files'
            },
            'add_pdb': {
                'cmd': f'current_dir=$(pwd); awk \'{{system("cd "$1"; echo "$1"; add_pdb_ref.py; cd \'"$current_dir"\'")}}\'  {filtered_triage}',
                'desc': 'Adding PDB references'
            },
            'add_abb': {
                'cmd': f'current_dir=$(pwd); awk \'{{system("cd "$1"; echo "$1"; add_abb_ref.py; cd \'"$current_dir"\'")}}\'  {filtered_triage}',
                'desc': 'Adding abbreviation references'
            },
            'add_paperblast': {
                'cmd': f'current_dir=$(pwd); awk \'{{system("cd "$1"; echo "$1"; query_paperblast.py; cd \'"$current_dir"\'")}}\'  {filtered_triage}',
                'desc': 'Adding PaperBLAST files'
            },
            'add_species': {
                'cmd': f'shuf {filtered_triage} | awk \'{{system("species_summary.pl "$1";")}}\' ',
                'desc': 'Adding species files'
            },
            'add_swiss': {
                'cmd': f'current_dir=$(pwd); awk \'$8>0 && $8<1000{{system("cd "$1"; echo "$1"; add_swiss_ref.pl; cd \'"$current_dir"\'")}}\'  {filtered_triage}',
                'desc': 'Adding Swiss-Prot references'
            },
            'add_ted': {
                'cmd': f'current_dir=$(pwd); awk \'{{system("cd "$1"; echo "$1"; ted_ali.pl SEED; cd \'"$current_dir"\'")}}\'  {filtered_triage}',
                'desc': 'Adding TED files'
            },
            'add_foldseek': {
                'cmd': f'current_dir=$(pwd); awk \'{{system("cd "$1"; echo "$1"; add_foldseek.sh; cd \'"$current_dir"\'")}}\'  {filtered_triage}',
                'desc': 'Adding Foldseek files'
            }
        }
        
        for step_name, script_info in scripts.items():
            if self.is_step_complete(step_name):
                print(f"⊳ {script_info['desc']} already complete, skipping...")
                continue
            
            print(f"\n{script_info['desc']}...")
            try:
                subprocess.run(
                    script_info['cmd'],
                    shell=True,
                    check=True,
                    cwd=self.working_dir,
                    executable='/bin/bash'
                )
                print(f"✓ {script_info['desc']} completed")
                self.mark_step_complete(step_name)
            except subprocess.CalledProcessError as e:
                print(f"⚠ Error running {step_name}: {e}", file=sys.stderr)
                print("  Continuing with remaining scripts...")
    
    def run_pipeline(self):
        """Run the complete curation pipeline"""
        print(f"\n{'='*60}")
        print("PFAM CURATION PIPELINE")
        print(f"{'='*60}")
        print(f"Working directory: {self.working_dir.absolute()}")
        print(f"Number of accessions: {len(self.accessions)}")
        print(f"{'='*60}\n")
        
        # Step 1: Run TED build
        self.run_ted_build()
        
        # Step 2-3: Run iterations
        self.run_iterations()
        
        # Step 4: Run annotation scripts
        self.run_annotation_scripts()
        
        print(f"\n{'='*60}")
        print("✓ PIPELINE COMPLETE")
        print(f"{'='*60}\n")
        
        # Print summary
        print("Summary of completed steps:")
        for step, progress_file in self.PROGRESS_FILES.items():
            status = "✓" if self.is_step_complete(step) else "✗"
            print(f"  {status} {step}")

def main():
    parser = argparse.ArgumentParser(
        description='Automate Pfam curation directory building',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --list accessions.txt
  %(prog)s --proteome UP000000625
  %(prog)s --list accessions.txt --resume
        """
    )
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '--list',
        help='File containing UniProt accessions (one per line)'
    )
    input_group.add_argument(
        '--proteome',
        help='UniProt Proteome ID (e.g., UP000000625)'
    )
    
    parser.add_argument(
        '--working-dir',
        default='.',
        help='Working directory for curation (default: current directory)'
    )
    
    parser.add_argument(
        '--resume',
        action='store_true',
        help='Resume from last completed step'
    )
    
    parser.add_argument(
        '--reset',
        action='store_true',
        help='Reset all progress markers and start from beginning'
    )
    
    args = parser.parse_args()
    
    # Initialize pipeline
    pipeline = CurationPipeline(args.working_dir)
    
    # Reset progress if requested
    if args.reset:
        print("Resetting all progress markers...")
        for progress_file in pipeline.PROGRESS_FILES.values():
            progress_path = pipeline.working_dir / progress_file
            if progress_path.exists():
                progress_path.unlink()
        print("✓ All progress markers removed")
    
    # Load accessions
    if args.list:
        pipeline.accessions = pipeline.load_accessions_from_file(args.list)
    else:
        pipeline.accessions = pipeline.get_proteome_accessions(args.proteome)
    
    # Run pipeline
    try:
        pipeline.run_pipeline()
    except KeyboardInterrupt:
        print("\n\n⚠ Pipeline interrupted by user")
        print("You can resume by running the same command with --resume")
        sys.exit(1)
    except Exception as e:
        print(f"\n⚠ Pipeline failed with error: {e}", file=sys.stderr)
        print("You can resume by running the same command with --resume")
        sys.exit(1)

if __name__ == '__main__':
    main()
