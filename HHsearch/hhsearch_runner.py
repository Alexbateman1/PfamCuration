#!/usr/bin/env python3
"""
HHsearch runner for all-against-all Pfam family comparisons.
Uses HHsearch for profile-profile comparisons (HMM vs HMM).
"""

import subprocess
import os
import logging
from pathlib import Path
import re

class HHsearchRunner:
    """Manages HHsearch searches and result parsing."""

    def __init__(self, seed_dir, results_dir, e_value_threshold=1.0):
        """
        Initialize HHsearch runner.

        Args:
            seed_dir: Directory containing SEED alignment files
            results_dir: Directory for HHsearch results
            e_value_threshold: E-value cutoff for significant hits (default: 1.0)
        """
        self.seed_dir = Path(seed_dir)
        self.results_dir = Path(results_dir)
        self.e_value_threshold = e_value_threshold

        # Create subdirectories
        self.a3m_dir = self.results_dir / "A3M"
        self.hhr_dir = self.results_dir / "HHR"
        self.parsed_dir = self.results_dir / "PARSED"

        for dir_path in [self.a3m_dir, self.hhr_dir, self.parsed_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)

        logging.info(f"HHsearch runner initialized")
        logging.info(f"Results directory: {self.results_dir}")
        logging.info(f"E-value threshold: {self.e_value_threshold}")

    def mul_to_a3m(self, seed_file, output_file):
        """
        Convert Pfam SEED alignment (Stockholm format) to aligned FASTA for hhmake.

        Note: We use simple conversion rather than reformat.pl because:
        - reformat.pl creates proper A3M with insertion encoding (lowercase)
        - This causes sequences to have different lengths
        - hhmake with -M first can't handle variable-length input

        This creates aligned FASTA (all sequences same length, gaps as '-')
        which hhmake can process with -M first flag.

        Args:
            seed_file: Path to SEED alignment file (Stockholm format)
            output_file: Path to output A3M file

        Returns:
            True if successful, False otherwise
        """
        try:
            sequences = {}
            with open(seed_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    # Skip Stockholm annotation lines (start with #)
                    if line.startswith('#'):
                        continue

                    parts = line.split(None, 1)
                    if len(parts) >= 2:
                        seq_id = parts[0]
                        seq = parts[1]
                        # Convert gap character from '.' to '-' for A3M format
                        seq = seq.replace('.', '-')
                        # Concatenate if this seq_id was seen before (multi-line sequences in Stockholm)
                        sequences[seq_id] = sequences.get(seq_id, '') + seq

            # Write aligned FASTA format (all sequences same length)
            with open(output_file, 'w') as f:
                for seq_id, seq in sequences.items():
                    f.write(f">{seq_id}\n{seq}\n")

            logging.debug(f"Converted {seed_file} to aligned FASTA format ({len(sequences)} sequences)")
            return True

        except Exception as e:
            logging.error(f"Failed to convert {seed_file} to A3M: {e}")
            return False

    def build_hhm_database(self):
        """
        Build HH-suite database from SEED alignments using hhmake.
        Creates three ffindex databases with consistent base name: pfam_a3m, pfam_hhm, pfam_cs219.

        Returns:
            Path to database base name (for use with hhsearch -d)
        """
        db_dir = self.results_dir / "hhsuite_db"
        db_dir.mkdir(exist_ok=True)

        try:
            seed_files = sorted(self.seed_dir.glob("PF*_SEED"))
            logging.info(f"Building HH-suite database from {len(seed_files)} SEED alignments...")
            logging.info("This may take some time (converting alignments to HH-suite profiles)...")

            # Create temporary directories for building databases
            a3m_build_dir = db_dir / "a3m_dir"
            hhm_build_dir = db_dir / "hhm_dir"
            a3m_build_dir.mkdir(exist_ok=True)
            hhm_build_dir.mkdir(exist_ok=True)

            failed = 0
            skipped = 0
            built = 0
            for i, seed_file in enumerate(seed_files, 1):
                if i % 100 == 0:
                    logging.info(f"Progress: {i}/{len(seed_files)} processed - {built} built, {skipped} skipped (already exist), {failed} failed")

                family_id = seed_file.stem.replace('_SEED', '')
                a3m_file = self.a3m_dir / f"{family_id}.a3m"
                hhm_file = hhm_build_dir / f"{family_id}.hhm"
                a3m_copy = a3m_build_dir / f"{family_id}.a3m"

                # Skip if already built
                if hhm_file.exists() and a3m_copy.exists():
                    skipped += 1
                    continue

                # Convert SEED to A3M using reformat.pl for proper A3M format
                if not a3m_file.exists():
                    cmd = f"reformat.pl sto a3m {seed_file} {a3m_file}"
                    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                    if result.returncode != 0:
                        logging.warning(f"reformat.pl failed for {family_id}: {result.stderr}")
                        failed += 1
                        continue

                # Copy A3M to build directory
                if not a3m_copy.exists():
                    import shutil
                    shutil.copy(a3m_file, a3m_copy)

                # Build HH-suite profile with hhmake
                if not hhm_file.exists():
                    cmd = f"hhmake -i {a3m_file} -o {hhm_file} -v 0"
                    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

                    if result.returncode != 0:
                        logging.warning(f"hhmake failed for {family_id}: {result.stderr}")
                        failed += 1
                    else:
                        built += 1
                else:
                    built += 1

            logging.info(f"Profile building complete: {built} built, {skipped} skipped (already exist), {failed} failed")

            # Build ffindex databases with consistent base name 'pfam'
            base_name = "pfam"

            logging.info("Building pfam_a3m ffindex database...")
            cmd = f"ffindex_build -as {db_dir}/{base_name}_a3m.ffdata {db_dir}/{base_name}_a3m.ffindex {a3m_build_dir}/"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if result.returncode != 0:
                logging.error(f"ffindex_build for A3M failed: {result.stderr}")
                raise RuntimeError("Failed to build A3M ffindex database")

            logging.info("Building pfam_hhm ffindex database...")
            cmd = f"ffindex_build -as {db_dir}/{base_name}_hhm.ffdata {db_dir}/{base_name}_hhm.ffindex {hhm_build_dir}/"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if result.returncode != 0:
                logging.error(f"ffindex_build for HHM failed: {result.stderr}")
                raise RuntimeError("Failed to build HHM ffindex database")

            logging.info("Building pfam_cs219 database with cstranslate...")
            cmd = f"cstranslate -i {db_dir}/{base_name}_a3m -o {db_dir}/{base_name}_cs219 -b -f -I a3m"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if result.returncode != 0:
                logging.warning(f"cstranslate failed: {result.stderr}")
                logging.warning("Proceeding without CS219 prefiltering (searches will be slower)")

            logging.info(f"HH-suite database created: {db_dir}/{base_name}")
            return str(db_dir / base_name)

        except Exception as e:
            logging.error(f"Failed to build HH-suite database: {e}")
            raise

    def _build_concatenated_hhm_database(self, hhm_dir, db_dir):
        """Fallback: concatenate HHM files into single file."""
        logging.info("Building concatenated HH-suite database...")

        db_file = db_dir / "pfam_db_concat"
        hhm_files = sorted(hhm_dir.glob("*.hhm"))

        with open(db_file, 'w') as outf:
            for hhm_file in hhm_files:
                with open(hhm_file, 'r') as inf:
                    outf.write(inf.read())
                    outf.write('\n')

        logging.info(f"Concatenated database created: {db_file}")
        return str(db_file)

    def clean_results(self, clean_a3m=True, clean_hhm=True, clean_hhr=True, clean_parsed=True, clean_db=True, clean_summary=True):
        """
        Clean up results from previous runs to ensure fresh start.

        Args:
            clean_a3m: Remove A3M alignment files
            clean_hhm: Remove HHM profile files
            clean_hhr: Remove HHR search result files
            clean_parsed: Remove parsed TSV files
            clean_db: Remove ffindex database files
            clean_summary: Remove summary files (aggregated results)
        """
        logging.info("=== Cleaning up old results ===")

        cleaned_counts = {
            'a3m': 0,
            'hhm': 0,
            'hhr': 0,
            'parsed': 0,
            'db_files': 0,
            'summary': 0
        }

        # Clean A3M files
        if clean_a3m and self.a3m_dir.exists():
            a3m_files = list(self.a3m_dir.glob("*.a3m"))
            cleaned_counts['a3m'] = len(a3m_files)
            for f in a3m_files:
                f.unlink()
            logging.info(f"Removed {cleaned_counts['a3m']} A3M files from {self.a3m_dir}")

        # Clean HHM files
        if clean_hhm:
            hhm_dir = self.results_dir / "hhsuite_db" / "hhm_files"
            if hhm_dir.exists():
                hhm_files = list(hhm_dir.glob("*.hhm"))
                cleaned_counts['hhm'] = len(hhm_files)
                for f in hhm_files:
                    f.unlink()
                logging.info(f"Removed {cleaned_counts['hhm']} HHM files from {hhm_dir}")

        # Clean HHR files
        if clean_hhr and self.hhr_dir.exists():
            hhr_files = list(self.hhr_dir.glob("*.hhr"))
            cleaned_counts['hhr'] = len(hhr_files)
            for f in hhr_files:
                f.unlink()
            logging.info(f"Removed {cleaned_counts['hhr']} HHR files from {self.hhr_dir}")

        # Clean parsed files
        if clean_parsed and self.parsed_dir.exists():
            parsed_files = list(self.parsed_dir.glob("*_hits.tsv"))
            cleaned_counts['parsed'] = len(parsed_files)
            for f in parsed_files:
                f.unlink()
            logging.info(f"Removed {cleaned_counts['parsed']} parsed TSV files from {self.parsed_dir}")

        # Clean database files
        if clean_db:
            db_dir = self.results_dir / "hhsuite_db"
            if db_dir.exists():
                db_files = list(db_dir.glob("pfam_db_hhm.*"))
                cleaned_counts['db_files'] = len(db_files)
                for f in db_files:
                    f.unlink()
                logging.info(f"Removed {cleaned_counts['db_files']} database files from {db_dir}")

        # Clean summary files
        if clean_summary:
            summary_dir = self.results_dir / "SUMMARY"
            if summary_dir.exists():
                summary_files = list(summary_dir.glob("*.tsv")) + list(summary_dir.glob("*.json"))
                cleaned_counts['summary'] = len(summary_files)
                for f in summary_files:
                    f.unlink()
                logging.info(f"Removed {cleaned_counts['summary']} summary files from {summary_dir}")

        total_cleaned = sum(cleaned_counts.values())
        logging.info(f"Cleanup complete: {total_cleaned} total files removed")
        logging.info(f"  - A3M: {cleaned_counts['a3m']}, HHM: {cleaned_counts['hhm']}, "
                    f"HHR: {cleaned_counts['hhr']}, Parsed: {cleaned_counts['parsed']}, "
                    f"DB files: {cleaned_counts['db_files']}, Summary: {cleaned_counts['summary']}")

        return cleaned_counts

    def run_hhsearch(self, query_a3m, database, output_hhr, num_threads=8):
        """
        Run HHsearch profile-profile comparison.

        Args:
            query_a3m: Path to query A3M alignment file
            database: Path to HHM database base name (without suffix)
            output_hhr: Path to output HHR file
            num_threads: Number of CPU threads

        Returns:
            True if successful, False otherwise
        """
        try:
            # HHsearch command for profile-profile comparison
            # Uses A3M alignment as query against HHM database
            cmd = (f"hhsearch -i {query_a3m} -d {database} -o {output_hhr} "
                   f"-cpu {num_threads} -e {self.e_value_threshold}")

            result = subprocess.run(cmd, shell=True, check=True,
                                  capture_output=True, text=True)
            logging.debug(f"HHsearch completed for {query_a3m}")
            return True

        except subprocess.CalledProcessError as e:
            logging.error(f"HHsearch failed for {query_a3m}: {e}")
            return False

    def parse_hhr(self, hhr_file):
        """
        Parse HHR output file to extract hits.

        Args:
            hhr_file: Path to HHR file

        Returns:
            List of hit dictionaries
        """
        hits = []

        try:
            with open(hhr_file, 'r') as f:
                lines = f.readlines()

            # Find the hits table section
            in_table = False
            for line in lines:
                if line.startswith(' No Hit'):
                    in_table = True
                    continue

                if in_table:
                    if line.strip() == '' or line.startswith('No '):
                        continue
                    if line.startswith('Done'):
                        break

                    # Parse hit line
                    # Format: No Hit Prob E-value P-value Score SS Cols Query HMM Template HMM
                    match = re.match(r'\s*(\d+)\s+(\S+)\s+([\d.]+)\s+([\d.E-]+)\s+([\d.E-]+)\s+([\d.]+)\s+([\d.]+)\s+(\d+)\s+([\d-]+)\s+([\d-]+)\s+\((\d+)\)', line)

                    if match:
                        hit = {
                            'rank': int(match.group(1)),
                            'target': match.group(2),
                            'probability': float(match.group(3)),
                            'e_value': float(match.group(4)),
                            'p_value': float(match.group(5)),
                            'score': float(match.group(6)),
                            'ss': float(match.group(7)),
                            'cols': int(match.group(8)),
                            'query_range': match.group(9),
                            'template_range': match.group(10),
                            'template_length': int(match.group(11))
                        }
                        hits.append(hit)

            logging.debug(f"Parsed {len(hits)} hits from {hhr_file}")
            return hits

        except Exception as e:
            logging.error(f"Failed to parse {hhr_file}: {e}")
            return []

    def run_search_for_family(self, family_id, database):
        """
        Run complete HHsearch for one family.

        Args:
            family_id: Pfam family ID (e.g., PF00001)
            database: Path to HHM database base name (e.g., /path/to/pfam)

        Returns:
            List of hits
        """
        # Get A3M query file (HHsearch uses A3M as query, not HHM)
        a3m_file = self.a3m_dir / f"{family_id}.a3m"

        if not a3m_file.exists():
            logging.warning(f"A3M file not found: {a3m_file}")
            return []

        # Run HHsearch
        hhr_file = self.hhr_dir / f"{family_id}.hhr"
        if not self.run_hhsearch(a3m_file, database, hhr_file):
            return []

        # Parse results
        hits = self.parse_hhr(hhr_file)

        # Save parsed hits
        parsed_file = self.parsed_dir / f"{family_id}_hits.tsv"
        self.save_hits(family_id, hits, parsed_file)

        return hits

    def save_hits(self, query_family, hits, output_file):
        """Save parsed hits to TSV file."""
        try:
            with open(output_file, 'w') as f:
                # Header
                f.write("query_family\ttarget_family\tprobability\te_value\tp_value\t"
                       "score\tss\tcols\tquery_range\ttemplate_range\ttemplate_length\n")

                # Hits
                for hit in hits:
                    f.write(f"{query_family}\t{hit['target']}\t{hit['probability']}\t"
                           f"{hit['e_value']}\t{hit['p_value']}\t{hit['score']}\t"
                           f"{hit['ss']}\t{hit['cols']}\t{hit['query_range']}\t"
                           f"{hit['template_range']}\t{hit['template_length']}\n")

            logging.debug(f"Saved {len(hits)} hits to {output_file}")

        except Exception as e:
            logging.error(f"Failed to save hits to {output_file}: {e}")

    def create_family_batches(self, families, batch_size=100):
        """
        Split families into batches for SLURM array jobs.

        Args:
            families: List of family IDs
            batch_size: Number of families per batch

        Returns:
            List of batches
        """
        batches = []
        for i in range(0, len(families), batch_size):
            batches.append(families[i:i + batch_size])

        logging.info(f"Created {len(batches)} batches of ~{batch_size} families each")
        return batches

    def generate_slurm_script(self, database, batch_dir, num_batches,
                             partition="standard", time_limit="02:00:00",
                             memory="16GB", cpus=8):
        """
        Generate SLURM array job script for batch processing.

        Args:
            database: Path to HHM database
            batch_dir: Directory containing batch files
            num_batches: Number of batches (for array size)
            partition: SLURM partition
            time_limit: Job time limit
            memory: Memory per job
            cpus: CPUs per job

        Returns:
            Path to generated SLURM script
        """
        script_path = self.results_dir / "run_hhsearch_batch.sh"

        # Get the directory containing the pipeline scripts
        pipeline_dir = Path(__file__).parent.absolute()

        script_content = f"""#!/bin/bash
#SBATCH --job-name=pfam_hhsearch
#SBATCH --array=1-{num_batches}
#SBATCH --partition={partition}
#SBATCH --time={time_limit}
#SBATCH --mem={memory}
#SBATCH --cpus-per-task={cpus}
#SBATCH --output={self.results_dir}/slurm_logs/hhsearch_%A_%a.out
#SBATCH --error={self.results_dir}/slurm_logs/hhsearch_%A_%a.err

# Add pipeline directory to PYTHONPATH
export PYTHONPATH={pipeline_dir}:$PYTHONPATH

# Load modules if needed
# module load hhsuite

# Get batch file for this array task
BATCH_FILE={batch_dir}/batch_${{SLURM_ARRAY_TASK_ID}}.txt

# Read families from batch file
while read FAMILY_ID; do
    echo "Processing $FAMILY_ID..."

    A3M_FILE={self.a3m_dir}/${{FAMILY_ID}}.a3m
    HHR_FILE={self.hhr_dir}/${{FAMILY_ID}}.hhr

    # Skip if already processed
    if [ -f "{self.parsed_dir}/${{FAMILY_ID}}_hits.tsv" ]; then
        echo "  Already processed, skipping..."
        continue
    fi

    # Run HHsearch (uses A3M alignment as query against HHM database)
    hhsearch -i $A3M_FILE -d {database} -o $HHR_FILE \\
        -cpu {cpus} -e {self.e_value_threshold}

    # Parse results
    python3 -c "
from hhsearch_runner import HHsearchRunner
runner = HHsearchRunner('{self.seed_dir}', '{self.results_dir}', {self.e_value_threshold})
hits = runner.parse_hhr('$HHR_FILE')
parsed_file = '{self.parsed_dir}/${{FAMILY_ID}}_hits.tsv'
runner.save_hits('$FAMILY_ID', hits, parsed_file)
"

done < $BATCH_FILE

echo "Batch ${{SLURM_ARRAY_TASK_ID}} complete"
"""

        with open(script_path, 'w') as f:
            f.write(script_content)

        # Make executable
        script_path.chmod(0o755)

        logging.info(f"SLURM script created: {script_path}")
        return script_path

    def setup_slurm_batches(self, families, database, batch_size=100):
        """
        Set up SLURM batch processing infrastructure.

        Args:
            families: List of family IDs
            database: Path to HHM database
            batch_size: Families per batch

        Returns:
            Tuple of (slurm_script_path, num_batches)
        """
        # Create batch directory
        batch_dir = self.results_dir / "batches"
        batch_dir.mkdir(exist_ok=True)

        # Create SLURM logs directory
        (self.results_dir / "slurm_logs").mkdir(exist_ok=True)

        # Create batches
        batches = self.create_family_batches(families, batch_size)

        # Write batch files
        for i, batch in enumerate(batches, 1):
            batch_file = batch_dir / f"batch_{i}.txt"
            with open(batch_file, 'w') as f:
                f.write('\n'.join(batch))

        # Generate SLURM script
        slurm_script = self.generate_slurm_script(database, batch_dir, len(batches))

        logging.info(f"SLURM batch setup complete: {len(batches)} batches, {len(families)} families")
        return slurm_script, len(batches)
