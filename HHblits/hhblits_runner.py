#!/usr/bin/env python3
"""
Simple HHblits runner for all-against-all Pfam family comparisons.
"""

import subprocess
import os
import logging
from pathlib import Path
import re

class HHblitsRunner:
    """Manages HHblits searches and result parsing."""

    def __init__(self, seed_dir, hmm_dir, results_dir, e_value_threshold=1e-10):
        """
        Initialize HHblits runner.

        Args:
            seed_dir: Directory containing SEED alignment files
            hmm_dir: Directory containing HMM files
            results_dir: Directory for HHblits results
            e_value_threshold: E-value cutoff for significant hits
        """
        self.seed_dir = Path(seed_dir)
        self.hmm_dir = Path(hmm_dir)
        self.results_dir = Path(results_dir)
        self.e_value_threshold = e_value_threshold

        # Create subdirectories
        self.a3m_dir = self.results_dir / "A3M"
        self.hhr_dir = self.results_dir / "HHR"
        self.parsed_dir = self.results_dir / "PARSED"

        for dir_path in [self.a3m_dir, self.hhr_dir, self.parsed_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)

        logging.info(f"HHblits runner initialized")
        logging.info(f"Results directory: {self.results_dir}")
        logging.info(f"E-value threshold: {self.e_value_threshold}")

    def mul_to_a3m(self, seed_file, output_file):
        """
        Convert mul format SEED to A3M format for HHblits.
        Mul format is: one line per sequence, format "seq_id<whitespace>sequence"

        Args:
            seed_file: Path to mul format SEED file
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
                    # Split on whitespace - first part is ID, rest is sequence
                    parts = line.split(None, 1)
                    if len(parts) >= 2:
                        seq_id = parts[0]
                        seq = parts[1]
                        sequences[seq_id] = seq

            if not sequences:
                logging.warning(f"No sequences found in {seed_file}")
                return False

            # Write A3M format (FASTA-like)
            with open(output_file, 'w') as f:
                for seq_id, seq in sequences.items():
                    f.write(f">{seq_id}\n{seq}\n")

            logging.debug(f"Converted {seed_file} to A3M format ({len(sequences)} sequences)")
            return True

        except Exception as e:
            logging.error(f"Failed to convert {seed_file} to A3M: {e}")
            return False

    def build_hhm_database(self):
        """
        Build HH-suite database from SEED alignments using hhmake.
        Note: Pfam HMM files are in HMMER format, but HHblits needs HH-suite format.
        We build HH-suite profiles from SEED alignments instead.

        Returns:
            Path to database directory (for use with hhblits -d)
        """
        db_dir = self.results_dir / "hhsuite_db"
        db_dir.mkdir(exist_ok=True)

        try:
            seed_files = sorted(self.seed_dir.glob("PF*_SEED"))
            logging.info(f"Building HH-suite database from {len(seed_files)} SEED alignments...")
            logging.info("This may take some time (converting alignments to HH-suite profiles)...")

            # Build HHsuite profiles from SEED alignments
            hhm_dir = db_dir / "hhm_files"
            hhm_dir.mkdir(exist_ok=True)

            failed = 0
            for i, seed_file in enumerate(seed_files, 1):
                if i % 100 == 0:
                    logging.info(f"Progress: {i}/{len(seed_files)} profiles built")

                family_id = seed_file.stem.replace('_SEED', '')
                a3m_file = self.a3m_dir / f"{family_id}.a3m"
                hhm_file = hhm_dir / f"{family_id}.hhm"

                # Skip if already built
                if hhm_file.exists():
                    continue

                # Convert SEED to A3M if not already done
                if not a3m_file.exists():
                    if not self.mul_to_a3m(seed_file, a3m_file):
                        failed += 1
                        continue

                # Build HH-suite profile with hhmake
                # Use -M first to handle alignment column variations more permissively
                cmd = f"hhmake -i {a3m_file} -o {hhm_file} -M first"
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

                if result.returncode != 0:
                    logging.warning(f"hhmake failed for {family_id}: {result.stderr}")
                    failed += 1

            logging.info(f"Profile building complete: {len(seed_files) - failed} succeeded, {failed} failed")

            # Build ffindex database from HHM files
            ffdata = db_dir / "pfam_db_hhm.ffdata"
            ffindex = db_dir / "pfam_db_hhm.ffindex"

            logging.info("Building ffindex database from HH-suite profiles...")
            cmd = f"ffindex_build -s {ffdata} {ffindex} {hhm_dir}/"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

            if result.returncode != 0:
                logging.error(f"ffindex_build failed: {result.stderr}")
                logging.info("Falling back to concatenated database...")
                return self._build_concatenated_hhm_database(hhm_dir, db_dir)

            logging.info(f"HH-suite database created: {db_dir}/pfam_db")
            return str(db_dir / "pfam_db")

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

    def run_hhblits(self, query_a3m, database, output_hhr, num_threads=8, max_iterations=2):
        """
        Run HHblits search.

        Args:
            query_a3m: Path to query alignment in A3M format
            database: Path to HHM database
            output_hhr: Path to output HHR file
            num_threads: Number of CPU threads
            max_iterations: Number of HHblits iterations

        Returns:
            True if successful, False otherwise
        """
        try:
            cmd = (f"hhblits -i {query_a3m} -d {database} -o {output_hhr} "
                   f"-cpu {num_threads} -n {max_iterations} -e {self.e_value_threshold}")

            result = subprocess.run(cmd, shell=True, check=True,
                                  capture_output=True, text=True)
            logging.debug(f"HHblits search completed for {query_a3m}")
            return True

        except subprocess.CalledProcessError as e:
            logging.error(f"HHblits failed for {query_a3m}: {e}")
            return False

    def parse_hhr(self, hhr_file):
        """
        Parse HHblits HHR output file.

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
        Run complete HHblits search for one family.

        Args:
            family_id: Pfam family ID (e.g., PF00001)
            database: Path to HHM database

        Returns:
            List of hits
        """
        # Convert SEED to A3M
        seed_file = self.seed_dir / f"{family_id}_SEED"
        a3m_file = self.a3m_dir / f"{family_id}.a3m"

        if not seed_file.exists():
            logging.warning(f"SEED file not found: {seed_file}")
            return []

        if not self.mul_to_a3m(seed_file, a3m_file):
            return []

        # Run HHblits
        hhr_file = self.hhr_dir / f"{family_id}.hhr"
        if not self.run_hhblits(a3m_file, database, hhr_file):
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
        script_path = self.results_dir / "run_hhblits_batch.sh"

        # Get the directory containing the pipeline scripts
        pipeline_dir = Path(__file__).parent.absolute()

        script_content = f"""#!/bin/bash
#SBATCH --job-name=pfam_hhblits
#SBATCH --array=1-{num_batches}
#SBATCH --partition={partition}
#SBATCH --time={time_limit}
#SBATCH --mem={memory}
#SBATCH --cpus-per-task={cpus}
#SBATCH --output={self.results_dir}/slurm_logs/hhblits_%A_%a.out
#SBATCH --error={self.results_dir}/slurm_logs/hhblits_%A_%a.err

# Add pipeline directory to PYTHONPATH
export PYTHONPATH={pipeline_dir}:$PYTHONPATH

# Load modules if needed
# module load hhsuite

# Get batch file for this array task
BATCH_FILE={batch_dir}/batch_${{SLURM_ARRAY_TASK_ID}}.txt

# Read families from batch file
while read FAMILY_ID; do
    echo "Processing $FAMILY_ID..."

    SEED_FILE={self.seed_dir}/${{FAMILY_ID}}_SEED
    A3M_FILE={self.a3m_dir}/${{FAMILY_ID}}.a3m
    HHR_FILE={self.hhr_dir}/${{FAMILY_ID}}.hhr

    # Skip if already processed
    if [ -f "{self.parsed_dir}/${{FAMILY_ID}}_hits.tsv" ]; then
        echo "  Already processed, skipping..."
        continue
    fi

    # Convert SEED to A3M
    python3 -c "
from hhblits_runner import HHblitsRunner
runner = HHblitsRunner('{self.seed_dir}', '{self.hmm_dir}', '{self.results_dir}')
runner.mul_to_a3m('$SEED_FILE', '$A3M_FILE')
"

    # Run HHblits (note: database path without extension)
    hhblits -i $A3M_FILE -d {database} -o $HHR_FILE \\
        -cpu {cpus} -n 2 -e {self.e_value_threshold}

    # Parse results
    python3 -c "
from hhblits_runner import HHblitsRunner
runner = HHblitsRunner('{self.seed_dir}', '{self.hmm_dir}', '{self.results_dir}')
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


if __name__ == "__main__":
    # Test the HHblits runner
    logging.basicConfig(level=logging.INFO,
                       format='%(asctime)s - %(levelname)s - %(message)s')

    seed_dir = "./test_output/SEED"
    hmm_dir = "./test_output/HMM"
    results_dir = "./test_results"

    runner = HHblitsRunner(seed_dir, hmm_dir, results_dir)

    # Test building database
    print("Building HHM database...")
    db_file = runner.build_hhm_database()
    print(f"Database: {db_file}")

    # Test running search for one family
    families = [f.stem.replace('_SEED', '') for f in Path(seed_dir).glob("PF*_SEED")]
    if families:
        test_family = families[0]
        print(f"\nTesting HHblits search for {test_family}...")
        hits = runner.run_search_for_family(test_family, db_file)
        print(f"Found {len(hits)} hits")
        if hits:
            print(f"Top hit: {hits[0]}")
