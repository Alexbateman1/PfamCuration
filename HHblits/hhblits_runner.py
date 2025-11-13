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

    def stockholm_to_a3m(self, seed_file, output_file):
        """
        Convert Stockholm format SEED to A3M format for HHblits.

        Args:
            seed_file: Path to Stockholm format SEED file
            output_file: Path to output A3M file

        Returns:
            True if successful, False otherwise
        """
        try:
            # Simple conversion: extract sequences, remove gaps in first sequence
            sequences = {}
            with open(seed_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#') or line.startswith('//'):
                        continue
                    parts = line.split()
                    if len(parts) >= 2:
                        seq_id = parts[0]
                        seq = parts[1]
                        sequences[seq_id] = seq

            if not sequences:
                logging.warning(f"No sequences found in {seed_file}")
                return False

            # Write A3M format
            with open(output_file, 'w') as f:
                for seq_id, seq in sequences.items():
                    f.write(f">{seq_id}\n{seq}\n")

            logging.debug(f"Converted {seed_file} to A3M format")
            return True

        except Exception as e:
            logging.error(f"Failed to convert {seed_file} to A3M: {e}")
            return False

    def build_hhm_database(self):
        """
        Build HHblits database from HMM files.
        This concatenates all HMMs into a single database.

        Returns:
            Path to database file
        """
        db_file = self.results_dir / "pfam_hhm_db"

        try:
            # Concatenate all HMM files
            hmm_files = sorted(self.hmm_dir.glob("PF*.hmm"))
            logging.info(f"Building HHM database from {len(hmm_files)} HMM files...")

            with open(db_file, 'w') as outf:
                for hmm_file in hmm_files:
                    with open(hmm_file, 'r') as inf:
                        outf.write(inf.read())
                        outf.write('\n')

            logging.info(f"HHM database created: {db_file}")
            return db_file

        except Exception as e:
            logging.error(f"Failed to build HHM database: {e}")
            raise

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

        if not self.stockholm_to_a3m(seed_file, a3m_file):
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
