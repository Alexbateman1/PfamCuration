#!/usr/bin/env python3
"""
Main pipeline for Pfam HHsearch all-against-all profile comparison.

This pipeline:
1. Extracts SEED files from Pfam SVN repository (can reuse from HHblits)
2. Builds HH-suite HMM database
3. Runs all-against-all HHsearch profile comparisons
4. Aggregates results into summary files
5. Supports incremental updates via SVN revision tracking
"""

import argparse
import logging
import sys
import subprocess
import time
from pathlib import Path
import json
from datetime import datetime

from svn_manager import SVNManager
from hhsearch_runner import HHsearchRunner


class PfamHHsearchPipeline:
    """Main pipeline orchestrator for HHsearch all-against-all."""

    def __init__(self, work_dir, svn_url, e_value_threshold=1.0, incremental=True,
                 use_slurm=True, batch_size=100, poll_interval=60):
        """
        Initialize pipeline.

        Args:
            work_dir: Working directory for all pipeline files
            svn_url: Pfam SVN repository URL
            e_value_threshold: E-value cutoff for HHsearch (default: 1.0)
            incremental: Enable incremental updates
            use_slurm: Use SLURM for parallel processing
            batch_size: Number of families per SLURM batch
            poll_interval: Seconds between SLURM job status checks (default: 60)
        """
        self.work_dir = Path(work_dir)
        self.svn_url = svn_url
        self.e_value_threshold = e_value_threshold
        self.incremental = incremental
        self.use_slurm = use_slurm
        self.batch_size = batch_size
        self.poll_interval = poll_interval

        # Track SLURM job ID
        self.slurm_job_file = None

        # Create directory structure
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.progress_dir = self.work_dir / ".progress"
        self.progress_dir.mkdir(exist_ok=True)

        self.data_dir = self.work_dir / "DATA"
        self.results_dir = self.work_dir / "RESULTS"
        self.summary_dir = self.results_dir / "SUMMARY"

        for d in [self.data_dir, self.results_dir, self.summary_dir]:
            d.mkdir(parents=True, exist_ok=True)

        # Progress markers
        self.markers = {
            'svn_discovery': self.progress_dir / '.svn_discovery_complete',
            'extraction': self.progress_dir / '.extraction_complete',
            'db_build': self.progress_dir / '.db_build_complete',
            'hhsearch': self.progress_dir / '.hhsearch_complete',
            'aggregation': self.progress_dir / '.aggregation_complete'
        }

        # Revision tracking
        self.revision_file = self.work_dir / '.last_revision'

        logging.info(f"Pipeline initialized in {self.work_dir}")

    def is_step_complete(self, step):
        """Check if a pipeline step is complete."""
        return self.markers[step].exists()

    def mark_step_complete(self, step):
        """Mark a pipeline step as complete."""
        self.markers[step].touch()
        logging.info(f"Step '{step}' marked complete")

    def reset_progress(self, from_step=None, clean_files=True):
        """
        Reset progress markers from a specific step onwards.

        Args:
            from_step: Step to reset from (if None, reset all)
            clean_files: Also clean up result files associated with the steps being reset
        """
        steps = list(self.markers.keys())
        if from_step:
            idx = steps.index(from_step)
            steps = steps[idx:]

        # Clean up files before resetting markers
        if clean_files and steps:
            logging.info("Cleaning up old result files...")
            hhsearch_runner = HHsearchRunner(
                seed_dir=self.data_dir / "SEED",
                results_dir=self.results_dir,
                e_value_threshold=self.e_value_threshold
            )

            # Determine what to clean based on which steps are being reset
            # db_build creates: A3M, HHM, DB files
            # hhsearch creates: HHR, parsed TSV files
            # aggregation creates: summary files (uses parsed TSV as input)
            clean_a3m = 'db_build' in steps
            clean_hhm = 'db_build' in steps
            clean_db = 'db_build' in steps
            clean_hhr = 'hhsearch' in steps
            clean_parsed = 'hhsearch' in steps
            clean_summary = 'aggregation' in steps

            hhsearch_runner.clean_results(
                clean_a3m=clean_a3m,
                clean_hhm=clean_hhm,
                clean_hhr=clean_hhr,
                clean_parsed=clean_parsed,
                clean_db=clean_db,
                clean_summary=clean_summary
            )

        for step in steps:
            if self.markers[step].exists():
                self.markers[step].unlink()
                logging.info(f"Reset progress marker for '{step}'")

    def get_last_revision(self):
        """Get the last processed SVN revision."""
        if self.revision_file.exists():
            with open(self.revision_file, 'r') as f:
                return f.read().strip()
        return None

    def save_current_revision(self, revision):
        """Save current SVN revision."""
        with open(self.revision_file, 'w') as f:
            f.write(str(revision))
        logging.info(f"Saved current revision: {revision}")

    def run_svn_discovery(self, svn_manager):
        """Step 1: Discover families in SVN and track revision."""
        if self.is_step_complete('svn_discovery'):
            logging.info("SVN discovery already complete, skipping...")
            return None

        logging.info("=== Step 1: SVN Discovery ===")

        # Get current revision
        current_revision = svn_manager.get_current_revision()
        last_revision = self.get_last_revision()

        # Determine which families to process
        if self.incremental and last_revision:
            logging.info(f"Incremental mode: checking changes since revision {last_revision}")
            families = svn_manager.get_changed_families(last_revision, current_revision)
            if not families:
                logging.info("No families changed since last run")
                return []
        else:
            logging.info("Full mode: processing all families")
            families = svn_manager.list_all_families()

        # Save family list
        family_list_file = self.work_dir / 'families_to_process.txt'
        with open(family_list_file, 'w') as f:
            f.write('\n'.join(families))

        logging.info(f"Families to process: {len(families)}")
        self.save_current_revision(current_revision)
        self.mark_step_complete('svn_discovery')

        return families

    def run_extraction(self, svn_manager, families):
        """Step 2: Extract SEED files."""
        if self.is_step_complete('extraction'):
            logging.info("Extraction already complete, skipping...")
            return

        logging.info("=== Step 2: SEED Extraction ===")

        results = svn_manager.extract_all_families(families)

        # Save extraction manifest
        manifest_file = self.work_dir / 'extraction_manifest.tsv'
        with open(manifest_file, 'w') as f:
            f.write("family_id\tseed_path\tseed_exists\n")
            for family_id, (seed_path, _) in results.items():
                seed_exists = seed_path is not None and Path(seed_path).exists()
                f.write(f"{family_id}\t{seed_path}\t{seed_exists}\n")

        logging.info(f"Extraction complete: {len(results)} families")
        self.mark_step_complete('extraction')

    def run_db_build(self, hhsearch_runner):
        """Step 3: Build HHsearch database."""
        if self.is_step_complete('db_build'):
            logging.info("Database build already complete, skipping...")
            db_file = self.results_dir / "hhsuite_db" / "pfam_db_hhm"
            return str(db_file)

        logging.info("=== Step 3: HH-suite Database Build ===")

        db_file = hhsearch_runner.build_hhm_database()

        self.mark_step_complete('db_build')
        return db_file

    def submit_slurm_job(self, slurm_script):
        """
        Submit SLURM batch job and return job ID.

        Args:
            slurm_script: Path to SLURM script

        Returns:
            Job ID as string
        """
        logging.info(f"Submitting SLURM job: {slurm_script}")

        result = subprocess.run(
            ['sbatch', str(slurm_script)],
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            logging.error(f"sbatch failed: {result.stderr}")
            raise RuntimeError(f"Failed to submit SLURM job: {result.stderr}")

        # Parse job ID from output like "Submitted batch job 12345678"
        output = result.stdout.strip()
        job_id = output.split()[-1]

        logging.info(f"SLURM job submitted: {job_id}")

        # Save job ID to file for recovery
        job_id_file = self.progress_dir / '.slurm_job_id'
        with open(job_id_file, 'w') as f:
            f.write(job_id)

        return job_id

    def get_slurm_job_status(self, job_id):
        """
        Check status of SLURM job array.

        Args:
            job_id: SLURM job ID

        Returns:
            Dict with 'running', 'pending', 'completed', 'failed' counts
        """
        # Use sacct to get job status (works even after jobs complete)
        result = subprocess.run(
            ['sacct', '-j', job_id, '--format=JobID,State', '--noheader', '--parsable2'],
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            # Fallback to squeue for running jobs
            result = subprocess.run(
                ['squeue', '-j', job_id, '--format=%T', '--noheader'],
                capture_output=True,
                text=True
            )

            if result.returncode != 0:
                logging.warning(f"Could not get job status for {job_id}")
                return None

            states = result.stdout.strip().split('\n')
            return {
                'running': sum(1 for s in states if s == 'RUNNING'),
                'pending': sum(1 for s in states if s == 'PENDING'),
                'completed': 0,
                'failed': 0,
                'total': len([s for s in states if s])
            }

        # Parse sacct output
        status = {
            'running': 0,
            'pending': 0,
            'completed': 0,
            'failed': 0,
            'total': 0
        }

        for line in result.stdout.strip().split('\n'):
            if not line or '_' in line.split('|')[0]:  # Skip sub-jobs (batch, extern)
                continue

            parts = line.split('|')
            if len(parts) >= 2:
                state = parts[1]
                status['total'] += 1

                if state in ('RUNNING',):
                    status['running'] += 1
                elif state in ('PENDING', 'REQUEUED'):
                    status['pending'] += 1
                elif state in ('COMPLETED',):
                    status['completed'] += 1
                elif state in ('FAILED', 'CANCELLED', 'TIMEOUT', 'OUT_OF_MEMORY', 'NODE_FAIL'):
                    status['failed'] += 1

        return status

    def wait_for_slurm_jobs(self, job_id, num_batches):
        """
        Wait for all SLURM jobs to complete.

        Args:
            job_id: SLURM job ID
            num_batches: Expected number of array tasks

        Returns:
            True if all jobs completed successfully, False if any failed
        """
        logging.info(f"Waiting for SLURM jobs to complete (job ID: {job_id}, {num_batches} tasks)...")
        logging.info(f"Polling every {self.poll_interval} seconds...")

        start_time = time.time()
        last_log_time = 0

        while True:
            status = self.get_slurm_job_status(job_id)

            if status is None:
                logging.warning("Could not get job status, will retry...")
                time.sleep(self.poll_interval)
                continue

            running = status['running']
            pending = status['pending']
            completed = status['completed']
            failed = status['failed']
            total = status['total']

            # Log progress every 5 minutes or on significant change
            current_time = time.time()
            if current_time - last_log_time >= 300:  # 5 minutes
                elapsed = (current_time - start_time) / 60
                logging.info(f"SLURM status after {elapsed:.1f} min: "
                           f"{completed} completed, {running} running, "
                           f"{pending} pending, {failed} failed (of {num_batches} batches)")
                last_log_time = current_time

            # Check if all jobs are done
            if running == 0 and pending == 0:
                elapsed = (time.time() - start_time) / 60
                logging.info(f"All SLURM jobs finished after {elapsed:.1f} minutes")
                logging.info(f"Final status: {completed} completed, {failed} failed")

                if failed > 0:
                    logging.warning(f"{failed} jobs failed - check SLURM logs in RESULTS/slurm_logs/")
                    return False

                return True

            time.sleep(self.poll_interval)

    def run_hhsearch_searches(self, hhsearch_runner, families, database, use_slurm=True, batch_size=100):
        """Step 4: Run all-against-all HHsearch comparisons."""
        if self.is_step_complete('hhsearch'):
            logging.info("HHsearch searches already complete, skipping...")
            return

        logging.info("=== Step 4: HHsearch All-Against-All Profile Comparisons ===")

        total = len(families)
        logging.info(f"Running HHsearch for {total} families...")

        if use_slurm:
            # Check if we have a pending job from a previous run
            job_id_file = self.progress_dir / '.slurm_job_id'
            num_batches_file = self.progress_dir / '.slurm_num_batches'

            if job_id_file.exists() and num_batches_file.exists():
                # Resume monitoring existing job
                with open(job_id_file, 'r') as f:
                    job_id = f.read().strip()
                with open(num_batches_file, 'r') as f:
                    num_batches = int(f.read().strip())

                logging.info(f"Found existing SLURM job: {job_id}")

                # Check if job is still running
                status = self.get_slurm_job_status(job_id)
                if status and (status['running'] > 0 or status['pending'] > 0):
                    logging.info("Job still running, resuming monitoring...")
                else:
                    logging.info("Previous job completed, will set up new batch")
                    job_id = None
            else:
                job_id = None

            if job_id is None:
                # Set up SLURM batch processing
                slurm_script, num_batches = hhsearch_runner.setup_slurm_batches(
                    families, database, batch_size=batch_size
                )

                logging.info(f"SLURM batch setup complete: {num_batches} batches")

                # Save num_batches for recovery
                with open(num_batches_file, 'w') as f:
                    f.write(str(num_batches))

                # Submit the job
                job_id = self.submit_slurm_job(slurm_script)

            # Wait for jobs to complete
            success = self.wait_for_slurm_jobs(job_id, num_batches)

            # Clean up job tracking files
            if job_id_file.exists():
                job_id_file.unlink()
            if num_batches_file.exists():
                num_batches_file.unlink()

            if not success:
                logging.error("Some SLURM jobs failed. Check logs and re-run pipeline.")
                raise RuntimeError("SLURM jobs failed")

            logging.info(f"HHsearch searches complete: {total} families processed")
            self.mark_step_complete('hhsearch')

        else:
            # Sequential processing (for testing or small runs)
            logging.info("Running in sequential mode (no SLURM)...")
            all_hits = {}
            for i, family_id in enumerate(families, 1):
                if i % 100 == 0:
                    logging.info(f"Progress: {i}/{total} searches completed")

                hits = hhsearch_runner.run_search_for_family(family_id, database)
                all_hits[family_id] = hits

            logging.info(f"HHsearch searches complete: {total} families processed")
            self.mark_step_complete('hhsearch')

            return all_hits

    def build_sequence_to_family_map(self):
        """
        Build mapping from sequence IDs to Pfam families by reading A3M files.

        Returns:
            Dictionary mapping sequence IDs to Pfam family IDs
        """
        seq_to_family = {}
        a3m_dir = self.results_dir / "A3M"

        a3m_files = list(a3m_dir.glob("PF*.a3m"))
        logging.info(f"Building sequence-to-family map from {len(a3m_files)} A3M files...")

        for i, a3m_file in enumerate(a3m_files, 1):
            if i % 2000 == 0:
                logging.info(f"  Progress: {i}/{len(a3m_files)} A3M files processed")

            family_id = a3m_file.stem  # PF00001.a3m -> PF00001

            try:
                with open(a3m_file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith('>'):
                            # Extract sequence ID (remove '>')
                            seq_id = line[1:]
                            seq_to_family[seq_id] = family_id
                            # Also map without /start-end if present
                            if '/' in seq_id:
                                base_id = seq_id.split('/')[0]
                                seq_to_family[base_id] = family_id
            except Exception as e:
                logging.warning(f"Failed to read {a3m_file}: {e}")

        logging.info(f"Mapped {len(seq_to_family)} sequences to {len(a3m_files)} families")
        return seq_to_family

    def run_aggregation(self):
        """Step 5: Aggregate results into summary files with Pfam family mapping."""
        if self.is_step_complete('aggregation'):
            logging.info("Aggregation already complete, skipping...")
            return

        logging.info("=== Step 5: Result Aggregation ===")

        # Build sequence-to-family mapping for Pfam lookups
        seq_to_family = self.build_sequence_to_family_map()

        # Aggregate all parsed hits into single file with Pfam mapping
        parsed_dir = self.results_dir / "PARSED"
        all_hits_file = self.summary_dir / "hhsearch_all_vs_all.tsv"

        hit_count = 0
        mapped_count = 0
        unmapped_count = 0

        with open(all_hits_file, 'w') as outf:
            # Write header with target_pfam_family column
            outf.write("query_family\ttarget_family\ttarget_pfam_family\tprobability\te_value\tp_value\t"
                      "score\tss\tcols\tquery_range\ttemplate_range\ttemplate_length\n")

            # Aggregate all hit files
            for hit_file in sorted(parsed_dir.glob("PF*_hits.tsv")):
                with open(hit_file, 'r') as inf:
                    lines = inf.readlines()[1:]  # Skip header

                    for line in lines:
                        parts = line.strip().split('\t')
                        if len(parts) < 2:
                            continue

                        target_seq = parts[1]

                        # Look up Pfam family for target sequence
                        seq_id = target_seq.split('/')[0] if '/' in target_seq else target_seq
                        target_pfam = seq_to_family.get(target_seq, seq_to_family.get(seq_id, 'UNKNOWN'))

                        if target_pfam != 'UNKNOWN':
                            mapped_count += 1
                        else:
                            unmapped_count += 1

                        # Insert target_pfam_family after target_family
                        new_parts = [parts[0], parts[1], target_pfam] + parts[2:]
                        outf.write('\t'.join(new_parts) + '\n')
                        hit_count += 1

        logging.info(f"Aggregated {hit_count} hits into {all_hits_file}")
        logging.info(f"Pfam mapping: {mapped_count} mapped, {unmapped_count} unmapped")

        # Generate statistics
        self.generate_statistics(all_hits_file)

        self.mark_step_complete('aggregation')

    def generate_statistics(self, all_hits_file):
        """Generate summary statistics from aggregated hits."""
        stats_file = self.summary_dir / "family_statistics.tsv"
        high_conf_file = self.summary_dir / "high_confidence_overlaps.tsv"

        family_stats = {}
        high_conf_hits = []

        # Parse all hits
        # Column indices (with target_pfam_family at index 2):
        # 0: query_family, 1: target_family, 2: target_pfam_family,
        # 3: probability, 4: e_value, 5: p_value, 6: score, ...
        with open(all_hits_file, 'r') as f:
            lines = f.readlines()[1:]  # Skip header

            for line in lines:
                parts = line.strip().split('\t')
                if len(parts) < 12:  # Now 12 columns with target_pfam_family
                    continue

                query = parts[0]
                target_pfam = parts[2]  # Use Pfam family ID, not sequence ID
                probability = float(parts[3])
                e_value = float(parts[4])
                score = float(parts[6])

                # Update family statistics
                if query not in family_stats:
                    family_stats[query] = {
                        'num_hits': 0,
                        'best_e_value': float('inf'),
                        'best_target': None,
                        'best_score': 0
                    }

                family_stats[query]['num_hits'] += 1
                if e_value < family_stats[query]['best_e_value']:
                    family_stats[query]['best_e_value'] = e_value
                    family_stats[query]['best_target'] = target_pfam
                    family_stats[query]['best_score'] = score

                # Collect high confidence hits (e.g., e-value < 1e-10)
                if e_value < 1e-10:
                    high_conf_hits.append(line)

        # Write family statistics
        with open(stats_file, 'w') as f:
            f.write("family\tnum_hits\tbest_e_value\tbest_target\tbest_score\n")
            for family, stats in sorted(family_stats.items()):
                f.write(f"{family}\t{stats['num_hits']}\t{stats['best_e_value']}\t"
                       f"{stats['best_target']}\t{stats['best_score']}\n")

        logging.info(f"Family statistics saved to {stats_file}")

        # Write high confidence overlaps
        with open(high_conf_file, 'w') as f:
            f.write("query_family\ttarget_family\ttarget_pfam_family\tprobability\te_value\tp_value\t"
                   "score\tss\tcols\tquery_range\ttemplate_range\ttemplate_length\n")
            f.writelines(high_conf_hits)

        logging.info(f"High confidence overlaps saved to {high_conf_file} ({len(high_conf_hits)} hits)")

        # Save pipeline metadata
        metadata = {
            'run_date': datetime.now().isoformat(),
            'svn_revision': self.get_last_revision(),
            'e_value_threshold': self.e_value_threshold,
            'total_families': len(family_stats),
            'total_hits': len(lines),
            'high_confidence_hits': len(high_conf_hits)
        }

        metadata_file = self.summary_dir / "pipeline_metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)

        logging.info(f"Pipeline metadata saved to {metadata_file}")

    def run(self):
        """Run the complete pipeline."""
        logging.info("=" * 60)
        logging.info("Pfam HHsearch All-Against-All Pipeline")
        logging.info("=" * 60)

        # Initialize managers
        svn_manager = SVNManager(self.svn_url, self.data_dir)
        hhsearch_runner = HHsearchRunner(
            seed_dir=self.data_dir / "SEED",
            results_dir=self.results_dir,
            e_value_threshold=self.e_value_threshold
        )

        # Run pipeline steps
        families = self.run_svn_discovery(svn_manager)

        if families is None:
            # Load from file
            family_list_file = self.work_dir / 'families_to_process.txt'
            with open(family_list_file, 'r') as f:
                families = [line.strip() for line in f if line.strip()]

        if not families:
            logging.info("No families to process. Pipeline complete.")
            return

        self.run_extraction(svn_manager, families)

        database = self.run_db_build(hhsearch_runner)

        self.run_hhsearch_searches(hhsearch_runner, families, database,
                                   use_slurm=self.use_slurm, batch_size=self.batch_size)

        self.run_aggregation()

        logging.info("=" * 60)
        logging.info("Pipeline complete!")
        logging.info(f"Results available in: {self.summary_dir}")
        logging.info("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="Pfam HHsearch all-against-all profile comparison pipeline"
    )
    parser.add_argument(
        '--work-dir', '-w',
        default='./HHSEARCH_PIPELINE',
        help='Working directory for pipeline (default: ./HHSEARCH_PIPELINE)'
    )
    parser.add_argument(
        '--svn-url',
        default='https://xfam-svn-hl.ebi.ac.uk/svn/pfam/',
        help='Pfam SVN repository URL'
    )
    parser.add_argument(
        '--e-value', '-e',
        type=float,
        default=1.0,
        help='E-value threshold for HHsearch (default: 1.0)'
    )
    parser.add_argument(
        '--full',
        action='store_true',
        help='Force full run (disable incremental mode)'
    )
    parser.add_argument(
        '--no-slurm',
        action='store_true',
        help='Run sequentially without SLURM (for testing)'
    )
    parser.add_argument(
        '--batch-size',
        type=int,
        default=100,
        help='Number of families per SLURM batch (default: 100)'
    )
    parser.add_argument(
        '--reset',
        help='Reset progress from specified step (svn_discovery, extraction, db_build, hhsearch, aggregation)'
    )
    parser.add_argument(
        '--no-clean',
        action='store_true',
        help='Skip cleaning result files when using --reset (default: clean files)'
    )
    parser.add_argument(
        '--poll-interval',
        type=int,
        default=60,
        help='Seconds between SLURM job status checks (default: 60)'
    )
    parser.add_argument(
        '--log-level',
        default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        help='Logging level (default: INFO)'
    )

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(Path(args.work_dir) / 'pipeline.log'),
            logging.StreamHandler(sys.stdout)
        ]
    )

    # Initialize pipeline
    pipeline = PfamHHsearchPipeline(
        work_dir=args.work_dir,
        svn_url=args.svn_url,
        e_value_threshold=args.e_value,
        incremental=not args.full,
        use_slurm=not args.no_slurm,
        batch_size=args.batch_size,
        poll_interval=args.poll_interval
    )

    # Reset if requested
    if args.reset:
        pipeline.reset_progress(from_step=args.reset, clean_files=not args.no_clean)

    # Run pipeline
    try:
        pipeline.run()
    except Exception as e:
        logging.error(f"Pipeline failed: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
