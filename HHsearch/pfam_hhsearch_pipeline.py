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
from pathlib import Path
import json
from datetime import datetime

from svn_manager import SVNManager
from hhsearch_runner import HHsearchRunner


class PfamHHsearchPipeline:
    """Main pipeline orchestrator for HHsearch all-against-all."""

    def __init__(self, work_dir, svn_url, e_value_threshold=1.0, incremental=True,
                 use_slurm=True, batch_size=100):
        """
        Initialize pipeline.

        Args:
            work_dir: Working directory for all pipeline files
            svn_url: Pfam SVN repository URL
            e_value_threshold: E-value cutoff for HHsearch (default: 1.0)
            incremental: Enable incremental updates
            use_slurm: Use SLURM for parallel processing
            batch_size: Number of families per SLURM batch
        """
        self.work_dir = Path(work_dir)
        self.svn_url = svn_url
        self.e_value_threshold = e_value_threshold
        self.incremental = incremental
        self.use_slurm = use_slurm
        self.batch_size = batch_size

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
                hmm_dir=self.data_dir / "HMM",
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
            f.write("family_id\tseed_path\thmm_path\tseed_exists\thmm_exists\n")
            for family_id, (seed_path, hmm_path) in results.items():
                seed_exists = seed_path is not None and seed_path.exists()
                hmm_exists = hmm_path is not None and hmm_path.exists()
                f.write(f"{family_id}\t{seed_path}\t{hmm_path}\t{seed_exists}\t{hmm_exists}\n")

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

    def run_hhsearch_searches(self, hhsearch_runner, families, database, use_slurm=True, batch_size=100):
        """Step 4: Run all-against-all HHsearch comparisons."""
        if self.is_step_complete('hhsearch'):
            logging.info("HHsearch searches already complete, skipping...")
            return

        logging.info("=== Step 4: HHsearch All-Against-All Profile Comparisons ===")

        total = len(families)
        logging.info(f"Running HHsearch for {total} families...")

        if use_slurm:
            # Set up SLURM batch processing
            slurm_script, num_batches = hhsearch_runner.setup_slurm_batches(
                families, database, batch_size=batch_size
            )

            logging.info(f"SLURM batch setup complete: {num_batches} batches")
            logging.info(f"To submit jobs, run:")
            logging.info(f"  sbatch {slurm_script}")
            logging.info(f"")
            logging.info(f"Monitor with: squeue -u $USER")
            logging.info(f"")
            logging.info(f"When complete, re-run pipeline to aggregate results")

            # Don't mark as complete - user needs to run SLURM jobs
            logging.info("NOTE: Pipeline paused. Submit SLURM jobs, then re-run to continue.")

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

    def run_aggregation(self):
        """Step 5: Aggregate results into summary files."""
        if self.is_step_complete('aggregation'):
            logging.info("Aggregation already complete, skipping...")
            return

        logging.info("=== Step 5: Result Aggregation ===")

        # Aggregate all parsed hits into single file
        parsed_dir = self.results_dir / "PARSED"
        all_hits_file = self.summary_dir / "hhsearch_all_vs_all.tsv"

        hit_count = 0
        with open(all_hits_file, 'w') as outf:
            # Write header
            outf.write("query_family\ttarget_family\tprobability\te_value\tp_value\t"
                      "score\tss\tcols\tquery_range\ttemplate_range\ttemplate_length\n")

            # Aggregate all hit files
            for hit_file in sorted(parsed_dir.glob("PF*_hits.tsv")):
                with open(hit_file, 'r') as inf:
                    lines = inf.readlines()[1:]  # Skip header
                    outf.writelines(lines)
                    hit_count += len(lines)

        logging.info(f"Aggregated {hit_count} hits into {all_hits_file}")

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
        with open(all_hits_file, 'r') as f:
            lines = f.readlines()[1:]  # Skip header

            for line in lines:
                parts = line.strip().split('\t')
                if len(parts) < 11:
                    continue

                query = parts[0]
                target = parts[1]
                probability = float(parts[2])
                e_value = float(parts[3])
                score = float(parts[5])

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
                    family_stats[query]['best_target'] = target
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
            f.write("query_family\ttarget_family\tprobability\te_value\tp_value\t"
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
            hmm_dir=self.data_dir / "HMM",
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
        batch_size=args.batch_size
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
