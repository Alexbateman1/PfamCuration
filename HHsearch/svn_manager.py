#!/usr/bin/env python3
"""
SVN manager for Pfam SEED alignments using checkout/update workflow.

This module handles SVN operations for extracting and updating SEED alignments
from the Pfam SVN repository. It uses svn checkout/update for efficient
incremental updates rather than individual file exports.
"""

import subprocess
import logging
import shutil
from pathlib import Path


class SVNManager:
    """Manages SVN operations for Pfam SEED extraction."""

    def __init__(self, repo_url, data_dir):
        """
        Initialize SVN manager.

        Args:
            repo_url: Base SVN repository URL (e.g., https://xfam-svn-hl.ebi.ac.uk/svn/pfam/)
            data_dir: Base directory for data files
        """
        self.repo_url = repo_url.rstrip('/')
        self.data_dir = Path(data_dir)
        self.families_url = f"{self.repo_url}/trunk/Data/Families"

        # SVN working copy directory (contains nested family directories)
        self.checkout_dir = self.data_dir / "Families"

        # Flat SEED directory for pipeline compatibility
        # Files are named PF00001_SEED, PF00002_SEED, etc.
        self.seed_dir = self.data_dir / "SEED"
        self.seed_dir.mkdir(parents=True, exist_ok=True)

        logging.info(f"SVN Manager initialized: {self.families_url}")
        logging.info(f"Checkout directory: {self.checkout_dir}")
        logging.info(f"SEED directory: {self.seed_dir}")

    def is_checked_out(self):
        """Check if SVN working copy exists."""
        return (self.checkout_dir / ".svn").exists()

    def checkout(self):
        """
        Perform initial SVN checkout of the Families directory.

        Note: This checks out the entire Families directory structure including
        files we don't use (HMM, DESC, etc.). This is a one-time cost and enables
        efficient incremental updates via 'svn update'.
        """
        logging.info(f"Performing initial checkout from {self.families_url}...")
        logging.info("This may take several hours for the initial checkout...")

        self.checkout_dir.parent.mkdir(parents=True, exist_ok=True)

        cmd = f"svn checkout {self.families_url} {self.checkout_dir}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode != 0:
            logging.error(f"Checkout failed: {result.stderr}")
            raise RuntimeError(f"SVN checkout failed: {result.stderr}")

        logging.info("Checkout complete")
        return result.stdout

    def update(self):
        """
        Update existing SVN working copy.

        Returns:
            SVN update output showing what changed
        """
        if not self.is_checked_out():
            raise RuntimeError("No checkout exists. Run checkout() first.")

        logging.info("Updating SVN working copy...")
        cmd = f"svn update {self.checkout_dir}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode != 0:
            logging.error(f"Update failed: {result.stderr}")
            raise RuntimeError(f"SVN update failed: {result.stderr}")

        logging.info("Update complete")

        # Log summary of changes
        if result.stdout:
            lines = result.stdout.strip().split('\n')
            updated = [l for l in lines if l.startswith('U ')]
            added = [l for l in lines if l.startswith('A ')]
            deleted = [l for l in lines if l.startswith('D ')]
            if updated or added or deleted:
                logging.info(f"SVN changes: {len(added)} added, {len(updated)} updated, {len(deleted)} deleted")

        return result.stdout

    def checkout_or_update(self):
        """
        Checkout if no working copy exists, otherwise update.

        Returns:
            SVN command output
        """
        if self.is_checked_out():
            return self.update()
        else:
            return self.checkout()

    def get_current_revision(self):
        """Get current SVN revision number."""
        if self.is_checked_out():
            cmd = f"svn info {self.checkout_dir} --show-item revision"
        else:
            cmd = f"svn info {self.families_url} --show-item revision"

        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode != 0:
            logging.error(f"Failed to get revision: {result.stderr}")
            raise RuntimeError(f"Failed to get SVN revision: {result.stderr}")

        revision = result.stdout.strip()
        logging.info(f"Current SVN revision: {revision}")
        return revision

    def list_all_families(self):
        """
        List all family IDs.

        Returns:
            Sorted list of family IDs (e.g., ['PF00001', 'PF00002', ...])
        """
        if self.is_checked_out():
            # List from local checkout (fast)
            families = [d.name for d in self.checkout_dir.iterdir()
                       if d.is_dir() and d.name.startswith('PF') and not d.name.startswith('.')]
        else:
            # List from SVN server (slower, requires network)
            cmd = f"svn list {self.families_url}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

            if result.returncode != 0:
                raise RuntimeError(f"Failed to list families: {result.stderr}")

            families = [line.strip('/') for line in result.stdout.strip().split('\n')
                       if line.strip() and line.startswith('PF')]

        logging.info(f"Found {len(families)} families")
        return sorted(families)

    def get_changed_families(self, last_revision, current_revision):
        """
        Get list of families that changed between revisions.

        Args:
            last_revision: Previous revision number
            current_revision: Current revision number

        Returns:
            List of changed family IDs
        """
        cmd = f"svn log -r {last_revision}:{current_revision} -v {self.families_url}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode != 0:
            logging.error(f"Failed to get SVN log: {result.stderr}")
            raise RuntimeError(f"Failed to get SVN log: {result.stderr}")

        changed = set()
        for line in result.stdout.split('\n'):
            if '/trunk/Data/Families/PF' in line:
                parts = line.split('/trunk/Data/Families/')
                if len(parts) > 1:
                    family = parts[1].split('/')[0]
                    if family.startswith('PF'):
                        changed.add(family)

        logging.info(f"Found {len(changed)} changed families between r{last_revision} and r{current_revision}")
        return sorted(changed)

    def sync_seeds(self, families=None):
        """
        Sync SEED files from SVN checkout to flat SEED directory.

        This copies SEED files from the nested SVN structure (Families/PF00001/SEED)
        to the flat structure expected by the pipeline (SEED/PF00001_SEED).

        Args:
            families: List of families to sync (if None, sync all)

        Returns:
            Dictionary mapping family_id to seed_path
        """
        if not self.is_checked_out():
            raise RuntimeError("No checkout exists. Run checkout_or_update() first.")

        if families is None:
            families = self.list_all_families()

        logging.info(f"Syncing {len(families)} SEED files to flat directory...")

        results = {}
        synced = 0
        skipped = 0
        missing = 0

        for i, family_id in enumerate(families, 1):
            src = self.checkout_dir / family_id / "SEED"
            dst = self.seed_dir / f"{family_id}_SEED"

            if not src.exists():
                missing += 1
                results[family_id] = None
                continue

            # Copy if destination doesn't exist or source is newer
            if not dst.exists():
                shutil.copy2(src, dst)
                synced += 1
            elif src.stat().st_mtime > dst.stat().st_mtime:
                shutil.copy2(src, dst)
                synced += 1
            else:
                skipped += 1

            results[family_id] = dst

            if i % 2000 == 0:
                logging.info(f"Progress: {i}/{len(families)} processed")

        logging.info(f"Sync complete: {synced} copied, {skipped} up-to-date, {missing} missing SEED files")
        return results

    def extract_all_families(self, families=None):
        """
        Main entry point: checkout/update SVN and sync SEED files.

        This method:
        1. Performs SVN checkout (first run) or update (subsequent runs)
        2. Syncs SEED files from nested SVN structure to flat directory

        Args:
            families: List of families to process (if None, process all)

        Returns:
            Dictionary mapping family_id to (seed_path, None)
            Note: Second element is None for backwards compatibility
                  (previously was hmm_path, now removed)
        """
        # Checkout or update SVN
        self.checkout_or_update()

        # Sync SEED files to flat directory
        seed_results = self.sync_seeds(families)

        # Return in format compatible with existing pipeline
        # (family_id -> (seed_path, hmm_path)) where hmm_path is now always None
        return {fam: (path, None) for fam, path in seed_results.items()}

    def get_seed_path(self, family_id):
        """
        Get path to SEED file for a family.

        Args:
            family_id: Pfam family ID (e.g., PF00001)

        Returns:
            Path to SEED file in flat directory, or None if not found
        """
        seed_path = self.seed_dir / f"{family_id}_SEED"
        if seed_path.exists():
            return seed_path
        return None


if __name__ == "__main__":
    # Test the SVN manager
    import sys

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    repo_url = "https://xfam-svn-hl.ebi.ac.uk/svn/pfam/"
    data_dir = "./test_svn_checkout"

    manager = SVNManager(repo_url, data_dir)

    # Test getting revision (doesn't require checkout)
    try:
        revision = manager.get_current_revision()
        print(f"Current SVN revision: {revision}")
    except Exception as e:
        print(f"Failed to get revision: {e}")
        sys.exit(1)

    # Test listing families (from server if no checkout)
    print("\nListing first 10 families from SVN server...")
    try:
        families = manager.list_all_families()[:10]
        print(f"First 10 families: {families}")
    except Exception as e:
        print(f"Failed to list families: {e}")

    print("\nTo test checkout/update, run:")
    print(f"  python3 {__file__} --checkout")
    print("\nNote: Full checkout takes several hours for ~20k families")
