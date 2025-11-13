#!/usr/bin/env python3
"""
Simple SVN manager for extracting Pfam SEED alignments and HMMs from SVN repository.
"""

import subprocess
import os
import logging
from pathlib import Path

class SVNManager:
    """Manages SVN operations for Pfam family extraction."""

    def __init__(self, repo_url, output_dir):
        """
        Initialize SVN manager.

        Args:
            repo_url: Base SVN repository URL
            output_dir: Base directory for extracted files
        """
        self.repo_url = repo_url.rstrip('/')
        self.output_dir = Path(output_dir)
        self.families_url = f"{self.repo_url}/trunk/Data/Families"

        # Create output directories
        self.seed_dir = self.output_dir / "SEED"
        self.hmm_dir = self.output_dir / "HMM"
        self.seed_dir.mkdir(parents=True, exist_ok=True)
        self.hmm_dir.mkdir(parents=True, exist_ok=True)

        logging.info(f"SVN Manager initialized: {self.families_url}")
        logging.info(f"Output directory: {self.output_dir}")

    def get_current_revision(self):
        """Get current SVN revision number."""
        try:
            cmd = f"svn info {self.families_url} --show-item revision"
            result = subprocess.run(cmd, shell=True, check=True,
                                  capture_output=True, text=True)
            revision = result.stdout.strip()
            logging.info(f"Current SVN revision: {revision}")
            return revision
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to get SVN revision: {e}")
            raise

    def list_all_families(self):
        """List all family IDs in the repository."""
        try:
            cmd = f"svn list {self.families_url}"
            result = subprocess.run(cmd, shell=True, check=True,
                                  capture_output=True, text=True)
            families = [line.strip('/') for line in result.stdout.strip().split('\n')
                       if line.strip() and line.startswith('PF')]
            logging.info(f"Found {len(families)} families in SVN")
            return sorted(families)
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to list families: {e}")
            raise

    def get_changed_families(self, last_revision, current_revision):
        """
        Get list of families that changed between revisions.

        Args:
            last_revision: Previous revision number
            current_revision: Current revision number

        Returns:
            List of changed family IDs
        """
        try:
            cmd = f"svn log -r {last_revision}:{current_revision} --verbose {self.families_url}"
            result = subprocess.run(cmd, shell=True, check=True,
                                  capture_output=True, text=True)

            # Parse log output to find changed families
            changed_families = set()
            for line in result.stdout.split('\n'):
                if 'trunk/Data/Families/PF' in line:
                    # Extract family ID from path like: /trunk/Data/Families/PF00001/SEED
                    parts = line.split('trunk/Data/Families/')
                    if len(parts) > 1:
                        family = parts[1].split('/')[0]
                        if family.startswith('PF'):
                            changed_families.add(family)

            logging.info(f"Found {len(changed_families)} changed families between r{last_revision} and r{current_revision}")
            return sorted(changed_families)
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to get changed families: {e}")
            raise

    def extract_seed(self, family_id):
        """
        Extract SEED alignment for a family.

        Args:
            family_id: Pfam family ID (e.g., PF00001)

        Returns:
            Path to extracted SEED file
        """
        seed_url = f"{self.families_url}/{family_id}/SEED"
        seed_path = self.seed_dir / f"{family_id}_SEED"

        try:
            cmd = f"svn export --force {seed_url} {seed_path}"
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
            logging.debug(f"Extracted SEED for {family_id}")
            return seed_path
        except subprocess.CalledProcessError as e:
            logging.warning(f"Failed to extract SEED for {family_id}: {e}")
            return None

    def extract_hmm(self, family_id):
        """
        Extract HMM profile for a family.

        Args:
            family_id: Pfam family ID (e.g., PF00001)

        Returns:
            Path to extracted HMM file
        """
        hmm_url = f"{self.families_url}/{family_id}/HMM"
        hmm_path = self.hmm_dir / f"{family_id}.hmm"

        try:
            cmd = f"svn export --force {hmm_url} {hmm_path}"
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
            logging.debug(f"Extracted HMM for {family_id}")
            return hmm_path
        except subprocess.CalledProcessError as e:
            logging.warning(f"Failed to extract HMM for {family_id}: {e}")
            return None

    def extract_family(self, family_id):
        """
        Extract both SEED and HMM for a family.

        Args:
            family_id: Pfam family ID (e.g., PF00001)

        Returns:
            Tuple of (seed_path, hmm_path)
        """
        seed_path = self.extract_seed(family_id)
        hmm_path = self.extract_hmm(family_id)
        return (seed_path, hmm_path)

    def bulk_export_all_families(self, skip_existing=True):
        """
        Bulk export entire Families directory (much faster than individual exports).
        Resumable: skips families that have already been extracted.

        Args:
            skip_existing: Skip families that already have SEED/HMM files

        Returns:
            Dictionary mapping family_id to (seed_path, hmm_path)
        """
        import shutil
        from tempfile import mkdtemp

        # Check what's already been extracted
        if skip_existing:
            existing_seeds = {f.stem.replace('_SEED', '') for f in self.seed_dir.glob("PF*_SEED")}
            existing_hmms = {f.stem for f in self.hmm_dir.glob("PF*.hmm")}
            already_extracted = existing_seeds & existing_hmms  # Both must exist
            logging.info(f"Found {len(already_extracted)} families already extracted, will skip them")
        else:
            already_extracted = set()

        temp_dir = Path(mkdtemp(prefix="pfam_bulk_"))
        logging.info(f"Performing bulk SVN export to {temp_dir}...")

        try:
            # Export entire Families directory
            cmd = f"svn export {self.families_url} {temp_dir}/Families"
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
            logging.info("Bulk export complete, copying files...")

            results = {}
            families_dir = temp_dir / "Families"
            copied = 0
            skipped = 0

            # Iterate through exported families
            for family_dir in sorted(families_dir.glob("PF*")):
                if not family_dir.is_dir():
                    continue

                family_id = family_dir.name

                # Skip if already extracted
                if skip_existing and family_id in already_extracted:
                    seed_path = self.seed_dir / f"{family_id}_SEED"
                    hmm_path = self.hmm_dir / f"{family_id}.hmm"
                    results[family_id] = (seed_path, hmm_path)
                    skipped += 1
                    continue

                seed_src = family_dir / "SEED"
                hmm_src = family_dir / "HMM"

                seed_path = None
                hmm_path = None

                # Copy SEED if exists
                if seed_src.exists():
                    seed_path = self.seed_dir / f"{family_id}_SEED"
                    shutil.copy2(seed_src, seed_path)

                # Copy HMM if exists
                if hmm_src.exists():
                    hmm_path = self.hmm_dir / f"{family_id}.hmm"
                    shutil.copy2(hmm_src, hmm_path)

                results[family_id] = (seed_path, hmm_path)
                copied += 1

                if (copied + skipped) % 1000 == 0:
                    logging.info(f"Progress: {copied} copied, {skipped} skipped, {copied + skipped} total")

            logging.info(f"Bulk extraction complete: {copied} families copied, {skipped} skipped")
            return results

        finally:
            # Clean up temp directory
            if temp_dir.exists():
                shutil.rmtree(temp_dir)
                logging.info(f"Cleaned up temp directory: {temp_dir}")

    def extract_all_families(self, families=None, use_bulk=True):
        """
        Extract SEED and HMM files for all families.

        Args:
            families: List of family IDs to extract (if None, extracts all)
            use_bulk: Use bulk export method (much faster for large sets)

        Returns:
            Dictionary mapping family_id to (seed_path, hmm_path)
        """
        if families is None:
            families = self.list_all_families()

        # Use bulk export if processing many families
        if use_bulk and len(families) > 100:
            logging.info(f"Using bulk export for {len(families)} families...")
            all_results = self.bulk_export_all_families()

            # Filter to requested families if specified
            if families != self.list_all_families():
                results = {fam: all_results.get(fam, (None, None)) for fam in families}
            else:
                results = all_results

            return results

        # Fall back to individual extraction for small sets
        results = {}
        total = len(families)

        logging.info(f"Extracting {total} families individually...")

        for i, family_id in enumerate(families, 1):
            if i % 100 == 0:
                logging.info(f"Progress: {i}/{total} families extracted")

            seed_path, hmm_path = self.extract_family(family_id)
            results[family_id] = (seed_path, hmm_path)

        logging.info(f"Extraction complete: {total} families processed")
        return results


if __name__ == "__main__":
    # Test the SVN manager
    logging.basicConfig(level=logging.INFO,
                       format='%(asctime)s - %(levelname)s - %(message)s')

    repo_url = "https://xfam-svn-hl.ebi.ac.uk/svn/pfam/"
    output_dir = "./test_output"

    manager = SVNManager(repo_url, output_dir)

    # Test getting revision
    revision = manager.get_current_revision()
    print(f"Current revision: {revision}")

    # Test listing families (first 10)
    families = manager.list_all_families()[:10]
    print(f"First 10 families: {families}")

    # Test extracting one family
    if families:
        test_family = families[0]
        print(f"\nTesting extraction of {test_family}...")
        seed_path, hmm_path = manager.extract_family(test_family)
        print(f"SEED: {seed_path}")
        print(f"HMM: {hmm_path}")
