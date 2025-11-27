#!/usr/bin/env python3
"""
Local TED domain data access module.

This module provides access to TED domain information from a local TSV file,
replacing the need for TED API calls. All scripts that need TED domain data
should use this module instead of calling the API directly.

Usage:
    from ted_local import TEDLocal

    ted = TEDLocal()  # Uses default data path
    domains = ted.get_domains('P12345')

    # Or with custom path:
    ted = TEDLocal('/path/to/ted_data.tsv.gz')

Data file format (tab-separated):
    ted_id                          chopping            consensus_level
    AF-A0A000-F1-model_v4_TED01     11-41_290-389       high
"""

import gzip
import os
import re
from typing import Dict, List, Optional
from collections import defaultdict


# Default path to the TED data file
DEFAULT_TED_DATA_PATH = '/nfs/production/agb/pfam/data/TED/ted_365m_domain_boundaries_consensus_level.tsv.gz'


class TEDLocal:
    """
    Class for accessing local TED domain data.

    Loads TED domain boundaries from a local TSV file and provides
    methods to query domains by UniProt accession.
    """

    def __init__(self, data_path: str = None):
        """
        Initialize TEDLocal.

        Args:
            data_path: Path to the TED data file (gzipped TSV).
                      If None, uses the default path.
        """
        self.data_path = data_path or DEFAULT_TED_DATA_PATH
        self._data: Dict[str, List[Dict]] = None
        self._loaded = False

    def _load_data(self):
        """Load and index the TED data file by UniProt accession."""
        if self._loaded:
            return

        self._data = defaultdict(list)

        if not os.path.exists(self.data_path):
            raise FileNotFoundError(f"TED data file not found: {self.data_path}")

        # Determine if file is gzipped
        if self.data_path.endswith('.gz'):
            open_func = gzip.open
            mode = 'rt'
        else:
            open_func = open
            mode = 'r'

        with open_func(self.data_path, mode) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) < 2:
                    continue

                ted_id = parts[0]
                chopping = parts[1]
                consensus_level = parts[2] if len(parts) > 2 else 'unknown'

                # Extract UniProt accession from ted_id
                # Format: AF-A0A000-F1-model_v4_TED01 -> A0A000
                uniprot_acc = self._extract_uniprot_acc(ted_id)

                if uniprot_acc:
                    # Parse chopping into segments
                    segments = self._parse_chopping(chopping)

                    # Extract TED domain number (e.g., TED01)
                    ted_suffix = ted_id.split('_')[-1] if '_' in ted_id else ted_id

                    self._data[uniprot_acc].append({
                        'ted_id': ted_id,
                        'ted_suffix': ted_suffix,
                        'chopping': chopping,
                        'segments': segments,
                        'consensus_level': consensus_level
                    })

        self._loaded = True

    @staticmethod
    def _extract_uniprot_acc(ted_id: str) -> Optional[str]:
        """
        Extract UniProt accession from TED ID.

        Format: AF-A0A000-F1-model_v4_TED01 -> A0A000

        Args:
            ted_id: The TED domain identifier

        Returns:
            UniProt accession or None if not found
        """
        # Pattern: AF-{UNIPROT_ACC}-F1-model_v4_TED##
        match = re.match(r'^AF-([A-Z0-9]+)-F\d+-', ted_id)
        if match:
            return match.group(1)
        return None

    @staticmethod
    def _parse_chopping(chopping: str) -> List[Dict[str, int]]:
        """
        Parse chopping string into list of segment dictionaries.

        Args:
            chopping: Chopping string like "11-41_290-389"

        Returns:
            List of dicts with 'start' and 'end' keys
        """
        segments = []
        for segment in chopping.split('_'):
            if '-' in segment:
                try:
                    start, end = segment.split('-')
                    segments.append({
                        'start': int(start),
                        'end': int(end)
                    })
                except ValueError:
                    continue
        return segments

    def get_domains(self, uniprot_acc: str) -> List[Dict]:
        """
        Get TED domains for a UniProt accession.

        Args:
            uniprot_acc: UniProt accession (with or without version, e.g., 'P12345' or 'P12345.2')

        Returns:
            List of domain dictionaries, each containing:
                - ted_id: Full TED identifier
                - ted_suffix: Short TED domain name (e.g., 'TED01')
                - chopping: Original chopping string
                - segments: List of {'start': int, 'end': int}
                - consensus_level: 'high', 'medium', etc.
        """
        self._load_data()

        # Strip version number if present (P12345.2 -> P12345)
        base_acc = uniprot_acc.split('.')[0]

        return self._data.get(base_acc, [])

    def get_domain_count(self, uniprot_acc: str) -> int:
        """
        Get count of TED domains for a UniProt accession.

        Args:
            uniprot_acc: UniProt accession

        Returns:
            Number of TED domains
        """
        return len(self.get_domains(uniprot_acc))

    def get_domain_boundaries(self, uniprot_acc: str) -> List[tuple]:
        """
        Get simplified domain boundaries (start, end) for a UniProt accession.

        For multi-segment domains, returns the overall start and end.

        Args:
            uniprot_acc: UniProt accession

        Returns:
            List of (start, end) tuples, one per domain
        """
        domains = self.get_domains(uniprot_acc)
        boundaries = []

        for domain in domains:
            segments = domain['segments']
            if segments:
                start = min(seg['start'] for seg in segments)
                end = max(seg['end'] for seg in segments)
                boundaries.append((start, end))

        return boundaries

    def has_domains(self, uniprot_acc: str) -> bool:
        """
        Check if a UniProt accession has any TED domains.

        Args:
            uniprot_acc: UniProt accession

        Returns:
            True if domains exist, False otherwise
        """
        return self.get_domain_count(uniprot_acc) > 0

    def get_all_accessions(self) -> List[str]:
        """
        Get list of all UniProt accessions in the dataset.

        Returns:
            List of UniProt accessions
        """
        self._load_data()
        return list(self._data.keys())


# Module-level convenience functions using a shared instance
_shared_instance: TEDLocal = None


def get_shared_instance(data_path: str = None) -> TEDLocal:
    """
    Get or create a shared TEDLocal instance.

    This is useful when multiple parts of a script need TED data
    without loading the file multiple times.

    Args:
        data_path: Path to TED data file (only used on first call)

    Returns:
        Shared TEDLocal instance
    """
    global _shared_instance
    if _shared_instance is None:
        _shared_instance = TEDLocal(data_path)
    return _shared_instance


def get_domains(uniprot_acc: str, data_path: str = None) -> List[Dict]:
    """
    Convenience function to get TED domains for a UniProt accession.

    Uses a shared TEDLocal instance.

    Args:
        uniprot_acc: UniProt accession
        data_path: Optional path to TED data file

    Returns:
        List of domain dictionaries
    """
    return get_shared_instance(data_path).get_domains(uniprot_acc)


def get_domain_count(uniprot_acc: str, data_path: str = None) -> int:
    """
    Convenience function to get TED domain count for a UniProt accession.

    Args:
        uniprot_acc: UniProt accession
        data_path: Optional path to TED data file

    Returns:
        Number of TED domains
    """
    return get_shared_instance(data_path).get_domain_count(uniprot_acc)


if __name__ == '__main__':
    # Simple test/demo when run directly
    import sys

    if len(sys.argv) < 2:
        print("Usage: ted_local.py <uniprot_accession>")
        print("\nExample:")
        print("  ted_local.py P12345")
        sys.exit(1)

    acc = sys.argv[1]
    ted = TEDLocal()

    print(f"Looking up TED domains for: {acc}")
    domains = ted.get_domains(acc)

    if not domains:
        print(f"No TED domains found for {acc}")
    else:
        print(f"Found {len(domains)} TED domain(s):")
        for domain in domains:
            print(f"  {domain['ted_suffix']}: {domain['chopping']} ({domain['consensus_level']})")
            for seg in domain['segments']:
                print(f"    Segment: {seg['start']}-{seg['end']}")
