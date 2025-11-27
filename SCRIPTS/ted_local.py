#!/usr/bin/env python3
"""
Local TED domain data access module with SQLite indexing for fast lookups.

This module provides access to TED domain information from a local SQLite
database (indexed from the original TSV file), replacing the need for TED
API calls. All scripts that need TED domain data should use this module.

First-time setup (creates the index - takes a few minutes):
    ted_local.py --create-index

Usage in scripts:
    from ted_local import TEDLocal

    ted = TEDLocal()  # Uses default database path
    domains = ted.get_domains('P12345')

    # Or with custom path:
    ted = TEDLocal('/path/to/ted_index.db')

Command-line usage:
    ted_local.py --create-index              # Create SQLite index from TSV
    ted_local.py P12345                      # Look up domains for accession
    ted_local.py --stats                     # Show database statistics
"""

import argparse
import gzip
import os
import re
import sqlite3
import sys
from typing import Dict, List, Optional


# Default paths
DEFAULT_TED_TSV_PATH = '/nfs/production/agb/pfam/data/TED/ted_365m_domain_boundaries_consensus_level.tsv.gz'
DEFAULT_TED_DB_PATH = '/nfs/production/agb/pfam/data/TED/ted_index.db'


class TEDLocal:
    """
    Class for accessing local TED domain data via SQLite index.

    Uses a SQLite database for fast indexed lookups by UniProt accession.
    """

    def __init__(self, db_path: str = None):
        """
        Initialize TEDLocal.

        Args:
            db_path: Path to the SQLite database file.
                    If None, uses the default path.
        """
        self.db_path = db_path or DEFAULT_TED_DB_PATH
        self._conn: sqlite3.Connection = None

    def _get_connection(self) -> sqlite3.Connection:
        """Get or create database connection."""
        if self._conn is None:
            if not os.path.exists(self.db_path):
                raise FileNotFoundError(
                    f"TED database not found: {self.db_path}\n"
                    f"Run 'ted_local.py --create-index' to create it."
                )
            self._conn = sqlite3.connect(self.db_path)
            self._conn.row_factory = sqlite3.Row
        return self._conn

    def close(self):
        """Close database connection."""
        if self._conn:
            self._conn.close()
            self._conn = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

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
            uniprot_acc: UniProt accession (with or without version)

        Returns:
            List of domain dictionaries, each containing:
                - ted_id: Full TED identifier
                - ted_suffix: Short TED domain name (e.g., 'TED01')
                - chopping: Original chopping string
                - segments: List of {'start': int, 'end': int}
                - consensus_level: 'high', 'medium', etc.
        """
        # Strip version number if present (P12345.2 -> P12345)
        base_acc = uniprot_acc.split('.')[0]

        conn = self._get_connection()
        cursor = conn.execute(
            "SELECT ted_id, ted_suffix, chopping, consensus_level "
            "FROM domains WHERE uniprot_acc = ?",
            (base_acc,)
        )

        domains = []
        for row in cursor:
            domains.append({
                'ted_id': row['ted_id'],
                'ted_suffix': row['ted_suffix'],
                'chopping': row['chopping'],
                'segments': self._parse_chopping(row['chopping']),
                'consensus_level': row['consensus_level']
            })

        return domains

    def get_domain_count(self, uniprot_acc: str) -> int:
        """Get count of TED domains for a UniProt accession."""
        return len(self.get_domains(uniprot_acc))

    def get_domain_boundaries(self, uniprot_acc: str) -> List[tuple]:
        """
        Get simplified domain boundaries (start, end) for a UniProt accession.

        For multi-segment domains, returns the overall start and end.
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
        """Check if a UniProt accession has any TED domains."""
        base_acc = uniprot_acc.split('.')[0]
        conn = self._get_connection()
        cursor = conn.execute(
            "SELECT 1 FROM domains WHERE uniprot_acc = ? LIMIT 1",
            (base_acc,)
        )
        return cursor.fetchone() is not None

    def get_stats(self) -> Dict:
        """Get database statistics."""
        conn = self._get_connection()

        stats = {}

        # Total domains
        cursor = conn.execute("SELECT COUNT(*) FROM domains")
        stats['total_domains'] = cursor.fetchone()[0]

        # Unique proteins
        cursor = conn.execute("SELECT COUNT(DISTINCT uniprot_acc) FROM domains")
        stats['unique_proteins'] = cursor.fetchone()[0]

        # Consensus level breakdown
        cursor = conn.execute(
            "SELECT consensus_level, COUNT(*) FROM domains GROUP BY consensus_level"
        )
        stats['by_consensus'] = dict(cursor.fetchall())

        return stats


def create_index(tsv_path: str = None, db_path: str = None, verbose: bool = True):
    """
    Create SQLite index from the TED TSV file.

    Memory-optimized: uses small batches and streaming to avoid OOM errors.

    Args:
        tsv_path: Path to the gzipped TSV file
        db_path: Path for the output SQLite database
        verbose: Print progress information
    """
    import gc

    tsv_path = tsv_path or DEFAULT_TED_TSV_PATH
    db_path = db_path or DEFAULT_TED_DB_PATH

    if not os.path.exists(tsv_path):
        print(f"Error: TSV file not found: {tsv_path}", file=sys.stderr)
        sys.exit(1)

    if verbose:
        print(f"Creating TED index database...")
        print(f"  Source: {tsv_path}")
        print(f"  Output: {db_path}")

    # Remove existing database if present
    if os.path.exists(db_path):
        os.remove(db_path)

    # Create database and table with optimized settings for bulk insert
    conn = sqlite3.connect(db_path)

    # Optimize SQLite for bulk loading (much faster, uses less memory)
    conn.execute("PRAGMA synchronous = OFF")
    conn.execute("PRAGMA journal_mode = OFF")
    conn.execute("PRAGMA cache_size = 10000")
    conn.execute("PRAGMA temp_store = MEMORY")

    conn.execute("""
        CREATE TABLE domains (
            uniprot_acc TEXT NOT NULL,
            ted_id TEXT NOT NULL,
            ted_suffix TEXT NOT NULL,
            chopping TEXT NOT NULL,
            consensus_level TEXT
        )
    """)

    # Read and insert data using a generator to minimize memory
    if verbose:
        print("  Reading TSV file and inserting records...")

    record_count = 0
    batch_size = 10000  # Smaller batches to reduce memory usage

    def record_generator():
        """Generator that yields records one at a time to minimize memory."""
        # Pattern to extract UniProt accession from ted_id
        acc_pattern = re.compile(r'^AF-([A-Z0-9]+)-F\d+-')

        if tsv_path.endswith('.gz'):
            open_func = gzip.open
            mode = 'rt'
        else:
            open_func = open
            mode = 'r'

        with open_func(tsv_path, mode) as f:
            for line in f:
                line = line.rstrip('\n\r')
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) < 2:
                    continue

                ted_id = parts[0]
                chopping = parts[1]
                consensus_level = parts[2] if len(parts) > 2 else 'unknown'

                # Extract UniProt accession
                match = acc_pattern.match(ted_id)
                if not match:
                    continue

                uniprot_acc = match.group(1)

                # Extract TED suffix (e.g., TED01)
                ted_suffix = ted_id.split('_')[-1] if '_' in ted_id else ted_id

                yield (uniprot_acc, ted_id, ted_suffix, chopping, consensus_level)

    # Process in small batches
    batch = []
    for record in record_generator():
        batch.append(record)
        record_count += 1

        if len(batch) >= batch_size:
            conn.executemany(
                "INSERT INTO domains (uniprot_acc, ted_id, ted_suffix, chopping, consensus_level) "
                "VALUES (?, ?, ?, ?, ?)",
                batch
            )
            batch.clear()  # Clear list in-place (more memory efficient than batch = [])

            if verbose and record_count % 1000000 == 0:
                print(f"    Processed {record_count:,} records...")
                gc.collect()  # Force garbage collection periodically

    # Insert remaining records
    if batch:
        conn.executemany(
            "INSERT INTO domains (uniprot_acc, ted_id, ted_suffix, chopping, consensus_level) "
            "VALUES (?, ?, ?, ?, ?)",
            batch
        )
        batch.clear()

    conn.commit()
    gc.collect()

    if verbose:
        print(f"  Processed {record_count:,} records total")
        print(f"  Creating index on uniprot_acc (this may take a minute)...")

    # Create index for fast lookups
    conn.execute("CREATE INDEX idx_uniprot_acc ON domains (uniprot_acc)")
    conn.commit()

    conn.close()

    # Get file size
    db_size = os.path.getsize(db_path) / (1024 * 1024)

    if verbose:
        print(f"\nIndex created successfully!")
        print(f"  Total records: {record_count:,}")
        print(f"  Database size: {db_size:.1f} MB")


# Module-level shared instance
_shared_instance: TEDLocal = None


def get_shared_instance(db_path: str = None) -> TEDLocal:
    """
    Get or create a shared TEDLocal instance.

    Args:
        db_path: Path to database file (only used on first call)

    Returns:
        Shared TEDLocal instance
    """
    global _shared_instance
    if _shared_instance is None:
        _shared_instance = TEDLocal(db_path)
    return _shared_instance


def get_domains(uniprot_acc: str, db_path: str = None) -> List[Dict]:
    """Convenience function to get TED domains for a UniProt accession."""
    return get_shared_instance(db_path).get_domains(uniprot_acc)


def get_domain_count(uniprot_acc: str, db_path: str = None) -> int:
    """Convenience function to get TED domain count for a UniProt accession."""
    return get_shared_instance(db_path).get_domain_count(uniprot_acc)


def main():
    parser = argparse.ArgumentParser(
        description='Local TED domain data access with SQLite indexing'
    )
    parser.add_argument(
        'accession',
        nargs='?',
        help='UniProt accession to look up'
    )
    parser.add_argument(
        '--create-index',
        action='store_true',
        help='Create SQLite index from TSV file'
    )
    parser.add_argument(
        '--tsv-path',
        default=DEFAULT_TED_TSV_PATH,
        help=f'Path to TSV file (default: {DEFAULT_TED_TSV_PATH})'
    )
    parser.add_argument(
        '--db-path',
        default=DEFAULT_TED_DB_PATH,
        help=f'Path to SQLite database (default: {DEFAULT_TED_DB_PATH})'
    )
    parser.add_argument(
        '--stats',
        action='store_true',
        help='Show database statistics'
    )

    args = parser.parse_args()

    if args.create_index:
        create_index(args.tsv_path, args.db_path)
        return

    if args.stats:
        ted = TEDLocal(args.db_path)
        stats = ted.get_stats()
        print("TED Database Statistics:")
        print(f"  Total domains: {stats['total_domains']:,}")
        print(f"  Unique proteins: {stats['unique_proteins']:,}")
        print(f"  By consensus level:")
        for level, count in sorted(stats['by_consensus'].items()):
            print(f"    {level}: {count:,}")
        return

    if not args.accession:
        parser.print_help()
        print("\nExamples:")
        print("  ted_local.py --create-index    # Create index (first-time setup)")
        print("  ted_local.py P12345            # Look up domains")
        print("  ted_local.py --stats           # Show statistics")
        sys.exit(1)

    # Look up accession
    ted = TEDLocal(args.db_path)
    domains = ted.get_domains(args.accession)

    if not domains:
        print(f"No TED domains found for {args.accession}")
    else:
        print(f"Found {len(domains)} TED domain(s) for {args.accession}:")
        for domain in domains:
            print(f"  {domain['ted_suffix']}: {domain['chopping']} ({domain['consensus_level']})")


if __name__ == '__main__':
    main()
