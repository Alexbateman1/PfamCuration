#!/usr/bin/env python3
"""
Analyze STRING protein networks for Pfam families to identify over-represented domains.

This script:
1. Extracts STRING identifiers from a Pfam family's seq_info file
2. Streams through global STRING network data files
3. Retrieves network partners for each protein (direct or multi-hop)
4. Maps proteins to UniProt accessions
5. Looks up Pfam domains for each UniProt accession
6. For proteins without domains, identifies UniRef_50 clusters
7. Calculates enrichment statistics for domains and clusters
8. Reports over-represented domains/clusters above frequency threshold

Usage:
    analyze_string_network.py <pfam_directory> [options]

Example:
    analyze_string_network.py /path/to/PF06267 --score-threshold 400 --min-frequency 0.2
"""

import argparse
import gzip
import os
import re
import subprocess
import sys
from collections import defaultdict
from typing import Dict, List, Set, Tuple, Optional
import urllib.request
import urllib.error


# Default paths
DEFAULT_STRING_DATA_DIR = "/nfs/production/agb/pfam/data/STRING"
DEFAULT_MYSQL_CONFIG = "~/.my.cnf"


def parse_seq_info_for_string_ids(seq_info_path: str) -> List[Tuple[str, str]]:
    """
    Parse seq_info file to extract STRING identifiers.

    Args:
        seq_info_path: Path to seq_info file

    Returns:
        List of tuples: (taxid, protein_id) e.g., ('588602', 'SAMN04487991_1451')
    """
    string_ids = []

    if not os.path.exists(seq_info_path):
        print(f"Warning: seq_info file not found at {seq_info_path}")
        return string_ids

    with open(seq_info_path, 'r') as f:
        for line in f:
            line = line.strip()
            # Match lines like: DR   STRING; 588602.SAMN04487991_1451; -.
            match = re.match(r'DR\s+STRING;\s+(\d+)\.(\S+);\s+', line)
            if match:
                taxid = match.group(1)
                protein_id = match.group(2)
                string_ids.append((taxid, protein_id))

    return string_ids


def run_species_summary(pfam_dir: str) -> bool:
    """
    Run species_summary_new.pl on the Pfam directory to generate seq_info.

    Args:
        pfam_dir: Path to Pfam family directory

    Returns:
        True if successful, False otherwise
    """
    try:
        print(f"Running species_summary_new.pl on {pfam_dir}...")
        result = subprocess.run(
            ['species_summary_new.pl', pfam_dir],
            capture_output=True,
            text=True,
            timeout=300
        )

        if result.returncode != 0:
            print(f"Warning: species_summary_new.pl returned non-zero exit code: {result.returncode}")
            if result.stderr:
                print(f"stderr: {result.stderr}")

        return True

    except subprocess.TimeoutExpired:
        print("Error: species_summary_new.pl timed out after 5 minutes")
        return False
    except FileNotFoundError:
        print("Error: species_summary_new.pl not found in PATH")
        print("Attempting to continue with existing seq_info file...")
        return True
    except Exception as e:
        print(f"Error running species_summary_new.pl: {e}")
        return False


def build_networks_from_links(links_file: str, query_proteins: Set[str],
                              score_threshold: int = 400, network_depth: int = 1) -> Dict[str, Set[str]]:
    """
    Stream through protein.links.full file to build networks for query proteins.

    Args:
        links_file: Path to protein.links.full.v12.0.txt.gz file
        query_proteins: Set of STRING protein IDs to find networks for
        score_threshold: Minimum combined score (0-999) to include link
        network_depth: Number of network hops

    Returns:
        Dictionary mapping protein -> set of connected proteins
    """
    print(f"  Building networks from links file (threshold={score_threshold}, depth={network_depth})...")

    # For depth > 1, we need to iteratively expand
    proteins_of_interest = query_proteins.copy()
    all_networks = defaultdict(set)

    for depth in range(network_depth):
        print(f"    Depth {depth + 1}: Looking for networks of {len(proteins_of_interest):,} proteins...")

        # Stream through file to find links for proteins of interest
        links_found = 0
        lines_processed = 0

        with gzip.open(links_file, 'rt') as f:
            # Read header
            header = f.readline().strip().split()

            # Find column indices
            try:
                protein1_idx = header.index('protein1')
                protein2_idx = header.index('protein2')
                combined_idx = header.index('combined_score')
            except ValueError as e:
                print(f"  Error: Expected column not found in header: {e}")
                return all_networks

            for line in f:
                lines_processed += 1
                parts = line.strip().split()

                if len(parts) > max(protein1_idx, protein2_idx, combined_idx):
                    protein1 = parts[protein1_idx]
                    protein2 = parts[protein2_idx]
                    score = int(parts[combined_idx])

                    if score >= score_threshold:
                        # Check if either protein is in our current set of interest
                        if protein1 in proteins_of_interest:
                            all_networks[protein1].add(protein2)
                            links_found += 1

                        if protein2 in proteins_of_interest:
                            all_networks[protein2].add(protein1)
                            links_found += 1

                if lines_processed % 10000000 == 0:
                    print(f"      Processed {lines_processed:,} lines, found {links_found:,} links...")

        print(f"    Depth {depth + 1}: Found {links_found:,} links from {lines_processed:,} total lines")

        # For next depth, expand to neighbors
        if depth + 1 < network_depth:
            new_proteins = set()
            for protein in proteins_of_interest:
                if protein in all_networks:
                    new_proteins.update(all_networks[protein])

            # Only add proteins we haven't seen before
            proteins_of_interest = new_proteins - proteins_of_interest - query_proteins

            if not proteins_of_interest:
                print(f"    No new proteins to expand at depth {depth + 2}, stopping early")
                break

    return all_networks


def get_uniprot_mappings(aliases_file: str, string_ids: Set[str]) -> Dict[str, Set[str]]:
    """
    Stream through protein.aliases file to get UniProt mappings for specific STRING IDs.

    Args:
        aliases_file: Path to protein.aliases.v12.0.txt.gz file
        string_ids: Set of STRING protein IDs to find mappings for

    Returns:
        Dictionary mapping STRING protein ID -> set of UniProt accessions
    """
    print(f"  Getting UniProt mappings for {len(string_ids):,} STRING proteins...")

    string_to_uniprot = defaultdict(set)
    lines_processed = 0
    mappings_found = 0

    with gzip.open(aliases_file, 'rt') as f:
        # Read header
        header = f.readline().strip().split('\t')

        try:
            string_id_idx = header.index('#string_protein_id')
            alias_idx = header.index('alias')
            source_idx = header.index('source')
        except ValueError as e:
            print(f"  Error: Expected column not found in header: {e}")
            return string_to_uniprot

        for line in f:
            lines_processed += 1
            parts = line.strip().split('\t')

            if len(parts) > max(string_id_idx, alias_idx, source_idx):
                string_id = parts[string_id_idx]

                # Only process if this is a STRING ID we care about
                if string_id in string_ids:
                    alias = parts[alias_idx]
                    source = parts[source_idx]

                    # Look specifically for UniProt_AC source
                    if source == 'UniProt_AC':
                        string_to_uniprot[string_id].add(alias)
                        mappings_found += 1

            if lines_processed % 10000000 == 0:
                print(f"    Processed {lines_processed:,} lines, found {mappings_found:,} UniProt mappings...")

    print(f"  Found {mappings_found:,} UniProt mappings from {lines_processed:,} total lines")
    print(f"  Mapped {len(string_to_uniprot):,} STRING IDs to UniProt")

    return string_to_uniprot


def query_pfam_domains(uniprot_acc: str, config_file: str = DEFAULT_MYSQL_CONFIG) -> List[Tuple[str, int, int, str]]:
    """
    Query pfam_live database for Pfam domains in a UniProt sequence.

    Args:
        uniprot_acc: UniProt accession
        config_file: Path to MySQL config file

    Returns:
        List of tuples: (pfamseq_acc, seq_start, seq_end, pfamA_acc)
    """
    import mysql.connector
    from configparser import ConfigParser

    config_file = os.path.expanduser(config_file)

    try:
        # Parse MySQL config file
        config = ConfigParser()
        config.read(config_file)

        # Get credentials
        db_config = {
            'host': config.get('client', 'host', fallback='localhost'),
            'user': config.get('client', 'user'),
            'password': config.get('client', 'password', fallback=''),
            'database': 'pfam_live'
        }

        # Connect and query
        conn = mysql.connector.connect(**db_config)
        cursor = conn.cursor()

        query = """
            SELECT pfamseq_acc, seq_start, seq_end, pfamA_acc
            FROM pfamA_reg_full_significant
            WHERE pfamseq_acc = %s
        """

        cursor.execute(query, (uniprot_acc,))
        results = cursor.fetchall()

        cursor.close()
        conn.close()

        return results

    except Exception as e:
        # Silently handle errors - not all proteins will be in database
        return []


def fetch_uniref50_from_uniprot(uniprot_acc: str) -> Optional[str]:
    """
    Fetch UniRef50 cluster from UniProt REST API.

    Args:
        uniprot_acc: UniProt accession

    Returns:
        UniRef50 cluster ID, or None if not found
    """
    try:
        # Use the search API to find UniRef50 cluster
        url = f"https://www.uniprot.org/uniprotkb/{uniprot_acc}.txt"

        with urllib.request.urlopen(url, timeout=10) as response:
            content = response.read().decode('utf-8')

            # Look for UniRef50 cross-reference
            # Format: DR   UniRef; UniRef50_P12345; -.
            match = re.search(r'DR\s+UniRef;\s+(UniRef50_\S+);', content)
            if match:
                return match.group(1)

    except Exception as e:
        pass

    return None


def analyze_networks(pfam_dir: str, string_data_dir: str, score_threshold: int = 400,
                     network_depth: int = 1, mysql_config: str = DEFAULT_MYSQL_CONFIG) -> Dict:
    """
    Main analysis function to process STRING networks for a Pfam family.

    Args:
        pfam_dir: Path to Pfam family directory
        string_data_dir: Directory for STRING data files
        score_threshold: Minimum STRING score to include interactions
        network_depth: Number of network hops (1=direct neighbors, 2=neighbors of neighbors)
        mysql_config: Path to MySQL config file

    Returns:
        Dictionary with analysis results
    """
    results = {
        'pfam_domains': defaultdict(lambda: {'networks': set(), 'proteins': set()}),
        'uniref_clusters': defaultdict(lambda: {'networks': set(), 'proteins': set()}),
        'total_networks': 0,
        'total_proteins': 0,
        'query_proteins': []
    }

    # Get seq_info path
    seq_info_path = os.path.join(pfam_dir, 'seq_info')

    # Parse STRING IDs from seq_info
    print("\n1. Extracting STRING identifiers from seq_info...")
    string_ids = parse_seq_info_for_string_ids(seq_info_path)

    if not string_ids:
        print("No STRING identifiers found in seq_info!")
        return results

    print(f"   Found {len(string_ids)} STRING identifiers")
    results['query_proteins'] = string_ids

    # Build full STRING IDs (taxid.protein_id)
    query_proteins = set()
    for taxid, protein_id in string_ids:
        full_string_id = f"{taxid}.{protein_id}"
        query_proteins.add(full_string_id)

    print(f"   Query proteins: {len(query_proteins)}")

    # Get paths to STRING files
    links_file = os.path.join(string_data_dir, 'protein.links.full.v12.0.txt.gz')
    aliases_file = os.path.join(string_data_dir, 'protein.aliases.v12.0.txt.gz')

    # Check files exist
    if not os.path.exists(links_file):
        print(f"ERROR: Links file not found: {links_file}")
        return results

    if not os.path.exists(aliases_file):
        print(f"ERROR: Aliases file not found: {aliases_file}")
        return results

    print(f"   Using links file: {links_file}")
    print(f"   Using aliases file: {aliases_file}")

    # Build networks by streaming through links file
    print("\n2. Building networks...")
    networks = build_networks_from_links(links_file, query_proteins, score_threshold, network_depth)

    if not networks:
        print("No networks found for query proteins!")
        return results

    # Collect all proteins we need UniProt mappings for
    all_network_proteins = set()
    for query_protein in query_proteins:
        if query_protein in networks:
            all_network_proteins.update(networks[query_protein])
            results['total_networks'] += 1

    print(f"\n   Found networks for {results['total_networks']} query proteins")
    print(f"   Total network partners: {len(all_network_proteins):,}")

    # Get UniProt mappings for all network proteins
    print("\n3. Mapping network partners to UniProt...")
    string_to_uniprot = get_uniprot_mappings(aliases_file, all_network_proteins)

    # Analyze each query protein's network
    print("\n4. Analyzing Pfam domains and UniRef clusters...")
    for query_protein in query_proteins:
        if query_protein not in networks:
            continue

        network_id = query_protein
        network_partners = networks[query_protein]

        print(f"\n   Analyzing {query_protein}: {len(network_partners)} partners")

        partners_with_uniprot = 0
        for partner_string_id in network_partners:
            results['total_proteins'] += 1

            # Get UniProt accessions for this partner
            uniprot_accs = string_to_uniprot.get(partner_string_id, set())

            if not uniprot_accs:
                continue

            partners_with_uniprot += 1

            # For each UniProt accession
            for uniprot_acc in uniprot_accs:
                # Query Pfam domains
                domains = query_pfam_domains(uniprot_acc, mysql_config)

                if domains:
                    # Has Pfam domains
                    for pfamseq_acc, seq_start, seq_end, pfamA_acc in domains:
                        results['pfam_domains'][pfamA_acc]['networks'].add(network_id)
                        results['pfam_domains'][pfamA_acc]['proteins'].add(partner_string_id)
                else:
                    # No Pfam domains - get UniRef50 cluster
                    uniref_cluster = fetch_uniref50_from_uniprot(uniprot_acc)

                    if uniref_cluster:
                        results['uniref_clusters'][uniref_cluster]['networks'].add(network_id)
                        results['uniref_clusters'][uniref_cluster]['proteins'].add(partner_string_id)

        print(f"     {partners_with_uniprot}/{len(network_partners)} partners have UniProt mappings")

    return results


def print_report(results: Dict, min_frequency: float = 0.2, output_file: Optional[str] = None):
    """
    Print analysis report.

    Args:
        results: Results dictionary from analyze_networks()
        min_frequency: Minimum frequency threshold (0-1) to report
        output_file: Optional file to write report to
    """
    output = []

    output.append("=" * 80)
    output.append("STRING Network Analysis Report")
    output.append("=" * 80)
    output.append("")

    output.append(f"Query proteins: {len(results['query_proteins'])}")
    output.append(f"Total networks analyzed: {results['total_networks']}")
    output.append(f"Total network proteins: {results['total_proteins']}")
    output.append(f"Reporting frequency threshold: {min_frequency:.1%}")
    output.append("")

    # Pfam domains
    output.append("=" * 80)
    output.append("Pfam Domain Enrichment")
    output.append("=" * 80)
    output.append("")

    # Filter by frequency threshold and sort by network count
    pfam_data = []
    for domain, data in results['pfam_domains'].items():
        network_freq = len(data['networks']) / results['total_networks'] if results['total_networks'] > 0 else 0
        if network_freq >= min_frequency:
            pfam_data.append((domain, data['networks'], data['proteins'], network_freq))

    pfam_data.sort(key=lambda x: len(x[1]), reverse=True)

    if pfam_data:
        output.append(f"{'Pfam Domain':<15} {'Networks':<12} {'Proteins':<12} {'Frequency':<12}")
        output.append("-" * 80)

        for domain, networks, proteins, freq in pfam_data:
            output.append(f"{domain:<15} {len(networks):<12} {len(proteins):<12} {freq:>6.1%}")
    else:
        output.append(f"No Pfam domains found above {min_frequency:.1%} frequency threshold.")

    output.append("")

    # UniRef clusters
    output.append("=" * 80)
    output.append("UniRef50 Cluster Enrichment (proteins without Pfam domains)")
    output.append("=" * 80)
    output.append("")

    # Filter by frequency threshold and sort by network count
    uniref_data = []
    for cluster, data in results['uniref_clusters'].items():
        network_freq = len(data['networks']) / results['total_networks'] if results['total_networks'] > 0 else 0
        if network_freq >= min_frequency:
            uniref_data.append((cluster, data['networks'], data['proteins'], network_freq))

    uniref_data.sort(key=lambda x: len(x[1]), reverse=True)

    if uniref_data:
        output.append(f"{'UniRef50 Cluster':<30} {'Networks':<12} {'Proteins':<12} {'Frequency':<12}")
        output.append("-" * 80)

        for cluster, networks, proteins, freq in uniref_data[:50]:  # Top 50
            output.append(f"{cluster:<30} {len(networks):<12} {len(proteins):<12} {freq:>6.1%}")

        if len(uniref_data) > 50:
            output.append(f"\n... and {len(uniref_data) - 50} more clusters above threshold")
    else:
        output.append(f"No UniRef50 clusters found above {min_frequency:.1%} frequency threshold.")

    output.append("")
    output.append("=" * 80)

    # Print to console
    report = "\n".join(output)
    print(report)

    # Write to file if requested
    if output_file:
        with open(output_file, 'w') as f:
            f.write(report)
        print(f"\nReport written to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze STRING protein networks for Pfam families",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument(
        'pfam_dir',
        help='Path to Pfam family directory'
    )

    parser.add_argument(
        '--string-data-dir',
        default=DEFAULT_STRING_DATA_DIR,
        help=f'Directory containing STRING data files (default: {DEFAULT_STRING_DATA_DIR})'
    )

    parser.add_argument(
        '--score-threshold',
        type=int,
        default=400,
        help='Minimum STRING combined score (0-999) to include interactions (default: 400 = medium confidence)'
    )

    parser.add_argument(
        '--network-depth',
        type=int,
        default=1,
        help='Network depth: 1=direct neighbors, 2=neighbors of neighbors, etc. (default: 1)'
    )

    parser.add_argument(
        '--min-frequency',
        type=float,
        default=0.2,
        help='Minimum frequency (0.0-1.0) to report domain/cluster (default: 0.2 = 20%%)'
    )

    parser.add_argument(
        '--mysql-config',
        default=DEFAULT_MYSQL_CONFIG,
        help=f'Path to MySQL config file (default: {DEFAULT_MYSQL_CONFIG})'
    )

    parser.add_argument(
        '--skip-species-summary',
        action='store_true',
        help='Skip running species_summary_new.pl (use existing seq_info)'
    )

    parser.add_argument(
        '--output',
        help='Output file for report (default: print to console only)'
    )

    args = parser.parse_args()

    # Validate Pfam directory
    if not os.path.isdir(args.pfam_dir):
        print(f"Error: Pfam directory not found: {args.pfam_dir}")
        sys.exit(1)

    # Validate STRING data directory
    if not os.path.isdir(args.string_data_dir):
        print(f"Error: STRING data directory not found: {args.string_data_dir}")
        print(f"Please ensure STRING data files are available at this location")
        sys.exit(1)

    # Validate parameters
    if not 0 <= args.min_frequency <= 1:
        print(f"Error: --min-frequency must be between 0.0 and 1.0")
        sys.exit(1)

    if args.network_depth < 1:
        print(f"Error: --network-depth must be >= 1")
        sys.exit(1)

    # Run species_summary_new.pl if requested
    if not args.skip_species_summary:
        if not run_species_summary(args.pfam_dir):
            print("Warning: species_summary_new.pl failed, continuing anyway...")

    # Run analysis
    print("\n" + "=" * 80)
    print(f"Analyzing STRING networks for: {args.pfam_dir}")
    print(f"STRING data directory: {args.string_data_dir}")
    print(f"Score threshold: {args.score_threshold}")
    print(f"Network depth: {args.network_depth}")
    print(f"Minimum frequency: {args.min_frequency:.1%}")
    print("=" * 80)

    results = analyze_networks(
        args.pfam_dir,
        args.string_data_dir,
        args.score_threshold,
        args.network_depth,
        args.mysql_config
    )

    # Print report
    print_report(results, args.min_frequency, args.output)


if __name__ == '__main__':
    main()
