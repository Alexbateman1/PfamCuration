#!/usr/bin/env python3
"""
Analyze STRING protein networks for Pfam families to identify over-represented domains.

This script:
1. Extracts STRING identifiers from a Pfam family's seq_info file
2. Downloads and caches STRING network data
3. Retrieves network partners for each protein
4. Maps proteins to UniProt accessions
5. Looks up Pfam domains for each UniProt accession
6. For proteins without domains, identifies UniRef_50 clusters
7. Calculates enrichment statistics for domains and clusters
8. Reports over-represented domains/clusters (excluding singletons)

Usage:
    analyze_string_network.py <pfam_directory> [options]

Example:
    analyze_string_network.py /path/to/PF06267 --score-threshold 400
"""

import argparse
import gzip
import os
import re
import subprocess
import sys
from collections import defaultdict, Counter
from pathlib import Path
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


def download_string_file(taxid: str, filename: str, output_dir: str,
                         base_url: str = "https://stringdb-downloads.org/download") -> Optional[str]:
    """
    Download a STRING data file for a specific taxon.

    Args:
        taxid: NCBI taxonomy ID
        filename: Name of the file to download (e.g., 'protein.links.v12.0.txt.gz')
        output_dir: Directory to save the downloaded file
        base_url: Base URL for STRING downloads

    Returns:
        Path to downloaded file, or None if download failed
    """
    # Create species-specific subdirectory
    species_dir = os.path.join(output_dir, taxid)
    os.makedirs(species_dir, exist_ok=True)

    # Full filename with taxid
    full_filename = f"{taxid}.{filename}"
    output_path = os.path.join(species_dir, full_filename)

    # Check if file already exists
    if os.path.exists(output_path):
        print(f"  Using cached file: {output_path}")
        return output_path

    # Construct download URL
    url = f"{base_url}/{full_filename}"

    print(f"  Downloading {full_filename}...")
    try:
        urllib.request.urlretrieve(url, output_path)
        print(f"  Downloaded to: {output_path}")
        return output_path
    except urllib.error.HTTPError as e:
        print(f"  Error downloading {url}: HTTP {e.code}")
        return None
    except Exception as e:
        print(f"  Error downloading {url}: {e}")
        return None


def parse_protein_links(links_file: str, score_threshold: int = 400) -> Dict[str, Set[str]]:
    """
    Parse STRING protein.links file to get network connections.

    Args:
        links_file: Path to protein.links.v12.0.txt.gz file
        score_threshold: Minimum combined score (0-999) to include link

    Returns:
        Dictionary mapping protein -> set of connected proteins
    """
    networks = defaultdict(set)

    print(f"  Parsing protein links (threshold={score_threshold})...")

    with gzip.open(links_file, 'rt') as f:
        # Skip header
        header = f.readline()

        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                protein1 = parts[0]
                protein2 = parts[1]
                score = int(parts[2])

                if score >= score_threshold:
                    networks[protein1].add(protein2)
                    networks[protein2].add(protein1)

    return networks


def parse_protein_aliases(aliases_file: str) -> Dict[str, Set[str]]:
    """
    Parse STRING protein.aliases file to map STRING IDs to UniProt accessions.

    Args:
        aliases_file: Path to protein.aliases.v12.0.txt.gz file

    Returns:
        Dictionary mapping STRING protein ID -> set of UniProt accessions
    """
    string_to_uniprot = defaultdict(set)

    print(f"  Parsing protein aliases...")

    with gzip.open(aliases_file, 'rt') as f:
        # Skip header
        header = f.readline()

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                string_id = parts[0]  # e.g., "588602.SAMN04487991_1451"
                alias = parts[1]
                source = parts[2]

                # Look for UniProt accessions (from various sources)
                if 'UniProt' in source or 'Ensembl' in source:
                    # Extract UniProt-like accessions (e.g., P12345, Q9NY99)
                    if re.match(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$', alias):
                        string_to_uniprot[string_id].add(alias)

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


def query_uniref50_cluster(uniprot_acc: str, config_file: str = DEFAULT_MYSQL_CONFIG) -> Optional[str]:
    """
    Query for UniRef50 cluster membership for a UniProt accession.

    Args:
        uniprot_acc: UniProt accession
        config_file: Path to MySQL config file

    Returns:
        UniRef50 cluster ID, or None if not found
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

        # Try to get UniRef50 from pfamseq table
        query = """
            SELECT auto_uniprot_reg_full
            FROM pfamseq
            WHERE pfamseq_acc = %s
        """

        cursor.execute(query, (uniprot_acc,))
        result = cursor.fetchone()

        cursor.close()
        conn.close()

        if result:
            return result[0]

    except Exception as e:
        pass

    # If not found in database, try to fetch from UniProt API
    return fetch_uniref50_from_uniprot(uniprot_acc)


def fetch_uniref50_from_uniprot(uniprot_acc: str) -> Optional[str]:
    """
    Fetch UniRef50 cluster from UniProt REST API.

    Args:
        uniprot_acc: UniProt accession

    Returns:
        UniRef50 cluster ID, or None if not found
    """
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_acc}.txt"

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
                     mysql_config: str = DEFAULT_MYSQL_CONFIG) -> Dict:
    """
    Main analysis function to process STRING networks for a Pfam family.

    Args:
        pfam_dir: Path to Pfam family directory
        string_data_dir: Directory for STRING data files
        score_threshold: Minimum STRING score to include interactions
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

    # Group by taxonomy ID
    taxids = defaultdict(list)
    for taxid, protein_id in string_ids:
        taxids[taxid].append(protein_id)

    print(f"   Covering {len(taxids)} species")

    # Process each species
    for taxid, protein_ids in taxids.items():
        print(f"\n2. Processing species {taxid} ({len(protein_ids)} proteins)...")

        # Download necessary files
        print("   Downloading STRING data files...")
        links_file = download_string_file(taxid, 'protein.links.v12.0.txt.gz', string_data_dir)
        aliases_file = download_string_file(taxid, 'protein.aliases.v12.0.txt.gz', string_data_dir)

        if not links_file or not aliases_file:
            print(f"   Skipping species {taxid} - failed to download data files")
            continue

        # Parse network data
        print("   Parsing network data...")
        networks = parse_protein_links(links_file, score_threshold)
        string_to_uniprot = parse_protein_aliases(aliases_file)

        # Analyze each query protein
        for protein_id in protein_ids:
            full_string_id = f"{taxid}.{protein_id}"

            print(f"\n3. Analyzing network for {full_string_id}...")

            if full_string_id not in networks:
                print(f"   No network found for {full_string_id}")
                continue

            network_partners = networks[full_string_id]
            print(f"   Found {len(network_partners)} network partners")

            results['total_networks'] += 1
            network_id = full_string_id

            # Process each network partner
            for partner_string_id in network_partners:
                results['total_proteins'] += 1

                # Get UniProt accessions for this partner
                uniprot_accs = string_to_uniprot.get(partner_string_id, set())

                if not uniprot_accs:
                    continue

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
                        uniref_cluster = query_uniref50_cluster(uniprot_acc, mysql_config)

                        if uniref_cluster:
                            results['uniref_clusters'][uniref_cluster]['networks'].add(network_id)
                            results['uniref_clusters'][uniref_cluster]['proteins'].add(partner_string_id)

    return results


def print_report(results: Dict, output_file: Optional[str] = None):
    """
    Print analysis report.

    Args:
        results: Results dictionary from analyze_networks()
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
    output.append("")

    # Pfam domains
    output.append("=" * 80)
    output.append("Pfam Domain Enrichment")
    output.append("=" * 80)
    output.append("")

    # Filter out singletons and sort by network count
    pfam_data = [(domain, data['networks'], data['proteins'])
                 for domain, data in results['pfam_domains'].items()
                 if len(data['networks']) > 1]
    pfam_data.sort(key=lambda x: len(x[1]), reverse=True)

    if pfam_data:
        output.append(f"{'Pfam Domain':<15} {'Networks':<12} {'Proteins':<12} {'Network %':<12}")
        output.append("-" * 80)

        for domain, networks, proteins in pfam_data:
            network_pct = (len(networks) / results['total_networks'] * 100) if results['total_networks'] > 0 else 0
            output.append(f"{domain:<15} {len(networks):<12} {len(proteins):<12} {network_pct:>6.1f}%")
    else:
        output.append("No Pfam domains found in multiple networks.")

    output.append("")

    # UniRef clusters
    output.append("=" * 80)
    output.append("UniRef50 Cluster Enrichment (proteins without Pfam domains)")
    output.append("=" * 80)
    output.append("")

    # Filter out singletons and sort by network count
    uniref_data = [(cluster, data['networks'], data['proteins'])
                   for cluster, data in results['uniref_clusters'].items()
                   if len(data['networks']) > 1]
    uniref_data.sort(key=lambda x: len(x[1]), reverse=True)

    if uniref_data:
        output.append(f"{'UniRef50 Cluster':<30} {'Networks':<12} {'Proteins':<12} {'Network %':<12}")
        output.append("-" * 80)

        for cluster, networks, proteins in uniref_data[:50]:  # Top 50
            network_pct = (len(networks) / results['total_networks'] * 100) if results['total_networks'] > 0 else 0
            output.append(f"{cluster:<30} {len(networks):<12} {len(proteins):<12} {network_pct:>6.1f}%")

        if len(uniref_data) > 50:
            output.append(f"\n... and {len(uniref_data) - 50} more clusters")
    else:
        output.append("No UniRef50 clusters found in multiple networks.")

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
        help=f'Directory for STRING data files (default: {DEFAULT_STRING_DATA_DIR})'
    )

    parser.add_argument(
        '--score-threshold',
        type=int,
        default=400,
        help='Minimum STRING combined score (0-999) to include interactions (default: 400 = medium confidence)'
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

    # Create STRING data directory if needed
    os.makedirs(args.string_data_dir, exist_ok=True)

    # Run species_summary_new.pl if requested
    if not args.skip_species_summary:
        if not run_species_summary(args.pfam_dir):
            print("Warning: species_summary_new.pl failed, continuing anyway...")

    # Run analysis
    print("\n" + "=" * 80)
    print(f"Analyzing STRING networks for: {args.pfam_dir}")
    print(f"STRING data directory: {args.string_data_dir}")
    print(f"Score threshold: {args.score_threshold}")
    print("=" * 80)

    results = analyze_networks(
        args.pfam_dir,
        args.string_data_dir,
        args.score_threshold,
        args.mysql_config
    )

    # Print report
    print_report(results, args.output)


if __name__ == '__main__':
    main()
