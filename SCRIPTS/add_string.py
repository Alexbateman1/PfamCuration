#!/usr/bin/env python3
"""
Analyze STRING protein networks for Pfam families to identify over-represented domains.

This script:
1. Extracts STRING identifiers from a Pfam family's seq_info file
2. Downloads missing STRING network data files per species (cached locally)
3. Retrieves network partners for each protein (direct or multi-hop)
4. Maps proteins to UniProt accessions
5. Looks up Pfam domains for each UniProt accession
6. For proteins without domains, records the UniProt accession
7. Calculates enrichment statistics for domains
8. Reports over-represented domains above frequency threshold

Usage:
    add_string.py <pfam_directory> [options]

Example:
    add_string.py /path/to/PF06267 --score-threshold 400 --min-frequency 0.3
"""

import argparse
import gzip
import os
import re
import subprocess
import sys
from collections import defaultdict
from typing import Dict, List, Set, Tuple, Optional
import time
import urllib.request
import urllib.error


# Default paths
DEFAULT_STRING_DATA_DIR = "/nfs/production/agb/pfam/data/STRING"
DEFAULT_MYSQL_CONFIG = "~/.my.cnf"
STRING_DOWNLOAD_BASE_URL = "https://stringdb-downloads.org/download"


def parse_desc_file(pfam_dir: str) -> Tuple[str, str]:
    """
    Parse DESC file to get Pfam accession and ID.

    Args:
        pfam_dir: Path to Pfam family directory

    Returns:
        Tuple of (pfamA_acc, pfamA_id), using placeholders if not found
    """
    desc_path = os.path.join(pfam_dir, 'DESC')
    pfam_acc = 'PFXXXXX'
    pfam_id = 'Unknown'

    if not os.path.exists(desc_path):
        return pfam_acc, pfam_id

    with open(desc_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('AC   '):
                pfam_acc = line[5:].rstrip('.')
            elif line.startswith('ID   '):
                pfam_id = line[5:]

    return pfam_acc, pfam_id


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


def download_string_file(taxid: str, filename: str, output_dir: str, max_retries: int = 3) -> Optional[str]:
    """
    Download a STRING data file for a specific taxon with retries.

    Args:
        taxid: NCBI taxonomy ID
        filename: Name of the file to download (e.g., 'protein.links.full.v12.0.txt.gz')
        output_dir: Directory to save the downloaded file
        max_retries: Maximum number of download attempts

    Returns:
        Path to downloaded file, or None if download failed
    """
    # Full filename with taxid
    full_filename = f"{taxid}.{filename}"
    output_path = os.path.join(output_dir, full_filename)

    # Check if file already exists
    if os.path.exists(output_path):
        return output_path

    # Construct download URL
    # Directory name is filename without .txt.gz extension
    # e.g., protein.links.full.v12.0.txt.gz -> protein.links.full.v12.0
    dir_name = filename.replace('.txt.gz', '')
    url = f"{STRING_DOWNLOAD_BASE_URL}/{dir_name}/{full_filename}"

    print(f"  Downloading {full_filename}...")

    for attempt in range(max_retries):
        try:
            urllib.request.urlretrieve(url, output_path)
            print(f"  Downloaded to: {output_path}")
            return output_path
        except urllib.error.HTTPError as e:
            if e.code == 404:
                print(f"  ERROR: File not available on STRING server (HTTP 404): {url}")
                return None
            else:
                print(f"  Download attempt {attempt + 1} failed: HTTP {e.code}")
                if attempt < max_retries - 1:
                    wait_time = 2 ** attempt
                    print(f"  Retrying in {wait_time} seconds...")
                    time.sleep(wait_time)
        except Exception as e:
            print(f"  Download attempt {attempt + 1} failed: {e}")
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt
                print(f"  Retrying in {wait_time} seconds...")
                time.sleep(wait_time)

    print(f"  ERROR: Failed to download after {max_retries} attempts")
    return None


def get_string_file_path(taxid: str, filename: str, string_data_dir: str, auto_download: bool = True) -> Optional[str]:
    """
    Get path to a STRING data file, downloading if necessary.

    Args:
        taxid: NCBI taxonomy ID
        filename: Name of the file (e.g., 'protein.links.full.v12.0.txt.gz')
        string_data_dir: Base directory for STRING data
        auto_download: Whether to auto-download missing files

    Returns:
        Path to file if available, None otherwise
    """
    full_filename = f"{taxid}.{filename}"
    file_path = os.path.join(string_data_dir, full_filename)

    # Check if file exists
    if os.path.exists(file_path):
        return file_path

    # Try to download if auto_download is enabled
    if auto_download:
        print(f"  File not found locally: {full_filename}")
        return download_string_file(taxid, filename, string_data_dir)

    return None


def parse_protein_links_full(links_file: str, query_proteins: Set[str],
                             score_threshold: int = 400) -> Dict[str, Set[str]]:
    """
    Parse STRING protein.links.full file to get network connections for query proteins.

    Args:
        links_file: Path to protein.links.full.v12.0.txt.gz file
        query_proteins: Set of STRING protein IDs to find networks for (within this species)
        score_threshold: Minimum combined score (0-999) to include link

    Returns:
        Dictionary mapping protein -> set of connected proteins
    """
    networks = defaultdict(set)

    print(f"    Parsing links file (threshold={score_threshold})...")

    with gzip.open(links_file, 'rt') as f:
        # Read header
        header = f.readline().strip().split()

        # Find column indices
        try:
            protein1_idx = header.index('protein1')
            protein2_idx = header.index('protein2')
            combined_idx = header.index('combined_score')
        except ValueError as e:
            print(f"    Error: Expected column not found in header: {e}")
            return networks

        line_count = 0
        link_count = 0

        for line in f:
            line_count += 1
            parts = line.strip().split()

            if len(parts) > max(protein1_idx, protein2_idx, combined_idx):
                protein1 = parts[protein1_idx]
                protein2 = parts[protein2_idx]
                score = int(parts[combined_idx])

                if score >= score_threshold:
                    # Only keep links involving our query proteins
                    if protein1 in query_proteins:
                        networks[protein1].add(protein2)
                        link_count += 1
                    if protein2 in query_proteins:
                        networks[protein2].add(protein1)
                        link_count += 1

            if line_count % 100000 == 0:
                print(f"      Processed {line_count:,} lines, {link_count:,} links found...")

    print(f"    Found {link_count:,} links for query proteins from {line_count:,} total lines")
    return networks


def expand_network_depth(networks: Dict[str, Set[str]], query_protein: str,
                         depth: int, all_links_file: str, score_threshold: int) -> Set[str]:
    """
    Expand network to additional depths by reading links file again.

    Args:
        networks: Existing network dictionary
        query_protein: Starting protein
        depth: Target depth (must be > 1)
        all_links_file: Path to links file
        score_threshold: Score threshold

    Returns:
        Set of all proteins at target depth
    """
    if depth <= 1 or query_protein not in networks:
        return networks.get(query_protein, set())

    current_level = {query_protein}
    all_proteins = {query_protein}

    for hop in range(depth):
        # Get all neighbors of current level
        next_level = set()
        for protein in current_level:
            if protein in networks:
                next_level.update(networks[protein])

        # Remove already visited
        next_level -= all_proteins

        if not next_level:
            break

        # For depths > 1, we need to find neighbors of next_level
        if hop + 1 < depth:
            # Build mini-network for next_level proteins
            temp_networks = parse_protein_links_full(all_links_file, next_level, score_threshold)
            # Merge into main networks
            for protein, neighbors in temp_networks.items():
                networks[protein] = neighbors

        all_proteins.update(next_level)
        current_level = next_level

    return all_proteins - {query_protein}


def parse_protein_aliases(aliases_file: str, string_ids: Set[str]) -> Dict[str, Set[str]]:
    """
    Parse STRING protein.aliases file to get UniProt mappings for specific STRING IDs.

    Args:
        aliases_file: Path to protein.aliases.v12.0.txt.gz file
        string_ids: Set of STRING protein IDs to find mappings for

    Returns:
        Dictionary mapping STRING protein ID -> set of UniProt accessions
    """
    string_to_uniprot = defaultdict(set)

    print(f"    Parsing aliases for {len(string_ids):,} proteins...")

    with gzip.open(aliases_file, 'rt') as f:
        # Read header
        header = f.readline().strip().split('\t')

        try:
            string_id_idx = header.index('#string_protein_id')
            alias_idx = header.index('alias')
            source_idx = header.index('source')
        except ValueError as e:
            print(f"    Error: Expected column not found in header: {e}")
            return string_to_uniprot

        for line in f:
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

    print(f"    Mapped {len(string_to_uniprot):,} STRING IDs to UniProt")
    return string_to_uniprot


def load_pfama_id_mapping(config_file: str = DEFAULT_MYSQL_CONFIG) -> Dict[str, str]:
    """
    Load pfamA_acc -> pfamA_id mapping into memory.

    Args:
        config_file: Path to MySQL config file

    Returns:
        Dictionary mapping pfamA_acc to pfamA_id
    """
    config_file = os.path.expanduser(config_file)

    try:
        query = "SELECT pfamA_acc, pfamA_id FROM pfamA"

        cmd = [
            'mysql',
            f'--defaults-file={config_file}',
            'pfam_live',
            '--quick',
            '--silent',
            '--skip-column-names',
            '-e',
            query
        ]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60
        )

        if result.returncode != 0:
            print(f"WARNING: Failed to load pfamA mapping, will use JOINs instead")
            return {}

        # Parse output
        mapping = {}
        for line in result.stdout.strip().split('\n'):
            if line:
                parts = line.split('\t')
                if len(parts) >= 2:
                    mapping[parts[0]] = parts[1]

        return mapping

    except Exception as e:
        print(f"WARNING: Failed to load pfamA mapping: {e}")
        return {}


def batch_query_pfam_domains(uniprot_accs: Set[str], config_file: str = DEFAULT_MYSQL_CONFIG,
                             pfama_mapping: Optional[Dict[str, str]] = None) -> Dict[str, List[Tuple[str, int, int, str, str]]]:
    """
    Query pfam_live database for Pfam domains for multiple UniProt sequences at once.

    Args:
        uniprot_accs: Set of UniProt accessions to query
        config_file: Path to MySQL config file
        pfama_mapping: Optional pre-loaded pfamA_acc -> pfamA_id mapping

    Returns:
        Dictionary mapping uniprot_acc to list of tuples: (pfamseq_acc, seq_start, seq_end, pfamA_acc, pfamA_id)
    """
    if not uniprot_accs:
        return {}

    config_file = os.path.expanduser(config_file)

    # Build IN clause with proper escaping
    acc_list = "','".join(uniprot_accs)

    try:
        # Query without JOIN if we have the mapping cached
        if pfama_mapping:
            query = f"SELECT r.pfamseq_acc, r.seq_start, r.seq_end, r.pfamA_acc FROM pfamA_reg_full_significant r WHERE r.pfamseq_acc IN ('{acc_list}') AND r.in_full = 1"
        else:
            query = f"SELECT r.pfamseq_acc, r.seq_start, r.seq_end, r.pfamA_acc, a.pfamA_id FROM pfamA_reg_full_significant r JOIN pfamA a ON r.pfamA_acc = a.pfamA_acc WHERE r.pfamseq_acc IN ('{acc_list}') AND r.in_full = 1"

        cmd = [
            'mysql',
            f'--defaults-file={config_file}',
            'pfam_live',
            '--quick',
            '--silent',
            '--skip-column-names',
            '-e',
            query
        ]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300  # Longer timeout for batch query
        )

        if result.returncode != 0:
            print(f"\n         ERROR: mysql returned code {result.returncode}")
            if result.stderr:
                print(f"         {result.stderr}")
            return {}

        # Parse output and group by uniprot_acc
        results = defaultdict(list)
        for line in result.stdout.strip().split('\n'):
            if line:
                parts = line.split('\t')
                if pfama_mapping and len(parts) >= 4:
                    pfamseq_acc = parts[0]
                    seq_start = int(parts[1])
                    seq_end = int(parts[2])
                    pfamA_acc = parts[3]
                    pfamA_id = pfama_mapping.get(pfamA_acc, pfamA_acc)  # Fallback to acc if not in mapping
                    results[pfamseq_acc].append((pfamseq_acc, seq_start, seq_end, pfamA_acc, pfamA_id))
                elif not pfama_mapping and len(parts) >= 5:
                    pfamseq_acc = parts[0]
                    seq_start = int(parts[1])
                    seq_end = int(parts[2])
                    pfamA_acc = parts[3]
                    pfamA_id = parts[4]
                    results[pfamseq_acc].append((pfamseq_acc, seq_start, seq_end, pfamA_acc, pfamA_id))

        return dict(results)

    except subprocess.TimeoutExpired:
        print(f"\n         ERROR: Batch query timeout")
        return {}
    except Exception as e:
        print(f"\n         ERROR in batch query: {e}")
        return {}


def should_skip_analysis(pfam_dir: str, force: bool = False) -> bool:
    """
    Check if analysis should be skipped because STRING file is already up-to-date.

    Args:
        pfam_dir: Path to Pfam family directory
        force: If True, never skip (always re-run)

    Returns:
        True if analysis should be skipped, False otherwise
    """
    if force:
        return False

    string_file = os.path.join(pfam_dir, 'STRING')
    seq_info_file = os.path.join(pfam_dir, 'seq_info')

    # Skip if STRING file exists and is newer than seq_info
    if os.path.exists(string_file) and os.path.exists(seq_info_file):
        string_mtime = os.path.getmtime(string_file)
        seq_info_mtime = os.path.getmtime(seq_info_file)

        if string_mtime > seq_info_mtime:
            return True

    return False


def analyze_networks(pfam_dir: str, string_data_dir: str, score_threshold: int = 400,
                     network_depth: int = 1, mysql_config: str = DEFAULT_MYSQL_CONFIG,
                     debug_limit: Optional[int] = None, auto_download: bool = True, force: bool = False) -> Dict:
    """
    Main analysis function to process STRING networks for a Pfam family.

    Args:
        pfam_dir: Path to Pfam family directory
        string_data_dir: Directory for STRING data files
        score_threshold: Minimum STRING score to include interactions
        network_depth: Number of network hops (1=direct neighbors, 2=neighbors of neighbors)
        mysql_config: Path to MySQL config file
        debug_limit: If set, only process first N networks
        auto_download: Whether to auto-download missing STRING files
        force: If True, re-run even if STRING file is up-to-date

    Returns:
        Dictionary with analysis results
    """
    # Check if we should skip this analysis
    if should_skip_analysis(pfam_dir, force):
        print(f"STRING file is up-to-date for {pfam_dir}, skipping analysis (use --force to override)")
        return None

    results = {
        'pfam_domains': defaultdict(lambda: {'networks': set(), 'proteins': defaultdict(set), 'network_proteins': defaultdict(set), 'pfamA_id': None}),
        'proteins_without_domains': set(),
        'total_networks': 0,
        'total_proteins': 0,
        'query_proteins': [],
        'query_protein_uniprot': {},  # network_id -> set of UniProt accessions
        'query_pfam_acc': None,
        'query_pfam_id': None
    }

    # Get query family info from DESC file
    results['query_pfam_acc'], results['query_pfam_id'] = parse_desc_file(pfam_dir)

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

    if debug_limit:
        print(f"   DEBUG MODE: Processing only first {debug_limit} networks")

    # Load pfamA_acc -> pfamA_id mapping once for efficiency
    print("\n   Loading pfamA mappings...")
    pfama_mapping = load_pfama_id_mapping(mysql_config)
    if pfama_mapping:
        print(f"   Loaded {len(pfama_mapping)} pfamA mappings")

    # Process each species
    networks_processed = 0

    for taxid, protein_ids in taxids.items():
        print(f"\n2. Processing species {taxid} ({len(protein_ids)} proteins)...")

        # Get or download STRING files
        print("   Getting STRING data files...")
        links_file = get_string_file_path(taxid, 'protein.links.full.v12.0.txt.gz',
                                          string_data_dir, auto_download)
        aliases_file = get_string_file_path(taxid, 'protein.aliases.v12.0.txt.gz',
                                           string_data_dir, auto_download)

        if not links_file:
            print(f"   Skipping species {taxid} - links file not available")
            continue

        if not aliases_file:
            print(f"   Skipping species {taxid} - aliases file not available")
            continue

        print(f"   Using: {os.path.basename(links_file)}")
        print(f"   Using: {os.path.basename(aliases_file)}")

        # Build full STRING IDs for this species
        query_proteins_this_species = set()
        for protein_id in protein_ids:
            full_string_id = f"{taxid}.{protein_id}"
            query_proteins_this_species.add(full_string_id)

        # Parse network data
        print("   Parsing network data...")
        networks = parse_protein_links_full(links_file, query_proteins_this_species, score_threshold)

        # Collect all network proteins for UniProt mapping
        all_network_proteins = set()
        for query_protein in query_proteins_this_species:
            if query_protein in networks:
                if network_depth == 1:
                    network_partners = networks[query_protein]
                else:
                    network_partners = expand_network_depth(networks, query_protein, network_depth,
                                                           links_file, score_threshold)
                all_network_proteins.update(network_partners)

        if not all_network_proteins:
            print(f"   No networks found for this species")
            continue

        print(f"   Total network partners across all query proteins: {len(all_network_proteins):,}")

        # Get UniProt mappings for network partners
        print("   Mapping to UniProt...")
        string_to_uniprot = parse_protein_aliases(aliases_file, all_network_proteins)

        # Also get UniProt mappings for query proteins themselves
        query_to_uniprot = parse_protein_aliases(aliases_file, query_proteins_this_species)
        for query_protein, uniprot_accs in query_to_uniprot.items():
            results['query_protein_uniprot'][query_protein] = uniprot_accs

        # Collect all unique UniProt accessions for this species
        all_uniprot_accs = set()
        for string_id_set in string_to_uniprot.values():
            all_uniprot_accs.update(string_id_set)

        # Batch query Pfam domains for all UniProt accessions at once
        print(f"   Querying Pfam domains for {len(all_uniprot_accs)} UniProt accessions...")
        domain_cache = batch_query_pfam_domains(all_uniprot_accs, mysql_config, pfama_mapping)
        print(f"   Found domains for {len(domain_cache)} proteins")

        # Analyze each query protein's network
        print(f"\n3. Processing networks for species {taxid}...")
        for query_protein in query_proteins_this_species:
            if query_protein not in networks:
                continue

            # Check debug limit
            if debug_limit and networks_processed >= debug_limit:
                print(f"\n   DEBUG LIMIT REACHED: Processed {debug_limit} networks, stopping")
                return results

            network_id = query_protein

            if network_depth == 1:
                network_partners = networks[query_protein]
            else:
                network_partners = expand_network_depth(networks, query_protein, network_depth,
                                                       links_file, score_threshold)

            results['total_networks'] += 1
            networks_processed += 1

            partners_with_uniprot = 0
            pfam_domains_found = 0
            proteins_without_domains = 0

            for partner_string_id in network_partners:
                results['total_proteins'] += 1

                # Get UniProt accessions
                uniprot_accs = string_to_uniprot.get(partner_string_id, set())

                if not uniprot_accs:
                    continue

                partners_with_uniprot += 1

                # For each UniProt accession
                for uniprot_acc in uniprot_accs:
                    # Look up domains from cache
                    domains = domain_cache.get(uniprot_acc, [])

                    if domains:
                        pfam_domains_found += len(domains)
                        for pfamseq_acc, seq_start, seq_end, pfamA_acc, pfamA_id in domains:
                            results['pfam_domains'][pfamA_acc]['networks'].add(network_id)
                            results['pfam_domains'][pfamA_acc]['proteins'][partner_string_id].add(uniprot_acc)
                            results['pfam_domains'][pfamA_acc]['network_proteins'][network_id].add(uniprot_acc)
                            results['pfam_domains'][pfamA_acc]['pfamA_id'] = pfamA_id
                    else:
                        # No Pfam domains - record this UniProt accession
                        results['proteins_without_domains'].add(uniprot_acc)
                        proteins_without_domains += 1

            # Print summary for this network
            if networks_processed % 10 == 0:
                print(f"   Processed {networks_processed} networks...")

        # Print species summary
        print(f"   Species {taxid} complete: {networks_processed} networks processed")

    return results


def write_csv_output(results: Dict, pfam_dir: str):
    """
    Write CSV file with Pfam accessions as rows and species/networks as columns.
    Cells contain UniProt accessions for proteins, making it easy to identify
    potential interacting protein pairs for AlphaFold modeling.

    The first data row contains the query family (the family being analyzed).

    Args:
        results: Results dictionary from analyze_networks()
        pfam_dir: Path to Pfam family directory
    """
    output_file = os.path.join(pfam_dir, 'STRING.csv')

    # Collect all unique networks and sort them
    all_networks = set()
    for domain_data in results['pfam_domains'].values():
        all_networks.update(domain_data['networks'])

    all_networks = sorted(all_networks)

    if not all_networks:
        print(f"  No networks found, skipping CSV output")
        return

    # Sort Pfam domains by network count (most frequent first)
    pfam_sorted = sorted(results['pfam_domains'].items(),
                         key=lambda x: len(x[1]['networks']),
                         reverse=True)

    with open(output_file, 'w') as f:
        # Write header row with network IDs as columns
        header = ['Pfam_Acc', 'Pfam_ID'] + all_networks
        f.write(','.join(header) + '\n')

        # Write query family as first data row
        query_row = [results['query_pfam_acc'], results['query_pfam_id']]
        for network_id in all_networks:
            uniprot_accs = results['query_protein_uniprot'].get(network_id, set())
            cell_value = ';'.join(sorted(uniprot_accs)) if uniprot_accs else ''
            query_row.append(cell_value)
        f.write(','.join(query_row) + '\n')

        # Write other Pfam domain rows
        for pfam_acc, domain_data in pfam_sorted:
            pfam_id = domain_data['pfamA_id'] or ''
            row = [pfam_acc, pfam_id]

            # For each network, get the UniProt accessions
            for network_id in all_networks:
                uniprot_accs = domain_data['network_proteins'].get(network_id, set())
                # Join multiple accessions with semicolon
                cell_value = ';'.join(sorted(uniprot_accs)) if uniprot_accs else ''
                row.append(cell_value)

            f.write(','.join(row) + '\n')

    print(f"\n  CSV output written to: {output_file}")


def write_detailed_results(results: Dict, pfam_dir: str):
    """
    Write detailed results file showing ALL domains and clusters found.

    Args:
        results: Results dictionary from analyze_networks()
        pfam_dir: Path to Pfam family directory
    """
    # Write to STRING file in pfam_dir
    output_file = os.path.join(pfam_dir, 'STRING')

    # Remove old STRING directory if it exists
    if os.path.isdir(output_file):
        import shutil
        shutil.rmtree(output_file)

    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("STRING Network Analysis - Detailed Results (ALL findings)\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"Query proteins: {len(results['query_proteins'])}\n")
        f.write(f"Total networks analyzed: {results['total_networks']}\n")
        f.write(f"Total network proteins: {results['total_proteins']}\n\n")

        # All Pfam domains (no threshold)
        f.write("=" * 80 + "\n")
        f.write("ALL Pfam Domains Found (no frequency filter)\n")
        f.write("=" * 80 + "\n\n")

        if results['pfam_domains']:
            # Sort by network count
            pfam_sorted = sorted(results['pfam_domains'].items(),
                               key=lambda x: len(x[1]['networks']),
                               reverse=True)

            f.write(f"{'Pfam Domain':<25} {'Networks':<12} {'Proteins':<12} {'Frequency':<12}\n")
            f.write("-" * 80 + "\n")

            for domain, data in pfam_sorted:
                network_freq = len(data['networks']) / results['total_networks'] if results['total_networks'] > 0 else 0
                domain_display = f"{domain} ({data['pfamA_id']})" if data['pfamA_id'] else domain
                f.write(f"{domain_display:<25} {len(data['networks']):<12} {len(data['proteins']):<12} {network_freq:>6.1%}\n")

            f.write(f"\nTotal unique Pfam domains: {len(results['pfam_domains'])}\n\n")

            # Detailed breakdown
            f.write("\nDetailed Pfam Domain Breakdown:\n")
            f.write("-" * 80 + "\n")
            for domain, data in pfam_sorted:
                domain_display = f"{domain} ({data['pfamA_id']})" if data['pfamA_id'] else domain
                f.write(f"\n{domain_display}:\n")
                f.write(f"  Found in {len(data['networks'])} network(s):\n")
                # List all proteins with their UniProt accessions
                for string_id in sorted(data['proteins'].keys()):
                    uniprot_accs = data['proteins'][string_id]
                    # Show each UniProt accession associated with this STRING ID
                    for uniprot_acc in sorted(uniprot_accs):
                        f.write(f"    - {uniprot_acc} ({string_id})\n")
                f.write(f"  Total proteins: {len(data['proteins'])}\n")
        else:
            f.write("No Pfam domains found in any network partners.\n")

        f.write("\n")

        # Proteins without Pfam domains
        f.write("=" * 80 + "\n")
        f.write("Proteins Without Pfam Domains\n")
        f.write("=" * 80 + "\n\n")

        if results['proteins_without_domains']:
            f.write(f"Total UniProt accessions without Pfam domains: {len(results['proteins_without_domains'])}\n\n")
            for uniprot_acc in sorted(results['proteins_without_domains']):
                f.write(f"  {uniprot_acc}\n")
        else:
            f.write("All proteins have Pfam domain assignments.\n")

        f.write("\n" + "=" * 80 + "\n")

    print(f"\n  Detailed results written to: {output_file}")


def print_report(results: Dict, min_frequency: float = 0.3, output_file: Optional[str] = None):
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
            pfam_data.append((domain, data['pfamA_id'], data['networks'], data['proteins'], network_freq))

    pfam_data.sort(key=lambda x: len(x[2]), reverse=True)

    if pfam_data:
        output.append(f"{'Pfam Domain':<25} {'Networks':<12} {'Proteins':<12} {'Frequency':<12}")
        output.append("-" * 80)

        for domain, pfam_id, networks, proteins, freq in pfam_data:
            domain_display = f"{domain} ({pfam_id})" if pfam_id else domain
            output.append(f"{domain_display:<25} {len(networks):<12} {len(proteins):<12} {freq:>6.1%}")
    else:
        output.append(f"No Pfam domains found above {min_frequency:.1%} frequency threshold.")

    output.append("")
    output.append("=" * 80)
    output.append(f"Proteins Without Pfam Domains: {len(results['proteins_without_domains'])}")
    output.append("=" * 80)
    output.append("")
    output.append("(See detailed results file for full list)")
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
        '--network-depth',
        type=int,
        default=1,
        help='Network depth: 1=direct neighbors, 2=neighbors of neighbors, etc. (default: 1)'
    )

    parser.add_argument(
        '--min-frequency',
        type=float,
        default=0.3,
        help='Minimum frequency (0.0-1.0) to report domain/cluster (default: 0.3 = 30%%)'
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
        '--no-download',
        action='store_true',
        help='Do not auto-download missing STRING files (skip species if files not found)'
    )

    parser.add_argument(
        '--debug-limit',
        type=int,
        help='Debug mode: only process first N networks (for testing)'
    )

    parser.add_argument(
        '--output',
        help='Output file for report (default: print to console only)'
    )

    parser.add_argument(
        '--force',
        action='store_true',
        help='Force re-analysis even if STRING file is up-to-date'
    )

    args = parser.parse_args()

    # Validate Pfam directory
    if not os.path.isdir(args.pfam_dir):
        print(f"Error: Pfam directory not found: {args.pfam_dir}")
        sys.exit(1)

    # Create STRING data directory if needed
    os.makedirs(args.string_data_dir, exist_ok=True)

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
    print(f"Auto-download: {'No' if args.no_download else 'Yes'}")
    if args.debug_limit:
        print(f"DEBUG MODE: Limit to first {args.debug_limit} networks")
    print("=" * 80)

    results = analyze_networks(
        args.pfam_dir,
        args.string_data_dir,
        args.score_threshold,
        args.network_depth,
        args.mysql_config,
        args.debug_limit,
        not args.no_download,
        args.force
    )

    # Write detailed results file and report (if analysis was not skipped)
    if results is not None:
        write_detailed_results(results, args.pfam_dir)
        write_csv_output(results, args.pfam_dir)
        print_report(results, args.min_frequency, args.output)


if __name__ == '__main__':
    main()
