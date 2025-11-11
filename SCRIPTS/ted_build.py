#!/usr/bin/env python3
"""
Script to create Pfam curation directories from TED domains
Features:
- Downloads TED domain PDB files
- Skips domains with >50% overlap with Pfam domains
"""

import sys
import os
import json
import argparse
import subprocess
from pathlib import Path
import requests
from typing import Dict, List, Tuple, Optional

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Fetches domain information from TED API, creates directories for each domain, '
                    'downloads PDB files.'
    )
    parser.add_argument('uniprot_acc', help='UniProt accession')
    parser.add_argument('--afdb-swissprot', action='store_true', 
                        help='Include afdb-swissprot in FoldSeek search')
    parser.add_argument('--afdb-proteome', action='store_true',
                        help='Include afdb-proteome in FoldSeek search')
    parser.add_argument('--afdb50', action='store_true',
                        help='Include afdb50 in FoldSeek search')
    parser.add_argument('--bfvd', action='store_true',
                        help='Include BFVD in FoldSeek search')
    
    return parser.parse_args()

def download_ted_pdb(ted_id: str, output_file: str) -> bool:
    """Download TED domain PDB file"""
    pdb_url = f"https://ted.cathdb.info/api/v1/files/{ted_id}.pdb"
    
    print(f"Downloading PDB file from {pdb_url}")
    
    try:
        response = requests.get(pdb_url, timeout=30)
        response.raise_for_status()
        
        with open(output_file, 'wb') as f:
            f.write(response.content)
        
        print(f"PDB file downloaded successfully to {output_file}")
        return True
    except requests.exceptions.RequestException as e:
        print(f"Error: Failed to download PDB file for {ted_id}: {e}")
        return False

def check_overlap_with_pfam(ted_start: int, ted_end: int, 
                           pfam_domains: List[Dict]) -> Dict:
    """Check for overlap with Pfam domains"""
    ted_start = int(ted_start)
    ted_end = int(ted_end)
    ted_length = ted_end - ted_start + 1
    
    print(f"  Checking TED domain {ted_start}-{ted_end} (length: {ted_length})")
    
    for pfam in pfam_domains:
        pfam_start = int(pfam['start'])
        pfam_end = int(pfam['end'])
        
        print(f"    Checking against Pfam {pfam['pfam_id']} at {pfam_start}-{pfam_end}")
        
        # Calculate overlap
        overlap_start = max(ted_start, pfam_start)
        overlap_end = min(ted_end, pfam_end)
        
        if overlap_start <= overlap_end:
            overlap_length = overlap_end - overlap_start + 1
            overlap_percentage = (overlap_length / ted_length) * 100
            
            print(f"      Overlap found: {overlap_start}-{overlap_end} "
                  f"(length: {overlap_length}, percentage: {overlap_percentage:.1f}%)")
            
            if overlap_percentage > 50:
                return {
                    'overlaps': True,
                    'pfam_id': pfam['pfam_id'],
                    'pfam_name': pfam['pfam_name'],
                    'overlap_percentage': overlap_percentage
                }
    
    return {'overlaps': False}

def fetch_ted_domains(uniprot_acc: str) -> Dict:
    """Fetch TED domains from API"""
    ted_api_url = f"https://ted.cathdb.info/api/v1/uniprot/summary/{uniprot_acc}?skip=0&limit=100"
    
    print("Fetching TED domains...")
    
    try:
        response = requests.get(ted_api_url, timeout=30)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error: Failed to fetch domain information from TED for {uniprot_acc}: {e}")
        sys.exit(1)
    except json.JSONDecodeError:
        print(f"Error: Failed to parse domain information for {uniprot_acc}.")
        sys.exit(1)

def fetch_pfam_domains(uniprot_acc: str) -> List[Dict]:
    """Fetch Pfam domains from InterPro API"""
    interpro_api_url = f"https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/{uniprot_acc}"
    
    print("Fetching Pfam domains from InterPro...")
    
    pfam_domains = []
    
    try:
        response = requests.get(interpro_api_url, timeout=30)
        response.raise_for_status()
        interpro_data = response.json()
        
        # Extract Pfam domains
        if 'results' in interpro_data and isinstance(interpro_data['results'], list):
            for result in interpro_data['results']:
                if ('metadata' in result and 
                    'member_databases' in result['metadata'] and 
                    'pfam' in result['metadata']['member_databases']):
                    
                    pfam_dict = result['metadata']['member_databases']['pfam']
                    pfam_id = list(pfam_dict.keys())[0]
                    pfam_name = pfam_dict[pfam_id]
                    
                    print(f"Found Pfam domain: {pfam_id} ({pfam_name})")
                    
                    # Extract positions
                    if 'proteins' in result and isinstance(result['proteins'], list):
                        for protein in result['proteins']:
                            if protein['accession'].lower() == uniprot_acc.lower():
                                if ('entry_protein_locations' in protein and 
                                    isinstance(protein['entry_protein_locations'], list)):
                                    
                                    for location in protein['entry_protein_locations']:
                                        if 'fragments' in location and isinstance(location['fragments'], list):
                                            for fragment in location['fragments']:
                                                if 'start' in fragment and 'end' in fragment:
                                                    pfam_domains.append({
                                                        'pfam_id': pfam_id,
                                                        'pfam_name': pfam_name,
                                                        'start': int(fragment['start']),
                                                        'end': int(fragment['end'])
                                                    })
                                                    print(f"  Position: {fragment['start']}-{fragment['end']}")
    
    except requests.exceptions.RequestException as e:
        print(f"Warning: Failed to fetch InterPro information for {uniprot_acc}: {e}. "
              "Continuing without Pfam checking.")
    except (json.JSONDecodeError, KeyError) as e:
        print(f"Warning: Failed to parse InterPro information for {uniprot_acc}: {e}. "
              "Continuing without Pfam checking.")
    
    return pfam_domains

def run_command(cmd: str, error_msg: str = None):
    """Run shell command and handle errors"""
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
        return result
    except subprocess.CalledProcessError as e:
        if error_msg:
            print(f"{error_msg}: {e}")
        raise

def main():
    args = parse_arguments()
    uniprot_acc = args.uniprot_acc
    
    print(f"Fetching domain information for UniProt accession: {uniprot_acc}")
    
    # Fetch TED domains
    domain_data = fetch_ted_domains(uniprot_acc)
    
    # Extract domain count
    domain_count = domain_data.get('count', 0)
    if domain_count == 0:
        print(f"No domains found for {uniprot_acc} in TED.")
        sys.exit(0)
    
    print(f"Found {domain_count} domains for {uniprot_acc} in TED.")
    
    # Fetch Pfam domains
    pfam_domains = fetch_pfam_domains(uniprot_acc)
    print(f"Found {len(pfam_domains)} Pfam domains for {uniprot_acc}.")
    
    # Process each TED domain
    for domain in domain_data.get('data', []):
        ted_id = domain['ted_id']
        chopping = domain['chopping']
        
        # Extract TED part from ted_id
        ted_part = ted_id.split('_')[-1] if '_' in ted_id else ted_id
        
        # Create shorter directory name
        dir_name = f"{uniprot_acc}_{ted_part}"
        
        print(f"Processing domain: {ted_id} ({chopping})")
        
        # Extract start and end points from chopping
        segments = chopping.split('_')
        
        # Get first start point
        start_residue = int(segments[0].split('-')[0])
        
        # Get last end point
        end_residue = int(segments[-1].split('-')[1])
        
        print(f"TED domain coordinates: {start_residue}-{end_residue}")
        print(f"Number of Pfam domains to check against: {len(pfam_domains)}")
        
        # Check if directory already exists in current dir or subdirectories
        subdirs_to_check = ['.', 'OVERLAP', 'IGNORE', 'DONE']
        dir_exists = False
        existing_location = None
        
        for subdir in subdirs_to_check:
            if subdir == '.':
                check_path = Path(dir_name)
            else:
                check_path = Path(subdir) / dir_name
            
            if check_path.exists() and check_path.is_dir():
                dir_exists = True
                existing_location = check_path
                break
        
        if dir_exists:
            print(f"Skipping domain {ted_id}: Directory already exists at {existing_location}")
            continue
        
        # Check for overlap with Pfam domains
        overlap_result = check_overlap_with_pfam(start_residue, end_residue, pfam_domains)
        
        if overlap_result['overlaps']:
            print(f"Skipping domain {ted_id}: Overlaps {overlap_result['overlap_percentage']:.2f}% "
                  f"with Pfam domain {overlap_result['pfam_id']} ({overlap_result['pfam_name']})")
            continue
        
        # Create directory
        Path(dir_name).mkdir(exist_ok=True)
        
        # Store current directory
        current_dir = os.getcwd()
        
        # Change to domain directory
        os.chdir(dir_name)
        
        try:
            # Build pfetch command
            pfetch_cmd = f"pfetch {uniprot_acc}"
            if start_residue:
                pfetch_cmd += f" -s {start_residue}"
            if end_residue:
                pfetch_cmd += f" -e {end_residue}"
            
            # Execute commands
            print(f"Running: {pfetch_cmd} (using region {start_residue}-{end_residue})")
            run_command(f"{pfetch_cmd} > FA", "Failed to run pfetch")
            
            print("Creating alignment...")
            run_command("create_alignment.pl -fasta FA -mu > SEED", "Failed to create alignment")
            
            print("Building profile...")
            run_command("pfbuild -withpfmake", "Failed to build profile")
            
            # Download TED domain PDB file
            pdb_file = f"{ted_id}.pdb"
            download_ted_pdb(ted_id, pdb_file)
            
            print(f"Processing of domain {ted_id} completed.")
            
        finally:
            # Return to original directory
            os.chdir(current_dir)
    
    print("All domains processed successfully.")

if __name__ == "__main__":
    main()
