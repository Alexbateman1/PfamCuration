#!/usr/bin/env python3
"""
This curation script identifies if any PDB structures match the family
by using the SIFTS mapping between UniProt and PDB. It then fetches PDB
entries and adds any associated papers to the DESC file.
"""

import sys
import os
import subprocess
import urllib.request
import re

SIFTS_FILE = "/homes/agb/Scripts/pdb_chain_uniprot.csv"

def main():
    print(f"Running {sys.argv[0]}", file=sys.stderr)
    
    # Check if already run
    if os.path.exists(".add_pdb_ref"):
        print("add_pdb_ref.py has already been run in this directory.", file=sys.stderr)
        sys.exit(0)
    
    # Create marker file
    open(".add_pdb_ref", 'w').close()
    
    # Check ALIGN file exists
    if not os.path.exists("ALIGN"):
        print("ALIGN file not found", file=sys.stderr)
        sys.exit(1)
    
    # Extract UniProt accessions and regions from ALIGN
    print("Extracting UniProt accessions from ALIGN", file=sys.stderr)
    align_entries = parse_align_file("ALIGN")
    
    if not align_entries:
        print("No entries found in ALIGN file", file=sys.stderr)
        sys.exit(0)
    
    # Write accessions to temporary file for grep
    with open("uniprot_list.tmp", 'w') as f:
        for acc in align_entries.keys():
            f.write(f"{acc}\n")
    
    # Grep SIFTS file for matching accessions
    print("Searching SIFTS mapping file", file=sys.stderr)
    pdb_matches = find_pdb_matches("uniprot_list.tmp", align_entries)
    
    # Clean up temp file
    os.remove("uniprot_list.tmp")
    
    if not pdb_matches:
        print("No PDB matches found", file=sys.stderr)
        sys.exit(0)
    
    print(f"Found {len(pdb_matches)} PDB structures with overlapping regions", file=sys.stderr)
    
    # Track PMIDs we've already added
    pubmed = set()
    total = 0
    
    # Process each PDB match
    for pdb_id in sorted(set(pdb_matches)):
        result = add_paper(pdb_id, pubmed)
        total += result
        if total > 4:
            print("Found enough structure papers", file=sys.stderr)
            break

def parse_align_file(filename):
    """Parse ALIGN file and extract accession, start, end positions"""
    entries = {}
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # Parse format: A0A0R2GW98.1/17-380
            match = re.match(r'^(\S+?)\.?\d*/(\d+)-(\d+)', line)
            if match:
                accession = match.group(1)
                start = int(match.group(2))
                end = int(match.group(3))
                midpoint = (start + end) / 2
                entries[accession] = {
                    'start': start,
                    'end': end,
                    'midpoint': midpoint
                }
    
    return entries

def find_pdb_matches(uniprot_file, align_entries):
    """Use grep to find matching entries in SIFTS file and check overlap"""
    pdb_matches = []
    
    try:
        # Grep for matching UniProt accessions
        result = subprocess.run(
            ['grep', '-f', uniprot_file, SIFTS_FILE],
            capture_output=True,
            text=True
        )
        
        for line in result.stdout.strip().split('\n'):
            if not line or line.startswith('#'):
                continue
            
            parts = line.split(',')
            if len(parts) < 9:
                continue
            
            pdb_id = parts[0].lower()
            uniprot_acc = parts[2]
            
            # Parse SP_BEG and SP_END (columns 8 and 9)
            try:
                sp_beg = int(parts[7])
                sp_end = int(parts[8])
            except (ValueError, IndexError):
                continue
            
            # Check if this UniProt is in our align entries
            if uniprot_acc in align_entries:
                midpoint = align_entries[uniprot_acc]['midpoint']
                
                # Check if midpoint overlaps with PDB region
                if sp_beg <= midpoint <= sp_end:
                    print(f"Match: {pdb_id} overlaps with {uniprot_acc} (midpoint {midpoint:.1f} in [{sp_beg}-{sp_end}])", file=sys.stderr)
                    pdb_matches.append(pdb_id)
    
    except subprocess.CalledProcessError as e:
        print(f"Error running grep: {e}", file=sys.stderr)
    
    return pdb_matches

def add_paper(pdb_id, pubmed):
    """Fetch PDB file and add paper if PMID found"""
    result = 0
    print(f"Attempting to add paper for {pdb_id}", file=sys.stderr)
    
    url = f"http://files.rcsb.org/view/{pdb_id}.pdb"
    
    try:
        with urllib.request.urlopen(url) as response:
            for line in response:
                line = line.decode('utf-8')
                
                # Look for PMID in JRNL records
                match = re.search(r'JRNL\s+PMID\s+(\d+)', line)
                if match:
                    pmid = match.group(1)
                    
                    if pmid in pubmed:
                        print(f"{pmid} inserted already. Skipping.", file=sys.stderr)
                        return 0
                    else:
                        pubmed.add(pmid)
                        comment = f"RC   Paper describing PDB structure {pdb_id}\n"
                        
                        # Append comment to DESC
                        with open("DESC", 'a') as desc:
                            desc.write(comment)
                        
                        # Add reference using existing Perl script
                        print(f"Adding {pmid} to DESC", file=sys.stderr)
                        subprocess.run(['add_ref.pl', pmid])
                        
                        # Touch .3d file
                        open('.3d', 'w').close()
                        
                        result = 1
                        break
    
    except Exception as e:
        print(f"Cannot fetch {pdb_id}: {e}", file=sys.stderr)
    
    return result

if __name__ == "__main__":
    main()
