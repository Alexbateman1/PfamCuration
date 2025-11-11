#!/usr/bin/env python3
"""
Query UniProt additional bibliography for proteins in Pfam ALIGN file

This script extracts bibliography entries from UniProt's additional bibliography
file for all proteins found in a Pfam ALIGN file.
"""

import argparse
import os
import sys
import re
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict


def extract_accessions_from_align(align_file):
    """Extract UniProt accessions from ALIGN file"""
    accessions = set()
    
    with open(align_file, 'r') as f:
        for line in f:
            # Stockholm format: accession/start-end or accession.version/start-end
            if line.startswith('#') or not line.strip():
                continue
            
            # Extract accession from first field
            match = re.match(r'^(\S+)', line)
            if match:
                full_id = match.group(1)
                # Extract just the accession (before . or /)
                accession = re.split(r'[./]', full_id)[0]
                accessions.add(accession)
    
    return accessions


def query_bibliography(bibl_file, accessions):
    """Query bibliography file for given accessions using grep for speed"""
    results = defaultdict(list)
    
    if not accessions:
        return results
    
    # Write accessions to temporary file for grep -f
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as tmp:
        tmp_file = tmp.name
        for acc in accessions:
            # Write just the accession for grep
            tmp.write(f"{acc}\n")
    
    try:
        # Use grep -F -f for fast filtering
        result = subprocess.run(
            ['grep', '-F', '-f', tmp_file, bibl_file],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='replace'
        )
        
        # Build set for filtering
        accession_set = set(accessions)
        
        # Parse grep output and filter to ensure accession is at start
        for line in result.stdout.split('\n'):
            if not line.strip():
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            
            accession = fields[0]
            
            # Only include if accession matches exactly
            if accession not in accession_set:
                continue
            
            source = fields[1] if len(fields) > 1 else ''
            pmid = fields[2] if len(fields) > 2 else ''
            gene_id = fields[3] if len(fields) > 3 else ''
            annotation = fields[4] if len(fields) > 4 else ''
            
            results[accession].append({
                'source': source,
                'pmid': pmid,
                'gene_id': gene_id,
                'annotation': annotation
            })
    
    finally:
        # Clean up temporary file
        if os.path.exists(tmp_file):
            os.remove(tmp_file)
    
    return results


def format_output(results, max_papers_per_protein=5):
    """Format results for output"""
    output_lines = []
    
    # Sort by accession
    for accession in sorted(results.keys()):
        entries = results[accession]
        
        # Group by PMID to avoid duplicates
        pmid_entries = {}
        for entry in entries:
            pmid = entry['pmid']
            if pmid not in pmid_entries:
                pmid_entries[pmid] = entry
        
        # Limit to max papers per protein
        unique_entries = list(pmid_entries.values())[:max_papers_per_protein]
        
        for entry in unique_entries:
            output_lines.append("=" * 80)
            output_lines.append(f"Protein: {accession}")
            output_lines.append(f"PMID: {entry['pmid']}")
            output_lines.append(f"Source: {entry['source']}")
            
            if entry['gene_id']:
                output_lines.append(f"Gene ID: {entry['gene_id']}")
            
            if entry['annotation']:
                # Clean up annotation
                annotation = entry['annotation'].strip()
                output_lines.append(f"Annotation: {annotation}")
            
            output_lines.append("")
    
    return "\n".join(output_lines)


def main():
    parser = argparse.ArgumentParser(
        description='Query UniProt additional bibliography for proteins in ALIGN file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  add_abb_ref.py
  add_abb_ref.py --align ALIGN --max-papers 10
  add_abb_ref.py --force

Environment Variable:
  UNIPROT_BIBL    Path to UniProt additional bibliography file
                  (default: /nfs/production/agb/pfam/users/agb/ABB/add_bibl_info.tb)

Output:
  Results are written to a file called 'uniprot_bibl' in the current directory
        """
    )
    
    parser.add_argument('-align', '--align',
                        default='ALIGN',
                        help='ALIGN file (default: ./ALIGN)')
    
    parser.add_argument('-bibl', '--bibliography',
                        default=os.environ.get('UNIPROT_BIBL',
                                             '/nfs/production/agb/pfam/users/agb/ABB/add_bibl_info.tb'),
                        help='Path to UniProt bibliography file')
    
    parser.add_argument('-max_papers', '--max-papers',
                        type=int, default=5,
                        help='Maximum papers per protein (default: 5)')
    
    parser.add_argument('-force', '--force',
                        action='store_true',
                        help='Force re-run even if results already exist')
    
    args = parser.parse_args()
    
    output_file = 'uniprot_bibl'
    
    # Check if already run (unless --force is specified)
    if os.path.exists(output_file) and not args.force:
        print(f"UniProt bibliography results already exist ({output_file})", file=sys.stderr)
        print("Use --force to re-run", file=sys.stderr)
        sys.exit(0)
    
    # Check required files exist
    if not os.path.exists(args.align):
        print(f"ERROR: ALIGN file {args.align} not found", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.bibliography):
        print(f"ERROR: UniProt bibliography file {args.bibliography} not found", file=sys.stderr)
        sys.exit(1)
    
    # Extract accessions from ALIGN
    print("Extracting accessions from ALIGN file...", file=sys.stderr)
    accessions = extract_accessions_from_align(args.align)
    print(f"Found {len(accessions)} unique accessions", file=sys.stderr)
    
    if not accessions:
        print("No accessions found in ALIGN file", file=sys.stderr)
        sys.exit(0)
    
    # Query bibliography
    print("Querying UniProt bibliography...", file=sys.stderr)
    results = query_bibliography(args.bibliography, accessions)
    
    # Count total entries
    total_entries = sum(len(entries) for entries in results.values())
    print(f"Found bibliography entries for {len(results)} proteins ({total_entries} total entries)", file=sys.stderr)
    
    if not results:
        print("No bibliography entries found for proteins in ALIGN", file=sys.stderr)
        # Create empty file to mark that we've run
        with open(output_file, 'w') as f:
            f.write("# No bibliography entries found\n")
        sys.exit(0)
    
    # Format and write output
    print(f"Writing results to {output_file}...", file=sys.stderr)
    output = format_output(results, args.max_papers)
    
    with open(output_file, 'w') as f:
        f.write(output)
    
    print("Done!", file=sys.stderr)
    print(f"Results written to: {output_file}", file=sys.stderr)
    print(f"To re-run, use: add_abb_ref.py --force", file=sys.stderr)


if __name__ == '__main__':
    main()
