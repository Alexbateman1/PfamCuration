#!/usr/bin/env python3
"""
Query PaperBLAST database using an HMM profile

This script searches a Pfam HMM against the PaperBLAST protein database
and retrieves relevant scientific papers for the top hits.
"""

import argparse
import os
import sqlite3
import subprocess
import sys
from pathlib import Path


def parse_tblout(tblout_file):
    """Parse hmmsearch tblout output file"""
    hits = []
    
    with open(tblout_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            fields = line.split()
            target = fields[0]
            evalue = float(fields[4])
            score = float(fields[5])
            bias = float(fields[6])
            
            hits.append({
                'target': target,
                'evalue': evalue,
                'score': score,
                'bias': bias
            })
    
    # Sort by E-value
    hits.sort(key=lambda x: x['evalue'])
    
    return hits


def get_gene_info(db_path, protein_id):
    """Get additional information about a gene/protein"""
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    
    info = {}
    
    # Get gene description
    cursor.execute("SELECT * FROM Gene WHERE geneId = ?", (protein_id,))
    gene_row = cursor.fetchone()
    if gene_row:
        info['description'] = gene_row['desc']
        info['organism'] = gene_row['organism']
    
    # Get curated information (often has better names/descriptions)
    cursor.execute("SELECT db, name, desc FROM CuratedGene WHERE protId = ?", (protein_id,))
    curated_rows = cursor.fetchall()
    if curated_rows:
        info['curated'] = [dict(row) for row in curated_rows]
    
    conn.close()
    return info


def query_papers(db_path, protein_id):
    """Query SQLite database for papers linked to a protein"""
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    
    # The geneId in PaperBLAST IS the protein identifier (UniProt accession, locus tag, etc.)
    gene_id = protein_id
    
    # Handle versioned identifiers (e.g., A0A0H3JPI5.1 -> A0A0H3JPI5)
    if '.' in gene_id and gene_id.split('.')[-1].isdigit():
        gene_id = gene_id.rsplit('.', 1)[0]
    
    results = []
    seen_pmids = set()
    
    # Query GenePaper for text-mined papers
    # queryTerm is important - it shows the locus tag or gene name that was actually mentioned in the paper
    try:
        cursor.execute("""
            SELECT gp.pmId as pubmedId, gp.title, gp.journal, gp.year, 
                   gp.queryTerm, s.snippet, 'EuropePMC' as source
            FROM GenePaper gp
            LEFT JOIN Snippet s ON gp.geneId = s.geneId AND gp.pmId = s.pmId
            WHERE gp.geneId = ?
            ORDER BY gp.year DESC
        """, (gene_id,))
        
        for row in cursor.fetchall():
            pmid = row['pubmedId']
            if pmid and pmid not in seen_pmids:
                seen_pmids.add(pmid)
                results.append(dict(row))
    except sqlite3.Error as e:
        print(f"Warning: Error querying GenePaper: {e}", file=sys.stderr)
    
    # Query CuratedPaper for curated annotations (using protId from CuratedGene)
    try:
        cursor.execute("""
            SELECT cp.pmId as pubmedId, cg.name as title, 
                   cg.comment as snippet, cg.db as source, cg.id2 as queryTerm
            FROM CuratedGene cg
            JOIN CuratedPaper cp ON cg.db = cp.db AND cg.protId = cp.protId
            WHERE cg.protId = ?
        """, (gene_id,))
        
        for row in cursor.fetchall():
            pmid = row['pubmedId']
            if pmid and pmid not in seen_pmids:
                seen_pmids.add(pmid)
                results.append(dict(row))
    except sqlite3.Error as e:
        print(f"Warning: Error querying CuratedPaper: {e}", file=sys.stderr)
    
    conn.close()
    return results


def get_identifiers(target):
    """Generate possible identifier variants for database lookup"""
    identifiers = [target]
    
    # Handle versioned identifiers (e.g., Q9X123.1 -> Q9X123)
    if '.' in target:
        base = target.rsplit('.', 1)[0]
        identifiers.append(base)
    
    return identifiers


def print_paper_info(protein, evalue, score, paper, gene_info=None, file_handle=None):
    """Print formatted paper information"""
    # Default to stdout if no file handle provided
    f = file_handle if file_handle else sys.stdout
    
    print("=" * 80, file=f)
    print(f"Protein: {protein} (E-value: {evalue}, Score: {score})", file=f)
    
    # Show the query term (locus tag/gene name mentioned in the paper)
    if paper.get('queryTerm'):
        print(f"Mentioned as: {paper['queryTerm']}", file=f)
    
    # Add gene/protein description if available
    if gene_info:
        if gene_info.get('description'):
            print(f"Description: {gene_info['description']}", file=f)
        if gene_info.get('organism'):
            print(f"Organism: {gene_info['organism']}", file=f)
        if gene_info.get('curated'):
            for curated in gene_info['curated']:
                curated_info = f"  {curated['db']}: {curated.get('name', '')}"
                if curated.get('desc'):
                    curated_info += f" - {curated['desc']}"
                print(curated_info, file=f)
    
    pmid = paper.get('pubmedId')
    if pmid:
        print(f"PMID: {pmid}", file=f)
    
    source = paper.get('source')
    if source:
        print(f"Source: {source}", file=f)
    
    title = paper.get('title')
    if title:
        title = ' '.join(title.split())
        print(f"Title: {title}", file=f)
    
    year = paper.get('year')
    journal = paper.get('journal')
    if journal and year:
        print(f"Journal: {journal} ({year})", file=f)
    
    if paper.get('snippet'):
        snippet = ' '.join(paper['snippet'].split())
        if len(snippet) > 500:
            snippet = snippet[:500] + "..."
        print(f"Snippet: {snippet}", file=f)
    
    print(file=f)


def main():
    parser = argparse.ArgumentParser(
        description='Query PaperBLAST database using an HMM profile',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Environment Variables:
  PAPERBLAST_DB       Path to PaperBLAST FASTA database (uniq.faa)
  PAPERBLAST_SQLITE   Path to PaperBLAST SQLite database (litsearch.db)

Example:
  query_paperblast.py
  query_paperblast.py -hmm HMM -max_hits 50 -max_papers 5
  
Output:
  Results are written to a file called 'paperblast' in the current directory
        """
    )
    
    parser.add_argument('-hmm', '--hmm', 
                        default='HMM',
                        help='HMM file (default: ./HMM)')
    
    parser.add_argument('-db', '--database',
                        default=os.environ.get('PAPERBLAST_DB', 
                                             '/nfs/production/agb/pfam/paperblast/uniq.faa'),
                        help='Path to PaperBLAST FASTA database')
    
    parser.add_argument('-sqlite', '--sqlite',
                        default=os.environ.get('PAPERBLAST_SQLITE',
                                             '/nfs/production/agb/pfam/paperblast/litsearch.db'),
                        help='Path to PaperBLAST SQLite database')
    
    parser.add_argument('-max_hits', '--max-hits',
                        type=int, default=20,
                        help='Maximum protein hits to check (default: 20)')
    
    parser.add_argument('-max_papers', '--max-papers',
                        type=int, default=3,
                        help='Maximum papers per protein (default: 3)')
    
    parser.add_argument('-evalue', '--evalue',
                        type=float, default=0.001,
                        help='E-value threshold for hmmsearch (default: 0.001)')
    
    parser.add_argument('-force', '--force',
                        action='store_true',
                        help='Force re-run even if results already exist')
    
    args = parser.parse_args()
    
    # Check if already run (unless --force is specified)
    tblout = 'paperblast_hmmsearch.tblout'
    if os.path.exists(tblout) and not args.force:
        print(f"PaperBLAST hmmsearch results already exist ({tblout})", file=sys.stderr)
        print("Use --force to re-run", file=sys.stderr)
        sys.exit(0)
    
    # Check required files exist
    if not os.path.exists(args.hmm):
        print(f"ERROR: HMM file {args.hmm} not found", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.database):
        print(f"ERROR: PaperBLAST FASTA database {args.database} not found", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.sqlite):
        print(f"ERROR: PaperBLAST SQLite database {args.sqlite} not found", file=sys.stderr)
        sys.exit(1)
    
    # Run hmmsearch
    tblout = 'paperblast_hmmsearch.tblout'
    hmmsearch_cmd = [
        'hmmsearch',
        '-E', str(args.evalue),
        '--tblout', tblout,
        '--cpu', '4',
        args.hmm,
        args.database
    ]
    
    print("Running hmmsearch against PaperBLAST database...", file=sys.stderr)
    
    try:
        result = subprocess.run(hmmsearch_cmd, 
                              stdout=subprocess.DEVNULL,
                              stderr=subprocess.DEVNULL,
                              check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: hmmsearch failed: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Parse results
    print("Parsing hmmsearch results...", file=sys.stderr)
    hits = parse_tblout(tblout)
    
    if not hits:
        print("No significant hits found in PaperBLAST database", file=sys.stderr)
        sys.exit(0)
    
    # Limit to top hits
    hits = hits[:args.max_hits]
    
    print(f"Found {len(hits)} significant hits", file=sys.stderr)
    
    # Open output file
    output_file = 'paperblast'
    print(f"Writing results to {output_file}...", file=sys.stderr)
    
    with open(output_file, 'w') as out_f:
        # Query for papers
        print("Querying PaperBLAST literature database...", file=sys.stderr)
        
        papers_seen = set()
        total_papers = 0
        
        for hit in hits:
            target = hit['target']
            evalue = hit['evalue']
            score = hit['score']
            
            # Get gene information for better context
            gene_info = get_gene_info(args.sqlite, target)
            
            papers_found = 0
            
            # Try different identifier formats
            for identifier in get_identifiers(target):
                papers = query_papers(args.sqlite, identifier)
                
                for paper in papers:
                    pmid = paper.get('pubmedId')
                    if not pmid or pmid in papers_seen:
                        continue
                    
                    papers_seen.add(pmid)
                    papers_found += 1
                    total_papers += 1
                    
                    # Write to file instead of stdout
                    print_paper_info(target, evalue, score, paper, gene_info, out_f)
                    
                    if papers_found >= args.max_papers:
                        break
                
                if papers_found > 0:
                    break
    
    print(f"Total papers found: {total_papers}", file=sys.stderr)
    
    # Clean up
    if os.path.exists(tblout):
        os.remove(tblout)
    
    print("Done!", file=sys.stderr)
    print(f"Results written to: paperblast", file=sys.stderr)


if __name__ == '__main__':
    main()
