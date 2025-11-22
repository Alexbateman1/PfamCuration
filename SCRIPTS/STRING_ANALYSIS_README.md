# STRING Network Analysis for Pfam Families

## Overview

This tool automates the analysis of STRING protein interaction networks for Pfam families, particularly useful for identifying potential functions of DUF (Domain of Unknown Function) families.

The pipeline identifies over-represented Pfam domains and UniRef50 clusters among network partners, which can provide clues about protein function through "guilt by association".

## How It Works

1. **Extract STRING IDs**: Parses the `seq_info` file from a Pfam family directory to find STRING protein identifiers
2. **Load Local Data**: Reads local STRING network data files (protein.links.full and protein.aliases) for relevant species
3. **Build Networks**: Constructs protein interaction networks based on STRING confidence scores
4. **Expand Networks**: Optionally expands to multi-hop neighbors (configurable depth)
5. **Map to UniProt**: Maps STRING protein IDs to UniProt accessions using protein.aliases
6. **Find Domains**: Queries the Pfam database to find domains in network partner proteins
7. **Find Clusters**: For proteins without Pfam domains, identifies UniRef50 cluster membership via UniProt API
8. **Calculate Statistics**: Computes enrichment frequencies showing which domains/clusters appear across networks
9. **Report Results**: Generates a report showing only entries above frequency threshold

## Installation

### Prerequisites

- Python 3.6+
- Required Python packages:
  ```bash
  pip install mysql-connector-python
  ```

- Access to:
  - Pfam MySQL database (pfam_live)
  - MySQL config file at `~/.my.cnf`
  - Local STRING data files (see Data Files section below)
  - Internet connection (for UniRef50 cluster lookups only)
  - `species_summary_new.pl` script in PATH (optional, but recommended)

### Setup

1. The script is located at: `SCRIPTS/analyze_string_network.py`

2. Download STRING data files for your species to `/nfs/production/agb/pfam/data/STRING/`

   For each species (taxonomy ID), you need:
   - `{taxid}.protein.links.full.v12.0.txt.gz` - Network links with scores
   - `{taxid}.protein.aliases.v12.0.txt.gz` - UniProt mappings

   Example for species 588602:
   ```bash
   cd /nfs/production/agb/pfam/data/STRING
   wget https://stringdb-downloads.org/download/protein.links.full.v12.0/588602.protein.links.full.v12.0.txt.gz
   wget https://stringdb-downloads.org/download/protein.aliases.v12.0/588602.protein.aliases.v12.0.txt.gz
   ```

3. Make sure your `~/.my.cnf` has proper database credentials:
   ```ini
   [client]
   host=your_db_host
   user=your_username
   password=your_password
   ```

## Usage

### Basic Usage

```bash
./SCRIPTS/analyze_string_network.py /path/to/PfamFamily
```

Example:
```bash
./SCRIPTS/analyze_string_network.py /path/to/PF06267
```

### With Custom Options

```bash
./SCRIPTS/analyze_string_network.py /path/to/PF06267 \
    --score-threshold 700 \
    --min-frequency 0.3 \
    --network-depth 2 \
    --output PF06267_string_analysis.txt
```

### Command Line Options

| Option | Default | Description |
|--------|---------|-------------|
| `pfam_dir` | (required) | Path to Pfam family directory |
| `--string-data-dir` | `/nfs/production/agb/pfam/data/STRING` | Directory containing local STRING data files |
| `--score-threshold` | 400 | Minimum STRING combined score (0-999)<br>400 = medium confidence<br>700 = high confidence |
| `--network-depth` | 1 | Network expansion depth<br>1 = direct neighbors only<br>2 = neighbors of neighbors, etc. |
| `--min-frequency` | 0.2 | Minimum frequency (0.0-1.0) to report<br>Only shows domains/clusters above this threshold |
| `--mysql-config` | `~/.my.cnf` | Path to MySQL config file |
| `--skip-species-summary` | False | Skip running species_summary_new.pl |
| `--output` | (none) | Output file for report |

## Understanding STRING Scores

STRING provides combined confidence scores from 0 to 999:

- **0-149**: Low confidence
- **150-399**: Low-medium confidence
- **400-699**: Medium confidence (default threshold)
- **700-899**: High confidence
- **900-999**: Highest confidence

Higher thresholds give more reliable interactions but fewer results.

## Output Format

The tool generates a report with two main sections:

### 1. Pfam Domain Enrichment

Shows Pfam domains found in network partner proteins, sorted by number of networks:

```
Pfam Domain     Networks     Proteins     Frequency
----------------------------------------------------------------
PF00001         15           45           75.0%
PF00002         12           38           60.0%
...
```

- **Networks**: Number of different query protein networks containing this domain
- **Proteins**: Total number of proteins with this domain across all networks
- **Frequency**: Fraction of query networks containing this domain (only shown if >= min-frequency)

### 2. UniRef50 Cluster Enrichment

Shows UniRef50 clusters for proteins without Pfam domains:

```
UniRef50 Cluster               Networks     Proteins     Frequency
----------------------------------------------------------------
UniRef50_P12345                8            23           40.0%
UniRef50_Q9NY99                5            15           25.0%
...
```

Only entries above the frequency threshold are reported (default 20%).

## STRING Data Files

The tool requires the following local files per species:

1. **protein.links.full.v12.0.txt.gz** (~8-15 MB per species)
   - Contains protein-protein interaction links with detailed subscores
   - Includes: neighborhood, fusion, cooccurrence, homology, coexpression, experiments, database, textmining channels
   - Format: `protein1 protein2 [channel_scores...] combined_score`
   - Example line:
     ```
     588602.SAMN04487991_0001 588602.SAMN04487991_0960 0 0 0 191 0 0 0 0 0 0 0 0 0 191
     ```

2. **protein.aliases.v12.0.txt.gz** (~200-500 KB per species)
   - Maps STRING IDs to external database identifiers (UniProt, RefSeq, etc.)
   - Format: `string_id alias source`
   - Example lines:
     ```
     588602.SAMN04487991_0001    A0A1I3IH28    UniProt_AC
     588602.SAMN04487991_0001    SFI47226.1    RefSeq
     ```

Files should be located in `/nfs/production/agb/pfam/data/STRING/` (directly or in `{taxid}/` subdirectories).

### Downloading STRING Data

To download data for a specific species:

```bash
TAXID=588602  # Replace with your taxonomy ID
cd /nfs/production/agb/pfam/data/STRING

# Download network links (full version with channel subscores)
wget https://stringdb-downloads.org/download/protein.links.full.v12.0/${TAXID}.protein.links.full.v12.0.txt.gz

# Download protein aliases (for UniProt mapping)
wget https://stringdb-downloads.org/download/protein.aliases.v12.0/${TAXID}.protein.aliases.v12.0.txt.gz
```

You can find available species and their taxonomy IDs at: https://string-db.org/cgi/download

## Workflow Example

```bash
# Basic analysis with medium confidence threshold (default)
./SCRIPTS/analyze_string_network.py /nfs/production/PF06267

# High confidence interactions only, higher frequency threshold
./SCRIPTS/analyze_string_network.py /nfs/production/PF06267 \
    --score-threshold 700 \
    --min-frequency 0.3

# Include neighbors of neighbors (2-hop network)
./SCRIPTS/analyze_string_network.py /nfs/production/PF06267 \
    --network-depth 2

# Save report to file
./SCRIPTS/analyze_string_network.py /nfs/production/PF06267 \
    --output PF06267_analysis.txt

# Use existing seq_info without regenerating
./SCRIPTS/analyze_string_network.py /nfs/production/PF06267 \
    --skip-species-summary

# Full custom analysis
./SCRIPTS/analyze_string_network.py /nfs/production/PF06267 \
    --score-threshold 700 \
    --network-depth 2 \
    --min-frequency 0.25 \
    --output PF06267_detailed.txt
```

## Troubleshooting

### "species_summary_new.pl not found"
- The script will continue with existing seq_info file
- Or use `--skip-species-summary` to suppress the warning

### "No STRING identifiers found in seq_info"
- Check that seq_info exists in the Pfam directory
- Verify it contains STRING cross-references (grep "STRING" seq_info)
- Run species_summary_new.pl manually to regenerate

### "ERROR: Links file not found for species {taxid}"
- Download the required STRING files for that species
- Files must be named: `{taxid}.protein.links.full.v12.0.txt.gz` and `{taxid}.protein.aliases.v12.0.txt.gz`
- Place files in `/nfs/production/agb/pfam/data/STRING/` or in a `{taxid}/` subdirectory
- Some species may not have STRING data available - check https://string-db.org/cgi/download

### "Database connection failed"
- Check `~/.my.cnf` credentials
- Verify access to pfam_live database
- Ensure mysql-connector-python is installed

### "No network found for protein"
- The protein may not have interactions in STRING
- Try lowering `--score-threshold`
- Check that the protein ID is correct in STRING

## Performance Notes

- STRING data files must be downloaded once per species (~8-15 MB for links, ~200-500 KB for aliases)
- File parsing is optimized with progress indicators for large networks
- Analysis speed depends on:
  - Number of query proteins
  - Number of species
  - Network size (number of interactions above threshold)
  - Network depth (2-hop networks are significantly larger)
  - Database query performance
  - UniRef50 API lookups (can be slow for many proteins without Pfam domains)

Typical runtime:
- Small network (1-10 query proteins, depth=1): 1-3 minutes
- Medium network (10-50 query proteins, depth=1): 3-10 minutes
- Large network (50+ query proteins, depth=1): 10-30 minutes
- 2-hop networks: 2-5x slower than direct neighbors

Tips for faster analysis:
- Use higher `--score-threshold` to reduce network size
- Use `--network-depth 1` for initial exploration
- Download STRING files in advance to avoid network delays

## Data Sources

- **STRING Database**: https://string-db.org/
  - Version: 12.0
  - License: CC BY 4.0

- **Pfam Database**: Via pfam_live MySQL database
  - Tables used: `pfamA_reg_full_significant`, `pfamseq`

- **UniProt**: REST API for UniRef50 cluster lookup
  - Used as fallback when database lookup fails

## Interpreting Results

### High Frequency Domains (>50% of networks)
Suggests strong functional association between query family and the enriched domain. These are prime candidates for functional annotation.

### Medium Frequency (20-50%)
Indicates possible functional relationship worth investigating. May represent:
- Context-specific interactions
- Regulatory relationships
- Pathway components

### Low Frequency (<20%, not shown by default)
May still be biologically relevant but less consistent across networks. Consider lowering `--min-frequency` to explore these.

### Multiple Related Domains
Look for co-enrichment of domains from:
- Same Pfam clan (structural relationship)
- Known pathways or protein complexes
- Similar GO terms or functions

### UniRef Clusters
Proteins without Pfam domains may represent:
- Uncharacterized protein families (potential new Pfam families)
- Species-specific proteins
- Novel domains not yet in Pfam
- Small proteins or fragments
- Disordered regions

High-frequency UniRef clusters are good candidates for:
- Literature search to find functional studies
- Creating new Pfam families
- Experimental characterization

### Network Depth Considerations
- **Depth=1**: Direct physical/functional partners - highest confidence
- **Depth=2**: Extended network - may include pathway members, regulatory relationships
- Higher depth increases sensitivity but may reduce specificity

## Advanced Usage

### Batch Processing Multiple Families

```bash
for dir in /path/to/pfam/families/PF*; do
    family=$(basename $dir)
    echo "Processing $family..."
    ./SCRIPTS/analyze_string_network.py $dir \
        --output ${family}_string_analysis.txt
done
```

### Comparing Different Parameters

```bash
# Compare different score thresholds
for threshold in 400 700 900; do
    ./SCRIPTS/analyze_string_network.py /path/to/PF06267 \
        --score-threshold $threshold \
        --output PF06267_score${threshold}.txt
done

# Compare different network depths
for depth in 1 2 3; do
    ./SCRIPTS/analyze_string_network.py /path/to/PF06267 \
        --network-depth $depth \
        --output PF06267_depth${depth}.txt
done

# Progressive frequency thresholds to explore results
for freq in 0.1 0.2 0.3 0.5; do
    ./SCRIPTS/analyze_string_network.py /path/to/PF06267 \
        --min-frequency $freq \
        --output PF06267_freq${freq}.txt
done
```

## Citation

If you use this tool, please cite:

- **STRING**: Szklarczyk D, et al. (2023) "The STRING database in 2023: protein-protein association networks and functional enrichment analyses for any sequenced genome of interest." Nucleic Acids Res. 51:D638-D646.

- **Pfam**: Mistry J, et al. (2021) "Pfam: The protein families database in 2021." Nucleic Acids Res. 49:D412-D419.

## Support

For issues or questions:
1. Check this README
2. Examine the verbose output for specific errors
3. Contact the Pfam curation team

## Author

Created for automated STRING network analysis of Pfam families.
