# STRING Network Analysis for Pfam Families

## Overview

This tool automates the analysis of STRING protein interaction networks for Pfam families, particularly useful for identifying potential functions of DUF (Domain of Unknown Function) families.

The pipeline identifies over-represented Pfam domains and UniRef50 clusters among network partners, which can provide clues about protein function through "guilt by association".

## How It Works

1. **Extract STRING IDs**: Parses the `seq_info` file from a Pfam family directory to find STRING protein identifiers
2. **Download Data**: Downloads and caches STRING network data files for relevant species
3. **Build Networks**: Constructs protein interaction networks based on STRING confidence scores
4. **Map to UniProt**: Maps STRING protein IDs to UniProt accessions
5. **Find Domains**: Queries the Pfam database to find domains in network partner proteins
6. **Find Clusters**: For proteins without Pfam domains, identifies UniRef50 cluster membership
7. **Calculate Statistics**: Computes enrichment statistics showing which domains/clusters appear in multiple networks
8. **Report Results**: Generates a report excluding singletons

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
  - Internet connection (for downloading STRING data)
  - `species_summary_new.pl` script in PATH (optional, but recommended)

### Setup

1. The script is located at: `SCRIPTS/analyze_string_network.py`

2. Ensure the STRING data directory exists:
   ```bash
   mkdir -p /nfs/production/agb/pfam/data/STRING
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
    --output PF06267_string_analysis.txt
```

### Command Line Options

| Option | Default | Description |
|--------|---------|-------------|
| `pfam_dir` | (required) | Path to Pfam family directory |
| `--string-data-dir` | `/nfs/production/agb/pfam/data/STRING` | Directory for caching STRING data files |
| `--score-threshold` | 400 | Minimum STRING combined score (0-999)<br>400 = medium confidence<br>700 = high confidence |
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
Pfam Domain     Networks     Proteins     Network %
----------------------------------------------------------------
PF00001         15           45           75.0%
PF00002         12           38           60.0%
...
```

- **Networks**: Number of different query protein networks containing this domain
- **Proteins**: Total number of proteins with this domain across all networks
- **Network %**: Percentage of query networks containing this domain

### 2. UniRef50 Cluster Enrichment

Shows UniRef50 clusters for proteins without Pfam domains:

```
UniRef50 Cluster               Networks     Proteins     Network %
----------------------------------------------------------------
UniRef50_P12345                8            23           40.0%
UniRef50_Q9NY99                5            15           25.0%
...
```

Only non-singleton entries (appearing in >1 network) are reported.

## STRING Data Files

The tool downloads and caches the following files per species:

1. **protein.links.v12.0.txt.gz** (~5-10 MB per species)
   - Contains protein-protein interaction links with confidence scores
   - Format: `protein1 protein2 combined_score`

2. **protein.aliases.v12.0.txt.gz** (~200-500 KB per species)
   - Maps STRING IDs to external database identifiers (UniProt, etc.)
   - Format: `string_id alias source`

Files are cached in `/nfs/production/agb/pfam/data/STRING/{taxid}/` and only downloaded once.

## Workflow Example

```bash
# Basic analysis with medium confidence threshold
./SCRIPTS/analyze_string_network.py /nfs/production/PF06267

# High confidence interactions only
./SCRIPTS/analyze_string_network.py /nfs/production/PF06267 \
    --score-threshold 700

# Save report to file
./SCRIPTS/analyze_string_network.py /nfs/production/PF06267 \
    --output PF06267_analysis.txt

# Use existing seq_info without regenerating
./SCRIPTS/analyze_string_network.py /nfs/production/PF06267 \
    --skip-species-summary
```

## Troubleshooting

### "species_summary_new.pl not found"
- The script will continue with existing seq_info file
- Or use `--skip-species-summary` to suppress the warning

### "No STRING identifiers found in seq_info"
- Check that seq_info exists in the Pfam directory
- Verify it contains STRING cross-references (grep "STRING" seq_info)
- Run species_summary_new.pl manually to regenerate

### "Error downloading STRING data"
- Check internet connectivity
- Verify the species taxonomy ID is valid
- Some species may not have STRING data available

### "Database connection failed"
- Check `~/.my.cnf` credentials
- Verify access to pfam_live database
- Ensure mysql-connector-python is installed

### "No network found for protein"
- The protein may not have interactions in STRING
- Try lowering `--score-threshold`
- Check that the protein ID is correct in STRING

## Performance Notes

- First run for a species downloads ~5-10 MB of data
- Subsequent runs use cached files (much faster)
- Analysis speed depends on:
  - Number of query proteins
  - Number of species
  - Network size (number of interactions)
  - Database query performance

Typical runtime: 1-5 minutes per species (first run), <1 minute (cached)

## Data Sources

- **STRING Database**: https://string-db.org/
  - Version: 12.0
  - License: CC BY 4.0

- **Pfam Database**: Via pfam_live MySQL database
  - Tables used: `pfamA_reg_full_significant`, `pfamseq`

- **UniProt**: REST API for UniRef50 cluster lookup
  - Used as fallback when database lookup fails

## Interpreting Results

### High Domain Enrichment (>50% of networks)
Suggests strong functional association between query family and the enriched domain.

### Medium Enrichment (20-50%)
Indicates possible functional relationship worth investigating.

### Multiple Related Domains
Look for co-enrichment of domains from the same clan or pathway.

### UniRef Clusters
Proteins without Pfam domains may represent:
- Uncharacterized protein families
- Species-specific proteins
- Novel domains not yet in Pfam

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

### Comparing Different Thresholds

```bash
for threshold in 400 700 900; do
    ./SCRIPTS/analyze_string_network.py /path/to/PF06267 \
        --score-threshold $threshold \
        --output PF06267_score${threshold}.txt
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
