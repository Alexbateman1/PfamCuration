# Pfam HHsearch All-Against-All Pipeline

Pipeline for running all-against-all HHsearch profile comparisons across Pfam families.

## Overview

**HHsearch** performs profile-profile comparisons (HMM vs HMM), making it ideal for comparing established Pfam family models to detect:
- Remote homology between families
- Potential overlaps and redundancies
- Clan relationships
- Domain architectures

### Key Differences: HHsearch vs HHblits

| Feature | HHsearch | HHblits |
|---------|----------|---------|
| **Purpose** | Profile-profile comparison | Database searching with iteration |
| **Input** | Pre-built HMM profiles | Sequence or MSA |
| **Iteration** | None (single comparison) | Multiple rounds to build profile |
| **Best for** | Comparing existing families | Finding remote homologs in databases |
| **Sensitivity** | High for profile comparison | High with iteration |

## Installation

### Prerequisites

- Python 3.7+
- HH-suite tools (`hhmake`, `hhsearch`, `reformat.pl`, `ffindex_build`)
- SVN client
- SLURM (for parallel processing)

### Setup

```bash
# Clone or navigate to the repository
cd /path/to/PfamCuration/HHsearch

# Ensure HH-suite is in your PATH
which hhsearch  # Should show path to hhsearch

# Check Python dependencies
python3 -c "import subprocess; from pathlib import Path"
```

## Quick Start

### Full Pipeline Run

```bash
python3 pfam_hhsearch_pipeline.py \
  --work-dir /path/to/work/directory \
  --e-value 1.0
```

This will:
1. Extract SEED files from Pfam SVN
2. Convert to A3M format using `reformat.pl`
3. Build HHM profiles using `hhmake`
4. Create FFindex database
5. Generate SLURM batch scripts
6. **Wait for you to submit SLURM jobs**
7. Aggregate results when jobs complete

### Submit SLURM Jobs

After the pipeline generates the batch script:

```bash
sbatch /path/to/work/directory/RESULTS/run_hhsearch_batch.sh
```

Monitor progress:
```bash
squeue -u $USER
```

### Aggregate Results

Once all jobs complete, re-run the pipeline to aggregate:

```bash
python3 pfam_hhsearch_pipeline.py \
  --work-dir /path/to/work/directory \
  --e-value 1.0
```

## Pipeline Steps

### 1. SVN Discovery
Discovers all Pfam families in the SVN repository.

**Output:** `families_to_process.txt`

### 2. SEED Extraction
Extracts SEED alignment files from SVN.

**Output:** `DATA/SEED/PF*_SEED`

### 3. Database Building
1. Converts SEED (Stockholm) → A3M using `reformat.pl`
2. Builds HHM profiles using `hhmake`
3. Creates FFindex database for fast searching

**Output:**
- `RESULTS/A3M/*.a3m`
- `RESULTS/hhsuite_db/hhm_files/*.hhm`
- `RESULTS/hhsuite_db/pfam_db_hhm.*`

### 4. HHsearch Searches
Runs all-against-all HHsearch comparisons using SLURM.

**Command per family:**
```bash
hhsearch -i query.hhm -d pfam_db_hhm -o output.hhr \
  -cpu 8 -e 1.0 -maxfilt 100000
```

**Output:**
- `RESULTS/HHR/*.hhr` (raw results)
- `RESULTS/PARSED/*_hits.tsv` (parsed TSV)

### 5. Aggregation
Combines all individual results into summary files.

**Output:**
- `RESULTS/SUMMARY/hhsearch_all_vs_all.tsv` (all hits)
- `RESULTS/SUMMARY/family_statistics.tsv` (per-family stats)
- `RESULTS/SUMMARY/high_confidence_overlaps.tsv` (E<1e-10)
- `RESULTS/SUMMARY/pipeline_metadata.json` (run info)

## Command-Line Options

```bash
python3 pfam_hhsearch_pipeline.py [OPTIONS]
```

### Main Options

- `--work-dir DIR` : Working directory (default: ./HHSEARCH_PIPELINE)
- `--e-value FLOAT` : E-value threshold (default: 1.0)
- `--batch-size N` : Families per SLURM batch (default: 100)

### Control Options

- `--full` : Force full run (ignore incremental mode)
- `--no-slurm` : Run sequentially without SLURM (for testing)
- `--reset STEP` : Reset from step: svn_discovery, extraction, db_build, hhsearch, aggregation
- `--no-clean` : Skip file cleanup when using --reset
- `--log-level LEVEL` : Logging level: DEBUG, INFO, WARNING, ERROR

## Advanced Usage

### Resetting Specific Steps

Reset and rebuild the database:
```bash
python3 pfam_hhsearch_pipeline.py \
  --work-dir /path/to/work \
  --reset db_build
```

This automatically cleans:
- A3M files
- HHM profiles
- Database files

Reset just the searches (keep database):
```bash
python3 pfam_hhsearch_pipeline.py \
  --work-dir /path/to/work \
  --reset hhsearch
```

Reset aggregation only:
```bash
python3 pfam_hhsearch_pipeline.py \
  --work-dir /path/to/work \
  --reset aggregation
```

### Reusing Data from HHblits

If you already ran the HHblits pipeline, you can reuse the SEED files:

```bash
# Create symlink to existing SEED directory
ln -s /path/to/HHblits/DATA/SEED /path/to/HHsearch/DATA/SEED

# Run pipeline (will skip extraction)
python3 pfam_hhsearch_pipeline.py --work-dir /path/to/HHsearch
```

You can also reuse A3M files:
```bash
ln -s /path/to/HHblits/RESULTS/A3M /path/to/HHsearch/RESULTS/A3M
```

### Testing on Small Subset

For testing, process only a few families:

```bash
# Edit families_to_process.txt to contain only test families
echo -e "PF00001\nPF00002\nPF00003" > /path/to/work/families_to_process.txt

# Reset extraction to force re-processing
python3 pfam_hhsearch_pipeline.py \
  --work-dir /path/to/work \
  --reset extraction \
  --no-slurm
```

## Output Format

### Main Results File: `hhsearch_all_vs_all.tsv`

| Column | Description |
|--------|-------------|
| query_family | Query Pfam family (e.g., PF00001) |
| target_family | Target Pfam family |
| probability | HHsearch probability (0-100) |
| e_value | E-value of the match |
| p_value | P-value of the match |
| score | HHsearch score |
| ss | Secondary structure score |
| cols | Number of aligned columns |
| query_range | Query alignment range |
| template_range | Target alignment range |
| template_length | Length of target profile |

### Adding Pfam Accessions

To map sequence IDs to Pfam families, use the helper script:

```bash
python3 add_pfam_mapping.py \
  /path/to/RESULTS/A3M \
  /path/to/RESULTS/SUMMARY/hhsearch_all_vs_all.tsv \
  /path/to/RESULTS/SUMMARY/hhsearch_all_vs_all_with_pfam.tsv
```

## Performance Expectations

For ~20,000 Pfam families:

- **Database building:** 2-4 hours (one-time)
- **HHsearch searches:** 20K × 20K comparisons
  - With SLURM (200 jobs): 4-8 hours
  - Sequential: Several weeks
- **Aggregation:** 5-10 minutes
- **Total disk space:** ~50-100 GB

## Troubleshooting

### Database Path Errors

If you see errors like "could not open pfam_db_hhm_cs219.ffdata":
- Check that database was built successfully
- Verify path in SLURM script matches actual database location
- Try `--reset db_build` to rebuild

### Missing HHM Files

If searches fail with "HHM file not found":
- Database building may have failed for some families
- Check `RESULTS/hhsuite_db/hhm_files/` for missing profiles
- Review build logs for conversion errors

### SLURM Jobs Immediately Finish

If jobs exit immediately:
- Check if parsed files already exist (script skips completed families)
- Use `--reset hhsearch` to force re-processing
- Review `.err` files in `RESULTS/slurm_logs/`

### Empty Aggregation Results

If aggregation finds 0 hits:
- Check that parsed files exist in `RESULTS/PARSED/`
- Verify SLURM jobs completed successfully
- Look for errors in `RESULTS/slurm_logs/*.err`

## Comparison with HHblits Pipeline

| Aspect | HHsearch | HHblits |
|--------|----------|---------|
| **Tool** | `hhsearch` | `hhblits` |
| **Comparison type** | Profile-profile | Sequence/MSA to database |
| **Iterations** | None | Multiple (-n 3) |
| **E-value default** | 1.0 | 1e-10 |
| **Best for** | Comparing families | Finding homologs |
| **Speed** | Faster (no iteration) | Slower (builds profile) |
| **Sensitivity** | High for profiles | Very high with iteration |

## Files and Directories

```
HHsearch/
├── pfam_hhsearch_pipeline.py  # Main pipeline
├── hhsearch_runner.py          # HHsearch execution
├── svn_manager.py              # SVN operations
├── add_pfam_mapping.py         # Add Pfam accessions
└── README.md                   # This file

Work directory:
HHSEARCH_PIPELINE/
├── .progress/                  # Progress markers
├── DATA/
│   ├── SEED/                  # SEED alignments
│   └── HMM/                   # HMMER HMMs (not used)
├── RESULTS/
│   ├── A3M/                   # A3M alignments
│   ├── hhsuite_db/            # HH-suite database
│   │   └── hhm_files/         # HMM profiles
│   ├── HHR/                   # HHsearch results
│   ├── PARSED/                # Parsed TSV files
│   ├── SUMMARY/               # Final results
│   ├── batches/               # SLURM batch files
│   └── slurm_logs/            # SLURM output
├── families_to_process.txt    # Family list
└── pipeline.log               # Pipeline log
```

## Citation

If you use this pipeline, please cite:
- **HH-suite:** Steinegger M. et al. (2019) Nucleic Acids Res.
- **Pfam:** Mistry J. et al. (2021) Nucleic Acids Res.

## Support

For issues or questions:
- Check troubleshooting section above
- Review log files in work directory
- Examine SLURM error logs in `RESULTS/slurm_logs/`
