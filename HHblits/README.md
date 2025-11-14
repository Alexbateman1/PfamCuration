# Pfam HHblits All-Against-All Pipeline

Simple pipeline for extracting Pfam SEED alignments and HMMs from SVN and performing all-against-all HHblits comparisons.

## Requirements

- Python 3.6+
- `svn` command-line client
- HHsuite (`hhblits`, `hhmake`, etc.)

## Quick Start

### First Run (Full)

```bash
python pfam_hhblits_pipeline.py --work-dir /path/to/work/dir
```

This will:
1. List all families in Pfam SVN
2. Extract SEED and HMM files for all families
3. Build HHblits database
4. Run all-against-all comparisons
5. Generate summary files

### Weekly Incremental Update

```bash
python pfam_hhblits_pipeline.py --work-dir /path/to/work/dir
```

The pipeline automatically runs in incremental mode. It will:
1. Check SVN for changes since last run
2. Extract only changed families
3. Re-run HHblits for changed families
4. Update results

### Force Full Run

```bash
python pfam_hhblits_pipeline.py --work-dir /path/to/work/dir --full
```

## Command-Line Options

```
--work-dir, -w      Working directory (default: ./HHBLITS_PIPELINE)
--svn-url           Pfam SVN URL (default: https://xfam-svn-hl.ebi.ac.uk/svn/pfam/)
--e-value, -e       E-value threshold (default: 1e-10)
--full              Force full run, disable incremental mode
--no-slurm          Run sequentially without SLURM (for testing)
--batch-size        Families per SLURM batch (default: 100)
--reset STEP        Reset progress from specified step
--log-level         Logging level: DEBUG, INFO, WARNING, ERROR
```

## Output Files

```
WORK_DIR/
├── .last_revision                          # SVN revision tracking
├── families_to_process.txt                 # List of families processed
├── extraction_manifest.tsv                 # Extraction status
├── pipeline.log                            # Pipeline log
├── DATA/
│   ├── SEED/                              # SEED alignments
│   └── HMM/                               # HMM profiles
└── RESULTS/
    ├── A3M/                               # Converted alignments
    ├── HHR/                               # Raw HHblits output
    ├── PARSED/                            # Parsed hit tables
    └── SUMMARY/
        ├── hhblits_all_vs_all.tsv        # Complete comparison matrix
        ├── family_statistics.tsv          # Per-family summary
        ├── high_confidence_overlaps.tsv   # Strong hits (E < 1e-20)
        └── pipeline_metadata.json         # Run metadata
```

## SLURM Batch Processing

The pipeline uses SLURM for efficient parallel processing:

### How it works:
1. **Families are batched** into groups of 100 (configurable)
2. **SLURM array job** submits ~288 jobs (for 28,761 families)
3. Each job processes 100 families sequentially
4. Much faster than 28,761 individual jobs

### Workflow:

```bash
# Step 1: Run pipeline (generates SLURM script)
python pfam_hhblits_pipeline.py -w /path/to/work

# Pipeline will pause after generating script, showing:
# "To submit jobs, run: sbatch /path/to/run_hhblits_batch.sh"

# Step 2: Submit SLURM jobs
sbatch /path/to/work/RESULTS/run_hhblits_batch.sh

# Step 3: Monitor jobs
squeue -u $USER

# Step 4: After jobs complete, re-run pipeline to aggregate
python pfam_hhblits_pipeline.py -w /path/to/work
```

### Batch Size Tuning:

```bash
# Larger batches = fewer jobs, longer per job
python pfam_hhblits_pipeline.py -w /path --batch-size 200  # ~144 jobs

# Smaller batches = more jobs, shorter per job
python pfam_hhblits_pipeline.py -w /path --batch-size 50   # ~575 jobs
```

### Sequential Mode (No SLURM):

For testing or small runs:

```bash
python pfam_hhblits_pipeline.py -w /path --no-slurm
```

## Progress Tracking

The pipeline uses hidden marker files to track completion:
- `.svn_discovery_complete`
- `.extraction_complete`
- `.db_build_complete`
- `.hhblits_complete`
- `.aggregation_complete`

If the pipeline is interrupted, it will resume from the last completed step.

### Reset Progress

To re-run from a specific step:

```bash
python pfam_hhblits_pipeline.py --reset extraction
```

## Incremental Updates

The pipeline tracks the SVN revision number in `.last_revision`. On subsequent runs:

1. Queries SVN for changes since last revision
2. Identifies modified families
3. Re-extracts only changed families
4. Re-runs HHblits for changed families
5. Updates aggregated results

## Example Usage

```bash
# First run - process all families
python pfam_hhblits_pipeline.py -w /data/pfam_hhblits

# Weekly update - only process changes
python pfam_hhblits_pipeline.py -w /data/pfam_hhblits

# Force full re-run with different E-value
python pfam_hhblits_pipeline.py -w /data/pfam_hhblits --full -e 1e-5

# Re-run just HHblits searches
python pfam_hhblits_pipeline.py -w /data/pfam_hhblits --reset hhblits
```

## Output Format

### hhblits_all_vs_all.tsv

Complete comparison matrix with columns:
- `query_family`: Query Pfam family
- `target_family`: Target Pfam family
- `probability`: HHblits probability score
- `e_value`: E-value
- `p_value`: P-value
- `score`: HHblits score
- `ss`: Secondary structure score
- `cols`: Aligned columns
- `query_range`: Query alignment range
- `template_range`: Target alignment range
- `template_length`: Target length

### family_statistics.tsv

Per-family summary:
- `family`: Pfam family ID
- `num_hits`: Number of significant hits
- `best_e_value`: Best E-value among hits
- `best_target`: Family with best hit
- `best_score`: Best HHblits score

### high_confidence_overlaps.tsv

High-confidence hits (E-value < 1e-20) suggesting potential overlaps or merges.

## Database Format

The pipeline builds an **ffindex database** for efficient HHM lookups:
- Uses `ffindex_build` to create indexed database
- Fast random access to 28,761+ HMM profiles
- Falls back to simple concatenation if ffindex is unavailable

The ffindex format creates two files:
- `pfam_db_hhm.ffdata` - Concatenated HMM data
- `pfam_db_hhm.ffindex` - Index mapping family IDs to byte offsets

## Notes

- The pipeline is designed to be simple and restartable
- All intermediate files are preserved for inspection
- SLURM batching provides efficient parallelization
- ffindex database enables fast lookups for incremental updates
