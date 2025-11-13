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

## Notes

- The pipeline is designed to be simple and restartable
- All intermediate files are preserved for inspection
- HHblits searches can be parallelized by running on SLURM (future enhancement)
- Database building concatenates all HMMs into a single file
