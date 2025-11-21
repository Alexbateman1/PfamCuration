# What to Do After Download Completes

## Step 1: Verify Download Completed

```bash
# Check that all families were downloaded
ls /nfs/production/agb/pfam/users/agb/HHblits/DATA/SEED/ | wc -l
# Should show: 28761

ls /nfs/production/agb/pfam/users/agb/HHblits/DATA/HMM/ | wc -l
# Should show: 28761
```

## Step 2: Re-run Pipeline (Builds Database & Generates SLURM Script)

The pipeline will automatically continue from where it left off:

```bash
python3 ~/Curation/PfamCuration/HHblits/pfam_hhblits_pipeline.py \
  --work-dir /nfs/production/agb/pfam/users/agb/HHblits \
  --e-value 1e-2
```

This will:
1. Skip the download (already complete)
2. Build ffindex database from HMMs
3. Create ~288 batch files
4. Generate SLURM script
5. Pause and show you the sbatch command

## Step 3: Submit SLURM Jobs

The pipeline will tell you exactly what to run, something like:

```bash
sbatch /nfs/production/agb/pfam/users/agb/HHblits/RESULTS/run_hhblits_batch.sh
```

This submits ~288 parallel jobs, each processing 100 families.

## Step 4: Monitor SLURM Jobs

```bash
# Check job status
squeue -u $USER

# Check specific job array
squeue -j <job_id>

# Watch for completion
watch -n 30 'squeue -u $USER | wc -l'
```

## Step 5: After SLURM Jobs Complete

Re-run the pipeline one more time to aggregate results:

```bash
python3 ~/Curation/PfamCuration/HHblits/pfam_hhblits_pipeline.py \
  --work-dir /nfs/production/agb/pfam/users/agb/HHblits \
  --e-value 1e-2
```

This will:
1. Aggregate all HHblits results
2. Generate summary files
3. Create statistics

## Step 6: Results

Your results will be in:

```
/nfs/production/agb/pfam/users/agb/HHblits/RESULTS/SUMMARY/
├── hhblits_all_vs_all.tsv           # Complete comparison matrix
├── family_statistics.tsv             # Per-family summary
├── high_confidence_overlaps.tsv      # Strong hits (E < 1e-2)
└── pipeline_metadata.json            # Run metadata
```

## Troubleshooting

### If some families failed to download:

```bash
# The pipeline will skip them and log warnings
# Check the log:
grep "Failed to extract" /nfs/production/agb/pfam/users/agb/HHblits/pipeline.log

# Re-run to try again:
python3 ~/Curation/PfamCuration/HHblits/pfam_hhblits_pipeline.py \
  --work-dir /nfs/production/agb/pfam/users/agb/HHblits \
  --e-value 1e-2 \
  --reset extraction
```

### If SLURM jobs fail:

```bash
# Check error logs
ls /nfs/production/agb/pfam/users/agb/HHblits/RESULTS/slurm_logs/

# Re-submit failed jobs (pipeline tracks completion per family)
sbatch /nfs/production/agb/pfam/users/agb/HHblits/RESULTS/run_hhblits_batch.sh
```

## Timeline

- **Download**: ~13-14 hours (currently running)
- **Database build**: ~5-10 minutes
- **HHblits searches**: ~2-4 hours (with ~288 parallel jobs)
- **Aggregation**: ~10-20 minutes

**Total**: ~15-18 hours end-to-end
