# Setting Up the Weekly HHsearch Scron Job

This guide explains how to set up the HHsearch pipeline to run automatically every Sunday at 4pm using SLURM's scron facility.

## Prerequisites

1. The initial SVN checkout must be complete (run `setup_svn_checkout.sh` first)
2. You must have access to scrontab on the SLURM cluster
3. The HHsearch scripts directory must be accessible from compute nodes

## Files

- `run_weekly_hhsearch.sh` - Main wrapper script that runs the pipeline
- `scrontab_entry.txt` - Template for the scrontab entry

## Setup Instructions

### Step 1: Update the script path

Edit `scrontab_entry.txt` and replace `/path/to/PfamCuration/HHsearch/` with the actual path:

```bash
# Example: if your scripts are in /homes/agb/Curation/PfamCuration/HHsearch/
0 16 * * 0 /homes/agb/Curation/PfamCuration/HHsearch/run_weekly_hhsearch.sh
```

### Step 2: Create the logs directory

```bash
mkdir -p /nfs/production/agb/pfam/data/all_vs_all/HHsearch/logs
```

### Step 3: View current scrontab (if any)

```bash
scrontab -l
```

### Step 4: Edit scrontab

```bash
scrontab -e
```

This opens an editor. Add the contents of `scrontab_entry.txt`:

```
# HHsearch Weekly Pipeline - runs every Sunday at 4pm (16:00)
#SCRON --partition=standard
#SCRON --time=24:00:00
#SCRON --mem=8G
#SCRON --job-name=hhsearch_weekly
#SCRON --output=/nfs/production/agb/pfam/data/all_vs_all/HHsearch/logs/scron_%j.out
#SCRON --error=/nfs/production/agb/pfam/data/all_vs_all/HHsearch/logs/scron_%j.err

0 16 * * 0 /homes/agb/Curation/PfamCuration/HHsearch/run_weekly_hhsearch.sh
```

Save and exit the editor.

### Step 5: Verify the scrontab was saved

```bash
scrontab -l
```

## Scrontab Format

```
# ┌───────────── minute (0 - 59)
# │ ┌───────────── hour (0 - 23)
# │ │ ┌───────────── day of month (1 - 31)
# │ │ │ ┌───────────── month (1 - 12)
# │ │ │ │ ┌───────────── day of week (0 - 6) (Sunday = 0)
# │ │ │ │ │
# * * * * * command
```

**Sunday at 4pm = `0 16 * * 0`**

## Monitoring

### Check if the scron job is scheduled

```bash
squeue -u $USER
```

### Check recent logs

```bash
ls -lrt /nfs/production/agb/pfam/data/all_vs_all/HHsearch/logs/
tail -f /nfs/production/agb/pfam/data/all_vs_all/HHsearch/logs/pipeline_*.log
```

### Check pipeline progress markers

```bash
ls -la /nfs/production/agb/pfam/data/all_vs_all/HHsearch/.progress/
```

## Manual Run

To run the pipeline manually (outside of scron):

```bash
cd /path/to/PfamCuration/HHsearch
./run_weekly_hhsearch.sh
```

Or run specific steps:

```bash
# Re-run just aggregation
python3 pfam_hhsearch_pipeline.py \
    --work-dir /nfs/production/agb/pfam/data/all_vs_all/HHsearch \
    --reset aggregation

# Run in incremental mode (only changed families)
python3 pfam_hhsearch_pipeline.py \
    --work-dir /nfs/production/agb/pfam/data/all_vs_all/HHsearch
```

## Troubleshooting

### Job not starting

1. Check scrontab is set correctly: `scrontab -l`
2. Check SLURM logs: `ls -la /nfs/production/agb/pfam/data/all_vs_all/HHsearch/logs/scron_*.err`

### Pipeline failing

1. Check the pipeline log: `tail /nfs/production/agb/pfam/data/all_vs_all/HHsearch/logs/pipeline_*.log`
2. Check SLURM job logs: `ls -la /nfs/production/agb/pfam/data/all_vs_all/HHsearch/RESULTS/slurm_logs/`
3. Check progress markers: `ls -la /nfs/production/agb/pfam/data/all_vs_all/HHsearch/.progress/`

### Reset and re-run

To reset from a specific step:

```bash
python3 pfam_hhsearch_pipeline.py \
    --work-dir /nfs/production/agb/pfam/data/all_vs_all/HHsearch \
    --reset <step> \
    --full

# Steps: svn_discovery, extraction, db_build, hhsearch, aggregation
```

## Removing the Scron Job

To remove the scheduled job:

```bash
scrontab -r
```

Or edit and remove the entry:

```bash
scrontab -e
# Delete the lines, save and exit
```
