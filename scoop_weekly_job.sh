#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem=4G
#SBATCH --job-name=scoop_weekly
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=your.email@ebi.ac.uk
#SBATCH --output=/nfs/production/agb/pfam/data/all_vs_all/SCOOP/scoop_job_%j.out
#SBATCH --error=/nfs/production/agb/pfam/data/all_vs_all/SCOOP/scoop_job_%j.err

# Change to the SCOOP directory
cd /nfs/production/agb/pfam/data/all_vs_all/SCOOP

# Run the SCOOP pipeline script
bash run_scoop.sh
