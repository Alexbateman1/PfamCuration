#!/bin/bash
#
# Weekly HHsearch Pipeline Runner
#
# This script is designed to be run via scron (SLURM cron) on a weekly basis.
# It performs the following:
# 1. Updates the SVN checkout with latest Pfam family data
# 2. Syncs SEED files to flat directory
# 3. Rebuilds the HH-suite database (for new/modified families)
# 4. Runs HHsearch all-against-all comparisons via SLURM
# 5. Aggregates results into summary files
#
# Usage:
#   ./run_weekly_hhsearch.sh
#
# Note: This script should be run from the HHsearch scripts directory
#       or the PYTHONPATH should include the scripts directory.
#

set -e  # Exit on error

# Configuration
WORK_DIR="/nfs/production/agb/pfam/data/all_vs_all/HHsearch"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${WORK_DIR}/logs"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="${LOG_DIR}/pipeline_${TIMESTAMP}.log"

# Create log directory if needed
mkdir -p "${LOG_DIR}"

echo "============================================================" | tee -a "${LOG_FILE}"
echo "HHsearch Weekly Pipeline Run" | tee -a "${LOG_FILE}"
echo "Started: $(date)" | tee -a "${LOG_FILE}"
echo "============================================================" | tee -a "${LOG_FILE}"
echo "" | tee -a "${LOG_FILE}"
echo "Configuration:" | tee -a "${LOG_FILE}"
echo "  WORK_DIR:   ${WORK_DIR}" | tee -a "${LOG_FILE}"
echo "  SCRIPT_DIR: ${SCRIPT_DIR}" | tee -a "${LOG_FILE}"
echo "  LOG_FILE:   ${LOG_FILE}" | tee -a "${LOG_FILE}"
echo "" | tee -a "${LOG_FILE}"

# Change to script directory (needed for imports)
cd "${SCRIPT_DIR}"

# Run the pipeline
# --full flag forces full reprocessing (recommended for weekly runs to catch all changes)
# Remove --full if you want incremental mode (only processes changed families)
python3 pfam_hhsearch_pipeline.py \
    --work-dir "${WORK_DIR}" \
    --full \
    --batch-size 100 \
    --poll-interval 60 \
    2>&1 | tee -a "${LOG_FILE}"

EXIT_CODE=${PIPESTATUS[0]}

echo "" | tee -a "${LOG_FILE}"
echo "============================================================" | tee -a "${LOG_FILE}"
echo "Pipeline finished: $(date)" | tee -a "${LOG_FILE}"
echo "Exit code: ${EXIT_CODE}" | tee -a "${LOG_FILE}"
echo "============================================================" | tee -a "${LOG_FILE}"

# Clean up old logs (keep last 12 weeks)
find "${LOG_DIR}" -name "pipeline_*.log" -mtime +84 -delete 2>/dev/null || true

exit ${EXIT_CODE}
