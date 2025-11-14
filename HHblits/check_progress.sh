#!/bin/bash
# Quick progress checker for Pfam download

WORK_DIR="${1:-./HHBLITS_PIPELINE}"

echo "=== Pfam Download Progress ==="
echo ""

# Count downloaded files
SEED_COUNT=$(ls "$WORK_DIR/DATA/SEED/" 2>/dev/null | wc -l)
HMM_COUNT=$(ls "$WORK_DIR/DATA/HMM/" 2>/dev/null | wc -l)
TOTAL=28761

echo "SEED files: $SEED_COUNT / $TOTAL"
echo "HMM files:  $HMM_COUNT / $TOTAL"

# Calculate percentage
if [ $SEED_COUNT -gt 0 ]; then
    PERCENT=$((SEED_COUNT * 100 / TOTAL))
    echo "Progress:   $PERCENT%"

    # Estimate time remaining (assumes ~35 families/minute)
    REMAINING=$((TOTAL - SEED_COUNT))
    MINUTES=$((REMAINING / 35))
    HOURS=$((MINUTES / 60))
    MINS=$((MINUTES % 60))
    echo "Estimated:  ~${HOURS}h ${MINS}m remaining"
fi

echo ""
echo "Latest files:"
ls -lt "$WORK_DIR/DATA/SEED/" 2>/dev/null | head -5

echo ""
echo "To check detailed logs:"
echo "  tail -f $WORK_DIR/pipeline.log"
echo ""
echo "To reattach to screen:"
echo "  screen -r pfam_download"
