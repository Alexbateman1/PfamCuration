#!/bin/bash
#
# Initial setup script for HHsearch pipeline SVN checkout.
#
# This script:
# 1. Cleans up any existing DATA directory (optional)
# 2. Performs a full SVN checkout of Pfam Families
# 3. Syncs SEED files to flat directory for pipeline use
#
# Usage:
#   ./setup_svn_checkout.sh [DATA_DIR]
#
# Example:
#   # Run in foreground (shows progress):
#   ./setup_svn_checkout.sh /nfs/production/agb/pfam/data/all_vs_all/HHsearch/DATA
#
#   # Run in background with nohup:
#   nohup ./setup_svn_checkout.sh /nfs/production/agb/pfam/data/all_vs_all/HHsearch/DATA > checkout.log 2>&1 &
#
#   # Or use screen:
#   screen -S svn_checkout
#   ./setup_svn_checkout.sh /nfs/production/agb/pfam/data/all_vs_all/HHsearch/DATA
#
# Note: Full checkout takes approximately 8 hours for ~20k families.
#

set -e  # Exit on error

# Configuration
SVN_URL="https://xfam-svn-hl.ebi.ac.uk/svn/pfam/trunk/Data/Families"
DATA_DIR="${1:-/nfs/production/agb/pfam/data/all_vs_all/HHsearch/DATA}"

FAMILIES_DIR="$DATA_DIR/Families"
SEED_DIR="$DATA_DIR/SEED"

echo "============================================================"
echo "HHsearch Pipeline - SVN Checkout Setup"
echo "============================================================"
echo ""
echo "Configuration:"
echo "  SVN URL:      $SVN_URL"
echo "  DATA_DIR:     $DATA_DIR"
echo "  Families dir: $FAMILIES_DIR"
echo "  SEED dir:     $SEED_DIR"
echo ""
echo "Start time: $(date)"
echo ""

# Check if DATA_DIR exists and has content
if [ -d "$DATA_DIR" ] && [ "$(ls -A $DATA_DIR 2>/dev/null)" ]; then
    echo "WARNING: DATA directory exists and is not empty: $DATA_DIR"
    echo ""
    ls -la "$DATA_DIR"
    echo ""
    read -p "Do you want to clean it and start fresh? (y/N): " confirm
    if [ "$confirm" = "y" ] || [ "$confirm" = "Y" ]; then
        echo ""
        echo "Cleaning up existing DATA directory..."
        rm -rf "$FAMILIES_DIR"
        rm -rf "$SEED_DIR"
        echo "Cleanup complete."
    else
        echo ""
        echo "Keeping existing data. Will update if SVN checkout exists."
    fi
fi

# Create directories
echo ""
echo "Creating directories..."
mkdir -p "$DATA_DIR"
mkdir -p "$SEED_DIR"

# Check if SVN checkout already exists
if [ -d "$FAMILIES_DIR/.svn" ]; then
    echo ""
    echo "SVN working copy already exists. Running svn update..."
    echo ""
    svn update "$FAMILIES_DIR"
else
    echo ""
    echo "============================================================"
    echo "Starting SVN checkout (this will take several hours)..."
    echo "============================================================"
    echo ""
    svn checkout "$SVN_URL" "$FAMILIES_DIR"
fi

echo ""
echo "SVN checkout/update complete!"
echo ""

# Get SVN revision
REVISION=$(svn info "$FAMILIES_DIR" --show-item revision)
echo "Current SVN revision: $REVISION"
echo ""

# Sync SEED files to flat directory
echo "============================================================"
echo "Syncing SEED files to flat directory..."
echo "============================================================"
echo ""

TOTAL=0
SYNCED=0
MISSING=0

for family_dir in "$FAMILIES_DIR"/PF*; do
    if [ -d "$family_dir" ]; then
        family_id=$(basename "$family_dir")
        src="$family_dir/SEED"
        dst="$SEED_DIR/${family_id}_SEED"

        TOTAL=$((TOTAL + 1))

        if [ -f "$src" ]; then
            cp "$src" "$dst"
            SYNCED=$((SYNCED + 1))
        else
            MISSING=$((MISSING + 1))
        fi

        # Progress every 1000 families
        if [ $((TOTAL % 1000)) -eq 0 ]; then
            echo "  Progress: $TOTAL families processed..."
        fi
    fi
done

echo ""
echo "Sync complete!"
echo "  Total families: $TOTAL"
echo "  SEED files synced: $SYNCED"
echo "  Missing SEED files: $MISSING"
echo ""

# Count files in SEED directory
SEED_COUNT=$(ls -1 "$SEED_DIR" | wc -l)
echo "Files in SEED directory: $SEED_COUNT"
echo ""

# Save revision to file
echo "$REVISION" > "$DATA_DIR/.svn_revision"
echo "SVN revision saved to: $DATA_DIR/.svn_revision"

echo ""
echo "============================================================"
echo "Setup complete!"
echo "============================================================"
echo ""
echo "End time: $(date)"
echo ""
echo "Next steps:"
echo "  1. The pipeline can now use SEED files from: $SEED_DIR"
echo "  2. For weekly updates, run: svn update $FAMILIES_DIR"
echo "  3. Then sync SEED files using the pipeline or this script"
echo ""
