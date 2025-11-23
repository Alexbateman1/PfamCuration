#!/bin/bash
#
# Run STRING analysis on Pfam directories listed in a triage file
#
# Usage: run_string_on_triage.sh <triage_file> [options]
#
# The triage file should contain one Pfam directory per line (e.g., PF00001)
# This script will shuffle the list and run add_string.py on each directory.
#
# Options are passed through to add_string.py (e.g., --score-threshold, --debug-limit)
#
# Examples:
#   run_string_on_triage.sh triage
#   run_string_on_triage.sh triage --score-threshold 500
#   run_string_on_triage.sh my_families.txt --debug-limit 10
#

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ADD_STRING="${SCRIPT_DIR}/add_string.py"

# Check if add_string.py exists
if [ ! -f "$ADD_STRING" ]; then
    echo "Error: add_string.py not found at $ADD_STRING"
    exit 1
fi

# Check for help flag or no arguments
if [ "$1" = "-h" ] || [ "$1" = "--help" ] || [ $# -eq 0 ]; then
    echo "Run STRING Analysis on Triage File"
    echo ""
    echo "Usage: $0 <triage_file> [options]"
    echo ""
    echo "The triage file should contain one Pfam directory per line."
    echo "Directories are processed in random order (shuffled)."
    echo ""
    echo "Options (passed to add_string.py):"
    echo "  --score-threshold N    Minimum STRING score (default: 400)"
    echo "  --min-frequency F      Minimum frequency threshold (default: 0.3)"
    echo "  --debug-limit N        Only process first N networks per family"
    echo ""
    echo "Examples:"
    echo "  $0 triage"
    echo "  $0 triage --score-threshold 500"
    echo "  $0 my_families.txt --debug-limit 10"
    echo ""
    exit 0
fi

# First argument is the triage file
TRIAGE_FILE="$1"
shift

# Check if triage file exists
if [ ! -f "$TRIAGE_FILE" ]; then
    echo "Error: Triage file not found: $TRIAGE_FILE"
    exit 1
fi

# Remaining arguments are options for add_string.py
OPTIONS="$@"

# Count total families
TOTAL=$(wc -l < "$TRIAGE_FILE")
echo "Processing $TOTAL families from $TRIAGE_FILE"
echo "Options: $OPTIONS"
echo ""

# Shuffle and process each directory
COUNTER=0
shuf "$TRIAGE_FILE" | while read -r PFAM_DIR; do
    # Skip empty lines
    if [ -z "$PFAM_DIR" ]; then
        continue
    fi

    COUNTER=$((COUNTER + 1))
    echo "================================================================================"
    echo "[$COUNTER/$TOTAL] Processing: $PFAM_DIR"
    echo "================================================================================"

    # Run add_string.py with the directory and any additional options
    "$ADD_STRING" "$PFAM_DIR" $OPTIONS

    EXIT_CODE=$?
    if [ $EXIT_CODE -ne 0 ]; then
        echo "WARNING: add_string.py failed for $PFAM_DIR with exit code $EXIT_CODE"
    fi

    echo ""
done

echo "Finished processing all families from $TRIAGE_FILE"
