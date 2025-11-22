#!/bin/bash
#
# Convenience wrapper for STRING network analysis
#
# Usage: analyze_string.sh <pfam_directory> [options]
#
# Examples:
#   analyze_string.sh PF06267
#   analyze_string.sh PF06267 --score-threshold 700
#   analyze_string.sh /path/to/PF06267 --output results.txt
#

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PYTHON_SCRIPT="${SCRIPT_DIR}/analyze_string_network.py"

# Check if Python script exists
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "Error: analyze_string_network.py not found at $PYTHON_SCRIPT"
    exit 1
fi

# Check for help flag
if [ "$1" = "-h" ] || [ "$1" = "--help" ] || [ $# -eq 0 ]; then
    echo "STRING Network Analysis for Pfam Families"
    echo ""
    echo "Usage: $0 <pfam_directory> [options]"
    echo ""
    echo "Options:"
    echo "  --score-threshold N    Minimum STRING score (default: 400)"
    echo "  --output FILE          Save report to file"
    echo "  --skip-species-summary Don't run species_summary_new.pl"
    echo "  -h, --help            Show this help"
    echo ""
    echo "Examples:"
    echo "  $0 PF06267"
    echo "  $0 PF06267 --score-threshold 700"
    echo "  $0 /full/path/to/PF06267 --output analysis.txt"
    echo ""
    exit 0
fi

# Run the Python script with all arguments
exec python3 "$PYTHON_SCRIPT" "$@"
