#!/bin/bash
# Fetch a single domain sequence and output with proper header
# Usage: fetch_domain_seq.sh <acc> <start> <end> <suffix>

ACC=$1
START=$2
END=$3
SUFFIX=$4

# Fetch sequence and replace header
pfetch -s "$START" -e "$END" "$ACC" 2>/dev/null | awk -v acc="$ACC" -v suf="$SUFFIX" -v s="$START" -v e="$END" '
NR==1 { print ">" acc "_" suf "/" s "-" e }
NR>1 { print }
'
