awk '$2>9 {print $1}' triage | shuf | while read -r family; do
    overlap_file="${family}/overlap"
    if [[ -f "$overlap_file" ]] && [[ $(find "$overlap_file" -mtime -1 2>/dev/null) ]]; then
        echo "Skipping $family - overlap file exists and is less than 1 day old"
    else
        pqc-overlap-rdb -no_sigP "$family"
    fi
done
