#!/bin/bash

triage_file="$1"
batch_size=25
max_jobs=40

if [[ -z "$triage_file" || ! -f "$triage_file" ]]; then
    echo "Usage: $0 <triage_file>"
    exit 1
fi

# Extract qualifying IDs
ids=($(awk '$5 > 0.8 && $4 > 9 {print $1}' "$triage_file"))

if [[ ${#ids[@]} -eq 0 ]]; then
    echo "No lines match criteria (\$5 > 0.8 and \$4 > 9)"
    exit 0
fi

echo "Found ${#ids[@]} qualifying entries"

# Process in batches
batch_num=0
for ((i=0; i<${#ids[@]}; i+=batch_size)); do
    batch=("${ids[@]:i:batch_size}")
    batch_num=$((batch_num + 1))
    
    # Wait if we have too many jobs running
    while [[ $(squeue -u "$USER" -h -n "annotate_batch" | wc -l) -ge $max_jobs ]]; do
        echo "Waiting for job slots (currently at max $max_jobs)..."
        sleep 30
    done
    
    # Create batch script
    batch_script=$(mktemp)
    cat > "$batch_script" << EOF
#!/bin/bash
#SBATCH --job-name=annotate_batch
#SBATCH --output=annotate_batch_${batch_num}_%j.out
#SBATCH --error=annotate_batch_${batch_num}_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=2G

for id in ${batch[@]}; do
    add_annotation.sh "\$id"
done
EOF
    
    sbatch "$batch_script"
    echo "Submitted batch $batch_num (${#batch[@]} IDs)"
    rm "$batch_script"
done

echo "All $batch_num batches submitted"
