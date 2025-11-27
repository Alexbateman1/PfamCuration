current_dir=$(pwd); shuf triage | awk '$5>0.8{system("cd "$1"; echo "$1"; ted_ali.pl SEED; cd '"$current_dir"'")}'
