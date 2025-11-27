current_dir=$(pwd); awk '$5>0.8{system("cd "$1"; echo "$1"; query_paperblast.py; cd '"$current_dir"'")}' triage
