current_dir=$(pwd); awk '$5<0.2{system("cd "$1"; echo "$1"; query_paperblast.py; cd '"$current_dir"'")}' triage
