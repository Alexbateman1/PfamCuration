current_dir=$(pwd); awk '{system("cd "$1"; echo "$1"; query_paperblast.py; cd '"$current_dir"'")}' triage
