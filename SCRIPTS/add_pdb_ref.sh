current_dir=$(pwd); awk '$5<0.2{system("cd "$1"; echo "$1"; add_pdb_ref.py; cd '"$current_dir"'")}' triage
