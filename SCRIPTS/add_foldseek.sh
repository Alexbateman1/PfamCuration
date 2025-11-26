shuf triage | awk '$5<0.2{system("add_foldseek.py "$1)}'
