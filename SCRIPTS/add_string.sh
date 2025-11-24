shuf triage | awk '$5<0.2{system("add_string.py "$1)}'
