shuf triage | grep Iterate | awk '$5>0.8 && $4>9{system("add_foldseek.py "$1)}'
