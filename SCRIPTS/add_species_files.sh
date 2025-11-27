shuf triage | awk '$5>0.8{system("species_summary_new.pl "$1";")}' 
