shuf triage | awk '$5<0.2{system("species_summary_new.pl "$1";")}' 
