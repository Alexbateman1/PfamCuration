shuf triage | grep -v "Iterate/Iterate/Iterate/Iterate/Iterate/Iterate" | awk '$5>0.99 && $4>1{system("/homes/agb/Scripts/iterate_inline.pl " $1)}'
