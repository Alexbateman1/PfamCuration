current_dir=$(pwd);
echo $1

species_summary_new.pl $1

cd $1; add_abb_ref.py; add_pdb_ref.py; add_swiss_ref.pl; query_paperblast.py; ted_ali.pl SEED; cd $current_dir

add_images.py $1
add_foldseek.py $1
add_string.py $1 --min-frequency 0.5 --debug-limit 50
