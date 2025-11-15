#! /usr/bin/perl -w

# Software to find matches between Pfam families.
#
# Basically the method compares two HMM output files and finds the
# number of common regions between the two. Two regions are counted as
# common if one of their midpoints lies within the other.  
#
# To use the software an input file must be provided. An example is
# shown below:
#
# O00050 PF06465  448     500     2.5e-31
# O00050 PF07453  19      52      1.8e+03
# O00050 PF07572  89      159     9.8e+02
# O00053 PF00012  37      642     0.0
# O00053 PF00416  507     576     3.9e+02
# O00053 PF00630  167     264     2.6e+02
# O00053 PF00929  37      180     8.8e+02
#
#
# The columns are as follows:
#  column 1 protein identifier
#  column 2 pfam family accession
#  column 3 start of match
#  column 4 end of match
#  cloumn 5 E-value of match
#
# NB The file MUST be sorted by the first column!
#
# An example of the output is shown here:
#
# 15.1651888357382 PF00001 23432 PF00002 8904 Both 3639
# 5.48080433630521 PF00001 23432 PF00003 23851 Both 3508
# 0.00261500319008754 PF00001 23432 PF00004 20335 Both 1
# 0.00396570101830259 PF00001 23432 PF00005 26835 Both 2
# 2.00941226589356 PF00001 23432 PF00007 1008 Both 40
# 0.00354743435012688 PF00001 23432 PF00009 14976 Both 1
# 0.0153235629090441 PF00001 23432 PF00010 3426 Both 1
# 0.127958318941515 PF00001 23432 PF00014 2030 Both 5
#
#
# The columns of the output are as follows:
#
# Column 1: SCOOP normalised score. Scores greater than 50 are very likely to be true matches.
# Column 2: Family 1 identifier (Pfam accession in this case).
# Column 3: Number of matches in total in family 1.
# Column 4: Family 2 identifier (Pfam accession in this case).
# Column 5: Number of matches in total in family 2.
# Column 6: Number of common matches between family 1 and family 2.

use strict;
use Region;
use RegionSet;
use Getopt::Long;

my $e_thresh=1000000;
&GetOptions('e=s'   => \$e_thresh); # Set to true to filter matches for redundancy using score

print STDERR "Only looking at regions with E-value less or equal than $e_thresh\n";


my $file = shift @ARGV;
if (! $file){
    die "Usage $0: <data file>";
}

print STDERR "Loading up regions\n";
my $total_region=0;
my $total_residue=0;
my %family_total; # Counts matches per family
my %uber_matrix; # Counts number of matches between regions.
my %sump_score_matrix; # Stores sum of probabilities of true matches

open (FH, $file) or die "Cannot open $file";
my ($old_id,$set,$new_set);
while(<FH>){
    my $line=$_;
    chomp $line;
    my @line=split(/\s+/,$line);
    my $id=$line[0];
    my $evalue=$line[4];

    if ($evalue>$e_thresh){
	next;
    }

    # Make new region object
    my $region = new Region('id'=>$id,
			    'start'=>$line[2],
			    'end'=>$line[3],
			    'family'=>$line[1],
			    'evalue'=>$evalue);

    $total_residue+=$line[3]-$line[2]+1;

    if ($id ne $old_id and $old_id){
	# Make new set
	$new_set=new RegionSet($id);
	$new_set->add($region);
	$total_region++; # This is counting number of regions! NB this is what original paper states.

	# Now add old set to uber_matrix
	add_set_to_matrix($set,\%family_total,\%uber_matrix,\%sump_score_matrix);
	# delete set
	undef $set;

	$set=$new_set;
    } elsif ($id eq $old_id) {
	$set->add($region);
	$total_region++;
    } elsif ($id ne $old_id and ! $old_id){ # Only reach here for first region in file
	$new_set=new RegionSet($id);
	$new_set->add($region);
	$total_region++;

	# No old set as this is the first protein
	# No old set to delete

	$set=$new_set;
    }
    $old_id=$id;
}
print STDERR "Finished processing regions\n";


# OK now have all regions loaded up into matrices.


# This routine takes a RegionSet and does comparison and adds data to  various matrices
sub add_set_to_matrix {
    my ($set,$family_total,$uber_matrix,$sump_score_matrix)=@_;

    my @regions=$set->each();
    my $n=@regions;

    # Deal with case of only a single region matching protein
    if ($n == 1){
	$family_total->{$regions[0]->family()}++; # Still add to family counts	$family_residue_total->{$regions[0]->family()}+=$regions[0]->{'end'}-$regions[0]->{'start'}+1; # Still add to family residue counts
	return;
    }
    
    # Now loop over diagonal half matrix
    for (my $i=0;$i<$n;$i++){
	# For every region add to family count
	if (! $regions[$i]){
	    warn "No region for $i\n"
	}

	$family_total->{$regions[$i]->family()}++;

	for (my $j=$i+1;$j<$n;$j++){
	    if (! $regions[$j]){
		warn "No region for $j\n";
	    }
	    # Are regions the same?
	    if ($regions[$i]->overlap($regions[$j])){
		# Test whether either region has counted towards a match between these families already.
		# To avoid any region being counted multiple times
		if (! $regions[$i]->match($regions[$j]->family()) and ! $regions[$j]->match($regions[$i]->family())){
			$uber_matrix->{$regions[$i]->family()}{$regions[$j]->family()}++;
			$sump_score_matrix->{$regions[$i]->family()}{$regions[$j]->family()}+=exp(-$regions[$i]->{'evalue'})*exp(-$regions[$j]->{'evalue'});
		    
		    # Add to list of outputs matched in region object
		    $regions[$i]->add_match($regions[$j]->family());
		    $regions[$j]->add_match($regions[$i]->family());
		}
	    }
	}
    }
    return();
}

# Now do final processing
print STDERR "Total regions:$total_region\n";


# Set up output file naming
my $file_suffix="E$e_thresh.output";

# Loop over every pair of families in sum of probability matrix
open (FH, "> sump.$file_suffix")or die "cannot write sump.$file_suffix";
foreach my $x (sort keys %sump_score_matrix){
    foreach my $y (sort keys %{$sump_score_matrix{$x}}){
	my $score=$sump_score_matrix{$x}{$y};
	my $observed=$uber_matrix{$x}{$y};
	# print out raw scores
	print FH "$score $x $family_total{$x} $y $family_total{$y} Both $observed\n";
    } 
}
close FH;

# Reorder file by score
system ("sort -S 1G -grk1 sump.$file_suffix > $$");
system ("mv $$ sump.$file_suffix");
