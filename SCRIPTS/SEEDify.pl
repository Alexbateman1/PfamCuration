#! /usr/bin/perl -w

use strict;

# Attempts to make a nice SEED from the ALIGN. Run it from within the family directory.

# This script will use the existing best scoring sequences in the scores file as the template sequence.
# It will after potentially subsampling the alignment extend it by some amount the use split_align to
# trim back to the coordinates of the template sequences. Then do a bunch of tidying steps.

# Ideally we would also like to make it so Swiss-Prot sequences are kept in preference to other sequences

# Here is a sample of the scores file:
# 353.7 D4GTU5.1/1-165 1-165 6.7e-103
# 172.1 A0A6G2G1W4.1/1-168 1-167 4.6e-47
# 171.2 A0A1I2RXG8.1/1-170 1-169 8.8e-47

##############################################
# Identify template sequence and coordinates #
##############################################

my $template_acc;
my $s;
my $e;
open (SCORES, "scores") or die "No scores file present";
while(<SCORES>){
    if (/\s+(\S+\.\d+)\/(\d+)-(\d+)\s/){
	$template_acc=$1;
	$s=$2;
	$e=$3;

	print STDERR "Template - $template_acc $s $e\n";
	last;
    }
}
close SCORES;

# Get template from ALIGN and remove it from ALIGN file
system("cp ALIGN ORIGALIGN");
system("grep '$template_acc' ALIGN > tmp.template");
system("grep -v '$template_acc' ALIGN > tmp.ALIGN");

# Check if sp file exists that contains any swissprot matches that we should keep
if (! -e "sp"){
    print STDERR "No sp file exists in this directory. Make it now\n";
    system("swissprot.pl");
}

unlink "tmp.sp"; # remove this file before starting as we are appending to it!
if (-z "sp"){
    print STDERR "No swissprot sequences in this family\n";



} elsif (-s "sp"){
    print STDERR "Swissprot sequences exist in this family\n";

    # get all ALIGN lines for swissprot entries
    open (SP, "sp") or die "Cannot open sp file";
    while(<SP>){
	if (/^(\S+)\./){
	    system("grep '$1' tmp.ALIGN >> tmp.sp");
	    system("grep -v '$1' tmp.ALIGN > tmp.ALIGN2");
	    system("mv tmp.ALIGN2 tmp.ALIGN");
	}
    }
    close SP;
}


# We must ensure that this accession is present in the SEED alignment!
if (-e "ALIGN"){
    # Find how many sequence in the ALIGN.
    my $size=num_seq("ALIGN");
    if ($size>400){

	print STDERR "Down sampling alignment of $size to 400 sequences.\n";
	system("shuf -n 400 tmp.ALIGN > tmp.rest");	
    } else {
	system("cp tmp.ALIGN tmp.rest");
    }

    # Now reassemble new SEED with template, swissprot seqs and all the rest.
    system("cat tmp.template tmp.sp tmp.rest > tmp");
    
    system("extend.pl -align tmp -n 50 -c 50 -m > SEED2");
    if (! -s "SEED2"){
	warn "In $0: Failed to build SEED2";
	exit;
    }
    system("split_align.pl -align SEED2 -acc $template_acc -s $s -e $e > SEED3");
    if (! -s "SEED3"){
	warn "In $0: Failed to build SEED3";
	exit;
    }
    # Remove partial sequences
    system("belvu -P -o mul SEED3 | grep -v '//' > SEED4");
    if (! -s "SEED4"){
	warn "In $0: Failed to build SEED4";
	exit;
    }

    # Make non-redundant at 80% identity
    my $nr_thresh = 80;
    system("belvu -n $nr_thresh SEED4 -o mul | grep -v '//' > NEWSEED");

    # Check if we ended up with only 1 sequence
    my $seed_size = num_seq("NEWSEED");

    if ($seed_size == 1) {
	print STDERR "Only 1 sequence in NEWSEED at ${nr_thresh}% identity, trying higher thresholds...\n";
	$nr_thresh = 85;  # Start at 85%

	while ($nr_thresh <= 100 && $seed_size == 1) {
	    system("belvu -n $nr_thresh SEED4 -o mul | grep -v '//' > NEWSEED");
	    $seed_size = num_seq("NEWSEED");

	    if ($seed_size > 1) {
		print STDERR "Found $seed_size sequences at ${nr_thresh}% identity\n";
		last;
	    }

	    $nr_thresh += 5;
	}

	if ($seed_size == 1) {
	    warn "In $0: Still only 1 sequence after trying all thresholds up to 100%. Cannot create multi-sequence SEED alignment.";
	    exit;
	}
    } else {
	print STDERR "Created NEWSEED with $seed_size sequences at ${nr_thresh}% identity\n";
    }
} else {
    die "Cannot find file ALIGN in current directory!";
}



# Get number of sequence in an alignment
#
sub num_seq {
    my $file = shift @_;
    my $nseq = 0; #Set to zero size
    
    unless (-e $file){warn("failed in num_seq $file not found [$!]");}
    open (TMP, "$file") or warn("failed to open $file [$!]");
    while (<TMP>) {
        $nseq++;
    }
    close(TMP);
    return $nseq;
}
