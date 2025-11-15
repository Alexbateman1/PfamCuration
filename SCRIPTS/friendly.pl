#! /usr/bin/perl -w

# Convert scoop output data to a more friendly format
# with pfam ids and false and true matches marked up.
# Can add multiple file names to command line.

use strict;
# Get live clan data
open(FH, "mysql --defaults-file=~/.my.cnf pfam_live -e 'select clan.clan_acc,pfamA.pfamA_acc,pfamA_id,type,model_length from clan join clan_membership join pfamA where clan.clan_acc = clan_membership.clan_acc and clan_membership.pfamA_acc = pfamA.pfamA_acc;' |");
my %clanmap;
while(<FH>){
    if (/^(CL\d{4})\s+(PF\d{5})\s+(\S+)\s+(\S+)\s+(\S+)/){
        my $clan_acc=$1;
        my $pfam_acc=$2;
        my $pfam_id=$3;
        my $type=$4;
        my $model_length=$5;
        $clanmap{$pfam_acc}=$clan_acc;
    }
}
close FH;


# Get live mapping of Pfam ids and accs data
open(FH, "mysql --defaults-file=~/.my.cnf pfam_live -e 'select pfamA_acc,pfamA_id,model_length,type from pfamA;' |");
my %accmap;
my %pfamid2acc;
my %lengthmap;
my %typemap;
while(<FH>){
    if (/^(PF\d{5})\s+(\S+)\s+(\d+)\s+(\S+)/){
        my $pfam_acc=$1;
        my $pfam_id=$2;
	my $model_length=$3;
	my $type=$4;
        $pfamid2acc{$pfam_id}=$pfam_acc;
	$accmap{$pfam_acc}=$pfam_id;
	$lengthmap{$pfam_acc}=$model_length;
	$typemap{$pfam_acc}=$type;
    }
}
close FH;

# Get a list of which pfam families are nested
print STDERR "Getting list of nested domains\n";
my %nestmap;
open(FH, "mysql --defaults-file=~/.my.cnf pfam_live -e 'select pfamA_acc,nests_pfamA_acc from nested_domains;' |");
while(<FH>){
    # Consider relationship in either direction
    if (/^(\S+)\s+(\S+)/){
	$nestmap{$1}{$2}=1;
	$nestmap{$2}{$1}=1;
    }
}
close FH;



foreach my $file (@ARGV){
    if (-s "$file.friendly"){
	print STDERR "$file.friendly already exists. Skipping.\n";
	next;
    }

    open (FH, "$file") or die "cannot open $file";
    # Write output file

    open (OUTPUT, "> $file.friendly") or die "Cannot write to output file $file.friendly";
    while(<FH>){
	if (/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+\s+\S+\s+\d+)/){
	print STDERR;
	    my $score=$1;
	    my $acc1=$2;
	    my $x1=$3;
	    my $acc2=$4;
	    my $x2=$5;
	    
	    my $value="";
	    if ($clanmap{$acc1} and $clanmap{$acc2} and $clanmap{$acc1} eq $clanmap{$acc2}){
		$value="TRUE";
	    } elsif ($nestmap{$acc1}{$acc2} or $nestmap{$acc2}{$acc1}){
		$value='NESTED';
	    } elsif ($clanmap{$acc1} and $clanmap{$acc2} and $clanmap{$acc1} ne $clanmap{$acc2}){
		$value="FALSE";
	    } elsif ($clanmap{$acc1} ){
		$value="LINKED_1";
	    } elsif ($clanmap{$acc2}){
		$value="LINKED_2";
	    } 
	    
	    my $clan1;
	    if ($clanmap{$acc1}){
		$clan1=$clanmap{$acc1};
	    } else {
		$clan1="------";
	    }
	    my $clan2;
	    if ($clanmap{$acc2}){
		$clan2=$clanmap{$acc2};
	    } else {
		$clan2="------";
	    }
	    print OUTPUT "$score $acc1 ",$accmap{$acc1}," ",$lengthmap{$acc1}," ",$typemap{$acc1}," $clan1 $x1 $acc2 ",$accmap{$acc2}," ",$lengthmap{$acc2}," ",$typemap{$acc2}," $clan2 $x2 $value\n";
	} else {
	    warn "Unrecognised line! $_";
	}
    }
    close FH;
    close OUTPUT;
}





