#!/usr/bin/env perl
use strict;
use warnings;
use LWP::UserAgent;
use JSON;
use Getopt::Long;

# Default to first positional argument for alignment file
my $alignment_file = $ARGV[0];
my $output_file = "TED";
my $mean_file = "TED_NUM";
my $verbose = 0;
my $force = 0;

GetOptions(
    "verbose" => \$verbose,
    "force"   => \$force,  # Add force option to overwrite existing files
    "help"    => \&usage
) or usage();

usage() unless $alignment_file;

# Check if output files already exist
if ((-e $output_file || -e $mean_file) && !$force) {
    print STDERR "Error: File '$output_file' or '$mean_file' already exists in this directory.\n";
    print STDERR "This indicates the script has already been run in this directory.\n";
    print STDERR "Use --force to overwrite the existing files if needed.\n";
    exit 1;
}

# Function to display usage
sub usage {
    print "Usage: $0 <alignment_file> [options]\n";
    print "Generates a report of TED domains for proteins in a Pfam alignment.\n\n";
    print "Options:\n";
    print "  --verbose    Show detailed progress information\n";
    print "  --force      Overwrite existing files if they exist\n";
    print "  --help       Display this help message\n";
    exit 1;
}

# Function to strip version from UniProt accession
sub strip_version {
    my $acc = shift;
    # Remove version number if present (e.g., O34341.2 -> O34341)
    $acc =~ s/\.\d+$//;
    return $acc;
}

# Parse alignment file to extract sequence IDs
print STDERR "Parsing alignment file...\n" if $verbose;
my %sequences;
open(my $fh, '<', $alignment_file) or die "Cannot open $alignment_file: $!";
while (<$fh>) {
    chomp;
    if (/^(\S+)\/(\d+)-(\d+)/) {
        my $id = $1;
        my $uniprot_acc = strip_version($id);
        $sequences{$id} = {
            uniprot => $id,
            uniprot_no_version => $uniprot_acc,
            pfam_start => $2,
            pfam_end => $3,
            ted_domains => []
        };
    }
}
close($fh);

my $seq_count = scalar(keys %sequences);
print STDERR "Found $seq_count sequences in alignment.\n" if $verbose;

# Initialize HTTP client
my $ua = LWP::UserAgent->new;
$ua->timeout(30);
$ua->agent('Mozilla/5.0');

# Open output file for writing
open(my $out_fh, '>', $output_file) or die "Cannot open output file $output_file for writing: $!";

# Print header to the file
print $out_fh "Protein_Accession\tTED_Domains\tDomain_Count\tDomain_IDs\tDomain_Ranges\tOverlap_with_Pfam\n";

# Process each sequence
my $processed = 0;
my $with_ted = 0;
my $total_domains = 0;
my $errors = 0;

foreach my $seq_id (sort keys %sequences) {
    $processed++;
    my $query_acc = $sequences{$seq_id}{uniprot_no_version};
    
    # Only print detailed progress in verbose mode
    if ($verbose) {
        print STDERR "Fetching TED domains for $seq_id (using $query_acc for TED query)...\n";
        print STDERR "Progress: $processed/$seq_count sequences processed\n";
    }
    
    my $ted_api_url = "https://ted.cathdb.info/api/v1/uniprot/summary/$query_acc?skip=0&limit=100";
    my $ted_response = $ua->get($ted_api_url);
    
    if (!$ted_response->is_success) {
        warn "Warning: Failed to fetch domain information from TED for $query_acc.\n" if $verbose;
        $errors++;
        next;
    }
    
    my $domain_data;
    eval {
        $domain_data = decode_json($ted_response->content);
    };
    
    if ($@ || !$domain_data) {
        warn "Warning: Failed to parse domain information for $query_acc.\n" if $verbose;
        $errors++;
        next;
    }
    
    my $domain_count = $domain_data->{count} || 0;
    if ($domain_count == 0) {
        # Skip proteins with no TED domains
        next;
    }
    
    $with_ted++;
    $total_domains += $domain_count;
    
    # Process each domain separately
    foreach my $domain (@{$domain_data->{data}}) {
        my $ted_id = $domain->{ted_id};
        my $chopping = $domain->{chopping};
        
        # Calculate overlap with Pfam domain
        my $pfam_start = $sequences{$seq_id}{pfam_start};
        my $pfam_end = $sequences{$seq_id}{pfam_end};
        
        # Parse domain boundaries
        my @segments = split(/_/, $chopping);
        my ($ted_start, $ted_end);
        
        if (@segments == 1) {
            # Single segment domain
            ($ted_start, $ted_end) = ($chopping =~ /^(\d+)-(\d+)$/);
        } else {
            # Multi-segment domain
            my $first_segment = $segments[0];
            my $last_segment = $segments[-1];
            ($ted_start) = ($first_segment =~ /^(\d+)/);
            ($ted_end) = ($last_segment =~ /(\d+)$/);
        }
        
        $ted_start = int($ted_start);
        $ted_end = int($ted_end);
        
        # Calculate overlap percentage with Pfam domain
        my $overlap_percent = 0;
        if ($ted_start <= $pfam_end && $ted_end >= $pfam_start) {
            my $pfam_length = $pfam_end - $pfam_start + 1;
            my $overlap_start = $pfam_start > $ted_start ? $pfam_start : $ted_start;
            my $overlap_end = $pfam_end < $ted_end ? $pfam_end : $ted_end;
            my $overlap_length = $overlap_end - $overlap_start + 1;
            $overlap_percent = sprintf("%.1f%%", ($overlap_length / $pfam_length) * 100);
        }
        
        # Write one line per domain to the file
        print $out_fh "$seq_id\t";
        print $out_fh "$ted_id:$chopping\t";
        print $out_fh "1\t";  # Each line reports exactly one domain
        print $out_fh "$ted_id\t";
        print $out_fh "$ted_start-$ted_end\t";
        print $out_fh "$overlap_percent\n";
    }
}

# Make sure to close the file handle
close($out_fh);

# Calculate the mean domains per protein
my $mean_domains = $with_ted > 0 ? $total_domains / $with_ted : 0;

# Write the mean domains per protein to TED_NUM file
open(my $num_fh, '>', $mean_file) or die "Cannot open file $mean_file for writing: $!";
printf $num_fh "%.2f\n", $mean_domains;
close($num_fh);

# Print summary to STDERR
print STDERR "\n=== TED DOMAIN REPORT SUMMARY ===\n";
print STDERR "Total sequences analyzed: $seq_count\n";
print STDERR "Sequences with TED domains: $with_ted\n";
print STDERR "Total TED domains found: $total_domains\n";
print STDERR "Mean domains per protein: ", sprintf("%.2f", $mean_domains), "\n";
print STDERR "API or parsing errors: $errors\n";
print STDERR "Results written to file: $output_file\n";
print STDERR "Mean domains per protein written to: $mean_file\n";
