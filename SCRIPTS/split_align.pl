#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

# Parse command line arguments
my $align_file;
my $accession;
my $start;
my $end;

GetOptions(
    "align=s" => \$align_file,
    "acc=s"   => \$accession,
    "s=i"     => \$start,
    "e=i"     => \$end
) or die "Usage: $0 -align <alignment_file> -acc <accession> -s <start> -e <end>\n";

# Check required parameters
die "Missing required parameter: -align <alignment_file>\n" unless $align_file;
die "Missing required parameter: -s <start>\n" unless defined $start;
die "Missing required parameter: -e <end>\n" unless defined $end;
die "Missing required parameter: -acc <accession>\n" unless defined $accession;
die "Start coordinate must be positive\n" if $start < 1;
die "End coordinate must be greater than or equal to start\n" if $end < $start;

# Open the alignment file
open my $fh, '<', $align_file or die "Could not open $align_file: $!\n";

# Read the alignment and build position mappings
my %seqs;
my %orig_ranges;
my %seq_to_aln_map;  # Maps sequence positions to alignment positions
my %aln_to_seq_map;  # Maps alignment positions to sequence positions
my @seq_order;       # Preserve original order of sequences

while (my $line = <$fh>) {
    chomp $line;
    next if $line =~ /^\s*$/ or $line =~ /^#/;
    
    if ($line =~ /^(\S+)\/(\d+)-(\d+)\s+(.+)$/) {
        my ($base_id, $orig_start, $orig_end, $seq) = ($1, $2, $3, $4);
        my $full_id = "$base_id/$orig_start-$orig_end";  # Keep the full identifier
        $seqs{$full_id} = $seq;
        $orig_ranges{$full_id}{start} = $orig_start;
        $orig_ranges{$full_id}{end} = $orig_end;
        $orig_ranges{$full_id}{base_id} = $base_id;  # Store base ID for matching
        
        # Store the order of sequences as they appear in the file
        push @seq_order, $full_id;
        
        # Create position mappings
        my $seq_pos = $orig_start;
        my @this_seq_to_aln;
        my @this_aln_to_seq;
        
        for (my $aln_pos = 0; $aln_pos < length($seq); $aln_pos++) {
            my $char = substr($seq, $aln_pos, 1);
            if ($char =~ /[A-Za-z]/) {  # It's a residue, not a gap
                $this_aln_to_seq[$aln_pos] = $seq_pos;
                $this_seq_to_aln[$seq_pos - $orig_start] = $aln_pos;
                $seq_pos++;
            } else {
                $this_aln_to_seq[$aln_pos] = undef;  # Gap
            }
        }
        
        $seq_to_aln_map{$full_id} = \@this_seq_to_aln;
        $aln_to_seq_map{$full_id} = \@this_aln_to_seq;
    }
}
close $fh;

# Find the target sequence - use partial matching if needed
my $target_id;
if (exists $seqs{$accession}) {
    $target_id = $accession;
} else {
    # Try to match by base accession ID
    my @matches = grep { $orig_ranges{$_}{base_id} =~ /^$accession/ } keys %seqs;
    if (@matches == 1) {
        $target_id = $matches[0];
    } elsif (@matches > 1) {
        die "Error: Multiple sequences match accession '$accession'\n";
    } else {
        die "Error: No sequences match accession '$accession'\n";
    }
}

# Get alignment positions corresponding to the sequence positions for the target sequence
my $target_orig_start = $orig_ranges{$target_id}{start};
my $target_orig_end = $orig_ranges{$target_id}{end};

# Validate that the requested coordinates are within the sequence range
if ($start < $target_orig_start) {
    die "Error: Start position $start is before the beginning of sequence $target_id (which starts at $target_orig_start)\n";
}
if ($end > $target_orig_end) {
    die "Error: End position $end is beyond the end of sequence $target_id (which ends at $target_orig_end)\n";
}
if ($start > $target_orig_end) {
    die "Error: Start position $start is beyond the end of sequence $target_id (which ends at $target_orig_end)\n";
}
if ($end < $target_orig_start) {
    die "Error: End position $end is before the beginning of sequence $target_id (which starts at $target_orig_start)\n";
}

my $target_start_pos = $start - $target_orig_start;
my $target_end_pos = $end - $target_orig_start;

my $aln_start_pos = $seq_to_aln_map{$target_id}[$target_start_pos];
my $aln_end_pos = $seq_to_aln_map{$target_id}[$target_end_pos];

die "Error: Start position $start is outside the range of sequence $target_id\n" 
    unless defined $aln_start_pos;
die "Error: End position $end is outside the range of sequence $target_id\n" 
    unless defined $aln_end_pos;

# Extract alignment region and calculate new coordinates for each sequence
my %results;
foreach my $full_id (keys %seqs) {
    my $seq = $seqs{$full_id};
    my $orig_start = $orig_ranges{$full_id}{start};
    my $base_id = $orig_ranges{$full_id}{base_id};
    
    # Extract the region from the alignment
    my $region = substr($seq, $aln_start_pos, $aln_end_pos - $aln_start_pos + 1);
    
    # Calculate new sequence range
    my $new_start = undef;
    my $new_end = undef;
    
    # Find first non-gap position in the extracted region
    for (my $i = $aln_start_pos; $i <= $aln_end_pos; $i++) {
        if (defined $aln_to_seq_map{$full_id}[$i]) {
            $new_start = $aln_to_seq_map{$full_id}[$i];
            last;
        }
    }
    
    # Find last non-gap position in the extracted region
    for (my $i = $aln_end_pos; $i >= $aln_start_pos; $i--) {
        if (defined $aln_to_seq_map{$full_id}[$i]) {
            $new_end = $aln_to_seq_map{$full_id}[$i];
            last;
        }
    }
    
    # Handle case of all gaps
    if (!defined $new_start || !defined $new_end) {
        # This could be improved with a better approach for all-gap regions
        $new_start = $orig_start;
        $new_end = $orig_start - 1;  # Empty range
    }
    
    # Store the results using base ID with new range
    $results{$full_id}{seq} = $region;
    $results{$full_id}{start} = $new_start;
    $results{$full_id}{end} = $new_end;
    $results{$full_id}{base_id} = $base_id;
}

# Output the results
my %seen_sequences;
my @final_output;

foreach my $full_id (@seq_order) {
    my $base_id = $results{$full_id}{base_id};
    my $seq = $results{$full_id}{seq};
    my $output_line = "$base_id/$results{$full_id}{start}-$results{$full_id}{end}";
    
    # Create a key combining the base_id and sequence content
    my $seq_key = "$base_id:$seq";
    
    if (exists $seen_sequences{$seq_key}) {
        # Report duplicate to STDERR
        warn "WARNING: Duplicate sequence found for $base_id - removing duplicate entry $output_line (kept: $seen_sequences{$seq_key})\n";
        next; # Skip this duplicate
    }
    
    # Store this sequence and its identifier
    $seen_sequences{$seq_key} = $output_line;
    push @final_output, [$output_line, $seq];
}

# Print the final deduplicated output
foreach my $entry (@final_output) {
    printf "%-25s %s\n", $entry->[0], $entry->[1];
}
