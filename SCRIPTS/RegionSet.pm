=head1 NAME

RegionSet

=head1 DESCRIPTION

B<RegionSet> contains a collection of region objects.

=head1 AUTHOR

B<Alex Bateman> Email agb@sanger.ac.uk

=cut

#
# Perl Module for RegionSet
#
# Cared for by Alex Bateman <agb@sanger.ac.uk>
#

package RegionSet;

use vars qw($AUTOLOAD @ISA @EXPORT_OK);
use Exporter;
use Carp;
use strict;
use Region;

#
# Place functions/variables you want to *export*, ie be visible from the caller package into @EXPORT_OK
#

@EXPORT_OK = qw();

#
# @ISA has our inheritance.
#

@ISA = ( 'Exporter' );


#######
# new #
#######

=head1 new

Creates RegionSet object.

$regionset = new RegionSet($id);

=cut

sub new {
    my ($class,$id)=@_;

    my $set = {
	"id"      => $id,
	"regions" => [],
    };

    bless $set, $class;
    return $set;
}

#######
# add #
#######

=head1 add

Add region to set object

    $regionset->add($region1..$regionN);

=cut
sub add {
  my $self = shift @_;
  push(@{$self->{'regions'}},@_);
}


########
# each #
########

=head1 each

Loops over and returns a list of objects

    @list = $self->each();

=cut

sub each {
  my $self= shift @_;
  my @list=@{$self->{'regions'}};
  return @list;
}

#########
# write #
#########

=head1 write

Loops over List and writes element

    $self->write();

=cut

sub write {
    my $self = shift @_;

    my @list=$self->each();
    foreach my $element (@list){
	$element->write();
    }
}




1;  # says use was ok
__END__

