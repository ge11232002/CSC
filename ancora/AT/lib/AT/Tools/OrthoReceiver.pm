############################################################################
# OrthoReceiver.pm - a module for receiving the genomic sequences of two 
# orthologous mapping objects.
#
############################################################################


=head1 NAME

    OrthoReceiver

=head1 SYNOPSIS

    use AT::Tools::OrthoReceiver;
    my $receiver = AT::Tools::OrthoReceiver->new(-hlocation   => $humanLocation,
						 -mlocation   => $mouseLocation);
    my $hseq = $receiver->hseq();
    my $mseq = $receiver->mseq();

=head1 DESCRIPTION

B<OrthoReceiver> is a module that enables extraction of the genomic sequences of two orthologous B<ChromosomalLocation> objects.

=head1 METHODS DESCRIPTION

=cut


############################################################################

package AT::Tools::OrthoReceiver;

use strict;
use warnings;

use lib '/home/malina/Projects/Phylofoot/programs/DEV/AT/lib';
use lib '/home/malina/Projects/Phylofoot/programs/DEV/GENELYNX/lib';

use GeneLynx::MySQLdb;
use GeneLynx::Search::Advanced;
use GeneLynx::Resource;
use GeneLynx::ResourceLinks;
use AT::DB::GenomeMapping;
use AT::DB::GenomeAssembly;
use AT::ChromosomalLocation;

=head2 B<new()>

I<Title:> new()

I<Usage:> $receiver = AT::Tools::OrthoReceiver->new(-hlocation   => $humanLocation,
                                                    -mlocation   => $mouseLocation);

I<Input:> A human B<ChromosomalLocation> object and the orthologous mouse B<ChromosomalLocaion> object.

I<Output:> Returns a new B<AT::Tools::OrthoReceiver> object.

I<Description:>

=cut

sub new {
  my ($caller, %args) = @_;
  my $self  = bless {}, ref $caller || $caller;
  $self->{orgs} = { human => {location  => ($args{-hlocation } or print "No -hlocation\n")},
		    mouse  => {location => ($args{-mlocation} or print "No -mlocation\n")}};
  $self->extract();
  return $self;
}




=head2 B<hseq()>

I<Title:> hseq()

I<Usage:> $receiver->hseq()

I<Input:>

I<Output:> Returns the genomic sequence corresponding to the human B<ChromosomalLocation> object.

=cut

sub hseq               {$_[0] -> {_hseq} }

=head2 B<mseq()>

I<Title:> mseq()

I<Usage:> $receiver->mseq()

<Input:>

I<Output:> Returns the genomic sequence corresponding to the B<ChromosomalLocation> object from mouse.

=cut

sub mseq               {$_[0] -> {_mseq} }

=head2 B<extract()>

I<Title:> extract()

I<Usage:> Used internally by constructor

I<Input:>

I<Output:>

I<Description:> Extracts the genomic sequences for the B<ChromolomalLocation> objects.

=cut

sub extract{
  my ($self, %args) = @_;
  my $hgendb = $self->{orgs}->{human}->{location}->db();
  my $mgendb = $self->{orgs}->{human}->{location}->db();
  my $hseq = $hgendb->get_genome_seq(start  => $self->{orgs}->{human}->{location}->start(),
				     end    => $self->{orgs}->{human}->{location}->end(),
				     chr    => $self->{orgs}->{human}->{location}->chr(),
				     strand => $self->{orgs}->{human}->{location}->strand()
					 );
  if ($self->{orgs}->{human}->{location}->strand() eq "-")  { $hseq = $hseq->revcom; }
  $self->{_hseq}=$hseq;
  my $mseq = $mgendb->get_genome_seq(start  => $self->{orgs}->{mouse}->{location}->start(),
				     end    => $self->{orgs}->{mouse}->{location}->end(),
				     chr    => $self->{orgs}->{mouse}->{location}->chr(),
				     strand => $self->{orgs}->{mouse}->{location}->strand()
				    );
  if ($self->{orgs}->{mouse}->{location}->strand() eq "-")  { $mseq = $mseq->revcom; }
  $self->{_mseq}=$mseq;

}



1;
