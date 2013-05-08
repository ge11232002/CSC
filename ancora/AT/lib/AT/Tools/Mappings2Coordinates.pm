############################################################################
# AT::Tools::Mappings2Coordiantes.pm - module for calculation of physical 
# coordinates for a Mapping object.
#
############################################################################


=head1 NAME

B<AT::Tools::Mappings2Coordinates>

=head1 SYNOPSIS

    use AT::Tools::Mappings2Coordinates;
    my $m2c = AT::Tools::Mappings2Coordinates -> new( -hmapping   => $hmapping,
						      -hdb        => $hgendb,
						      -mmapping   => $mmapping,
						      -mdb        => $mgendb,
						      -upstream   => $upstream,
						      -downstream => $downstream);
    my $humanLocation = $m2c -> humanLocation();
    my $mouseLocation = $m2c -> mouseLocation();

=head1 DESCRIPTION

B<Mappings2Coordinates> takes two B<Mapping> objects (from orthologous genes in human and mouse) as input and creates corresponding B<ChromosomalLocation> objects.

=head1 METHODS DESCRIPTION

=cut


############################################################################

package AT::Tools::Mappings2Coordinates;

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


=head1 B<new()>

I<Title:> new()

I<Usage:> $m2c = AT::Tools::Mappings2Coordinates -> new(-hmapping   => $hmapping,
						        -hdb        => $hgendb,
						        -mmapping   => $mmapping,
						        -mdb        => $mgendb,
						        -upstream   => $upstream,
						        -downstream => $downstream);

I<Input:> A human mapping object, a mouse mappign object, database objects for the human and mouse genome assembly databases, relative positions for the calculation of physical coordiantes (i.e. number of bases to include upstream and downstream of the coding part of the gene).

I<Output:> Returns a new B<Mapping2Coordiantes> object.

I<Description:> Constructor

=cut

sub new {
  my ($caller, %args) = @_;
  my $self  = bless {}, ref $caller || $caller;
  my $upstream=0;
  my $downstream=0;
  if (defined $args{-upstream}){
    $upstream=$args{-upstream};
    delete $args{-upstream};
  }
  if (defined $args{-downstream}){
    $downstream=$args{-downstream};
    delete $args{-downstream};
  }
  $self->{orgs} = { human => {mapping => ($args{-hmapping} or print "No -hmapping\n"),
			      db      => ($args{-hdb} or print "No -hdb\n")},
		    mouse  => {mapping => ($args{-mmapping} or print "No -mmapping\n"),
			       db      => ($args{-mdb} or print "No -mdb\n")}};
  $self->{relPos} = { upstream   =>  ($upstream),
		      downstream =>  ($downstream)};
  $self->chromosomalLocation();
  return $self;
}

=head2 B<chromosomalLocation()>

I<Title:> chromosomalLocation()

I<Usage:> $m2c = -> chromosomaLocation();

I<Input:> 

I<Output:> 

I<Description:> Used internally in constructor to create ChromosomalLocation objects for each mapping. The ChormosomalLocation objects are objects containing the physical coordinates (base pair positions and chromosome number) for a chromosomal location.

=cut

sub chromosomalLocation {
  my ($self, %args) = @_;
  my $hmapping = $self->{orgs}->{human}->{mapping};
  my $mmapping = $self->{orgs}->{mouse}->{mapping};

  my $upstream    = $self->{relPos}->{upstream};
  my $downstream  = $self->{relPos}->{downstream};

  my $hupstream   = $upstream;
  my $hdownstream = $downstream;
  my $mupstream   = $upstream;
  my $mdownstream = $downstream;

  if ($hmapping -> strand eq "-"){
    $hupstream = $downstream;
    $hdownstream = $upstream;}

  if ($mmapping -> strand eq "-"){
    $mupstream = $downstream;
    $mdownstream = $upstream;}

  my $humanLocation = AT::ChromosomalLocation -> new(-start   => $hmapping -> tStart-$hupstream,
						     -end     => $hmapping -> tEnd+$hdownstream,
						     -chr     => $hmapping -> tName,
						     -strand  => $hmapping -> strand,
						     -db      => $self->{orgs}->{human}->{db});
  my $mouseLocation = AT::ChromosomalLocation -> new(-start   => $mmapping -> tStart-$mupstream,
						     -end     => $mmapping -> tEnd+$mdownstream,
						     -chr     => $mmapping -> tName,
						     -strand  => $mmapping -> strand,
						     -db      => $self->{orgs}->{mouse}->{db});
  $self -> {humanLocation} = $humanLocation;
  $self -> {mouseLocation} = $mouseLocation;
}

=head2 B<humanLocation()>

I<Title:> humanLocation()

I<Usage:> $m2c = -> humanLocation();

I<Input:> 

I<Output:> returns the human B<ChromosomalLocation> objecet.

=cut

sub humanLocation             {$_[0] -> {humanLocation}}

=head2 B<mouseLocation()>

I<Title:> mouseLocation()

I<Usage:> $m2c = -> mouseLocation();

I<Input:> 

I<Output:> returns the mouse B<ChromosomalLocation> objecet.

=cut

sub mouseLocation             {$_[0] -> {mouseLocation}}

1;
