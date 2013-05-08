############################################################################
# SequenceMapper.pm - a module for creating mappings of exons for a gene.
#
############################################################################


=head1 NAME

SequenceMapper - module for creating exon mappings for a gene.

=head1 SYNOPSIS

    use AT::SequenceMapper;
    $mapper = AT::SequenceMapper -> new(-extres  => $exrtes, 
				        -id      => $id,
					-gldb    => $gldb,
					-species => "mouse", 
					-hgldb   => $hgldb, 
					-mgldb   => $mgldb, 
					-mmapdb  => $mmapdb)
                                        );

    $humanmapping = $mapper -> get_mapping();

=head1 DESCRIPTION

B<SequenceMapper> is a module for receiving a B<mapping> object for a gene id.

=head1 METHODS DESCRIPTION

=cut


############################################################################

package AT::Tools::SequenceMapper;

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


=head2 B<new()>

I<Title:> new()

I<Usage:> $mapper = AT::Tools::SequenceMapper->new(-extres => word_index, 
					      -id          => "pai1", 
                                              -species     => "mouse", 
                                              -hgldb       => $hgldb, 
                                              -mgldb       => $mgldb, 
                                              -mmapdb      => $mmapdb)
					     );

I<Input:> Resourse, id, species, human glid database, mouse glid database(for mouse sequences) human map database (for human sequences) and mouse map database (for mouse database).

I<Output:> Returns a new B<AT::Tools::SequenceMapper> object.

=cut
  
sub new{
  my ($caller, %args) = @_;
  my $self = bless {}, ref $caller || $caller;
  my ($extres, $id, $species);
  if (defined $args{-extres}){
    $self->{_extres} = $args{-extres};
    delete $args{-extres};
  }
  else{
    print "Must provide extres\n";
  }
  if (defined $args{-id}){
    $self->{_id} = $args{-id};
    delete $args{-id};
  }
  else{
    print "Must provide id\n";
  }
  if (defined $args{-hgldb}){
    $self->{_hgldb} = $args{-hgldb};
    delete $args{-hgldb};
  }
  else{
    print "Must provide human glid database \n";
  }
  if (defined $args{-mgldb}){
    $self->{_mgldb} = $args{-mgldb};
    delete $args{-mgldb};
  }
  elsif($args{-species} eq "mouse"){
    print "must provide mouse glid database\n";
  }

  if (defined $args{-hmapdb}){
    $self->{_hmapdb} = $args{-hmapdb};
    delete $args{-hmapdb};
  } 
  elsif( $args{-species} eq "human"){
    print "Must provide a human  map database\n";
  }
  if (defined $args{-mmapdb}){
    $self->{_mmapdb} = $args{-mmapdb};
    delete $args{-mmapdb};
  }
  elsif ($args{-species} eq "mouse") {
    print "Must provide a mouse map database\n";
  }
  if (defined $args{-species}){
    $self->{_species} = $args{-species};
    delete $args{-species};
  }
  else{
    print "Must provide species\n";
  }
  
  $self->find_glids();
  $self->find_mapping();
  return $self;  
}


=head2 B<get_mapping()>

I<Title:> get_mapping()

I<Usage:> $mapper->get_mapping()

I<Input:>

I<Output:> Returns the mapping object.

=cut

sub get_mapping               {$_[0] -> {_mapping} }


=head2 B<find_mapping()>

I<Title:> find_mapping()

I<Usage:> used internally

I<Input:>

I<Output:>

I<Description:> Performs the mapping

=cut


sub find_mapping{
  my ($self, %args) = @_;
  my $species = $self->{_species};
  my ($gl, $gldb, $mapdb);
  if ($species eq "human"){
    $gl=$self->{_glhuman};
    $gldb=$self->{_hgldb};
    $mapdb=$self->{_hmapdb};
  }
  elsif($species eq "mouse"){
    $gl=$self->{_glmouse};
    $gldb=$self->{_mgldb};
    $mapdb=$self->{_mmapdb};
  }
  my $resource = GeneLynx::Resource->new("ucsc_golden_path");
  my $mapping;
  # Get coordinates and acc 
  my $linksobj = GeneLynx::ResourceLinks->new(-resourceobj => $resource,
					      -db          => $gldb,
					      -glid        => $gl
					     );
  foreach my $raw_id  ($linksobj->_raw_id_list)  {
    my ($chr, $start, $end, $strand, $acc) = split ":", $raw_id;
    ($mapping) = $mapdb->get_mappings_for_acc($acc);
    $self->{_mapping} = $mapping;
  }
}



=head2 B<find_glids()>

I<Title:> get_glid()

I<Usage:> used internally

I<Input:>

I<Output:>

I<Description:> Finds the human and mouse glids

=cut


sub find_glids{
  my ($self, %args) = @_;
  my $extres=$self->{_extres};
  my $id=$self->{_id};
 

  my $advsearch = GeneLynx::Search::Advanced->new
      (-query  => ['AND', $extres, $id],
       -db     => $self->{_hgldb});

  my $resource2 = GeneLynx::Resource->new("glmouse");
  my @my_mapping;
  foreach my $hit ($advsearch->hit_list)  {
    # get human and mouse glid
    my $linksobj2 = GeneLynx::ResourceLinks->new(-resourceobj => $resource2,
						 -db =>   $self->{_hgldb},
						 -glid => $hit->glid
						);
    my ($raw_id2) =  $linksobj2->_raw_id_list;
    my ($glhuman, $glmouse) = split ":", $raw_id2;
    $self->{_glhuman}=$glhuman;
    $self->{_glmouse}=$glmouse;
  }
}

1;
