############################################################################
# AT::Tools::ChromosomalLocation.pm - module for representaion of the 
# physical coordiantes for a genomic location.
#
############################################################################


=head1 NAME

B<AT::ChromosomalLocation>

=head1 SYNOPSIS

    use AT::ChromosomalLocation;
    my $CL = AT::ChromosomalLocation -> new(-chr    => $chr,
					    -start  => $start,
					    -end    => $end,
					    -db     => $db);

=head1 DESCRIPTION

B<ChromosomalLocation> contains the physical cordinates for a genomic region, and a database object that enables extraction of the sequence corresponding to the region. 

=head1 METHODS DESCRIPTION

=cut


############################################################################

package AT::ChromosomalLocation;

use strict;
use warnings;

use lib '/home/malina/Projects/Phylofoot/programs/DEV/AT/lib';
use lib '/home/malina/Projects/Phylofoot/programs/DEV/GENELYNX/lib';

#use GeneLynx::MySQLdb;
#use GeneLynx::Search::Advanced;
#use GeneLynx::Resource;
#use GeneLynx::ResourceLinks;
#use AT::DB::GenomeMapping;
#use AT::DB::GenomeAssembly;


=head1 B<new()>

I<Title:> new()

I<Usage:> $CL = AT::ChromosomalLocation -> new(-start   => $start,
					       -end     => $end,
					       -chr     => $chr,
					       -strand  => $strand
					       -db      => $db);
					       
I<Input:> All the values needed to represent a chromosomal location:
chromosome name, startposition and  endposition, the strand of the sequence and a database object for extracting the sequence.

I<Output:> Returns a new B<ChromosomalLocation> object.

I<Description:> Constructor

=cut

sub new {
  my ($class, %args) = @_;
  my ($chr, $start, $end, $strand, $db);
  if (defined $args{-chr}){
    $chr = $args{-chr};
    delete $args{-chr};
  }
  else {print "No -chr\n"}
  if (defined $args{-start}){
    $start=$args{-start};
    delete $args{-start};
  }
  else {print "No -start\n"}
  if (defined $args{-end}){
    $end = $args{-end};
    delete $args{-end};
  }
  else {print "No -end\n"}
  if (defined $args{-strand}){
    $strand=$args{-strand};
    delete $args{-strand};
  } 
  else {print "No -strand\n"}
  if (defined $args{-db}){
    $db=$args{-db};
    delete $args{-db};
  } 
  else {print "No -db\n"}
  
  return bless {chr   => $chr,
		start => $start,
		end   => $end,
		strand=> $strand,
		db    => $db
	       },ref $class || $class;
}

=head2 B<start()>

I<Title:> start()

I<Usage:> start = $CL-> start();

I<Input:> 

I<Output:> Returns the start coordinate for the object.

=cut

sub start             {$_[0] -> {start} }

=head2 B<end()>

I<Title:> end()

I<Usage:> end = $CL-> end();

I<Input:> 

I<Output:> Returns the end coordinate for the object.

=cut

sub end            {$_[0] -> {end} }

=head2 B<chr()>

I<Title:> chr()

I<Usage:> chr = $CL-> chr();

I<Input:> 

I<Output:> Returns the chromosome name for the object.

=cut

sub chr             {$_[0] -> {chr} }

=head2 B<strand()>

I<Title:> strand()

I<Usage:> strand = $CL-> strand();

I<Input:> 

I<Output:> Returns the strand of the object.

=cut

sub strand            {$_[0] -> {strand} }

=head2 B<db()>

I<Title:> db()

I<Usage:> db = $CL-> db();

I<Input:> 

I<Output:> Returns the database object for the ChromosomalLocation object. 

=cut

sub db            {$_[0] -> {db} }

1;
