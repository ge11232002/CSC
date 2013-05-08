# AT (Alternative Transcripts) module for AT::GenomeAssemblyTwoBit
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::DB::GenomeAssemblyTwoBit - interface to a .2bit file for a
genome assembly

=head1 SYNOPSIS

=over 4

=item * creating a GenomeAssemblyTwobit object

my $db = AT::DB::GenomeAssemblyTwoBit->new (
    file => '/data/goldenpath/DroVir_JUL04/scaffolds.2bit'
    id => 'DroVir_JUL04' );

=item * retrieving a genome segment as a Bio::LocatableSeq object

my $genomeseq = $db->get_genome_seq (
    chr => 'chr1',
    start => 1000,
    end => 2000 );

=back

=head1 DESCRIPTION

This class provides read access to .2bit files containing genome
assemblies. It is a subclass of AT::DB::TwoBit. The only addition in this
subclass is to provide the interface of AT::DB::GenomeAssembly.
For additional documentation see AT::DB::TwoBit.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::DB::GenomeAssemblyTwoBit;

use strict;
use vars '@ISA';
use Bio::Seq;
use Bio::LocatableSeq;
use Carp;
use AT::DB::TwoBit;

@ISA = qw(AT::DB::TwoBit);


sub new  {
    my ($caller, %args) = @_;
    my $self = $caller->SUPER::new(%args);

    if($args{assembly_name}) {
	carp("Argument assembly_name is deprecated, use argument id instead");
	$self->id($args{assembly_name});
    }

    return $self;
}

=head2 get_genome_seq

 Title     : get_genome_seq
 Usage     : my $genomeseq = $db->get_genome_seq (chr => "chr1",
                                                  start => 1000,
                                                  end => 2000);
 Function  : Retrieves a genome segment as an object
 Returns   : Bio::LocatableSeq
 Args      : chr - id of sequence (chromosome/scaffold/contig etc) to get
             start, end - specifies the location of the segment (optional)
             strand - this optional argument, which can be "-" or "+",
             sets the strand variable of the returned sequence, but has
             no other effect (e.g. the sequence is not reversed).
             Defaults to "+".

=cut

sub get_genome_seq  {
    my ($self, %args) = @_;
    return $self->get_seq(%args, id => $args{chr});
}


=head2 get_genome_seq_str

 Title     : get_genome_seq_str
 Usage     : my $seq_str = $db->get_genome_seq_str (chr => "chr1",
                                                    start => 1000,
                                                    end => 2000);
 Function  : Retrieves a genome segment as a string
 Returns   : String
 Args      : chr - id of sequence (chromosome/scaffold/contig etc) to get
             start, end - specifies the location of the segment (optional)

=cut


sub get_genome_seq_str
{
    my ($self, %args) = @_;
    return $self->get_seq_str(%args, id => $args{chr});
}


=head2 get_chr_size

 Title     : get_chr_size
 Usage     : my $size = $db->get_chr_size('chr1');
 Function  : Retrieves the size of a sequence in basepairs
 Returns   : Scalar
 Args      : Sequence id, e.g. 'scaffold_1' or 'chr1'

=cut


sub get_chr_size { shift->get_seq_size(@_) }


=head2 get_all_chr_sizes

 Title     : get_all_chr_sizes
 Usage     : my $sizes = $db->get_all_chr_sizes();
 Function  : Retrieves the sizes of all sequences in basepairs
 Returns   : Hash reference (key: sequence id, value: size)
 Args      : None

=cut


sub get_all_chr_sizes { shift->get_all_seq_sizes(@_) }


=head2 get_chr_names

 Title     : get_chr_names
 Usage     : my @list = $db->get_chr_names();
 Function  : Returns list ids for sequences in the .2bit fil
 Returns   : Array
 Args      : -

=cut

sub get_chr_names { shift->get_seq_ids(@_) }


1;


