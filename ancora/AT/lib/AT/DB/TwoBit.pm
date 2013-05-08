# AT (Alternative Transcripts) module for AT::DB::TwoBit
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::DB::TwoBit - interface to .2bit file of DNA sequences

=head1 SYNOPSIS

=over 4

=item * creating a TwoBit object

my $db = $AT::DB::TwoBit->new
    (file => 'sequences.2bit',
     id => 'MySeqs')

=item * retrieving a sequence as a Bio::LocatableSeq object

my $seq = $db->get_seq (
    id => 'NM_123456'',
    start => 100,
    end => 200 );

=back

=head1 DESCRIPTION

This module provides read access to sequence files in .2bit format.
It is a wrapper for the programs twoBitToFa and twoBitInfo. These two
binaries must exist in a directory given in the in the command
path ($PATH environment variable). Alternatively, full paths to the
executables can be passed to the AT::DB::TwoBit constructor.
The binaries are part of the blatSuite which can be obtained from
Jim Kent's web page (http://www.cse.ucsc.edu/~kent/).

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::DB::TwoBit;

use strict;
use vars '@ISA';
use Bio::Seq;
use Bio::LocatableSeq;
use Carp;
use AT::Root;

@ISA = qw(AT::Root);


=head2 new

 Title     : new
 Usage     : my $db = $AT::DB::TwoBit->new
               (file => 'sequences.2bit',
                id => 'MySeqs')
 Function  : Constructor
 Returns   : AT::DB::TwoBit object
 Args      : Required:
             file - .2bit file
             Optional:
	     id - name that will be incorporated into
		identifiers for retrieved sequences
             dummy_mode - return sequences of Ns
             twoBitToFa_binary - path to twoBitToFa
              (default 'twoBitToFa')
             twoBitInfo_binary - path to twoBitInfo
              (default 'twoBitInfo')
             seq_size_cache_size - max nr of sequences
              to cache size for (default 100)

=cut


sub new
{
    my ($caller, %args) = @_;
  
    my $self = bless {
	id => ($args{'id'} || ''),
	twoBitToFa_binary => ($args{twoBitToFa_binary} || 'twoBitToFa'),
	twoBitInfo_binary => ($args{twoBitInfo_binary} || 'twoBitInfo'),
        seq_size_cache_size => ($args{seq_size_cache_size} || 100),
	dummy_mode => $args{dummy_mode},
	_seq_sizes => {},
        _nr_cached_seq_sizes => 0
    }, ref($caller) || $caller;

    unless($args{'file'}) {
	croak "Missing file argument";
    }
    $self->{file} = $args{'file'};

    return $self;
}


sub file { shift->{'file'} }


=head2 get_seq

 Title     : get_seq
 Usage     : my $seq = $db->get_seq (id => "NM_123456",
				     start => 100,
				     end => 200);
 Function  : Retrieves a sequence segment as an object
 Returns   : Bio::LocatableSeq
 Args      : id - id of sequence to retrieve
             start - start coordinate in sequence (optional)
             end - end coordinate in sequence (optional)
             By default the entire sequence is retrieved.
             strand - this optional argument, which can be "-" or "+",
             sets the strand variable of the returned sequence, but has
             no other effect (e.g. the sequence is not reversed).
             Defaults to "+".

=cut

sub get_seq  {
    my ($self, %args) = @_;
    my $id   = $args{'id'} || croak "No id for get_seq";
    my $start = $args{'start'};
    my $end   = $args{'end'};

    my $prefix = $self->id;
    my $full_id = $prefix ? $prefix.'_'.$id : $id;
    my $seq_str = $self->_read_seq_str($id,$start,$end);
    $start = 1 unless($start);
    $end = $start + length($seq_str) - 1;

    my $seqobj  = Bio::LocatableSeq->new
	(-id => $full_id,
	 -seq => $seq_str,
	 -start => $start,
	 -end => $end,
	 -strand =>
	    (($args{'strand'} and ($args{'strand'} eq "-1" or
				   $args{'strand'} eq "-")) ? -1 : 1)
	);

    return $seqobj;
}


=head2 get_seq_str

 Title     : get_seq_str
 Usage     : my $seq_str = $db->get_seq_str (id => "NM_123456",
					     start => 100,
					     end => 200);
 Function  : Retrieves a sequence as a string
 Returns   : String
 Args      : id - id of sequence to retrieve
             start - start coordinate in sequence (optional)
             end - end coordinate in sequence (optional)
             By default the entire sequence is retrieved.

=cut


sub get_seq_str
{
    my ($self, %args) = @_;
    my $id   = $args{'id'} || croak "No id for get_seq";
    my $start = $args{'start'};
    my $end   = $args{'end'};

    my $seq_str = $self->_read_seq_str($id,$start,$end);

    return $seq_str;
}


=head2 get_seq_size

 Title     : get_seq_size
 Usage     : my $size = $db->get_seq_size('NM_123456');
 Function  : Retrieves the size of a sequence in basepairs
 Returns   : Integer
 Args      : Sequence id

=cut


sub get_seq_size
{
    my ($self, $id) = @_;
    my $seq_size;
    unless($seq_size = $self->{_seq_sizes}{$id}) {
	$seq_size = $self->_read_seq_size($id);
	if($self->{_nr_cached_seq_sizes} < $self->{seq_size_cache_size}) {
	    $self->{_seq_sizes}{$id} = $seq_size;
	    $self->{_nr_cached_seq_sizes}++;
	}
    }
    return $seq_size;
}


=head2 get_all_seq_sizes

 Title     : get_all_seq_sizes
 Usage     : my $sizes = $db->get_all_seq_sizes();
 Function  : Retrieves the sizes of all sequences in basepairs
 Returns   : Hash reference (key: sequence id, value: size)
 Args      : None

=cut


sub get_all_seq_sizes
{
    return shift->_read_all_seq_sizes();
}


=head2 get_seq_ids

 Title     : get_seq_ids
 Usage     : my @id_list = $db->get_seq_ids();
 Function  : Returns a list of all sequence ids in the .2bit file
 Returns   : Array
 Args      : -

=cut


sub get_seq_ids
{
    my ($self) = @_;
    my @names;

    my $binary = $self->{twoBitInfo_binary};
    my $file = $self->{file};
    my $cmd = "$binary $file /dev/stdout |";
    open IN, $cmd or croak "failed to open pipe [$cmd]";
    while(my $line = <IN>) {
	chomp $line;
	my ($name, $size) = split /\s+/, $line;
	push @names, $name;
	if(!($self->{_seq_sizes}{$name}) and
	   $self->{_nr_cached_seq_sizes} < $self->{seq_size_cache_size}) {
	    $self->{_seq_sizes}{$name} = $size;
	    $self->{_nr_cached_seq_sizes}++;
	}
    }
    close IN;

    return @names;
}


=head2 get_seq_names

 Synonym for get_seq_ids

=cut


sub get_seq_names { get_seq_ids(@_) } 


sub _read_seq_str
{
    my ($self, $id, $start, $end) = @_;

    if($self->dummy_mode) {
	croak "start and/or end undefined in dummy mode" unless(defined $start and defined $end);
	return 'N'x($end-$start+1);
    }
    
    my $binary = $self->{twoBitToFa_binary};
    my $file = $self->{file};
    my $options = "-seq=$id";
    if($start) {
	$options = $options.' -start='.($start-1);
    }
    if($end) {
	$options = $options.' -end='.$end;
    }
    my $cmd = "$binary $options $file /dev/stdout |";
    open IN, $cmd or croak "failed to open pipe [$cmd]";
    my $header = <IN>;
    croak "failed to read sequence with command [$cmd]" unless($header);
    my @seq_lines;
    while (my $line = <IN>) {
	chomp $line;
	push @seq_lines, $line;
    }
    close IN;
    my $seq_str = join('', @seq_lines);

    return $seq_str;
}


sub _read_seq_size
{
    my ($self, $id) = @_;
    my @names;

    my $binary = $self->{twoBitInfo_binary};
    my $file = $self->{file};
    my $cmd = "$binary $file:$id /dev/stdout |";
    open IN, $cmd or croak "failed to open pipe [$cmd]";
    my $line = <IN>;
    close IN;
    chomp $line;
    my ($read_id, $size) = split /\s+/, $line;
    unless(defined $read_id and defined $size) { croak "Error executing [$cmd]" } 
    unless($id eq $read_id) { croak "Requested id $id does not match retrieved id $read_id" }

    return $size;
}


sub _read_all_seq_sizes
{
    my ($self, $id) = @_;
    my @names;

    my %sizes;
    my $binary = $self->{twoBitInfo_binary};
    my $file = $self->{file};
    my $cmd = "$binary $file /dev/stdout |";
    open IN, $cmd or croak "failed to open pipe [$cmd]";
    while(my $line = <IN>) {
	chomp $line;
	my ($read_id, $size) = split /\s+/, $line;
	unless(defined $read_id and defined $size) { croak "Error executing [$cmd]" } 
	$sizes{$read_id} = $size;
    }
    close IN;
    return \%sizes;
}


1;


