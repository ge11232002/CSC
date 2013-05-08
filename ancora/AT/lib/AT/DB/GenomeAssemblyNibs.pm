# AT (Alternative Transcripts) module for AT::GenomeAssemblyibs
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::DB::GenomeAssemblyNibs - interface to a collection of .nib files for a
genome assembly

=head1 SYNOPSIS

=over 4

=item * creating a nib object

my $db = AT::DB::GenomeAssemblyNibs->new (
    id => 'HS_APR03'
    dir => '/data/goldenpath/HS_APR03/nib' );

=item * retrieving a genome segment as a Bio::LocatableSeq object

my $genomeseq = $db->get_genome_seq (
    chr => 1,
    start => 4222301,
    end => 4282300 );

=back

=head1 DESCRIPTION

This module is an alternative to AT::DB::GenomeAssembly for quicker access
to sequence data when .nib-files are available locally.
Let's try to maintain a common interface for the two modules!

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::DB::GenomeAssemblyNibs;

use strict;
use vars '@ISA';
use Bio::Seq;
use Bio::LocatableSeq;
use Carp;
use AT::Root;

@ISA = qw(AT::Root);


use constant NIB_SIG => scalar 0x6BE93D3A;


=head2 new

 Title     : new
 Usage     : my $nib_catalog = $AT::DB::GenomeAssemblyNibs->new
		(dir => '/data/goldenpath/HS_NOV02/nib/',
		 id => 'HS_NOV02')
 Function  : Constructor
 Returns   : AT::DB::GenomeAssemblyNibs object
 Args      : dir - where the nib files are
	     id - name that will be incorporated into
		identifiers for retrieved sequences;
		optional, but recommended for applications
		that use several assemblies
             dummy_mode - if true, sequences of Ns will be returned

=cut


sub new
{
    my ($caller, %args) = @_;
  
    my $self = bless {
	id => ($args{'id'} || $args{'assembly_name'} || ''),
	dummy_mode => ($args{dummy_mode} || 0),
	_chr_sizes => {}
    }, ref($caller) || $caller;

    carp("Argument assembly_name is deprecated, use argument id instead") if($args{'assembly_name'});

    unless($args{'dir'}) {
	die "Missing dir argument";
    }
    $self->dir($args{'dir'});

    return $self;
}


=head2 dir

 Title     : dir
 Usage     : $dir = $nib_catalog->dir();
	     $nib_catalog->dir($new_dir);
 Function  : Get/set dir
 Returns   : dir as string
 Args      : dir as string (optional)

=cut


sub dir
{
    my ($self, $dir) = @_;
    if(defined $dir) {
	$dir .= '/' unless($dir =~ /\/$/);
	$self->{'dir'} = $dir;
	#$self->_set_assembly_name_from_dir();
    }
    return $self->{'dir'};
}


sub _set_id_from_dir
{
    my ($self);
    my ($name) = $self->dir =~/(\w+)\/$/;
    $self->{'id'} = $name || '';
}


=head2 get_genome_seq

 Title     : get_genome_seq
 Usage     : my $genomeseq = $db->get_genome_seq (chr => "Y",
                                                  start => 4222301,
                                                  end => 4282300);
 Function  : Retrieves a genome segment as an object
 Returns   : Bio::LocatableSeq
 Args      : chr, start, end - specifies the location of the segment
             strand - this optional argument, which can be "-" or "+",
             sets the strand variable of the returned sequence, but has
             no other effect (e.g. the sequence is not reversed).
             Defaults to "+".

=cut

sub get_genome_seq  {
    my ($self, %args) = @_;
    my $start = $args{'start'}; croak "No start for get_genome_seq" unless defined $start;
    my $end   = $args{'end'};   croak "No end for get_genome_seq" unless defined $end;
    my $chr   = $args{'chr'} || croak "No chr for get_genome_seq";
    $chr = $self->_chr_name($chr);

    unless($self->id) {
	warn "Id is not set for ".ref($self)." object with dir = ".
	    $self->dir;
    }
  
    my $seqobj  = Bio::LocatableSeq->new
	(-id => $self->id."_".$chr,
	 -seq => $self->get_genome_seq_str(%args),
	 -start => $start,
	 -end => $end,
	 -strand =>
	    (($args{'strand'} and ($args{'strand'} eq "-1" or
				   $args{'strand'} eq "-")) ? -1 : 1)
	);

    return $seqobj;
}


=head2 get_genome_seq_str

 Title     : get_genome_seq_str
 Usage     : my $seq_str = $db->get_genome_seq_str (chr => "Y",
                                                    start => 4222301,
                                                    end => 4282300);
 Function  : Retrieves a genome segment as a string
 Returns   : String
 Args      : chr, start, end - specifies the location of the segment

=cut


sub get_genome_seq_str
{
    my ($self, %args) = @_;
    my $start = $args{'start'}; croak "No start for get_genome_seq" unless defined $start;
    my $end   = $args{'end'};   croak "No end for get_genome_seq" unless defined $end;
    my $chr   = $args{'chr'} || croak "No chr for get_genome_seq";

    return 'N'x($end-$start+1) if($self->dummy_mode);

    $chr = $self->_chr_name($chr);

    my ($fh, $chr_size) = $self->_nib_open_verify($chr);
    if($start < 1 or $start > $end or $end > $chr_size) {
	carp "Invalid start and end ($start, $end) for $chr of size $chr_size bp";
	return "";
    }

    my $str = $self->_nib_input($fh, $start, $end, 0);
    close $fh;  # should keep a few nibs open
    return $str;
}


=head2 get_chr_size

 Title     : get_chr_size
 Usage     : my $size = $db->get_chr_size('chr21');
 Function  : Retrieves the size of a chromosome in basepairs
 Returns   : Scalar
 Args      : Chromosome name, e.g. 'chr21'

=cut


sub get_chr_size
{
    my ($self, $chr) = @_;
    $chr = $self->_chr_name($chr);
    my $chr_size;
    unless($chr_size = $self->{_chr_sizes}{$chr}) {
	my $fh;
        ($fh, $chr_size) = $self->_nib_open_verify($chr);
	$self->{_chr_sizes}{$chr} = $chr_size;
	close $fh;  # should keep a few nibs open
    }
    return $chr_size;
}


=head2 get_chr_names

 Title     : get_chr_names
 Usage     : my @list = $db->get_chr_names();
 Function  : Returns a list of chromosome names, one for each
             .nib file in the directory.
 Returns   : Array
 Args      : -

=cut


sub get_chr_names
{
    my ($self) = @_;
    my @chr_list;
    foreach my $fn (glob($self->dir.'*.nib')) {
	my ($chr_name) = $fn =~ /(\w+).nib$/;
	push @chr_list, $chr_name;
    }
    return @chr_list;
}


sub _nib_open_verify
{
    my ($self, $chr) = @_;

    my $binary_input;
    my $fn = $self->_filename_from_chr($chr);
    open(NIB_IN, '<:raw', $fn) or die "Could not open .nib file $fn";
    unless(read(NIB_IN, $binary_input, 8) == 8)
	{ croak "Error reading header from .nib file $fn"; }
    my ($sig, $size) = unpack('N N', $binary_input);
    unless($sig == NIB_SIG) {
	($sig, $size) = unpack('V V', $binary_input);
	unless($sig == NIB_SIG)
	    { croak "$fn is not a good .nib file"; }
    }

    return(\*NIB_IN, $size);    
}


sub _nib_input
{
    my ($self, $fh, $start, $end) = @_;

    # test
    # - quicker to use getc than read?
    # - quicker to use seqstring that seqarray? can i set req string size?

    my @valToNtTbl = ('T','C','A','G','N',undef,undef,undef,
		      't','c','a','g','n');

    $start--;
    my $size = $end - $start;

    my @seq;
    $#seq = $size - 1;
    my $seqpos = 0;

    my $bytePos = $start >> 1;
    seek($fh, $bytePos + 8, 0);

    my $bVal;
    if($start & 1) {
	unless(read($fh, $bVal, 1) == 1)
	    { croak "Read error in .nib file"; }
	$bVal = ord($bVal);
	$seq[$seqpos++] = $valToNtTbl[$bVal&0xf];
	$size--;
    }
    my $byteSize = $size >> 1;
    while(--$byteSize >= 0) {
	unless(read($fh, $bVal, 1) == 1)
	    { croak "Read error in .nib file"; }
	$bVal = ord($bVal);
	$seq[$seqpos++] = $valToNtTbl[$bVal>>4];
	$seq[$seqpos++] = $valToNtTbl[$bVal&0xf];
    }
    if($size & 1) {
	unless(read($fh, $bVal, 1) == 1)
	    { croak "Read error in .nib file"; }
	$bVal = ord($bVal);
	$seq[$seqpos++] = $valToNtTbl[$bVal>>4];
    }

    return join('', @seq);
}


sub _filename_from_chr
{
    my ($self, $chr) = @_;
    return $self->{'dir'}.$chr.'.nib';
}


sub _chr_name
{
    my ($self, $chr) = @_;
    croak "No chr argument" unless($chr);
    $chr = "chr$chr" unless $chr =~ /^chr/;
    return $chr;
}


1;


