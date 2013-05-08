# AT (Alternative Transcripts) module for AT::MySQLdb
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::DB::GenomeAssembly - interface to genome assembly database

=head1 SYNOPSIS

=over 4

=item * creating a database object by connecting to an existing AT-type assembly database

my $db = AT::DB::GenomeAssembly->connect (
    -dbname => "HS_APR02",
    -dbhost => "myhost.mydomain",
    -dbuser => "myusername",
    -dbpass => "mypassword" );

=item * retrieving a genome segment as a Bio::LocatableSeq object

my $genomeseq = $db->get_genome_seq (
    chr => 1,
    start => 4222301,
    end => 4282300 );

=back

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=head2 connect

 Title     : connect
 Usage     : my $db = $AT::DB::GenomeAssembly->connect(%args);
 Function  : Creates a new object of this class.
 Note      : For details, see docs for superclass AT::DB::MySQLdb.

=cut

package AT::DB::GenomeAssembly;

use strict;
use vars '@ISA';
use Bio::Seq;
use Bio::LocatableSeq;
use DBI;
use Carp;
use AT::Root;
use AT::DB::MySQLdb;


@ISA = qw(AT::DB::MySQLdb);


use constant DEFAULT_CHUNKSIZE => scalar 50000;


sub _init {
    my ($self, %args) = @_;
    $self->{chunksize} = $args{chunksize} || DEFAULT_CHUNKSIZE;
}


sub id {
    return $_[0]->dbname;
}

=head2 get_genome_seq

 Title     : get_genome_seq
 Usage     : my $genomeseq = $db->get_genome_seq (chr => "Y",
                                                  start => 4222301,
                                                  end => 4282300);
 Function  : Retrieves a genome segment
 Returns   : Bio::LocatableSeq
 Args      : chr, start, end - specifies the location of the segment
             strand - this optional argument, which can be "-" or "+",
             sets the strand variable of the returned sequence, but has
             no other effect (e.g. the sequence is not reversed).
             Defaults to "+".

=cut

sub get_genome_seq  {
    my ($self, %args) = @_;
    #my $genome = $args{'genome'} || croak "No genome specified";
    #$genome = uc $genome;
    my $start = $args{'start'} || croak "No start for get_genome_seq";
    my $end   = $args{'end'} || croak "No end for get_genome_seq";
    my $chr   = $args{'chr'} || croak "No chr for get_genome_seq";
    $chr = $self->_chr_name($chr);
    my $table = $self->_chr_table($chr);
    my $sth = $self->dbh->prepare(qq!SELECT start, end, sequence from $table
			    WHERE chr = '$chr' and start >= ? and start <= ?
			    order by start!);
    $sth->execute($start-$self->{chunksize}+1, $end);
    my $seqstring = "";
    while (my ($frag_start, $frag_end, $seqfragment) = $sth->fetchrow_array)  {
	my $ss_offset = ($frag_start < $start) ? $start - $frag_start : 0;
	my $ss_length = ($frag_end > $end) ?
	    $end - $frag_start - $ss_offset + 1 : $frag_end - $frag_start + 1;
        $seqstring .= substr($seqfragment, $ss_offset, $ss_length);
    }
    unless(length($seqstring) == ($end - $start + 1)) {
	warn "Could not get sequence $chr:$start-$end\n";
	return undef;
    }
    my $requested_seqobj  = Bio::LocatableSeq->new
	(-id => $self->id."_".$chr,
	 -seq => $seqstring,
	 -start => $start,
	 -end => $end,
	 -strand => 1);
    if ($args{'strand'} and 
	($args{'strand'} eq "-1" or $args{'strand'} eq "-")) 
    {
        $requested_seqobj->strand(-1);
    }
    return $requested_seqobj;
}


sub get_genome_seq_str {
    my ($self, %args) = @_; 
    my $start = $args{'start'} || croak "No start for get_genome_seq";
    my $end   = $args{'end'} || croak "No end for get_genome_seq";
    my $chr   = $args{'chr'} || croak "No chr for get_genome_seq";
    $chr = $self->_chr_name($chr);
    my $table = $self->_chr_table($chr);
    my $sth = $self->dbh->prepare(qq!SELECT start, end, sequence from $table
			    WHERE chr = '$chr' and start >= ? and start <= ?
			    order by start!);
    $sth->execute($start-$self->{chunksize}+1, $end);
    my $seqstring = "";
    while (my ($frag_start, $frag_end, $seqfragment) = $sth->fetchrow_array)  {
	my $ss_offset = ($frag_start < $start) ? $start - $frag_start : 0;
	my $ss_length = ($frag_end > $end) ?
	    $end - $frag_start - $ss_offset + 1 : $frag_end - $frag_start + 1;
        $seqstring .= substr($seqfragment, $ss_offset, $ss_length);
    }
    unless(length($seqstring) == ($end - $start + 1)) {
	warn "Could not get sequence $chr:$start-$end\n";
	return undef;
    }
    return $seqstring;
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
    my $table = $self->_chr_table($self->_chr_name($chr));
    my $sth = $self->dbh->prepare("SELECT max(end) from $table");
    $sth->execute();
    my ($size) = $sth->fetchrow_array;
    croak "Chromosome size query failed for chr $chr, table $table"
	unless ($size);
    return $size;
}


=head2 store_chromosome

 Title     : store_chromosome
 Usage     : $db->store_chromosome( file => "chr1.fa" )
 Function  : Stores the sequence of a whole chromosome in the database
 Returns   : 1 if successful
 Args      : file - the name of a fasta sequence file
             fh - a filehandle to a fasta sequece file 
             Either a file or fh is required. Only the first sequence
             in the file/stream will be considered (but the header of a
             following sequence will be read and discarded, making the
             method unsuitable for successive calls on the same fh).
             The fasta header should contain the chromosome name on
             the form chrN.
 Note      : This method does not read the whole chromosome sequence as
             an object using Bio::SeqIO, since handling chromosome-sized
             sequence objects is slow and memory-consuming. A simple
             fasta-parser is implemented in the method and only a part
             of the sequence is kept in memory at any instant.

=cut


sub store_chromosome {
    my ($self, %args) = @_;

    # Check args and open file if required
    my ($fh, $filename);
    if($args{'fh'}) {
	$fh = $filename = $args{'fh'}
    }
    elsif($args{'file'}) {
	$filename = $args{'file'};
	unless(open(CHR_FILE, $filename)) {
	    warn "Could not open file $filename";
	    return undef;
	}
	$fh = \*CHR_FILE;
    }
    else {
	warn "No file or fh argument";
	return undef;
    }
    my $chunksize = $self->{chunksize};

    # Read header
    my $header;
    while(defined($header = <$fh>) and not ($header =~ /^>/)) {};
    unless(defined($header)) {
	warn "No fasta entry in file $filename";
	return undef;
    }
    chomp $header;
    my ($chr) = ($header =~ />\s*(chr\w+)/);
    unless(defined($chr)) {
	warn ("Fasta header should be of the form \">chrN\", not $header ". 
	      "in file $filename");
	return undef;
    }
    $chr = $self->_chr_name($chr);

    # Create table for chromosome
    my $table_name = $self->_chr_table($chr);
    $self->dbh->do("CREATE TABLE $table_name (".
		   "fragment_id mediumint unsigned not null auto_increment, ".
		   "chr varchar(16), ".
		   "start bigint unsigned, ".
		   "end bigint unsigned, ".
		   "sequence text, ".
		   "primary key (fragment_id), ".
		   "key (start), ".
		   "key (end))")
	or do {
	    warn "Could not create chromosome table $table_name";
	    return undef;
	};

    # Prepare insert query
    my $sth = $self->dbh->prepare("INSERT INTO $table_name ".
				  "(chr, start, end, sequence) ".
				  "VALUES(?,?,?,?)");

    # Read the sequence in chunks of $chunksize bases and store the chunks
    my $seq = "";
    my $n = 1;
    while(my $line = <$fh>) {
	last if ($line =~ /^>/);
	$line =~ s/\s//g;
	$seq .= $line;
	if (length($seq) >= $chunksize) {
	    $sth->execute($chr, $n, $n+$chunksize-1,
			  substr($seq, 0, $chunksize)); 
	    substr($seq, 0, $chunksize) = "";
	    $n += $chunksize;
	}
    }
    unless($seq eq "") {
	    $sth->execute($chr, $n, $n+length($seq)-1, $seq);
	    $n += length($seq);
    }

    # Close file if we opened it and return number of bases stored
    close($fh) if ($args{file});
    return $n;
}


# AT::DB::Importer is not used for importing chromosomes
# since it does not seem to be more effective than successive INSERTs.
# The explanation is that the chromosome chunk records are few in relation
# to their size.
# Anyway, below is a _create_chunk_importer method in case anyone
# wants to use the importer for chromosome import
#
# Here's the code to make use of it:
# Create importer object (for batch import)
# my $importer = $self->_create_chunk_importer($table_name);
# Store a chunk (in loop)
# my $chunk = Bio::LocatableSeq->new( -seq => $seq->subseq($i+1, $end),
#				    -start => $i+1,
#				    -end => $end);
# $importer->store_object($chunk);
# Finish
# $importer->finish;

sub _create_chunk_importer {
    my ($self, $table) = @_;

    require AT::DB::Importer;

    # define wrapper for storing one chunk
    my $storing_wrapper = sub {
	my ($importer, $seqobj) = @_;
	unless (ref($seqobj) and $seqobj->isa("Bio::LocatableSeq")) {
	    croak "Wrong sequence object type";
	}
	$importer->store_record($table,
				[ $seqobj->id, $seqobj->start,
				  $seqobj->end, $seqobj->seq ]);
   };

    # prepare table info needed by importer
    my $tables_columns = {
	$table => [ 'chr', 'start', 'end', 'sequence' ],
    };

    # create importer object
    my $importer = AT::DB::Importer->new( db => $self,
					  tables_columns => $tables_columns,
					  storing_wrapper => $storing_wrapper,
					  lock => 0); # no lock required

    # return the creation
    return $importer;
}


sub _chr_name
{
    my ($self, $chr) = @_;
    croak "No chr argument" unless($chr);
    $chr = "chr$chr" unless $chr =~ /^chr/;
    return $chr;
}


sub _chr_table
{
    my ($self, $chr) = @_;
    croak "No chr argument" unless($chr);
    return "GENOME_$chr";
}


1;


