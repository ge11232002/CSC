# AT module for AT::AxtIO
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::AxtIO - read alignments in axt formats

=head1 SYNOPSIS

use AT::AxtIO;

 my $aln_in = AT::AxtIO->new(-fh => \*STDIN);

 while(my $a = $aln_in->next_aln) {
    # Do something with aligment
    # It's a Bio::SimpleAlign so BioPerl applies
 }


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut



package AT::AxtIO;

use strict;
use vars '@ISA';
use Carp;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use AT::Root;

@ISA = qw(AT::Root);

=head2 new

 Title     : new
 Usage     : my $aln_in = $AT::AxtIO->new(-file => 'aln.axt');
 Function  : Constructor.
 Returns   : AT::AxtIO
 Args      : -fh	  	filehandle
	     -file		filename
	     Give either -fh or -file. If -file is given,
	     the method will try to open the file with that name.

=cut

sub new {
    my ($caller, %args) = @_;

    my $fh = $args{'-fh'};
    unless ($fh) {
	unless (exists $args{'-file'})  {
	    $fh = \*STDIN;
	}
	else  {
	    open FILE, $args{'-file'} 
	    or croak "Could not open ".$args{'-file'};
            $fh = \*FILE;
	}
    }

    my $self = bless { _fh => $fh,
		       }, ref $caller || $caller;
    return $self;
}

=head2 next_aln

 Title     : next_aln
 Usage     : my $aln = $axt_in->next_aln();
 Function  : Gets the next alignment in the stream
 Returns   : Bio::SimpleAlign
 Args      : -

=cut

sub next_aln {
    my ($self) = @_;
    my $fh = $self->{'_fh'};

    # Read data from file
    my $info_line = <$fh>;
    while($info_line and $info_line =~ /^\#/) {
	$info_line = <$fh>;
    }
    return unless($info_line);
    my $seqstr1 = <$fh>;
    my $seqstr2 = <$fh>;
    unless($seqstr2) {
	warn "Early EOF in axt file";
	return;
    }
    <$fh>;  # each alignment should be followed by an empty line
    
    # Process data
    chomp $info_line;
    my ($nr, $id1, $start1, $end1, $id2, $start2, $end2, $strand, $score) =
	split (/\s+/, $info_line);
    chomp $seqstr1;
    chomp $seqstr2;
    if($strand eq '-') {
	$strand = -1;
	$id2 = $id2.'_rev';
    }
    else {
	$strand = 1;
    }

    # Construct Bio::SimpleAlign
    my $seq1 = Bio::LocatableSeq->new(-id => $id1,
				      -start => $start1,
				      -end => $end1,
				      -strand => 1,
				      -seq => $seqstr1);
    my $seq2 = Bio::LocatableSeq->new(-id => $id2,
				      -start => $start2,
				      -end => $end2,
				      -strand => $strand,
				      -seq => $seqstr2);
    my $aln = Bio::SimpleAlign->new();
    $aln->score($score);
    $aln->add_seq($seq1);
    $aln->add_seq($seq2);

    # Return alignment
    return $aln;
}



1;
