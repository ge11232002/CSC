# AT::Prediction::ExonHolder module
#
# Copyright Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::Prediction::ExonHolder - abstract class for a genomic prediction
that that can be annotated with exons and introns

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::Prediction::ExonHolder;

use strict;
use vars '@ISA';
use AT::Prediction::Genomic;
use Carp;

@ISA = qw(AT::Prediction::Genomic);


sub new
{
    my ($caller, %args) = @_;

    my $self = bless {
	_exonL => [],
	_intronL => [],
	_undef_partL => [],
	_minor_featureL => []
    }, ref $caller || $caller;

    $self->{_feature_lists} = [$self->{_exonL},
			       $self->{_intronL},
			       $self->{_minor_featureL}];

    return $self;
}

=head2 start, end, strand

 Title     : start, end, strand
 Usage     : print  "Prediction on ", $pred->chr,
		    " strand ", $pred->strand,
		    " from ", $pred->start,
		    " to ", $pred->end,
		    "\n";
 Function  : Get absolute location of prediction.
	     To avoid confusion when in the same context
	     handling features like exons and introns,
	     it may be better to use the methods abs_start,
	     abs_end and abs_strand. (The problem is that, for
	     features, the values of start, end and strand will
	     be relative to a parent sequence, e.g. the gene
	     sequence for introns).
 Returns   : Scalars.
	     strand returns -1 or 1.
 Args      : -

=cut

sub start  { $_[0]->seq->start; }
sub end    { $_[0]->seq->end; }
sub strand { $_[0]->seq->strand; }


=head2 set_exons_introns

 Title     : set_exons_introns
 Usage     : $pred->set_exons_introns(\@exons, \@introns);
 Function  : Set exon and intron annotation
 Returns   : 1
 Args      : References to arrays of exon and intron
	     objects, as exemplified.

=cut


sub set_exons_introns
{
    my ($self, $exons, $introns, $undef) = @_;
    $self->{_exonL} = [ sort {$a->start*$a->strand <=> $b->start*$b->strand}
			@$exons ];
    $self->{_intronL} = [ sort {$a->start*$a->strand <=> $b->start*$b->strand}
			@$introns ];
    $self->{_undef_partL} = [ sort {$a->start*$a->strand <=> $b->start*$b->strand}
			@$undef ] if($undef);
    return 1;
}


=head2 exon_list, intron_list

 Title     : exon_list, intron_list
 Usage     : my @exons = $pred->exon_list();
	     my @introns = $pred->intron_list();
 Function  : 
 Returns   : Array of exon and intron objects, respectively.
 Args      : -

=cut


=head2 exon, intron

 Title     : exon, intron
 Usage     : my $exon2 = $pred->exon(2);
	     my $intron1 = $pred->intron(1);
 Function  : Gets a specific exon or intron.
 Returns   : An exon or intron object.
 Args      : A number between 1 and the number of exons/introns.

=cut


sub exon
{
    my ($self, $nr) = @_;
    my $exon = ($self->exon_list)[$nr-1];
    return $exon;
}


sub intron
{
    my ($self, $nr) = @_;
    my $intron = ($self->intron_list)[$nr-1];
    return $intron;
}

=head2 nr_exons, nr_introns

 Title     : nr_exons, nr_introns
 Usage     : print  "The gene has ", $gene->nr_exons, " exons ",
		    " and ", $gene->nr_introns, " introns.\n";
 Function  : 
 Returns   : Scalar
 Args      : -

=cut

sub nr_exons { scalar($_[0]->exon_list); }

sub nr_introns { scalar($_[0]->intron_list); }


sub gene_structure_string
{
    my ($self) = @_;

    my @se;
  
    for(my $i = 1; $i <= $self->nr_exons; $i++) {
	if(my $exon = $self->exon($i)) {
	    push @se, $exon->abs_start.'-'.$exon->abs_end;
	}
    }

    return join(',',@se);
}


sub print_exons_introns
{
    my ($self, $stream) = @_;
    $stream = \*STDOUT unless (defined $stream);
    
    print $stream "Gene structure for ", $self->loc_str, "\n";
    for(my $i = 1; $i <= $self->nr_exons or $i <= $self->nr_introns; $i++) {
	if(my $exon = $self->exon($i)) {
	    my ($s, $e) = $self->abs2rel($exon->abs_start, $exon->abs_end);
	    print $stream "E ",$s, "-", $e, "\t";
	    print $stream $exon->abs_start, "-",$exon->abs_end, "\t";
	    print $stream "(", $exon->abs_end-$exon->abs_start+1, ")";
	}
	if(my $intron = $self->intron($i)) {
	    my ($s, $e) = $self->abs2rel($intron->abs_start, $intron->abs_end);
	    print $stream "\tI ",$s, "-",$e, "\t";
	    print $stream $intron->abs_start, "-",$intron->abs_end, "\t";
	    print $stream "(", $intron->abs_end-$intron->abs_start+1, ")";
	}
	print $stream "\n";
    }
}


# Return exons in sorted order
#sub exons
#{
#    my ($self) = @_;
#    my @exons;
#    foreach my $f ($self->get_seqFeatures) {
#	push @exons, $f if ($f->primary_tag eq 'exon');
#    }
#    @exons = sort {$a->start*$a->strand <=> $b->start*$b->strand} @exons;
#    return @exons;
#}


=head2 spliced_seq

 Title     : spliced_seq
 Usage     : my $seq = $gene->spliced_seq();
 Function  : Returns a concatenation of the exon sequences
	     * THIS METHOD HAS NOT BEEN TESTED *
 Returns   : Bio::Seq
 Args      : -

=cut


sub spliced_seq
{
    my ($self) = @_;

    unless($self->_spliced_seq) {
	my $spliced_seqstr;
	foreach my $e ($self->exon_list) {
	    $spliced_seqstr .= $e->seq->seq;
	}
	$self->_spliced_seq = Bio::Seq->new
	    ( -display_id => $self->seq->display_name,
	      -seq => $spliced_seqstr);
    }
    
    return $self->_spliced_seq;
}


=head2 mapping_alignment

 Title     : mapping_alignment
 Usage     : my @exons = $pred->exon_list();
	     my @introns = $pred->intron_list();
 Function  : Reconstruct cDNA-genomic alignments for the
	     mappings supporting this prediction.
	     * THIS METHOD HAS NOT BEEN TESTED *
 Returns   : An AT::MapAlignment object
 Args      : -

=cut


sub mapping_alignment
{
    my ($self) = @_;
    
    unless($self->_mapping_aln) {
	my $aln_factory = AT::AlignmentFactory->new;
	$aln_factory->run( mappings => [ $self->mappings ],
			   genomic_seq => $self->seq);
	($self->_mapping_aln) = $aln_factory->get_alignments;
    }

    return $self->_mapping_aln;
}


sub query_seqs
{
    # gather cDNA/EST seqs from $self->mappings; use $self->_nr to make nr
}


# Feature handling methods (need to be updated!)

sub get_SeqFeatures  {
    my ($self) = @_;

    return @{$self->{_features}};
}


sub add_SeqFeature {
    my ($self,@feat) = @_;

    foreach my $feat (@feat) {

	if ( !$feat->isa('Bio::SeqFeatureI') ) {
	    warn("$feat does not implement Bio::SeqFeatureI.");
	    return undef;
	}

	if ( !$self->contains($feat) ) {
	    warn("$feat is not contained within parent feature.");
	    return undef;
	}

	$feat->attach_seq($self->seq);
	push(@{$self->{_features}},$feat);
    }

    return 1;
}


sub remove_SeqFeatures {
   my ($self) = @_;

   my @feats = @{$self->{_features}};
   $self->{_features} = []; # zap the array implicitly.
   return @feats;
}
