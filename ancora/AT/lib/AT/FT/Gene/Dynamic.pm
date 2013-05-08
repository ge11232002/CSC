# AT::FT::Gene::Dynamic module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::Prediction::Gene - represents a predicted gene

=head1 SYNOPSIS



=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::FT::Gene::Dynamic;

use strict;
use vars '@ISA';
use vars '$AUTOLOAD';
use AT::FT::Gene::Primary;
use AT::FT::PAS;
use AT::FT::TSS;
use AT::Prediction::Exon;
use AT::Prediction::Intron;
use AT::Prediction::UndefGenePart;
use Carp;

@ISA = qw/AT::Root/;

# should interit/conform to Gene::PrimaryI

=head2 new

 Title     : new
 Usage     : my $gene = AT::Prediction::Gene->new
		(seq => $seq,
		 spliceform_list => [ \@spliceforms ] );				
 Function  : Constructor.
	     This method is not intended to be used directly.
	     Use AT::Prediction::GenePredictor to create gene
	     objects.
 Returns   : An AT::Prediction::Gene object
 Args      : 

=cut

sub new
{
    my ($caller, %args) = @_;
    
    my ($self) = bless{ }, ref($caller) || $caller;

    my $seq = $args{seq};
    unless ($seq) {
	die "No sequence provided. Cannot make gene object.";
    }

    $self->{gene_id} = $args{gene_id} || 0;
    $self->{seq} = $seq;
    $self->{strand} = $args{strand} || croak "No strand arg";
    $self->{is_spliced} = $args{is_spliced}; # does not fit in a dynamic gene
    $self->{_peps} = $args{peps} || croak "No peps arg";
    $self->{_stringency} = $args{stringency} || 0;
    $self->{_gfmappingL} = $args{gfmapping_list} || [];
    #$self->{_unique_mappingL} = $args{mapping_list} || $args{unique_mapping_list} || []; # remove?
    #$self->{_nonunique_mappingL} = $args{nonunique_mapping_list} || [];  # remove?
    #$self->{_excluded_mappingL} = $args{excluded_mapping_list} || [];
    #$self->{_commentL} = $args{comment_list} || [];
    $self->{_use_score} = 1;

    $self->{_primary_gene} = undef;
    $self->{_primary_genes_by_stringency} = {};

    return $self;
}


sub pep_list { @{shift->{_peps}} }
sub gfmapping_list { @{shift->{_gfmappingL}} }

sub min_start { shift->{_peps}->[0]->start }
sub max_end { shift->{_peps}->[-1]->end }

sub gene_id {
    my ($self, $value) = @_;
    $self->{gene_id} = $value if(defined $value);
    return $self->{gene_id};
}

sub chr
{
    my ($self) = @_;
    unless ($self->{'chr'}) {
	($self->{'chr'}) = $self->{seq}->id =~ /(chr.+)/;
    }
    return $self->{'chr'};
}


sub overall_loc_str
{
    my ($self) = @_;

    my $str = $self->chr.":".$self->{seq}->start."-".$self->{seq}->end.":".
	($self->{strand} == 1 ? '+' : '-');
    return $str;
}

sub stringency
{
    my ($self, $s) = @_;
    if(defined $s) {
	$self->{_stringency} = $s;
	$self->{_primary_gene} = undef;
    }
    return $self->{_stringency};
}


sub use_score
{
    my ($self, $value) = @_;
    if(defined $value) {
        $self->{_use_score} = $value;
	$self->{_primary_gene} = undef;
    }
    return $self->{_use_score};
}


sub primary_gene {
    my ($self) = @_;
    my $primary = $self->{_primary_gene};
    unless(defined $primary) {
        my $stringency = $self->stringency();
	my $use_score = $self->use_score();
	my $key = "${stringency}_${use_score}";
        $primary = $self->{_primary_genes_by_stringency}{$key};
        unless(defined($primary)) {
	    $primary = $self->{_primary_genes_by_stringency}{$key} =
		$self->_create_primary_gene();
	}
        $self->{_primary_gene} = $primary;
    }
    return ($primary || undef);
}


# The gene will contain a set of Gene::Primary objects, one for each stringency
# This class is a Gene::Adjustable (will change this)
# Can rewrite it so that this class functions as a true relay...?
#
# Problems: what if a gene is attached as a member of an antisense pair (eg)
#  and then altered to a different stringency -- it means that the antisense pair
#  will also change. This can give strange effects if eg the antisense pair
#  caches data based on calculation on the members
#  More down-to-earth to be able to emit a number of genes at different stringencies
#  Tell genemachine what kind of gene we want to create - an adjustable or a primary
#  The adjustable will generate primary genes based on its stringency
#  A primary gene is directly accessibly through $gene->primary_gene. This is
#  good for saving in other structures (e.g. antisensepairs)
#  However an adjustable gene should also have the interface of a primary gene
#  - relaying everything to the primary gene (don't need to implement this now)

sub _create_primary_gene {
    my ($self) = @_;

    my $stringency = $self->stringency();
    my $all_peps = $self->{_peps};
    my $target_seq = $self->{seq};
    my $strand = $self->{strand};

    # select those peps that are above the set stringency
    my $score_method = "score".($self->use_score);
    #print STDERR "score_method: $score_method\n";
    my @peps = 	grep {$_->$score_method >= $stringency} @$all_peps;
    #print STDERR "nr peps: ",scalar(@peps),"\n";
    #print STDERR "score2: ",join(' ',map {$_->score2} @peps), "\n";
    return 0 unless (@peps);
    
    my $gene_start = $peps[0]->start;
    my $gene_end = $peps[-1]->end;
    my $seqstr = $target_seq->subseq($gene_start - $target_seq->start + 1,
				    $gene_end - $target_seq->start + 1);
    my $seq = Bio::LocatableSeq->new(-id => $target_seq->id,
				    -start => $gene_start,
				    -end => $gene_end,
				    -strand => $strand,
				    -seq => $seqstr);
    my $gene = AT::FT::Gene::Primary->new
	(#unique_mapping_list => $self->{_unique_mappingL},
	 #nonunique_mapping_list => $self->{_nonunique_mappingL},
	 gfmapping_list => $self->{_gfmappingL},
	 is_spliced => $self->{is_spliced},  # not correct; should be derived from gfmappings
	 seq => $seq);
	# note: now all mappings are attached to the new gene; they should be filtered

    # Create exons and introns
    my (@exons, @introns, @undef);
    my $i = 0;
    while ($i < @peps) {
	my $j = $i+1;
	while($j < @peps and $peps[$j]->start == $peps[$j-1]->end + 1) {
	    $j++;
	}
	# make exon from peps $i..$j-1
	#print STDERR "E: ",$peps[$i]->start,"-",$peps[$j-1]->end,"\n";
	push @exons, AT::Prediction::Exon->new
	    ( -start => $peps[$i]->start - $gene_start + 1,
	      -end => $peps[$j-1]->end - $gene_start + 1,
	      -strand => $strand,
	      -abs_start => $peps[$i]->start,
	      -abs_end => $peps[$j-1]->end
	      );
	if($j < @peps) {
	    if($self->_gap_looks_like_intron($peps[$j-1]->end+1, $peps[$j]->start-1)) {
		#print STDERR "I: ",$peps[$j-1]->end+1,"-",$peps[$j]->start-1,"\n";
		push @introns, AT::Prediction::Intron->new
		    ( -start => $peps[$j-1]->end - $gene_start + 2,
		    -end => $peps[$j]->start - $gene_start,
		    -strand => $strand,
		    -abs_start => $peps[$j-1]->end + 1,
		    -abs_end => $peps[$j]->start - 1
		    )
	    }
	    else {
		push @undef, AT::Prediction::UndefGenePart->new
		    (-start => $peps[$j-1]->end - $gene_start + 2,
		    -end => $peps[$j]->start - $gene_start,
		    -strand => $strand,
		    -abs_start => $peps[$j-1]->end + 1,
		    -abs_end => $peps[$j]->start - 1
		    )
	    }
	}
	$i = $j;
    }
    $gene->set_exons_introns(\@exons, \@introns, \@undef);

    my @open_pos =
	map { [$_->start, $_->open_score] }
	grep { $_->open_score >= 20 } @peps;
    my @close_pos =
	map { [$_->end, $_->close_score] }
	grep { $_->close_score >= 20 } @peps;
    my ($putative_TSS_pos, $putative_PAS_pos) =
	($strand == 1) ?
	(\@open_pos, \@close_pos) :
	(\@close_pos, \@open_pos);
    my $TSS = $self->_create_borders($putative_TSS_pos,
				     $gene_start,
				     $strand,
				     'AT::FT::TSS');
    my $PAS = $self->_create_borders($putative_PAS_pos,
				     $gene_start,
				     $strand,
				     'AT::FT::PAS');
    $gene->TSS_list($TSS);
    $gene->PAS_list($PAS);

    if($self->debug) {
        print STDERR "Created gene: ", $gene->loc_str, "\n";
        $gene->print_exons_introns(\*STDERR);
        $gene->print_borders(\*STDERR);
        print STDERR "\n";
    }
    
    return $gene;
}


sub _gap_looks_like_intron
{
    my ($self, $start, $end) = @_;
    return 0 if ($end - $start + 1 < 12);
    my $seq = $self->{seq};
    $start = $start - $seq->start + 1;
    $end = $end - $seq->start + 1;
    my $spljnc = ($seq->subseq($start,$start+1)) . ($seq->subseq($end-1,$end));
    $spljnc = lc $spljnc;
    #print STDERR $spljnc, "\n";
    if($self->{strand} == 1) {
	return 1 if($spljnc eq "gtag" or $spljnc eq "gcag" or $spljnc eq "atac")
    }
    else {
	return 1 if($spljnc eq "ctac" or $spljnc eq "ctgc" or $spljnc eq "gtat")
    }
    return 0;
}


sub _create_borders
{
    my ($self, $a, $gene_start, $strand, $border_class) = @_;
    my $max_d = 10;
    my @borders;
    for(my $i = 0; $i < @$a;) {
	my ($start, $score) = @{$a->[$i]};
	if ($score >= 25 or ($i+1 < @$a and $a->[$i+1][0] == $start+1))
	{
	    while($i+1 < @$a and $a->[$i+1][0] <= $a->[$i][0] + $max_d) {
		$i++;
		$score += $a->[$i][1];
	    }
	    my $end = $a->[$i][0];
	    my ($rel_start, $rel_end) = ($start - $gene_start + 1,
					 $end - $gene_start + 1);
	    push @borders, $border_class->new
		(-start => $rel_start,
		 -end => $rel_end,
		 -strand => $strand,
		 -abs_start => $start,
		 -abs_end => $end,
		 -score => $score);
	}
	$i++;
    }
    return \@borders;
}


sub AUTOLOAD  {
    my $self = shift;
    my ($method) = $AUTOLOAD =~ /:([^:]+)$/;
    my $primary = $self->primary_gene();
    if($primary) {
	return $primary->$method(@_);
    }
    else {
	warn "Attempt to call method $method on undefined primary gene for dynamic gene ", $self->overall_loc_str, "\n";
	return;
    }   
}


## methods in gene
#sub mapping_list { shift->{_primary_gene}->mapping_list(@_) }
#sub mRNA_mapping_list { shift->{_primary_gene}->mRNA_mapping_list(@_) }
#sub EST_mapping_list { shift->{_primary_gene}->EST_mapping_list(@_) }
#sub print_borders { shift->{_primary_gene}->print_borders(@_) }
#
## methods in exonholder
#sub start { shift->{_primary_gene}->start(@_) }
#sub end { shift->{_primary_gene}->end(@_) }
#sub strand { shift->{_primary_gene}->strand(@_) }
#sub nr_exons { shift->{_primary_gene}->nr_exons(@_) }
#sub nr_introns { shift->{_primary_gene}->nr_introns(@_) }
#sub { shift->{_primary_gene}->(@_) }
#


1;
