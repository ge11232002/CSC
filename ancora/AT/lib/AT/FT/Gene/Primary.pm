# AT::FT::Gene::Primary module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::FT::Gene::Primary

=head1 SYNOPSIS



=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::FT::Gene::Primary;

use strict;
use vars '@ISA';
use AT::Prediction::ExonHolder;

@ISA = qw/AT::Prediction::ExonHolder/;


sub new
{
    my ($caller, %args) = @_;
    
    my ($self) = $caller->SUPER::new(%args);

    my $seq = $args{seq};
    unless ($seq) {   # really require seq?
	die "No sequence provided. Cannot make gene object.";
    }

    $self->{seq} = $seq;
    $self->{_TSSL} = $args{TSS_list} || [];
    $self->{_PASL} = $args{PAS_list} || [];
    $self->{_gfmappingL} = $args{gfmapping_list} || [];
    $self->{is_spliced} = $args{is_spliced} || 0;
    #$self->{_unique_mappingL} = $args{mapping_list} || $args{unique_mapping_list} || [];
    #$self->{_nonunique_mappingL} = $args{nonunique_mapping_list} || [];
    #$self->{_excluded_mappingL} = $args{excluded_mapping_list} || [];
    #$self->{_commentL} = $args{comment_list} || [];

    return $self;
}
			 

sub nr_EST_gfmappings
{
    my ($self) = @_;
    my $nr = 0;
    foreach my $m ($self->gfmapping_list) {
	$nr++ if($m->EST_acc_list);
    }
    return $nr;
}


sub gfmapping_score_sum
{
    my ($self) = @_;
    my $score = 0;
    foreach my $m ($self->gfmapping_list) {
	$score += $m->score1;
    }
    return $score;
}


sub representative_mRNA_acc
# we only get rep mrna accs, since getting rep. est accs is messy for merged mappings
{
    my ($self) = @_;

    my (@mrna_maps, @est_maps);
    foreach my $m ($self->gfmapping_list) {
	if($m->mRNA_acc_list) {
	    push @mrna_maps, $m;
	}
	#else {
	#    push @est_maps, $m;
	#}
    }

    ## use mrnas mappings if any, otherwise use est mappings
    #my $maps = @mrna_maps ? \@mrna_maps: \@est_maps;

    # find mapping(s) that span most of the gene
    my @rep_maps; 
    my $rep_maps_span = 0;
    foreach my $m (@mrna_maps) {
	my $start = $self->start > $m->start ? $self->start : $m->start;
	my $end = $self->end < $m->end ? $self->end : $m->end;
	my $span = $end - $start;
	if($span > $rep_maps_span) {
	    @rep_maps = ($m);
	    $rep_maps_span = $span;
	}
	elsif($span == $rep_maps_span) {
	    push @rep_maps, $m;
	}
    }

    # select one of those mappings.
    # for mrnas: if there is a refseq, take that, otw take the mapping with highest score
    # for ests: take one rep by only one est; avoid broken mappings
    my ($rep_map, $rep_acc);
    if(@mrna_maps) { # mrnas
	my $hiscore = 0;
	foreach my $map (@rep_maps) {
	    my ($refseq) = $map->refSeq_acc_list;
	    if($refseq) {
		$rep_map = $map;
		$rep_acc = $refseq;
		last;
	    }
	    elsif($map->score1 > $hiscore) {
		$rep_map = $map;
		($rep_acc) = $map->mRNA_acc_list;
		$hiscore = $rep_map->score1;
	    }
	}
    }
#    else {  # ests
#	my $nr_ests = 100;
#	foreach my $map (@rep_maps) {
#	    my @accs = $map->EST_acc_list;
#	    if(@accs < $nr_ests) {
#		$rep_map = $map;
#		$nr_ests = @accs;
#		last if($nr_ests == 1);
#	    }
#	}
#    }

    return ($rep_acc);
}


sub print_borders
{
    my ($self, $stream) = @_;
    $stream = \*STDOUT unless (defined $stream);
    print $stream "Borders for ", $self->loc_str, "\n";
    foreach my $b (sort {$a->abs_start <=> $b->abs_start}
		   ($self->TSS_list, $self->PAS_list)) {
        my ($rel_start, $rel_end) = $self->abs2rel($b->abs_start, $b->abs_end);
        print $stream (ref($b) =~ /TSS/) ? 'TSS' : 'PAS';
	print $stream " ", $rel_start, "-", $rel_end, "\t";
	print $stream $b->abs_start, "-", $b->abs_end, " ", $b->abs_pos, "\n";
    }
}


=head2 want more methods?

This class inherits AT::Prediction::ExonHolder, which in turn inherits
AT::Prediction::Genomic. See those modules for more method descriptions.

=cut

sub _nr
{
    my ($self, @redundant) = @_;
    my %non_redundant;
    foreach my $item (@redundant) {
	$non_redundant{$item} = $item;
    }
    return values(%non_redundant);
}
