# AT::FT::Factory::MappingMerger module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME

AT::FT::Factory::MappingMerger module

=head1 SYNOPSIS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::FT::Factory::MappingMerger;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;
use AT::FT::GFMappingHSP;

use Class::Struct '_mapping_element' => [
    start => '$',
    end => '$',
    start_type => '$',
    end_type => '$',
];


@ISA = qw/AT::Root/;


=head2 new

 Title     : new
 Usage     : 
 Function  : Constructor
 Returns   : AT::Prediction::NaiveMapPartitioner
 Args      :
 
=cut


sub new
{
    my ($caller, %args) = @_;

    my $self = bless {
	verbose => ($args{'verbose'} || 0)
    }, ref $caller || $caller;
    
    return $self;
}


sub run
{
    my ($self, %args) = @_;
    my $in_mappings = $args{mappings} || croak "No mappings arg";
    my @out_mappings;
    my $groups = $self->_group_by_clone($in_mappings);
    foreach my $group (@$groups) {
   	my @merged = $self->_merge_mappings(@$group);
	if(@merged) {
	    push @out_mappings, @merged;
	}
	else {
	    push @out_mappings, @$group;
	}
    }
    return \@out_mappings;
}


sub _group_by_clone
{
    my ($self, $mappings) = @_;
    my %clones;
    my @ungrouped;
    foreach my $m (@$mappings) {
	#print STDERR join("\t", ($m->qName, $m->qType,
	#			 $m->library, $m->mrnaClone)), "\n";
	my $lib = $m->library;
	my $clone = $m->mrnaClone;
	if($m->refSeq_acc_list == 0 and $lib and $clone and
	   $clone ne 'n/a' and $clone ne 'Not applicable' and $clone ne 'IMAGE:' and $clone ne 'CA') {
	    # we don't merge refseqs with other transcripts since they can differ
	    push @{$clones{$lib.'__'.$clone}}, $m;
	}
	else {
	    push @ungrouped, [$m];
	}
    }
    return [values(%clones), @ungrouped];
}


sub _merge_mappings
{
    my ($self, @maps) = @_;
    return unless(@maps);
    return if($self->_duplicate_acc_check(@maps));
    @maps = sort {$a->start <=> $b->start} @maps;
    if(@maps > 1) {
	$self->_debug_msg("-- going to merge mappings --\n");
	map { $_->print_HSPs } @maps if($self->debug);
    }
    my @all_merged;
    my $merged = shift @maps;
    foreach my $m (@maps) {
	if($m->start <= $merged->end + 100000) {
	    $merged = $merged->merge($m);
	    if(ref($merged)) {
		$self->_debug_msg("OK\n");
		$merged->print_HSPs if($self->debug);
	    }
	    else {
		$self->_debug_msg($merged, "\n");
		return;
	    }
	}
	else {
	    push @all_merged, $merged;
	    $merged = $m;
	}
    }
    push @all_merged, $merged;
    return @all_merged;
}


sub _merge_mappings__old
{
    my ($self, @maps) = @_;
    return () unless(@maps);
    return undef if($self->_duplicate_acc_check(@maps));
    my $cons;
    @maps = sort {$a->tStart <=> $b->tStart} @maps;
    if(@maps > 1) {
	#print STDERR "-- mappings2clonemap --\n";
	#print STDERR $maps[0]->qName, "\t", join(' ', map {$_->start.'-'.$_->end} @$cons), "\n";
    }
    foreach my $i (1..@maps-1) {
	my $new = $self->_mapping_structure($maps[$i]);  # one mapping structure
	print STDERR $maps[$i]->qName, "\t", join(' ', map {$_->start.'-'.$_->end} @$new), "\n";
	my $err = $self->_merge_mapping_structures($cons, $new);
	if($err) {
	    #print STDERR $err, "\t", join ("\t", map {$_->qName} @maps), "\n";
	    return undef;
	}
	else {
	    #print STDERR "OK\n";
	}
    }
    #print STDERR "consensus\t", join(' ',
	#    map {($_->start_type eq 'break' ? '*':'').
  	#          $_->start.'-'.$_->end.
	#	  ($_->end_type eq 'break' ? '*':'')} @$cons), "\n"
	#if(@maps > 1);
    my $cm = AT::FT::CloneMapping->new(mappings => \@maps);
    foreach my $mapping_element (@$cons) {
	$cm->add_part(start => $mapping_element->start,
		      end => $mapping_element->end,
		      start_type => $mapping_element->start_type,
		      end_type => $mapping_element->end_type);
    }
    return $cm;
}


# check whether an array of mappings contains several mappings for some
# query sequence
sub _duplicate_acc_check
{
    my ($self, @mappings) = @_;
    my %acc;
    foreach my $raw_mapping (map { $_->primary_mapping_list} @mappings) {
	return 1 if ($acc{$raw_mapping->qName});
	$acc{$raw_mapping->qName} = 1;
    }
}


# merge b into a
sub _merge_mapping_structures
{
    my ($self, $a, $b) = @_;

    my $glitch = 5;
    my($i, $j) = (0, 0);
  
    # Find the first element in a that does not end before b's first element
    while($i < @$a and $a->[$i]->end < $b->[0]->start) { $i++; }  
    
    if($i == @$a) {
	# No overlap between $a and $b. Mark the break in the mapping
	push @$a, @$b;
	$a->[$i-1]->end_type('break');
	$a->[$i]->start_type('break');
    }
    else {
	# Some overlap between a and b
	# From this first overlap, all internal borders should agree
	for(; $i < @$a and $j < @$b; $i++, $j++) {

	    my ($a_start, $a_end) = ($a->[$i]->start, $a->[$i]->end);
	    my ($b_start, $b_end) = ($b->[$j]->start, $b->[$j]->end);

	    if($i and $a_start > $b_start + $glitch) {
		return "Err 1";
	    }
	    if($j and $a_start != $b_start) {
		if($b_start > $a_start + $glitch) {
		    return "Err 2";
		}
		# pick the most reliable of $a_start and $b_start
		$a->[$i]->start($b_start)
		    if ($self->_gap_score(@$a[$i-1,$i]) <
			$self->_gap_score(@$b[$j-1,$j]));
	    }
    
	    if($i != @$a-1 and $a_end < $b_end - $glitch) {
		return "Err 3";
	    }
	    if($j != @$b-1 and $a_end != $b_end) {
		if($b_end < $a_end - $glitch) {
		    return "Err 4";
		}
		# pick the most reliable of $a_end and $b_end
		$a->[$i]->end($b_end)
		    if($i == @$a-1 or
		       $self->_gap_score(@$a[$i,$i+1]) <
		       $self->_gap_score(@$b[$j,$j+1]));
	    }
	}

	if($i and $j and $i == @$a) {
	    $i--; $j--;
	    $a->[$i]->end($b->[$j]->end) if($a->[$i]->end < $b->[$j]->end);
	    push @$a, @$b[$j+1..@$b-1];
	}
    }

    return 0;
}


sub _gap_score
{
    my ($self, $match1, $match2) = @_;
    my $ceiling = 10;
    my $a = $match1->end - $match1->start + 1;
    $a = $ceiling if ($a > $ceiling);
    my $b = $match2->end - $match2->start + 1;
    $b = $ceiling if ($b > $ceiling);
    return $a * $b;
}


sub _mapping_structure
{
    my ($self, $mapping) = @_;
    my $min_gap = 12;
    my @structure;
    my @hsps = $mapping->all_HSPs;

    my $i = 0;
    while ($i < @hsps) {
	my $j = $i+1;
	while($j < @hsps and $hsps[$j]->tStart <= $hsps[$j-1]->tEnd + $min_gap) {
	    $j++;
	}
	# make 'match' from HSPs $i..$j-1
	push @structure, _mapping_element->new
	    ( start => $hsps[$i]->tStart,
	      end => $hsps[$j-1]->tEnd,
	      start_type => '',
	      end_type => '',
	      );
	$i = $j;
    }
    
    return \@structure;
}


1;
