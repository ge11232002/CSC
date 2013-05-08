# AT::FT::PEPGraph module
#
# Copyright Par Engstrom and Boris Lenhard
# 
# You may distribute this module under the same terms as perl itself
#

# POD

=head1 NAME


=head1 SYNOPSIS



=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package AT::FT::PEPGraph;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;
use AT::FT::GFMapping;


#use Class::Struct '_pep_string' => [
#    left_pep => '_pep',
#    right_pep => '_pep',
#    score => '$'
#];
#
#use Class::Struct '_pep_support' => [
#    mapping => 'AT::FT::GFMapping',
#    HSP_nr => '$'
#];
#
#use Class::Struct '_pep' => [
#    start => '$',
#    end => '$',
#    score => '$',
#    prev => '_pep',
#    next => '_pep',
#    left_strings => '%',
#    right_strings => '%',
#    open_support => '@',
#    close_support => '@',
#    _visited => '$'
#];


@ISA = qw/AT::Root/;


sub new
{
    my ($caller, %args) = @_;

    my $self = bless {
	head => ($args{head} || croak "No PEPGraph head")
    }, ref $caller || $caller;
   
    return $self;
}


sub visit
{
    my ($self, %args) = @_;
    my @visited;
    my $pep_score_thr = $args{min_pep_score} || 0;
    my $string_score_thr = $args{min_string_score} || 0;
    $self->_unset_visited_flags();
    for(my $pep = $self->head->next; $pep; $pep = $pep->next) {
        next if($pep->visited);
        push @visited, $self->_visit_connected_peps($pep,
						    $pep_score_thr,
						    $string_score_thr,
						    2);
    }
    return \@visited;
}


sub _visit_connected_peps
{
    my ($self, $pep, $pep_score_thr, $string_score_thr,
	$rel_incl_factor) = @_;

    return [] if($pep->score < $pep_score_thr);
    $pep->visited(1);

    #print STDERR "---\n";
    my @q = ($pep);
    for(my $i = 0; $i < @q; $i++) {
	my $pep = $q[$i];
        my $my_nr_maps = $pep->nr_mappings;
	#print STDERR "pep: ",$pep->start,"-",$pep->end," ",$pep->score," ", $my_nr_maps, "\n";
	my $next = $pep->next;
	if($next and !($next->visited) and $next->start == $pep->end+1) {
	    #print STDERR " next: ", $next->start, " ", ($next->visited || 0), "\n";
	    if(my $next_nr_maps = $next->nr_mappings) {
		my $nr_maps_ratio = $my_nr_maps > $next_nr_maps ? $my_nr_maps/$next_nr_maps : $next_nr_maps/$my_nr_maps;
		if($next->score >= $pep_score_thr or $nr_maps_ratio <= $rel_incl_factor) {
		    $next->visited(1);
		    push @q, $next;
		    #print STDERR "    followed next\n";
		}
	    }
	}
	my $prev = $pep->prev;
	if($prev and !($prev->visited) and $prev->end == $pep->start-1) {
	    if(my $prev_nr_maps = $prev->nr_mappings) {
		my $nr_maps_ratio = $my_nr_maps > $prev_nr_maps ? $my_nr_maps/$prev_nr_maps : $prev_nr_maps/$my_nr_maps;
		#print STDERR " prev: ", $prev->end, " ", ($prev->visited || 0), " ", $prev->score, " ", $prev->nr_mappings, " ", $nr_maps_ratio, "\n";
		if($prev->score >= $pep_score_thr or $nr_maps_ratio <= $rel_incl_factor) {
		    $prev->visited(1);
		    push @q, $prev;
		    #print STDERR "    followed prev\n";
		}
	    }
	}
	foreach my $string (values %{$pep->left_strings}) {
	    #print STDERR " left string: ",
		#$string->right_pep->start, "->", $string->left_pep->end, " ",
		#$string->score," ",
		#($string->left_pep->visited || 0),"\n";
	    my $left = $string->left_pep;
	    if(!($left->visited) and $string->score >= $string_score_thr) {
		#my $left_nr_maps = $left->nr_mappings;
		my $nr_maps_ratio = 1000; #$my_nr_maps > $left_nr_maps ? $my_nr_maps/$left_nr_maps : $left_nr_maps/$my_nr_maps;
		#print STDERR " left: ", $left->end, " ", ($left->visited || 0), " ", $left->score, " ", $left->nr_mappings, " ", $nr_maps_ratio, "\n";
		if($left->score >= $pep_score_thr or $nr_maps_ratio <= $rel_incl_factor) {
		    $string->left_pep->visited(1);
		    push @q, $string->left_pep;
		    #print STDERR "    followed string\n";
		}
	    }
	}
	foreach my $string (values %{$pep->right_strings}) {
	    #print STDERR " right string: ",
		#$string->left_pep->end, "->", $string->right_pep->start, " ",
		#$string->score, " ",
		#($string->right_pep->visited || 0),"\n";
	    my $right = $string->right_pep;
	    if(!($right->visited) and $string->score >= $string_score_thr) {
		#my $right_nr_maps = $right->nr_mappings;
		my $nr_maps_ratio = 1000; #$my_nr_maps > $right_nr_maps ? $my_nr_maps/$right_nr_maps : $right_nr_maps/$my_nr_maps;
		if($right->score >= $pep_score_thr or $nr_maps_ratio <= $rel_incl_factor) {
		    $string->right_pep->visited(1);
		    push @q, $string->right_pep;
		    #print STDERR "    followed string\n";
		}
	    }
	}

    }

    #print STDERR "---\n";

    return \@q;
}


sub _unset_visited_flags
{
    my ($self) = @_;
    for(my $pep = $self->head->next; $pep; $pep = $pep->next) {
        $pep->visited(0);
    }
}


sub prune
{
    my ($self, %args) = @_;
    my $min_pep_score = $args{min_pep_score} || 0;
    my $min_string_score = $args{min_string_score} || 0;
    for(my $pep = $self->head->next; $pep; $pep = $pep->next) {
	if($pep->score < $min_pep_score) {
	    $pep->delete;
	}
	else {
	    foreach my $string (values %{$pep->left_strings},
				values %{$pep->right_strings}) {
		if($string->score < $min_string_score) {
		    $string->delete;
		}
	    }
	}
    }
}


sub delete_all_peps
{
    my ($self) = @_;
    for(my $pep = $self->head->next; $pep; $pep = $pep->next) {
	$pep->prev(0);
	$pep->left_strings(0);
	$pep->right_strings(0);
    }
    #$self->head->next(0);
}


sub DESTROY
{
    my ($self) = @_;
    $self->delete_all_peps;
}



##
## set pep->visited to 1  (gray)
## foreach left pep and string
##   if visited { left_conn = 1; }
##   elsif(my @p = visit()) { left_conn = 1; push @conn, @p }
## unless(left_conn or left_mapping_end) {
##   purge this pep
##   return empty
## }
## foreach right pep and string
## ...(analogous)...
## set pep->visited to 2 (black)
## return (self, @conn)
#
## relies on previous purging (we could do it here, but it's cleaner
## without
## also easier to implement relative scores in purging (?)
#
## PROBLEM: this approach will remove self-suficcient subparts
## can such exist? if found, they must be left unvisited so that they
## can later be visited directly from the outer loop
## maybe better to develop a separate algorithm for the end-trimming...?
#
#sub _visit2
#{
#    my ($self, $pep, $pep_score_thr, $string_score_thr) = @_;
#
#    $pep->visited(1); # gray
#    #return [] if($pep->score < $pep_score_thr); # if keep this, we need to purge
#
#    my @connections;
#
#    my $has_left_conn;
#    foreach my $n ($pep->next, map {$_->left_pep} values %{$pep->left_strings}) {
#        if($n and $n->start == $pep->end+1) {
#	    if($n->visited) {
#	        $has_left_conn = 1;
#	    }
#	    elsif(my @left_conn = $self->_visit2($n)) {
#	        $has_left_conn = 1;
#	        push @connections, @left_conn;
#	    }
#	}
#    }
#
#    #unless $has_left_conn or...
#
#    my $has_right_conn;
#    foreach my $n ($pep->prev, map {$_->right_pep} values %{$pep->right_strings}) {
#        if($n and $n->start == $pep->end+1) {
#	    if($n->visited) {
#	        $has_right_conn = 1;
#	    }
#	    elsif(my @left_conn = $self->_visit2($n)) {
#	        $has_right_conn = 1;
#	        push @connections, @left_conn;
#	    }
#	}
#    }
#
#	my $prev = $pep->prev;
#	if($prev and !($prev->visited) and $prev->end == $pep->start-1) {
#	    #print STDERR " prev: ", $prev->end, " ", ($prev->visited || 0), "\n";
#	    $prev->visited(1);
#	    if($prev->score >= $pep_score_thr) {
#		push @q, $prev;
#	    }
#	}
#
#
#	foreach my $string (values %{$pep->left_strings}) {
#	    #print STDERR " left string: ",
#		#$string->right_pep->start, "->", $string->left_pep->end, " ",
#		#$string->score," ",
#		#($string->left_pep->visited || 0),"\n";
#	    if(!($string->left_pep->visited) and
#	       $string->score >= $string_score_thr) {
#		$string->left_pep->visited(1);
#		if($string->left_pep->score >= $pep_score_thr) {
#		    push @q, $string->left_pep;
#		    #print STDERR "    followed string\n";
#		}
#	    }
#	}
#	foreach my $string (values %{$pep->right_strings}) {
#	    #print STDERR " right string: ",
#		#$string->left_pep->end, "->", $string->right_pep->start, " ",
#		#$string->score, " ",
#		#($string->right_pep->visited || 0),"\n";
#	    if(!($string->right_pep->visited) and
#	       $string->score >= $string_score_thr) {
#		$string->right_pep->visited(1);
#		if($string->right_pep->score >= $pep_score_thr) {
#		    push @q, $string->right_pep;
#		    #print STDERR "    followed string\n";
#		}
#	    }
#	}
#
#    }
#
#    #print STDERR "---\n";
#
#    return \@q;
#}



1;
