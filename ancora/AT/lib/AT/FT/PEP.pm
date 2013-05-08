# AT::FT::PEP module
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

package AT::FT::PEP;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;


@ISA = qw/AT::Root/;


sub new
{
    my ($caller, %args) = @_;

    my $self = bless {
	pep_id => ($args{pep_id} || 0),
	start => ($args{start} || 0),
	end => ($args{end} || 0),
	score1 => ($args{score1} || 0),
	score2 => ($args{score2} || 0),  # secondary score
	open_score => ($args{open_score} || 0),
	close_score => ($args{close_score} || 0),
	_open_qNames => ($args{_open_qNames} || []), # for dbg
	_close_qNames => ($args{_close_qNames} || []), # for dbg
	prev => ($args{prev} || undef),
	next => ($args{next} || undef),
	left_strings => ($args{left_strings} || {}),
	right_strings => ($args{right_strings} || {}),
	_open_supportL => ($args{open_support_list} || []),
	_close_supportL => ($args{close_support_list} || []),
	_mappingL => ($args{mapping_list} || []),
	#_modified => 1,
	visited => 0
    }, ref $caller || $caller;

    #print STDERR "PEP: ",$self->score1," ",$self->score2, "\n";
   
    return $self;
}

sub score { shift->score1(@_); }

sub add_open_support {
    push @{$_[0]->{_open_supportL}}, $_[1];
    #$_[0]->{_modified} = 1;
}
sub add_close_support {
    push @{$_[0]->{_close_supportL}}, $_[1];
    #$_[0]->{_modified} = 1;   
}
sub set_open_support {
    $_[0]->{_open_supportL} = $_[1];
    #$_[0]->{_modified} = 1;
}
sub set_close_support {
    $_[0]->{_close_supportL} = $_[1];
    #$_[0]->{_modified} = 1;
}
sub set_mapping_list {
    $_[0]->{_mappingL} = $_[1];
    #print STDERR "list set to: ", join(',',map{$_->qName_list}@{$_[1]}), "\n";
    #print STDERR $_[0]->nr_mappings, "\n";
}

sub nr_mappings { scalar(@{shift->{_mappingL}}) }


sub left_string_list { return values %{$_[0]->left_strings}; }
sub right_string_list { return values %{$_[0]->right_strings}; }

sub delete
{
    my ($pep) = @_;
    $pep->prev->next($pep->next);
    $pep->next->prev($pep->prev) if($pep->next);
    foreach my $string ($pep->left_string_list,
			$pep->right_string_list) {
	$string->delete();
    }
}


#sub _update
#{
#    my ($self) = @_;
#    my ($open_score, $close_score);
#    foreach my $support ($self->open_support_list) {
#	if($support->mapping->start == $self->start and
#	   $support->mapping->trust_start) {
#	    $open_score += $support->score;
#        }
#    }
#    $self->{_open_score} = $open_score;
#    foreach my $support ($self->close_support_list) {
#	if($support->mapping->end == $self->end and
#	   $support->mapping->trust_end) {
#	    $close_score += $support->score;
#        }
#    }
#    $self->{_close_score} = $close_score;  
#}
#
#
#sub open_score {
#    my ($self) = @_;
#    $self->_update() if($self->{_modified});
#    return $self->{_open_score};
#}
#
#sub close_score {
#    my ($self) = @_;
#    $self->_update() if($self->{_modified});
#    return $self->{_open_score};
#}


1;
