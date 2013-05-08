# AT::MapAlignment module
#
# Copyright Boris Lenhard
# You may distrubute this code under the same terms as perl itself
#
# POD documentation:

=head1 NAME

AT::MapAlignment - subclass of Bio::SimpleAlign with added methods for
intron and exon retrieval

=head2 SYNOPSIS

 # Create an alignment factory and use it to make an alignment object
 # from a mapping object

 my $af = AT::MapAlignmentFactory->new(target_db => $target_db);
 $af->run( mapping => $mapping );
 my ($aln) = $af->get_alignments();

 # Output the alignment

 my $outstream = Bio::AlignIO->new ( -fh => \*STDOUT,
		   		     -format => "clustalw" );
 $outstream->write_aln ($aln);

 # Print the locations of predicted exons

 my $tSeq = $aln->get_seq_by_pos(1);
 my $tLength = $tSeq->end - $tSeq->start + 1;
 for(my $i = 1; $i <= $aln->nr_exons; $i++) {
    print "Exon $i\t";
    my $te = $aln->target_exon($i);
    my $qe = $aln->query_exon($i);

    # Print absolute genomic coords
    my ($tAbsStart) = $te->each_tag_value('abs_start');
    my ($tAbsEnd) = $te->each_tag_value('abs_end');
    print $tAbsStart, '-', $tAbsEnd, "\t";

    # Print relative coords in the attached genomic sequence 
    my ($tStart, $tEnd);
    if($te->strand == -1) {
	$tStart = $tLength - $te->end + 1;
	$tEnd = $tLength - $te->start + 1;
    }
    else {
	$tStart = $te->start;
	$tEnd = $te->end;
    }
    print $tStart, '-', $tEnd, "\t";

    # Print coords in the attached query sequence
    print $qe->start, '-', $qe->end, "\n";
 }

=head2 DESCRIPTION

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


package AT::MapAlignment;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;
use Bio::SimpleAlign;

@ISA = qw(Bio::SimpleAlign AT::Root);

sub new  {
    my ($caller, %args) = @_;
    my $self = $caller->SUPER::new(%args);
    $self->_init(%args);
    return $self;
}

sub _init {
    my ($self, %args) = @_;
    $self->{_query_exonL}      = ($args{'query_exonL' }     or []);
    $self->{_target_exonL}     = ($args{'target_exonL'}     or []);
    $self->{_intronL}          = ($args{'intronL' }         or []);
    $self->{_target_mismatchL} = ($args{'target_mismatchL'} or []);
    $self->{_query_mismatchL}  = ($args{'query_mismatchL'}  or []);
    $self->{_target_gapL}      = ($args{'target_gapL'}      or []);
    
}

# additional nifty features go here

sub rel_to_abs_target {
    my ($self, $rel_start, $rel_end) = @_;
    my ($abs_start, $abs_end);
    my $tSeqObj =($self->each_seq)[0];
    my ($strand) = (($self->each_seq)[0]->strand 
		    or ($self->target_exon_list)
		    ? ($self->target_exon_list)[0]->strand : 1);

    # start

    if ($rel_start =~ /^[\+\-]{0,1}\d+/)  {
	$rel_start =~ s/\+//;
	$abs_start = _feature_start( $self->target_exon(1), $strand, $rel_start);
    }
    elsif ($rel_start =~ /^e(\d+)$/)  { #fixme
	$abs_start = _feature_start( $self->target_exon($1)
				    or $self->target_exon($self->nr_exons) , 
				     $strand);
	
    }
    elsif ($rel_start =~ /^i(\d+)$/)  {
	$abs_start = _feature_start(($self->intron($1)  
				     or $self->intron($self->nr_introns) 
				     or $self->target_exon(1)) , $strand);
    }
    elsif ($rel_start =~ /^(\d+)$/)  {
	$abs_start = $rel_start;
    }
    else  {
	croak "Illegal syntax of the first argument of rel_to_abs_target ";
    }

    # end 
    if ($rel_end =~ /^[\+\-]{0,1}\d+/)  {
	$rel_end =~ s/\+//;
	$abs_end = _feature_end( $self->target_exon($self->nr_exons), 
				 $strand, $rel_end);
    }
    elsif ($rel_end =~ /^e(\d+)$/)  {
	$abs_end = _feature_end( ($self->target_exon($1) or
				  $self->target_exon($self->nr_exons)), 
				 $strand);
	
    }
    elsif ($rel_end =~ /^i(\d+)$/)  {
	$abs_start = _feature_end( ($self->intron($1)
				   or $self->intron($self->nr_introns) 
				   or $self->target_exon(1)), $strand);
    }
    elsif ($rel_end =~ /^(\d+)$/)  {
	$abs_end = $rel_end;
    }
    else  {
	croak "Illegal syntax of the second argument of rel_to_abs_target ";
    }
    
    return sort {$a <=> $b} ($abs_start, $abs_end);
    
}


sub _feature_start {

    # a utility function
    my ($feature, $strand, $offset) = (@_,0);
    my ($loc) = $feature->each_tag_value("abs_location");
    my $start = $loc->start;
    my $end = $loc->end;
#    my ($start) = ($feature->each_tag_value("abs_start"));
#    my ($end) = ($feature->each_tag_value("abs_end"));
   return $strand == 1 ?
	$start+$offset
	: $end-$offset;
}


sub _feature_end {

    # a utility function
    my ($feature, $strand, $offset) = (@_,0);
    my ($loc) = $feature->each_tag_value("abs_location");
    my $start = $loc->start;
    my $end = $loc->end;
#    my ($start) = ($feature->each_tag_value("abs_start"));
#    my ($end) = ($feature->each_tag_value("abs_end"));

    return $strand == 1 ?
	$end+$offset
	: $start-$offset;
}



sub get_all_target_features {
    my ($self) = @_;
    return sort {$a->start*$a->strand <=> $b->start*$b->strand}
    ($self->target_exon_list,
     $self->intron_list,
     $self->target_mismatch_list);
}


sub get_all_query_features {
    my ($self) = @_;
    return sort {$a->start <=> $b->start}
    ($self->query_exon_list,
     $self->target_gap_list,
     $self->query_mismatch_list);
}


sub query_exon {
    my ($self, $exon_nr) = @_;
    croak "Illegal exon number" unless $exon_nr;
    return ($self->query_exon_list)[$exon_nr-1];
}


sub target_exon {
    my ($self, $exon_nr) = @_;
    croak "Illegal exon number" unless $exon_nr;
    my @exons = 
	sort {$a->start*$a->strand <=> $b->start*$b->strand} 
    ($self->target_exon_list);
    return $exons[$exon_nr-1];

}
    

sub intron {
    my ($self, $intron_nr) = @_;
    croak "Illegal intron number" unless $intron_nr;
    my @introns = 
	sort {$a->start*$a->strand <=> $b->start*$b->strand} 
    ($self->intron_list);
    return $introns[$intron_nr-1];

}
   
sub nr_exons  {
    return scalar @{$_[0]->{'_query_exonL'}};
}


sub nr_introns  {
    return scalar @{$_[0]->{'_intronL'}};
}


1;
