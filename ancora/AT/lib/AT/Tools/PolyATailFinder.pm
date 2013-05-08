package AT::Tools::PolyATailFinder;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;

@ISA = qw(AT::Root);


sub _get_args {
    my ($self, %args1) = @_;
    my %args2 = (
	win_size => ($args1{win_size} || 15),
	min_content => ($args1{min_content} || 10),
	min_consecutive => ($args1{min_consecutive} || 6),
    );
	#^^ could make these into attributes if obj is to be instantiated
    return %args2;
}


sub find_tail_in_seq {
    my $T_tail_len = find_T_tail_in_seq(@_);
    my $A_tail_len = find_A_tail_in_seq(@_);
    return($T_tail_len, $A_tail_len);
}


sub find_T_tail_in_seq {
    my ($self, $seq, %args) = @_;
     %args = $self->_get_args(%args);
    return $self->_search_tail($seq->seq, 1, 'T', %args);
}


sub find_A_tail_in_seq {
    my ($self, $seq, %args) = @_;
     %args = $self->_get_args(%args);
    return $self->_search_tail($seq->seq, 0, 'A', %args);
}


sub find_tail_in_mapping {
    my ($self, $mapping, $seq_db) = @_;
    
    my $gap_allowed = 2;
    my @hsps = $mapping->all_HSPs;

    # Look for T
    my $i = 0;
    while($i < $#hsps and
	  $hsps[$i+1]->tStart <= $hsps[$i]->tEnd + $gap_allowed and
	  $hsps[$i+1]->qStart <= $hsps[$i]->qEnd + $gap_allowed)
    {
	$i++;
    }    
    my $beg_seq = $seq_db->get_genome_seq_str(chr => $mapping->tName,
					  start => $mapping->tStart,
					  end => $hsps[$i]->tEnd);
    my $T_tail_len = $self->_search_tail($beg_seq, 1, 'T');
    
    # Look for A
    $i = $#hsps;
    while($i > 0 and
	  $hsps[$i-1]->tEnd >= $hsps[$i]->tStart - $gap_allowed and
	  $hsps[$i-1]->qEnd >= $hsps[$i]->qStart - $gap_allowed)
    {
	$i++;
    }    
    my $end_seq = $seq_db->get_genome_seq_str(chr => $mapping->tName,
					  start => $mapping->tStart,
					  end => $hsps[$i]->tEnd);
    my $A_tail_len = $self->_search_tail($end_seq, 0, 'A');

    return($T_tail_len, $A_tail_len);

}


sub tailless_mapping {
    my ($self, $mapping, $seq_db) = @_;

    my ($T_tail_len, $A_tail_len) = $self->find_polyT_polyA($mapping, $seq_db);

    # Recalculate score
    if($T_tail_len > $A_tail_len) {
	return $self->_trunc_start($mapping, $seq_db, $T_tail_len);
    }
    elsif($A_tail_len) {
	die "a tailless not implemented\n";
    }
}


sub _trunc_start {
    my ($self, $m, $seq_db, $tail_len) = @_;
        
    # find hsp where the mapping should start
    my $tNewStart = $m->tStart + $tail_len;
    my @hsps = $m->all_HSPs;
    my $i = 0;
    for(; $i < @hsps; $i++) {
	last if ($tNewStart <= $hsps[$i]->tStart or
		 $tNewStart <= $hsps[$i]->tEnd);
    }
    die "tail longer than entire mapping" if ($i == @hsps);
    
    # need to truncate the 1st hsp?
    my $qNewStart;
    if($tNewStart > $hsps[$i]->tStart) {
	$qNewStart = $hsps[$i]->qStart+ $tNewStart - $hsps[$i]->tStart;
    }
    else {
	$qNewStart = $hsps[$i]->qStart;
    }

    # get query and target seq for the region
    my $tSeq = seq_db->get_genome_seq_str(chr => $m->tName,
					  start => $m->tStart,
					  end => $tNewStart-1);
    my @tSeq = (split //, $tSeq);
    my $qSeq = $m->query_seq;
    $qSeq = $qSeq->revcom if ($m->strand eq '-');
    my @qSeq = (split //, $qSeq);

    # calc nr of matches, misMatches etc to subtract
    my ($matches, $repMatches, $misMatches, $nCount) = (0,0,0,0);
    my ($tNumInsert, $tBaseInsert, $qNumInsert, $qBaseInsert) = (0,0,0,0);
    my ($qPos, $tPos) = ($hsps[0]->qStart-1, 0);
    my $j = 0;
    while($tPos <= $tail_len) {
	my ($qNt, $tNt) = ($qSeq[$qPos], $tSeq[$tPos]);
	$qNt = lc $qNt;
	$tNt = lc $tNt;
	my $in_repeat = ($tNt eq $tSeq[$tPos]);
	if($qNt eq 'n' or $tNt eq 'n') {
	    $nCount++;
	}
	elsif($qNt eq $tNt) {
	    if($in_repeat) { $repMatches++; }
	    else { $matches++; }
	}
	else {
	    $misMatches++;
	}
	if($tPos+$m->tStart > $hsps[$j]->tEnd) {
	    $tPos += $hsps[$j+1]->tStart - $hsps[$j]->tStart;
	    $qPos += $hsps[$j+1]->qStart - $hsps[$j]->qStart;
	    $j++;
	}
    }

    # create new hsps
    my @new_hsps;
    push @new_hsps, AT::HSP->new
	    (qStart    => $qNewStart,
	     qEnd      => $hsps[$i]->qEnd,
	     tStart    => $tNewStart,
	     tEnd      => $hsps[$i]->tEnd,
	     blockSize => $hsps[$i]->tEnd - $tNewStart + 1
	     );
    for my $j ($i+1 .. $#new_hsps)
    {
	push @new_hsps, AT::HSP->new( %{$hsps[$j]} );
    }

    # create new mapping and return
    my $new_m =  AT::Mapping->new
	(matches      => $m->$matches - $matches,
	 misMatches   => $m->misMatches - $misMatches,
	 repMatches   => $m->repMatches - $repMatches,
	 nCount       => $m->nCount - $nCount,
	 qNumInsert   => $m->qNumInsert - $qNumInsert,  # Note: 
	 qBaseInsert  => $m->qBaseInsert - $qBaseInsert,  # these 4 are not
	 tNumInsert   => $m->tNumInsert - $tNumInsert,   # modified as of now
	 tBaseInsert  => $m->tBaseInsert - $tBaseInsert,
	 strand       => $m->strand,
	 qName        => $m->qName,
	 qSize        => $m->qSize,
	 qStart       => $m->qNewStart,
	 qEnd         => $m->qEnd,
	 tName        => $m->tName,
	 tSize        => $m->tSize,
	 tStart       => $m->tNewStart,
	 tEnd         => $m->tEnd,
	 blockCount   => scalar @new_hsps,
	 HSPs         => [@new_hsps]);    

    return $new_m;
}


sub _search_tail
{
    my ($self, $seq, $fwd, $nt, %args) = @_;

    my $win_size = $args{win_size};
    my $min_content = $args{min_content};
    my $min_consec = $args{min_consecutive};

    #print STDERR "searching with $win_size, $min_content, $min_consec\n";

    $nt = lc $nt;
    $seq = lc $seq;
    my @seq = (split //, $seq);
    return 0 if(@seq < 8);
    @seq = reverse @seq unless ($fwd);

    # Loop variables
    my @cumul;
    $#cumul = @seq;
    $cumul[0] = 0;  # score base
    my $i = 1;

    # Step 1 nt until we're in a window with too few target nts
    for(; $i <= $win_size and $i <= @seq; $i++) {   # Step through the first window
	my $nt_score = $seq[$i-1] eq $nt ? 1 : 0;
	$cumul[$i] = $cumul[$i-1] + $nt_score;
    }
    return 0 if($cumul[$i-1] < $min_content); # Return if too few in the first window
    for(; $i <= @seq; $i++) {		      # Otherwise keep stepping
	my $nt_score = $seq[$i-1] eq $nt ? 1 : 0;
	$cumul[$i] = $cumul[$i-1] + $nt_score;
	my $window_score = $cumul[$i] - $cumul[$i-$win_size];
	last if($window_score < $min_content);
    }

    # Backtrack to find first $min_consec consecutive target nt's
    $i--;  # step back 1 so we're inside the last window with >= $min_content target nts
    while($cumul[$i] - $cumul[$i-$min_consec] != $min_consec) {
	$i--;
        return 0 if($i < $min_consec);
    }

    return $i;  # return number of bases that the tail covers
}
