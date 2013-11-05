package CONSITE::Analysis;


use vars qw($AUTOLOAD);
use CGI qw(:standard);
use TFBS::Matrix::PWM;
use TFBS::Matrix::PFM;
use TFBS::Matrix;
use Bio::Root::Root;
use Bio::Seq;
use Bio::SimpleAlign;
use DBI;
use GD;
use PDL;
use IO::String;
use Class::MethodMaker
    get_set => [qw(cutoff window
		   alignseq1 alignseq2
		   fsiteset1 fsiteset2
		   seq1name seq2name
		   seq1length seq2length
		   conservation1 conservation2
		   start_at end_at)],
    list    => [qw(alignseq1array alignseq2array)],
    hash    => [qw(seq1siteset seq2siteset)],
    new_with_init => ['new'];


use strict;


use constant IMAGE_WIDTH =>600;
use constant IMAGE_MARGIN =>50;

my @seq1RGB   = (48, 176, 200);
my @seq2RGB   = (31, 225,   0);
my $seq1hex   = "#30B0C8";
my $seq2hex   = "#1FD100";
my $seq1bghex = "#CCFFFF";
my $seq2bghex = "#CCFFCC";


sub init  {

    my ($self, %args) = @_;
    $self->throw("Analysis method should not be instantiated  ".
		 "- use one of its subclasses instead"
		 );
    return undef;

}

sub _set_start_end  {

    # start and end are set relative to the alignment !!!
    my ($self, %args) = @_;

    # calculate overlap - temporary
    if (defined($args{-start}) and defined($args{-end}))  {
	no strict 'refs';
	my $rangeref = ($self->{job}->rangeref() or 1);
	my $seqstrlen ="seq".$rangeref."length";
	my $reflength = $self->$seqstrlen;
	use strict 'refs';
	print STDERR "RANGEREF REFLENGTH ".$rangeref." ".$reflength."\n";

	# flip start, end if reversed

	if ($args{-start} > $args{-end})  {
	    ($args{-start}, $args{-end}) = ($args{-end}, $args{-start});
	}

	# try to set start_at

	my $start_at =
	   ( # if legal
	    $self->pdlindex($args{-start}, $rangeref=>0)
	    # else if negative, try to slide upstream
	    or ($args{-start}<1 ?
		($self->pdlindex(1,$rangeref=>0) - $args{-start}):0)
	    #else if bigger than length, try to slide downstream
		or ($args{-start}>$reflength
		    ?
		    $self->pdlindex($reflength,$rangeref=>0)
		     +$args{-start} - $reflength
		    : 1));
	# still overflows?
	if ($args{-start} > $reflength)  {
	    $start_at = $self->pdlindex()->getdim(0)
		- IMAGE_WIDTH
		    + 2*IMAGE_MARGIN;
	}
	$start_at = 1 if $start_at < 1;

	# try to set end_at

	my $end_at =
	   ( # if legal
	    $self->pdlindex($args{-end}, $rangeref=>0)
	    # else if negative, try to slide upstream
	    or ($args{-end}<1 ?
		($self->pdlindex(1,$rangeref=>0) - $args{-end}):0)
	    #else if bigger than length, try to slide downstream
		or ($args{-end}>$reflength
		    ?
		    $self->pdlindex($reflength,$rangeref=>0)
		     +$args{-end} - $reflength
		    :
		    $self->pdlindex()->getdim(0)));
	$end_at = IMAGE_WIDTH - 2*IMAGE_MARGIN  if $end_at < $start_at;
	$end_at =  $self->pdlindex->getdim(0)
	    if $end_at >$self->pdlindex()->getdim(0);

	$self->start_at($start_at); $self->{job}->start_at($start_at);
	$self->end_at($end_at);     $self->{job}->end_at($end_at);
    }
    elsif (not defined $self->{job}->alignstring())  { # single seq analysis
	$self->start_at(1);
	$self->end_at($self->{job}->seq1->length);
    }
    else {
	my $overlap_slice1 =
	    $self->pdlindex->slice(':,1')->where(($self->pdlindex->slice(':,2')>0));
	$overlap_slice1 = $overlap_slice1->where($overlap_slice1>0);

	my ($overlap_start, $overlap_end) =
	    ($self->pdlindex(list ($overlap_slice1->slice(0)),  1=>0),
	     $self->pdlindex(list ($overlap_slice1->slice(-1)), 1=>0));

	if ($overlap_end -  $overlap_start < (IMAGE_WIDTH-2*IMAGE_MARGIN)) {
	    $overlap_start =
		int ($overlap_start + ($overlap_start - $overlap_end)/2
		     -(IMAGE_WIDTH-2*IMAGE_MARGIN)/2);
	    $overlap_end =
		int ($overlap_start + ($overlap_start - $overlap_end)/2
		     +(IMAGE_WIDTH-2*IMAGE_MARGIN)/2);
	    $overlap_start=1 if $overlap_start<1;
	    $overlap_end = $self->pdlindex->getdim(0)
		if $overlap_end > $self->pdlindex->getdim(0);

	}
	$self->start_at($overlap_start);
	$self->end_at($overlap_end);
    }
    print STDERR "START_AT END AT ".$self->start_at()." ".$self->end_at()."\n";

}



sub DESTROY {
    my $self = shift;
    delete $self->{job};
}


sub pdlindex {
    my ($self, $input, $p1, $p2) = @_ ;
    # print ("PARAMS ", join(":", @_), "\n");
    if (ref($input) eq "PDL")  {
	$self->{pdlindex} = $input;
    }

    unless (defined $p2)  {
	return $self->{pdlindex};
    }
    else {
	my @results = list
	    $self->{pdlindex}->xchg(0,1)->slice($p2)->where(
	                             $self->{pdlindex}->xchg(0,1)->slice($p1)==$input
						 );
	wantarray ? return @results : return $results[0];
    }
}

sub lower_pdlindex {
    my ($self, $input, $p1, $p2) = @_;
    unless (defined $p2)  {
	$self->throw("Wrong number of parameters passed to lower_pdlindex");
    }
    my $result;
    my $i = $input;
    #my $max = $self->pdlindex->getdim(0);
    until ($result = $self->pdlindex($i, $p1 => $p2))  {
	$i--;
	last if $i==0;
    }
    return $result;
}

sub higher_pdlindex {
    my ($self, $input, $p1, $p2) = @_;
    unless (defined $p2)  {
	$self->throw("Wrong number of parameters passed to lower_pdlindex");
    }
    my $result;
    my $i = $input;
    until ($result = $self->pdlindex($i, $p1 => $p2))  {
	$i++;
	last unless ($self->pdlindex($i, $p1=>0) > 0);
    }
    return $result;
}




sub _strip_gaps {
    # not OO
    my $seq = shift;
    $seq =~ s/\-|\.//g;
    return $seq;
}



sub _abs_pos  {
    # this is ugly, but I have no time to make it better

    my ($self, $pos, $seqnr) = @_;
    no strict 'refs';
    my ($strand, $seqobj) = ("strand$seqnr", "seq$seqnr");
    if ($self->{job}->$strand eq "-1")  {
	return $self->{job}->$seqobj->length() - $pos  +1;
    }
    else  {
	return $pos;
    }

}

sub _abs_sign  {
    # this is also ugly, but I have no time to make it better

    my ($self, $sign, $seqnr) = @_;
    no strict 'refs';
    my $strand = "strand$seqnr";

    if ($self->{job}->$strand eq "-")  {
	return ($sign eq '-') ? "+" : "-";
    }
    else  {
	return $sign;
    }

}


1;


