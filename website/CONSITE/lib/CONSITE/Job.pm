package CONSITE::Job;
use vars qw(@ISA);
use Class::MethodMaker
    new_with_init => 'new',
    get_set       => [qw(db jobID seq1 seq2 seq3 min_ic ID 
			 strand1 strand2 strand3 strandref rangeref
			 alignstring TFlistref tmpdir exclude_orf
			 alignobject threshold window MatrixSet 
			 user_matrix
			 cutoff)],
    list          => IDlist;
    
use Bio::Root::Root;
use Persistence::Object::Simple;
use CONSITE::Alignment;
use Bio::Seq;
use strict;

@ISA = qw(Bio::Root::Root Persistence::Object::Simple);

sub init  {
    my $self = shift;
    my %args = ( -jobID      => undef,  #time.sprintf("%04d", rand(10000)),
		 -alignstring  => undef,
		 -seq1       => undef,
		 -seq2       => undef,
		 -seq3       => undef,
		 -db         => undef,           # compulsory,
		 -IDlist     => [],
		 -min_ic     => 0,
		 -MatrixSet  => undef,
		 -tmpdir     => '/tmp',
		 -window     => 50,
		 -cutoff     => undef,
		 -threshold  => "80%",
		 @_);

    $self->{__Fn} = $args{-tmpdir}."/".$self->jobID();
    $self->tmpdir($args{-tmpdir});
    $self->jobID($args{-jobID}     or $self->throw ("No jobID provided."));

    $self->seq1($args{-seq1} 
		or $self->throw ("No sequence object seq1 provided.") );
    $self->seq2($args{-seq2}  or undef);
		#or $self->throw ("No sequence object seq2 provided.") );
    $self->seq3($args{-seq3}); 

    $self->TFlistref($args{-TFlist});
    $self->min_ic($args{-min_ic});
    $self->alignstring($args{-alignstring} or undef);
    $self->MatrixSet($args{-MatrixSet});
    $self->window($args{-window});
    $self->cutoff($args{-cutoff});
    $self->threshold($args{-threshold});
   

	      
	      
		    
    # check range : set to (1..length(seq1)) if not defined or out of range

    $self->start_at($args{-start_at} or $self->start_at() or 1);
    if ($self->start_at() < 1) { $self->start_at(1); }

    $self->end_at($args{-end_at} 
		   or $self->end_at() 
		   or $self->seq1->length());
    if ($self->end_at() > $self->seq1->length() ) { 
	$self->end_at($self->seq1->length()); 
    }

    return $self;
}
    
sub start_at {
    my ($self, $new_value) = @_;
    if (defined $new_value)  { 
	$new_value = 1 if $new_value<1 ;
	$self->{start_at} = $new_value;
    }
    return $self->{start_at};

}
sub end_at {
    my ($self, $new_value) = @_;
    if (defined $new_value)  { 
	$new_value = $self->seq1->length() 
	    if $new_value > $self->seq1->length() ;
	$self->{end_at} = $new_value;
    }
    return $self->{end_at};

}
    
    
sub load  {
    # class method
    my ($class, $filename) = @_;
    my $job = Persistence::Object::Simple->new(__Fn => $filename)
	or do {warn("Error loading $filename as Job object."); return undef;};
    return bless $job, ref($class) || $class;
}

1;






