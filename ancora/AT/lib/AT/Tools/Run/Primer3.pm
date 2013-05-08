package AT::Tools::Run::Primer3;

use strict;
use Boulder::Stream;
use File::Temp qw/ tempfile tempdir /;


sub new {
    my ($caller, %args) =@_;
    my $executable = "primer3_core";
    if (defined $args{-executable}) {
	$executable = $args{-executable};
	delete $args{-executable};
    }
    return bless {
		    _inputstone => Stone->new(%args),
		    _executable => $executable,
		    _tempdir    => tempdir("Primer3XXXXXXXX", TMPDIR => 1,
						CLEANUP => 1)
		}, ref $caller || $caller;
}


sub set_param {
    my ($self, %args) = @_;
    $self->{_inputstone}->replace(%args);
}

sub inputstone {
    $_[0]->{_inputstone};
}

sub resultstone {
    $_[0]->{_resultstone};
}

sub run {
    my ($self, %args) = @_;
    if (defined $args{'-sequence'}) {
	$self->set_sequence($args{'-sequence'});
    }
    
    unless($self->inputstone->get("SEQUENCE")) {
	die "No sequence provided";
    }
    my $temp_file = File::Temp::tempnam($self->{_tempdir}, "stone");
    open (OUTFH, "| " . $self->{_executable} . " > $temp_file") or die;
    my $bs = Boulder::Stream->new(
				  -out=>\*OUTFH);
    $bs->put($self->{_inputstone});
    #$bs->close;
    close OUTFH;
    open(INFH, "$temp_file") or die;
    my $is = Boulder::Stream->new(
				  -in =>\*INFH );

    return $self->{_resultstone} = $is->get;
}    


sub set_sequence {
    my ($self, $seq) = @_;
    if (ref($seq) and $seq->isa("Bio::SeqI")) {
	$self->set_param(SEQUENCE => $seq->seq,
			 PRIMER_SEQUENCE_ID => $seq->id);
    }
    else {
	$self->set_param(SEQUENCE => $seq);
	
    }
}

sub get_primer_pairs {
    my ($self) = @_;
    return undef unless $self->resultstone;
    my @tags = $self->resultstone->tags;
    my @primer_pairs;
    foreach my $tag (@tags)  {
	my ($primer_no, $newtag);
	if ($tag =~ /PRIMER_(\w+)_(\d+)(\w*)/) {
	    $primer_no = $2 +1;
	    $newtag = $1.$3;
	}
	elsif ($tag =~ /PRIMER_(\w+)/)  {
	    $primer_no = 1;
	    $newtag = $1;
	    
	}
	else {
	    next;
	}

	unless ($primer_pairs[$primer_no-1])  {
	    $primer_pairs[$primer_no-1] = Stone->new(RANK=>$primer_no);
	}
	#print "INSERTING $primer_no $newtag\n";
	$primer_pairs[$primer_no-1]->insert($newtag 
					   => $self->resultstone->get($tag)
					   );
    }
    #my @found;
    #foreach my $pp (@primer_pairs) { push (@found, $pp) if $pp ; }
    return @primer_pairs;
    
}    

1;
