package CONSITE::UserMatrix;
use vars qw(@ISA);

use Class::MethodMaker
    new_with_init => 'new',
    get_set       => [qw(ID name pmatrix type Error)];
    
use Bio::Root::Root;
use strict;

@ISA = qw(Bio::Root::Root);

sub init  {
    my ($self, %args) = @_;
    $self->ID   ($args{'-ID'} or 'UserDefined');
    $self->name ($args{'-name'} or 'MyMatrix');
    
    $self->type(uc($args{'-type'} or 'PFM'));

    unless (defined $args{-matrixstring})  { 
	$self->throw("No -matrix data provided for UserMatrix.");
    }
    print STDERR "MATRIXSTRING:".$args{-matrixstring};
    $self->_parse_matrix($args{-matrixstring}); # FIXME: error handling here
    
    return $self;

}



sub to_PWM  {
    
}



sub to_PFM  {

}

sub get_Matrix  {
    my ($self, $mt) = (@_, "PWM");
    my $matrixobj;
    eval("\$matrixobj= TFBS::Matrix::".$self->type()."->new".' 
            ( -ID    => $self->ID()."",
              -name  => $self->name()."",
              -class => "N.D.",
              -matrix=> $self->pmatrix()   # FIXME - temporary
              );');
    if ($@) {$self->throw($@); }
    if ($mt eq "PWM" and $self->type() ne "PWM") { 
	return $matrixobj->to_PWM();
    }
    elsif ($mt eq "PFM" and $self->type() eq "PWM")  {
	$self->throw("Cannot convert PWM to PFM");
    }
    else  {
	return $matrixobj;
    }
}


sub get_MatrixSet {
    # produce matrixset object from user matrix 
    my ($self, %args) = @_;
    my $mt = ($args{-matrixtype} or "PWM");
    my $matrixset = TFBS::MatrixSet->new();
    $matrixset->add_Matrix($self->get_Matrix(uc $mt));
    
    return $matrixset;
}

sub _parse_matrix  {
    my ($self, $matrixstring) = @_;
    $matrixstring =~ s/\r\n/\n/g;     # deal with DOS strings
    $matrixstring =~ s/\r/\n/g;       # deal with Mac strings
    $matrixstring =~ s/^[\n\s]//;     # remove any leading spaces and newlines
    $matrixstring =~ s/[\n\s]+$//;    # remove any trailing spaces and newlines
    $matrixstring =~ s/\s*\n\s*/\n/g; # strip spaces from rows
    print STDERR "PROCESSED MATRIXSTRING:".$matrixstring;
    my @matrix_rows = split ("\n", $matrixstring);
    if (@matrix_rows != 4) {
	$self->Error("Matrix has ".scalar @matrix_rows.
		     "instead of the required 4 rows");
	return undef;
    }
    my (@matrix, $nr_columns);
    foreach my $row (@matrix_rows)   {
	my @row_elements = split(/\s+/,$row);
	if (defined $nr_columns)  {
	    if ($nr_columns != scalar(@row_elements))  {
		$self->Error("Matrix rows have unequal numbers of elements.");
		return undef;
	    }
	}
	else  {
	    $nr_columns = scalar(@row_elements);
	}
	push @matrix, [split(/\s+/,$row)];
    }
    $self->pmatrix(\@matrix);
    return 1; # success
}
    
