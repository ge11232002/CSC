package JasparSubDB;
use strict;

# this is a class describing a sub-database,  filelds to show in browsing ans popup window. Mostly for holding data, but should also have a connect statement, which gives back a database handler object as usual


sub new {
 	my ($class, %args) = @_; 
    my $self =bless {
	    %args}, $class;
 	return $self;

}



# accessor functions
sub COLLECTION{ return $_[0]->{"COLLECTION"};  }
sub BROWSE_KEYS{ return $_[0]->{"BROWSE_KEYS"};}  
sub POP_UP_KEYS{ return $_[0]->{"POP_UP_KEYS"};}  
sub DESCRIPTION{ return $_[0]->{"DESCRIPTION"};}  

 sub _show{
    #sanity check: shows all information the object holds
    my $self=shift;
    require Data::Dumper;
    return Data::Dumper->Dump([$self]);


}

1;
