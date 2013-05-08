package AT::Root;

use strict;
use Carp;
use vars '$AUTOLOAD';

sub verbose {
    my ($self, $newvalue) = @_;
    if(defined($newvalue)) { $self->{verbose} = $newvalue; }
    elsif(!defined($self->{verbose})) { $self->{verbose} = 0; }
    return $self->{verbose};
}

sub debug {
    if(defined($_[1])) { carp("Method debug() does not take arguments. Use verbose().\n"); }
    return 1 if($_[0]->verbose() > 1);
}

sub _verbose_msg { print STDERR @_[1..$#_] if($_[0]->verbose()); }

sub _debug_msg { print STDERR @_[1..$#_] if ($_[0]->debug()); }


sub AUTOLOAD  {
    my ($self, $newvalue) = @_;
    my ($arg) = $AUTOLOAD =~ /:([^:]+)$/;
    if (exists $self->{$arg})  {
	$self->{$arg} = $newvalue if defined $newvalue;
	return $self->{$arg};
    }
    else {
	if ($arg =~ /(\w+)_list$/)  {
	    if (exists $self->{"_".$1."L"})  {
		$arg = '_'.$1.'L';
		$self->{$arg} = $newvalue if defined $newvalue;
		return @{$self->{$arg}};
	    }
	}
	if ($arg =~ /(\w+)_listref$/)  {
	    if (exists $self->{"_".$1."L"})  {
		$arg = '_'.$1.'L';
		$self->{$arg} = $newvalue if defined $newvalue;
		return $self->{$arg};
	    }
	}
    }

    croak "No such attribute: ".(ref($self)||$self)."::".$arg;

}

sub DESTROY {}

1;
