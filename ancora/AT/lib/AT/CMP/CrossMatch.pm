package AT::CMP::CrossMatch;

use strict;
use vars '@ISA';
use Carp;
use AT::Root;


@ISA = qw(AT::Root);


sub new  {
    my ($caller, %args) = @_;
    my $self = bless {
	mapping1   => $args{'mapping1'   },
	mapping2   => $args{'mapping2'   },
	alignment1 => $args{'alignment1' },
	alignment2 => $args{'alignment2' },
	genomic1   => $args{'genomic1'   },
	genomic2   => $args{'genomic2'   }
	}, ref $caller || $caller;
    

}


sub mappings  {

}




sub cross_align {

}


