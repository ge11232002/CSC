package InitDB;
use strict;
use DBI;
use TFBS::DB::JASPAR5;
use JasparSubDB;
use WEBSERVER_OPT;


# this will read in a file that says the overall connection, and lists the actual collections and their descriptions, and the keys with which to use
# wil give back an object with those thinsg (a glorified hash with a connected databse)

# this will read in the actual sub-databases used, their settings etc AND what filelds thwy will be using in browinsg and popups
# it will read a pre-defined file called SubDataBases.init, and give back an array of JasparSubDB objects
my %global;
my @sub_dbs;
my %data;


# This is a temporary hack. Use WEBSERVER_OPT
use constant BASE_DIR =>  "/opt/www/jaspar_2010/cgi-bin/";
use constant CGI_BASE_DIR => BASE_DIR."cgi-bin/";

my $dir = CGI_BASE_DIR;

sub init{
    open (INIT_FILE,"/opt/www/jaspar_2010/cgi-bin/SubDataBases.init")|| die "Cannot open ".$dir."SubDataBases.init";
    while (<INIT_FILE>){ # parse out the values to variosu dbs and make 
	next if /^#/; 
	chomp;
	
	next unless $_;
	if ($_ eq "START_ENTRY"){ # make a new entry
	    undef %data;
	    next;
	}


	if ($_ eq "START_ENTRY"){ # make a new entry
	    undef %data;
	    next;
	}
	if ($_ eq "END_ENTRY"){ # finish entry
	    my $db= JasparSubDB->new(%data); # this is now just a collection of keys, descriptions etc, not the actual connection
	    push (@sub_dbs,$db); 
	    
	    next;
	}
	
	
	my ($key, $val)= split(/\t/, $_);
	

	if ( $key eq "BROWSE_KEYS"){
	    $data{$key}= [split(",", $val)];
	}
	elsif ( $key eq "COLLECTION"){
	    $data{$key}=  $val;
	}
	


	elsif ( $key eq "POP_UP_KEYS"){
	    $data{$key}= [split(",", $val)];
	}
	elsif ( $key eq "DESCRIPTION"){
	    $data{$key}=  $val;
	}

	else{
	    $global{$key}= $val;
	}
       
	
	
    }
# establish connection
 
    
    my $string= "dbi:mysql:".$global{'DB_NAME'} .":".$global{'MYSQL_SERVER'};
    warn "$global{'MYSQL_USER'} $global{'MYSQL_PASSWD'} ";
    my $jaspar5=    TFBS::DB::JASPAR5->connect( $string,
					    $global{'MYSQL_USER'}, 
					    $global{'MYSQL_PASSWD'}, 
					    die_on_bad_params=>1
					);

    warn ref ($jaspar5);
  


#make an object that has 1) the connection, 2) the collections which you can query for browsekeys, etc
   
    my $self  = {};
    $self->{jaspar5} =$jaspar5;
    # remap the dbs into a hash keyed by collection name
    my %h;
    foreach my $subdb (@sub_dbs){
#	print $subdb->_show;
	$h{$subdb->COLLECTION}=$subdb;
    }
    
   
    $self->{sub_dbs}  = \%h;
       

    bless($self), "JASPAR_DB_INFO";           # but see below

  #my $self =bless {
	    #%args}, "JASPAR_DB_INFO";
 	return $self;


# make a hash from the subDBs - should really have amde a container class...
  #  return \@sub_dbs;
}

sub jaspar5{ # the jaspar 5 (a connected database object)
    return $_[0]->{'jaspar5'}
}


sub collections{ #gives an array of strings - the collection names
    my $self= shift @_;
    return keys %{ $self->{sub_dbs}};

}

sub description{ #given a collection, give the description
    my ($self, $collection)= @_;
    return $self->{sub_dbs}{$collection}->DESCRIPTION;
 
}



sub browse_keys{ # give an array of keys to use for borwser visualization. Input is the collection name
 my ($self, $collection)= @_;
    return $self->{sub_dbs}{$collection}->BROWSE_KEYS;

}

sub pop_up_keys{ # as above but keys for the actual popup page
 my ($self, $collection)= @_;
    return $self->{sub_dbs}{$collection}->POP_UP_KEYS;

}



1;
