package InitVersion;
use strict;
use WEBSERVER_OPT;
use JasparSubDB;
# this will read in teh version info and what server you are at 
# it will read a pre-defined file called DatabaseVersion.init, and give back a hash 

my %data;


use constant BASE_DIR =>  "/opt/www/jaspar_2010/";
use constant CGI_BASE_DIR => BASE_DIR."cgi-bin/";
my $dir = CGI_BASE_DIR;
sub  init{
    open (INIT_FILE, $dir. "DatabaseVersion.init")|| die "cannot open  DatabaseVersion.init";
    while (<INIT_FILE>){ # parse out the values to variosu dbs and make 
	next if /^#/; 
	chomp;
	next unless $_;
	my ($key, $val)= split(/\t/, $_);
	
	$data{$key}=$val;	
    }
# return the array
    return \%data;
close INIT_FILE
}









1;
