#! /usr/bin/perl -w
use DBI;

my $dbh=DBI->connect("dbi:mysql:JASPAR4:localhost",
				      "root",
				      "a2580s");

# open files

foreach my $TABLE( @ARGV){
my $sth=$dbh->prepare("SELECT * FROM $TABLE ORDER BY ID");
$sth->execute();
while (my @res= $sth->fetchrow_array()){
    
    print join ("\t", @res), "\n";
    
}
}