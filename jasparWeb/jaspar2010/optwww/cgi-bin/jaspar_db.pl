#! /usr/bin/perl -w
use lib "/opt/www/jaspar_2010/cgi-bin/";

use JasparDB;

my $web= JasparDB->new();
$web->run;
