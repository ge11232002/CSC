#! /usr/bin/perl -w
use lib "/opt/www/jaspar/cgi-bin/";

use JasparDB;

my $web= JasparDB->new();
$web->run;
