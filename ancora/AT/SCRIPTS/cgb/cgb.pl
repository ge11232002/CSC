#!/usr/bin/perl -w
use strict;
use lib '/Home/bccs/engstrom/DEVEL/AT/lib';
use AT::WWW::ComparativeGenomeBrowser;
use CGI;

$ENV{'PATH'} = $ENV{'PATH'}.':/net/bccs/cbu/linux/blatSuite';

my %config = ( html_dir => '/Home/bccs/engstrom/html/cgb',
               tmp_dir => '/Home/bccs/engstrom/html/cgb/tmp',
               tmp_dir_url => '/html/cgb/tmp',
	       drawing_config_file => 'draw_loci.conf');

my $interface = AT::WWW::ComparativeGenomeBrowser->new(config => \%config);
$interface->handle_request(CGI->new);

