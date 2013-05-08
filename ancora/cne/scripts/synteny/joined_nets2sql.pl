#!/usr/bin/perl
use strict;
use warnings;

# This script simply prints SQL commands instead of connecting via DBI, 
# so that we can see warnings more easily. 

#use DBI;
#my $DB_HOST = 'localhost';
#my $DB_NAME = 'AT_DR_MAR06';
#my $dbh = DBI->connect("dbi:mysql:host=$DB_HOST;database=$DB_NAME;port=$DB_PORT","at_read","tittut")
#    or die "could not connect to db $UCSC_DB_NAME @ $DB_HOST:$DB_PORT";

my ($MAX_GAP1, $MAX_GAP2) = @ARGV;

die "usage: perl $0 <max_gap1> <max_gap2>\n" unless($MAX_GAP2);

my $table1 = "joinedNets_hg18_gap".int($MAX_GAP1/1000)."k".int($MAX_GAP2/1000)."k";
my $table2 = "joinedNets_hg18_gap".int($MAX_GAP1/1000)."k".int($MAX_GAP2/1000)."k_detail";

my $file1 = "joinedNets_danRer4_hg18_gap".int($MAX_GAP1/1000)."k".int($MAX_GAP2/1000)."k.txt";
my $file2 = "joinedNets_danRer4_hg18_gap".int($MAX_GAP1/1000)."k".int($MAX_GAP2/1000)."k.danRer4_detail.bed";

# create tables

print "create table $table1 (
  `chr1` varchar(40) NOT NULL default '',
  `start1` int(10) unsigned NOT NULL default '0',
  `end1` int(10) unsigned NOT NULL default '0',
  `chr2` varchar(40) NOT NULL default '',
  `start2` int(10) unsigned NOT NULL default '0',
  `end2` int(10) unsigned NOT NULL default '0',
  `nr_of_nets` smallint(5) unsigned NOT NULL default '0',
  `status` varchar(20) NOT NULL default '',
  `ali` int(10) unsigned NOT NULL default '0',
  KEY `chr1` (`chr1`,`start1`),
  KEY `chr1_2` (`chr1`,`end1`),
  KEY `chr2` (`chr2`,`start2`),
  KEY `chr2_2` (`chr2`,`end2`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;\n";

print "create table $table2 (
  `chrom` varchar(255) NOT NULL default '',
  `chromStart` int(10) unsigned NOT NULL default '0',
  `chromEnd` int(10) unsigned NOT NULL default '0',
  `name` varchar(255) NOT NULL default '',
  `score` smallint(5) unsigned NOT NULL default '0',
  `strand` char(1) NOT NULL default '',
  `thickStart` int(10) unsigned NOT NULL default '0',
  `thickEnd` int(10) unsigned NOT NULL default '0',
  `itemRgb` varchar(255) default NULL,
  `blockCount` int(10) unsigned NOT NULL default '0',
  `blockSizes` longblob NOT NULL,
  `blockStarts` longblob NOT NULL,
  KEY `chrom` (`chrom`,`chromStart`),
  KEY `chrom_2` (`chrom`,`chromEnd`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;\n";

# load data

print "load data local infile '$file1' into table $table1;\n";
print "load data local infile '$file2' into table $table2;\n";
