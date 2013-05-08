#!/usr/bin/env perl

use strict;
use warnings;

use AT::Tools::SequenceMapper;
use AT::Tools::OrthoReceiver;



use constant HUMAN_GL_DBNAME => "GLYNX_1_2";
use constant MOUSE_GL_DBNAME => "GLMOUSE_1_RC1";
use constant GL_HOST   => "sql.cgb.ki.se";
use constant GL_USER   => "genelynx";
use constant GL_PASS   => "bezimeni";

use Test;
plan(tests => 6);

my $hgldb = GeneLynx::MySQLdb->connect(-dbname =>  HUMAN_GL_DBNAME,
					 -dbhost => GL_HOST,
					 -dbuser => GL_USER,
					 -dbpass => GL_PASS);

my $hmapdb = AT::DB::GenomeMapping->connect(-dbname => "AT_HS_JUN02",
					       -dbhost => "nautilus.cgb.ki.se",
					       -dbuser => "at_read",
					       -dbpass => "tittut");

my $hmapper  = AT::Tools::SequenceMapper->new(-extres => "word_index", 
					      -id => "pai1", 
					      -species => "human", 
					      -hgldb => $hgldb, 
					      -hmapdb => $hmapdb);
my $hmapping = $hmapper->get_mapping();

ok($hmapping->tStart ,99254010);
ok($hmapping->strand ,"+");

my $mgldb = GeneLynx::MySQLdb->connect(-dbname => MOUSE_GL_DBNAME,
					  -dbhost => GL_HOST,
					  -dbuser => GL_USER,
					  -dbpass => GL_PASS);
my $mmapdb = AT::DB::GenomeMapping->connect(-dbname => "AT_MM_FEB02_2",
					       -dbhost => "nautilus.cgb.ki.se",
					       -dbuser => "at_read",
					       -dbpass => "tittut");
my $mmapper  = AT::Tools::SequenceMapper->new(-extres => "word_index", 
					      -id => "pai1", 
					      -species => "mouse", 
					      -hgldb => $hgldb, 
					      -mgldb => $mgldb, 
					      -mmapdb => $mmapdb);
my $mmapping = $mmapper->get_mapping();

ok($mmapping->tStart ,135683497);
ok($mmapping->strand ,"-");



my $hgendb = AT::DB::GenomeAssembly->connect(-dbname => "HS_JUN02",
					     -dbhost => "nautilus.cgb.ki.se",
					     -dbuser => "at_read",
					     -dbpass => "tittut");

my $mgendb = AT::DB::GenomeAssembly->connect(-dbname => "MM_FEB02",
					     -dbhost => "nautilus.cgb.ki.se",
					     -dbuser => "at_read",
					     -dbpass => "tittut");

my $receiver= AT::Tools::OrthoReceiver->new(-hmapping => $hmapping,
					    -mmapping => $mmapping,
					    -refspec  => "human",
					    -hgendb   => $hgendb,
					    -mgendb   => $mgendb);

my $hseq = $receiver->get_hseq();
my $mseq = $receiver->get_mseq();

my $hsubseq = $hseq->subseq(1, 10);

ok($hsubseq, "GAATTCCTGC" );

my $msubseq = $mseq->subseq(1, 10);

ok($msubseq, "TCAGGCTGCA");


