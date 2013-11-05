#!/usr/bin/env perl -w
use strict;
use TFBS::DB::JASPAR2;
use TFBS::DB::FlatFileDir;

my $mysqldb = TFBS::DB::JASPAR2->connect
    ("dbi:mysql:JASPAR2:forkhead.cgb.ki.se", 
     "consite",	
     "bezimeni");
my $flatdb = TFBS::DB::FlatFileDir->create("/tmp/JASPAR2");

my $matrixset = $mysqldb->get_MatrixSet(-matrixtype=>"PFM");

my $iter = $matrixset->Iterator;

while (my $pfm = $iter->next) {
    $flatdb->store_Matrix($pfm);
}

system("tar -cz -C /tmp JASPAR2 > JASPAR2_FlatFileDir.tar.gz");

unlink </tmp/JASPAR2/*>;
rmdir "/tmp/JASPAR2";

print STDERR "File JASPAR2_FlatFileDir.tar.gz written to current directory\n";
