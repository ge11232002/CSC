#! /usr/ bin/perl -w
# prepare text dumps of two kinds for each database
# locations are hard-coded, change this if you are on another machine

use lib "/home/albin/public_html/jaspar_2010/";



use TFBS::DB::FlatFileDir;
use TFBS::DB::JASPAR5;
use TFBS::Matrix::PFM;

use InitDB;
my ($db_info)=InitDB->init;

my $dbh= $db_info->jaspar5->dbh;


# make flatfiles: one per table

# jaspar core
my $an= $core_dbh->prepare (qq!select * from MATRIX_ANNOTATION!);
open AN, ">/home/sandelin/DEVEL/JASPAR/html/DOWNLOAD/mySQL/JASPAR_CORE/MATRIX_ANNOTATION.txt";

$an->execute();
while (my (@ary)= $an->fetchrow_array()){
    print AN join ("\t", @ary), "\n";
}

close AN;
my $da=$core_dbh->prepare (qq!select * from MATRIX_DATA!);
open DA, ">/home/sandelin/DEVEL/JASPAR/html/DOWNLOAD/mySQL/JASPAR_CORE/MATRIX_DATA.txt";

$da->execute();
while (my (@ary)= $da->fetchrow_array()){
    print DA join ("\t", @ary), "\n";
}

close DA;


my $in=$core_dbh->prepare (qq!select * from MATRIX_INFO!);
open IN, ">/home/sandelin/DEVEL/JASPAR/html/DOWNLOAD/mySQL/JASPAR_CORE/MATRIX_INFO.txt";

$in->execute();
while (my (@ary)= $in->fetchrow_array()){
    print IN join ("\t", @ary), "\n";
}

close IN;






# jaspar fam
my $an= $fam_dbh->prepare (qq!select * from MATRIX_ANNOTATION!);
open AN, ">/home/sandelin/DEVEL/JASPAR/html/DOWNLOAD/mySQL/JASPAR_FAM/MATRIX_ANNOTATION.txt";

$an->execute();
while (my (@ary)= $an->fetchrow_array()){
    print AN join ("\t", @ary), "\n";
}

close AN;
my $da=$fam_dbh->prepare (qq!select * from MATRIX_DATA!);
open DA, ">/home/sandelin/DEVEL/JASPAR/html/DOWNLOAD/mySQL/JASPAR_FAM/MATRIX_DATA.txt";

$da->execute();
while (my (@ary)= $da->fetchrow_array()){
    print DA join ("\t", @ary), "\n";
}

close DA;


my $in=$fam_dbh->prepare (qq!select * from MATRIX_INFO!);
open IN, ">/home/sandelin/DEVEL/JASPAR/html/DOWNLOAD/mySQL/JASPAR_FAM/MATRIX_INFO.txt";

$in->execute();
while (my (@ary)= $in->fetchrow_array()){
    print IN join ("\t", @ary), "\n";
}

close IN;

# jaspar phylof

my $an= $phylo_dbh->prepare (qq!select * from MATRIX_ANNOTATION!);
open AN, ">/home/sandelin/DEVEL/JASPAR/html/DOWNLOAD/mySQL/JASPAR_PHYLOFACTS/MATRIX_ANNOTATION.txt";

$an->execute();
while (my (@ary)= $an->fetchrow_array()){
    print AN join ("\t", @ary), "\n";
}

close AN;
my $da=$phylo_dbh->prepare (qq!select * from MATRIX_DATA!);
open DA, ">/home/sandelin/DEVEL/JASPAR/html/DOWNLOAD/mySQL/JASPAR_PHYLOFACTS/MATRIX_DATA.txt";

$da->execute();
while (my (@ary)= $da->fetchrow_array()){
    print DA join ("\t", @ary), "\n";
}

close DA;


my $in=$phylo_dbh->prepare (qq!select * from MATRIX_INFO!);
open (IN, ">/home/sandelin/DEVEL/JASPAR/html/DOWNLOAD/mySQL/JASPAR_PHYLOFACTS/MATRIX_INFO.txt")||die;

$in->execute();
while (my (@ary)= $in->fetchrow_array()){
    print IN join ("\t", @ary), "\n";
}

close IN;


my $phylofacts_db=TFBS::DB::JASPAR4->connect("dbi:mysql:JASPAR_PHYLOFACTS4:anx095",
				      "albin",
				      "as4488", die_on_bad_params=>0);


my $matrixset = $phylofacts_db->get_MatrixSet();
system "rm -r /home/sandelin/DEVEL/JASPAR/html/DOWNLOAD/MatrixDir/JASPAR_PHYLOFACTS/";
my $flatfile_db = TFBS::DB::FlatFileDir->create("/home/sandelin/DEVEL/JASPAR/html/DOWNLOAD/MatrixDir/JASPAR_PHYLOFACTS/");
my $matrix_iterator = $matrixset->Iterator(-sort_by =>'ID');

while (my $matrix_object = $matrix_iterator->next) {
    
    
    $flatfile_db->store_Matrix($matrix_object);
    # print $matrix_object->name(), "\n";
               # do whatever you want with individual matrix objects
}





my $core_db=TFBS::DB::JASPAR4->connect("dbi:mysql:JASPAR4:anx095",
				      "albin",
				      "as4488", die_on_bad_params=>0);

my $matrixset = $core_db->get_MatrixSet();
system "rm -r /home/sandelin/DEVEL/JASPAR/html/DOWNLOAD/MatrixDir/JASPAR_CORE/";
my $flatfile_db = TFBS::DB::FlatFileDir->create("/home/sandelin/DEVEL/JASPAR/html/DOWNLOAD/MatrixDir/JASPAR_CORE/");
my $matrix_iterator = $matrixset->Iterator(-sort_by =>'ID');

while (my $matrix_object = $matrix_iterator->next) {
    
    
    $flatfile_db->store_Matrix($matrix_object);
    # print $matrix_object->name(), "\n";
               # do whatever you want with individual matrix objects
}



 


my $familial_db=TFBS::DB::JASPAR4->connect("dbi:mysql:JASPAR_FAM4:anx095",
				      "albin",
				      "as4488", die_on_bad_params=>0);
my $matrixset = $familial_db->get_MatrixSet();
system "rm -r /home/sandelin/DEVEL/JASPAR/html/DOWNLOAD/MatrixDir/JASPAR_FAM/";
my $flatfile_db = TFBS::DB::FlatFileDir->create("/home/sandelin/DEVEL/JASPAR/html/DOWNLOAD/MatrixDir/JASPAR_FAM/");
my $matrix_iterator = $matrixset->Iterator(-sort_by =>'ID');

while (my $matrix_object = $matrix_iterator->next) {
    
    
    $flatfile_db->store_Matrix($matrix_object);
    # print $matrix_object->name(), "\n";
               # do whatever you want with individual matrix objects
}
