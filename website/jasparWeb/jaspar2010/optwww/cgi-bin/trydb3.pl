#! /usr(bin/perl -w
use lib "/home_local/albin/DEVEL/TFBS/TFBS";
use DB::JASPAR4;
use DB::JASPAR2;


my $db2=TFBS::DB::JASPAR4->connect("dbi:mysql:JASPAR4:localhost",
				      "root",
				      "a2580s");


 
my $db1 =
	    TFBS::DB::JASPAR2->connect("dbi:mysql:JASPAR2:forkhead.cgb.ki.se",
					"albin",
					"as4488");
	    

my $matrixset = $db2->get_MatrixSet(
                                    -ID=>[MA0084], 
                            );


my $matrix_iterator = 
                   $matrixset->Iterator(-sort_by =>'name');
           while (my $matrix_object = $matrix_iterator->next) {
               
               print $matrix_object->name(), "\n";
               # do whatever you want with individual matrix objects
           }
#my $pfm3= $db2->get_Matrix_by_name('GAMYB');
#print $pfm3->prettyprint();
#$db2->delete_Matrix_having_ID('MA0034');