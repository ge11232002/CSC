#! /usr/bin/perl -w


use TFBS::Matrix::PFM;
use TFBS::Matrix::PWM;
use TFBS::Matrix::ICM;
use lib "/opt/www/jaspar/lib/TFBS";
use lib "/opt/www/jaspar/lib/";
use TFBS::DB::FlatFileDir;
use TFBS::DB::JASPAR4;

 my $db=TFBS::DB::JASPAR4->connect("dbi:mysql:JASPAR4:mordor.cgb.ki.se",
                                      "albin",
                                      "as4488", die_on_bad_params=>0);


my $selected = $db->get_MatrixSet(
                                    -min_ic=>0.0001
                                    );
my $matrixset_iterator =
                   $selected->Iterator();
                   


#
#

my $dirdb = TFBS::DB::FlatFileDir->create("/opt/www/jaspar/html/DOWNLOAD/MatrixDir2");




while (my $pfm = $matrixset_iterator->next) {

    $dirdb->store_Matrix($pfm);

}
