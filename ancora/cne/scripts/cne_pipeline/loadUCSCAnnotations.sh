#!/bin/bash

#hg19
hg19()
{
  cd /export/data/CNEs/hg19/annotation/
  gunzip *.txt.gz
  # The net file
  cat netDanRer7.sql | mysql -u root -p UCSC_hg19
  mysqlimport -u root -Lp UCSC_hg19 netDanRer7.txt
  cat netMm10.sql    | mysql -u root -p UCSC_hg19
  mysqlimport -u root -Lp UCSC_hg19 netMm10.txt
  cat netCanFam3.sql | mysql -u root -p UCSC_hg19
  mysqlimport -u root -Lp UCSC_hg19 netCanFam3.txt
  cat netEquCab2.sql | mysql -u root -p UCSC_hg19
  mysqlimport -u root -Lp UCSC_hg19 netEquCab2.txt
  cat netFr2.sql     | mysql -u root -p UCSC_hg19
  mysqlimport -u root -Lp UCSC_hg19 netFr2.txt
  cat netGalGal3.sql | mysql -u root -p UCSC_hg19
  mysqlimport -u root -Lp UCSC_hg19 netGalGal3.txt
  cat netGasAcu1.sql | mysql -u root -p UCSC_hg19
  mysqlimport -u root -Lp UCSC_hg19 netGasAcu1.txt
}


#danRer7
  cd /export/data/CNEs/danRer7/annotation/
  gunzip *.txt.gz
  # The net file
  cat netHg19.sql | mysql -u root -p UCSC_danRer7
  mysqlimport -u root -Lp UCSC_danRer7 netHg19.txt
  cat netTetNig2.sql | mysql -u root -p UCSC_danRer7
  mysqlimport -u root -Lp UCSC_danRer7 netTetNig2.txt
  # The rmsk
  cat rmsk.sql | mysql -u root -p UCSC_danRer7
  mysqlimport -u root -Lp UCSC_danRer7 rmsk.txt

##ce10
mysql -u root -p -e 'create database UCSC_ce10'
mysql -u root -p -e 'grant select on UCSC_ce10.* to nobody@localhost'
cd /export/data/CNEs/ce10/annotation/
gunzip *.txt.gz
# The rmsk
cat rmsk.sql | mysql -u root -p UCSC_ce10
mysqlimport -u root -Lp UCSC_ce10 rmsk.txt



