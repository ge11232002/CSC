#!/bin/bash

#hg19
hg19()
{
  cd /export/data/CNEs/hg19/annotation/
  gunzip *.txt.gz
  # The net file
  cat netDanRer7.sql | mysql -u root -pgenome UCSC_hg19
  mysqlimport -u root -Lpgenome UCSC_hg19 netDanRer7.txt

}


#danRer7
danRer7()
{
  cd /export/data/CNEs/danRer7/annotation/
  gunzip *.txt.gz
  # The net file
  cat netHg19.sql | mysql -u root -pgenome UCSC_danRer7
  mysqlimport -u root -Lpgenome UCSC_danRer7 netHg19.txt

}

