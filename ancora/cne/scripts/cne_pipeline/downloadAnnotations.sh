#!/bin/bash


#hg19
hg19()
{
  cd /export/data/CNEs/hg19/annotation/
  # The net file
  rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/netDanRer7* .




}
