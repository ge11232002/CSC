#!/bin/bash


#hg19
hg19()
{
  # The net file
  rsync -avzP rsync:://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/netDanRer7* \
    /export/data/CNEs/hg19/annotation/



}
