#!/bin/bash


#hg19
hg19()
{
  cd /export/data/CNEs/hg19/annotation/
  # The net file
  rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/netDanRer7.* .


}
#mm10
mm10()
{
  cd /export/data/CNEs/mm10/annotation/

}

#canFam3
canFam3()
{
  cd /export/data/CNEs/canFam3/annotation/

}




#danRer7
danRer7()
{
  cd /export/data/CNEs/danRer7/annotation/
  # The net file
  rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/danRer7/database/netHg19.* .


}
