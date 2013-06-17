#!/bin/bash


#hg19
hg19()
{
  cd /export/data/CNEs/hg19/annotation/
  # The net file
  rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/netDanRer7.* .
  rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/netCanFam3.* .
  rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/netEquCab2.* .
  rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/netFr2.* .
  rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/netGalGal3.* .
  rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/netGasAcu1.* .
  rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/netMm10.* .
  rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/netMonDom5.* .
  rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/netOryLat2.* .
  rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/netRn4.* .
  rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/netTetNig2.* .
  rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/netXenTro3.* .


}

#mm10
mm10()
{
  cd /export/data/CNEs/mm10/annotation/
  # The net file
  rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/netHg19.* .

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
  rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/danRer7/database/net* .

}



##ce10
cd /export/data/CNEs/ce10/annotation/
# The rmsk
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/ce10/database/rmsk.* .
# RefGene
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/ce10/database/refGene.* .
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/ce10/database/refLink.* .
# Gap
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/ce10/database/gap.* .
# 

##caeRem3
cd /export/data/CNEs/caeRem3/annotation/
# Gap
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/caeRem3/database/*_gap.* .

##caePb2
cd /export/data/CNEs/caePb2/annotation/
# Gap
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/caePb2/database/*_gap.* .

##petMar2
cd /export/data/CNEs/petMar2/annotation
# The rmsk
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/petMar2/database/rmsk.* .
# The gap
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/petMar2/database/gap.* .
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/petMar2/database/cpgIslandExt.* .

##ornAna1
cd /export/data/CNEs/ornAna1/annotation
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/ornAna1/database/rmsk.* .
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/ornAna1/database/gap.* .
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/ornAna1/database/refGene.* .
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/ornAna1/database/refLink.* .
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/ornAna1/database/cpgIslandExt.*



