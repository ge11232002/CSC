#!/bin/bash

#hg19
hg19()
{
  cd /export/data/CNEs/synteny
  perl /opt/www/cne/scripts/synteny/join_nets.pl \
    danRer7 hg19 UCSC_danRer7 100000 300000

}

#mm10
mm10()
{
  cd /export/data/CNEs/synteny
  perl /opt/www/cne/scripts/synteny/join_nets.pl \
    hg19 mm10 UCSC_hg19 300000 300000

}


#danRer7
danRer7()
{
  cd /export/data/CNEs/synteny
  perl /opt/www/cne/scripts/synteny/join_nets.pl \
    hg19 danRer7 UCSC_hg19 300000 100000

}







