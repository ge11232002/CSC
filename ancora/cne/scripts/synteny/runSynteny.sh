#!/bin/bash

#hg19
hg19()
{
  cd /export/data/CNEs/synteny
  perl /opt/www/cne/scripts/synteny/join_nets.pl \
    danRer7 hg19 UCSC_danRer7 100000 300000
  perl /opt/www/cne/scripts/synteny/join_nets.pl \
    galGal4 hg19 UCSC_galGal4 200000 300000
  perl /opt/www/cne/scripts/synteny/join_nets.pl \
    xenTro3 hg19 UCSC_xenTro3 150000 300000
  perl /opt/www/cne/scripts/synteny/join_nets.pl \
    tetNig2 hg19 UCSC_tetNig2 100000 300000
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

#canFam3
canFam3()
{
  cd /export/data/CNEs/synteny
  perl /opt/www/cne/scripts/synteny/join_nets.pl \
    hg19 canFam3 UCSC_hg19 300000 300000


}

#equCab2
equCab2()
{
  cd /export/data/CNEs/synteny
  perl /opt/www/cne/scripts/synteny/join_nets.pl \
    hg19 equCab2 UCSC_hg19 300000 300000
}

#tetNig2
tetNig2()
{
  cd /export/data/CNEs/synteny
  perl /opt/www/cne/scripts/synteny/join_nets.pl \
    danRer7 tetNig2 UCSC_danRer7 100000 100000
}

#dm3
dm3()
{
 cd /export/data/CNEs/synteny
 perl /opt/www/cne/scripts/synteny/join_nets.pl \


}

