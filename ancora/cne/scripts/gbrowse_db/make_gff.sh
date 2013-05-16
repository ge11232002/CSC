#!/bin/bash

#hg19
hg19()
{
  perl /opt/www/cne/scripts/gbrowse_db/ucsc2gff.pl \
    -a /export/data/goldenpath/hg19/assembly.2bit \
    -d UCSC_hg19 \
    assembly refGene knownGene oreganno cpgIsland rmsk gap \
    >/export/data/CNEs/gff/hg19.gff
  wc -l /export/data/CNEs/gff/hg19.gff
  perl /opt/www/cne/scripts/gbrowse_db/ens2gff.pl \
    /export/data/goldenpath/hg19/assembly.2bit \
    /export/data/CNEs/hg19/annotation/ensembl_genes.txt \
    >>/export/data/CNEs/gff/hg19.gff
  wc -l /export/data/CNEs/gff/hg19.gff
  perl /opt/www/cne/scripts/gbrowse_db/mirbase2gff.pl \
    /export/data/goldenpath/hg19/assembly.2bit \
    /export/data/CNEs/hg19/annotation/hsa.gff3 \
    >>/export/data/CNEs/gff/hg19.gff
  wc -l /export/data/CNEs/gff/hg19.gff
  perl /opt/www/cne/scripts/gbrowse_db/synteny2gff.pl jn danRer7 \
    /export/data/CNEs/synteny/joinedNets_danRer7_hg19_gap100k300k.txt \
    /export/data/CNEs/synteny/joinedNets_danRer7_hg19_gap100k300k.hg19_detail.bed \
    >>/export/data/CNEs/gff/hg19.gff
  wc -l /export/data/CNEs/gff/hg19.gff
}

#mm10
mm10()
{
  perl /opt/www/cne/scripts/gbrowse_db/ucsc2gff.pl \
    -a /export/data/goldenpath/mm10/assembly.2bit \
    -d UCSC_mm10 \
    assembly refGene knownGene cpgIsland rmsk gap \
    >/export/data/CNEs/gff/mm10.gff
  wc -l /export/data/CNEs/gff/mm10.gff
  perl /opt/www/cne/scripts/gbrowse_db/ens2gff.pl \
    /export/data/goldenpath/mm10/assembly.2bit \
    /export/data/CNEs/mm10/annotation/ensembl_genes.txt \
    >>/export/data/CNEs/gff/mm10.gff
  wc -l /export/data/CNEs/gff/mm10.gff
  perl /opt/www/cne/scripts/gbrowse_db/mirbase2gff.pl \
    /export/data/goldenpath/mm10/assembly.2bit \
    /export/data/CNEs/mm10/annotation/mmu.gff3 \
    >>/export/data/CNEs/gff/mm10.gff
  wc -l /export/data/CNEs/gff/mm10.gff
  perl /opt/www/cne/scripts/gbrowse_db/mgi2gff.pl \
    /export/data/goldenpath/mm10/assembly.2bit \
    /export/data/CNEs/mm10/annotation/MGI_Gene_Model_Coord.rpt \
    >>/export/data/CNEs/gff/mm10.gff
  wc -l /export/data/CNEs/gff/mm10.gff
  perl /opt/www/cne/scripts/gbrowse_db/synteny2gff.pl jn hg19 \
    /export/data/CNEs/synteny/joinedNets_hg19_mm10_gap300k300k.txt \
    /export/data/CNEs/synteny/joinedNets_hg19_mm10_gap300k300k.mm10_detail.bed \
    >>/export/data/CNEs/gff/mm10.gff
  wc -l /export/data/CNEs/gff/mm10.gff
}

#canFam3
canFam3()
{
  perl /opt/www/cne/scripts/gbrowse_db/ucsc2gff.pl \
    -a /export/data/goldenpath/canFam3/assembly.2bit \
    -d UCSC_canFam3 \
    assembly refGene cpgIsland rmsk gap \
    >/export/data/CNEs/gff/canFam3.gff
  wc -l /export/data/CNEs/gff/canFam3.gff
  perl /opt/www/cne/scripts/gbrowse_db/ens2gff.pl \
    /export/data/goldenpath/canFam3/assembly.2bit \
    /export/data/CNEs/canFam3/annotation/ensembl_genes.txt \
    >>/export/data/CNEs/gff/canFam3.gff

}



# hg18
#hg18()
#{
#    perl ucsc2gff.pl \
#	-a /export/data/goldenpath/hg18/assembly.2bit \
#	-d UCSC_HS_MAR06 \
#	assembly refGene knownGene oreganno cpgIsland rmsk gap \
#	>hg18.gff
#    perl ens2gff.pl \
#	/export/data/goldenpath/hg18/assembly.2bit \
#	ens46-human-hg18-mart.txt \
#	>>hg18.gff
#    perl mirbase2gff.pl \
#	/export/data/goldenpath/hg18/assembly.2bit mirbase/mb10.1/hsa.gff \
#	>>hg18.gff
#    perl synteny2gff.pl jn danRer4 \
#	/home/davidf/projects/cne/dr/synteny/joinedNets_danRer4_hg18_gap100k300k.txt \
#	/home/davidf/projects/cne/dr/synteny/joinedNets_danRer4_hg18_gap100k300k.hg18_detail.bed \
#	>>hg18.gff
#    perl synteny2gff.pl jn danRer5 \
#	/home/davidf/projects/cne/dr/synteny/joinedNets_danRer5_hg18_gap100k300k.txt \
#	/home/davidf/projects/cne/dr/synteny/joinedNets_danRer5_hg18_gap100k300k.hg18_detail.bed \
#	>>hg18.gff
#}
#
## mm8
#mm8() 
#{
#    perl ucsc2gff.pl \
#	-a /export/data/goldenpath/mm8/assembly.2bit \
#	-d UCSC_MM_MAR06 \
#	assembly refGene knownGene cpgIsland rmsk gap \
#	>mm8.gff
#    perl ens2gff.pl \
#	/export/data/goldenpath/mm8/assembly.2bit \
#	ens46-mouse-mm8-mart.txt \
#	>>mm8.gff   
#    perl mgi2gff.pl \
#	/export/data/goldenpath/mm8/assembly.2bit \
#	mgi/mgi3.54/MGI_Coordinate.rpt \
#        >>mm8.gff
#    perl mirbase2gff.pl \
#	/export/data/goldenpath/mm8/assembly.2bit mirbase/mb10.0/mmu.gff \
#	>>mm8.gff
#}
#
#
## mm9
#mm9() 
#{
#    perl ucsc2gff.pl \
#	-a /export/data/goldenpath/mm9/assembly.2bit \
#	-d UCSC_mm9 \
#	assembly refGene knownGene cpgIsland rmsk gap \
#	>mm9.gff
#    perl ens2gff.pl \
#	/export/data/goldenpath/mm9/assembly.2bit \
#	ens48-mouse-mm9-mart.txt \
#	>>mm9.gff   
#    perl mgi2gff.pl \
#	/export/data/goldenpath/mm9/assembly.2bit \
#	mgi/mgi3.54/MGI_Coordinate_mm9.rpt \
#        >>mm9.gff
#    perl mirbase2gff.pl \
#	/export/data/goldenpath/mm9/assembly.2bit mirbase/mb10.1/mmu.gff \
#	>>mm9.gff
#}
#
#
## danRer4
#danRer4()
#{
#    perl ucsc2gff.pl \
#	-a /export/data/goldenpath/danRer4/assembly.2bit \
#	-d UCSC_DR_MAR06 \
#	assembly refGene "blastHg18KG:UCSC_HS_MAR06:Aligned human protein" zfin:zfin/entrezgene.txt rmsk gap \
#	>danRer4.gff
#    perl ens2gff.pl \
#	/export/data/goldenpath/danRer4/assembly.2bit \
#	ens45-zebrafish-danRer4-mart.txt \
#	>>danRer4.gff
#    perl mirbase2gff.pl \
#	/export/data/goldenpath/danRer4/assembly.2bit mirbase/mb10.0/dre.gff \
#	>>danRer4.gff
#    perl synteny2gff.pl jn hg18 \
#        ../../dr/synteny/joinedNets_danRer4_hg18_gap100k300k.txt \
#	../../dr/synteny/joinedNets_danRer4_hg18_gap100k300k.danRer4_detail.bed \
#	>>danRer4.gff
#}
#
## danRer5
#danRer5()
#{
#    perl ucsc2gff.pl \
#	-a /export/data/goldenpath/danRer5/assembly.2bit \
#	-d UCSC_danRer5 \
#	assembly refGene "blastHg18KG:UCSC_HS_MAR06:Aligned human protein" zfin:zfin/entrezgene.txt rmsk gap \
#	>danRer5.gff
#    perl ens2gff.pl \
#	/export/data/goldenpath/danRer5/assembly.2bit \
#	ens47-zebrafish-danRer5-mart.txt \
#	>>danRer5.gff
#    perl mirbase2gff.pl \
#	/export/data/goldenpath/danRer5/assembly.2bit mirbase/mb10.1/dre.gff \
#	>>danRer5.gff
#    perl synteny2gff.pl jn hg18 \
#        /home/davidf/projects/cne/dr/synteny/joinedNets_danRer5_hg18_gap100k300k.txt \
#	/home/davidf/projects/cne/dr/synteny/joinedNets_danRer5_hg18_gap100k300k.danRer5_detail.bed \
#	>>danRer5.gff
#}
#
## tetNig1
#tetNig1()
#{
#    perl ucsc2gff.pl \
#	-a /export/data/goldenpath/tetNig1/assembly.2bit \
#	-d UCSC_TN_FEB04 \
#	assembly rmsk gap \
#	>tetNig1.gff
#    perl ens2gff.pl \
#	/export/data/goldenpath/tetNig1/assembly.2bit \
#	ens48-tetraodon-tetNig1-mart.txt \
#	>>tetNig1.gff
#    perl mirbase2gff.pl \
#	/export/data/goldenpath/tetNig1/assembly.2bit mirbase/mb10.1/tni.gff \
#	>>tetNig1.gff
#}
#
#
## dm3
#dm3()
#{
#    perl ucsc2gff.pl \
#	-a /export/data/goldenpath/dm3/assembly.2bit \
#	-d UCSC_dm3 \
#	assembly refGene oreganno rmsk gap \
#	>dm3.gff
#    perl extract_genes_from_gff3.pl flybase flybase/fb5.2/dmel-all-r5.2-noSeq.gff >>dm3.gff
#    perl redfly2gff.pl /export/data/goldenpath/dm3/assembly.2bit redfly/redfly_crm.gff CRM >>dm3.gff
#    perl redfly2gff.pl /export/data/goldenpath/dm3/assembly.2bit redfly/redfly_tfbs.gff TFBS >>dm3.gff
#    perl synteny2gff.pl -g 100 cns droAna3  ../../dm3/synteny/vsDrosophilas/droAna3/dm3.droAna3.synteny.bed >>dm3.gff
#    perl synteny2gff.pl -g 100 cns dp4  ../../dm3/synteny/vsDrosophilas/dp4/dm3.dp4.synteny.bed  >>dm3.gff
#    perl synteny2gff.pl -g 100 cns droVir3  ../../dm3/synteny/vsDrosophilas/droVir3/dm3.droVir3.synteny.bed  >>dm3.gff
#    perl synteny2gff.pl -g 100 cns droMoj3  ../../dm3/synteny/vsDrosophilas/droMoj3/dm3.droMoj3.synteny.bed  >>dm3.gff
#}
#
#
## ce4 (WS170)
#ce4()
#{
#    perl ucsc2gff.pl \
#	-a /export/data/goldenpath/ce4/assembly.2bit \
#	-d UCSC_ce4 \
#	assembly refGene "blastHg18KG:UCSC_HS_MAR06:Aligned human protein" rmsk gap \
#	>ce4.gff
#    # ce4 = Wormbase Release WS170. WS176 is the latest release without changes to the assembly since WS170.
#    perl extract_genes_from_gff3.pl wormbase \
#	wormbase/WS176/elegansWS176.gff3 \
#	wormbase/WS176/geneIDs.WS176 \
#	>>ce4.gff
#}
#
#
#danRer4
#danRer5
#tetNig1
#hg18
#mm8
#mm9
#dm3
#ce4
