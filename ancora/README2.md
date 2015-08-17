# GFF3 for ancora
This explans how to use GFF3 for gbrowse which can provide much better performance as claimed.

## GFF3 from UCSC
```{sh}
genePredToGtf danRer10 refGene refGene.danRer10.gtf
gtf2gff3 refGene.danRer10.gtf > refGene.danRer10.gff3
```

Then we need to manually add the chromosoms lines for GFF3 file which is required by gbrowse.
Use the script ```addChrLines.R``` to generate the chroms lines ```chroms.gff3```.

## GFF3 from Ensembl
Download the GFF3 file from Ensembl [ftp](http://www.ensembl.org/info/data/ftp/index.html).
Use the script ```processEnsemblGFF3.R``` to generate the ```Danio_rerio.GRCz10.81.gbrowse.gff3```.

## Load the GFF3 file into mysql database
```{sh}
mysql -u root -p -e 'create database gbrowse_gff_danRer10'
mysql -u root -p -e 'grant select on gbrowse_gff_danRer10.* to nobody@localhost'
```

Merge the GFF3 files into one
```{sh}
cat refGene.danRer10.gff3 chroms.gff3 Danio_rerio.GRCz10.81.gbrowse.gff3 > gbrowse_gff_danRer10.gff3
```

Then you can load the GFF3 file to mysql by
```{sh}
bp_seqfeature_load.pl -u $USER -p $PASS -d gbrowse_gff_danRer10 -c gbrowse_gff_danRer10.gff3 danRer10.fa
```
