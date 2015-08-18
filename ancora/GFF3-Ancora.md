# GFF3 for Ancora Browser
This explans how to use GFF3 for gbrowse which can provide much better performance as claimed.

## refGene from UCSC
```{sh}
genePredToGtf danRer10 refGene refGene.danRer10.gtf
gtf2gff3 refGene.danRer10.gtf > refGene.danRer10.gff3
```

Then we need to manually add the chromosoms lines for GFF3 file which is required by gbrowse.
Use the script __addChrLines.R__ to generate the chroms lines __chroms.gff3__.

## knownGene from UCSC
```{sh}
genePredToGtf hg38 knownGene knownGene.hg38.gtf
```

## ensGene from Ensembl
Download the GFF3 file from Ensembl [ftp](http://www.ensembl.org/info/data/ftp/index.html).
Use the script __processEnsemblGFF3.R__ to generate the __Danio_rerio.GRCz10.81.gbrowse.gff3__.


## CpG from UCSC
Use the script __CpG.R__ to generate the file __CpG.gff3__.

## rmsk from UCSC
Use the script __rmsk.R__ to generate the file __rmsk.gff3__.

## gap from UCSC
Use the script __gap.R__ to generate the file __gap.gff3__.


## Load the GFF3 file into mysql database
```{sh}
mysql -u root -p -e 'create database gbrowse_gff_danRer10'
mysql -u root -p -e 'grant select on gbrowse_gff_danRer10.* to nobody@localhost'
```

Merge the GFF3 files into one

```{sh}
cat refGene.danRer10.gff3 chroms.gff3 Danio_rerio.GRCz10.81.gbrowse.gff3 CpG.gff3 rmsk.gff3 gap.gff3 > gbrowse_gff_danRer10.gff3
```

One buggy thing from Bioconductor package __rtracklayer__ is that the attributes of GFF3 file from ```export.gff3``` have space before ";", which causes problem in GBrowse.
To fix it,

```{sh}
perl -i -pe's/;\ /;/g' gbrowse_gff_danRer10.gff3
```

Then you can load the GFF3 file to mysql by

```{sh}
bp_seqfeature_load.pl -u $USER -p $PASS -d gbrowse_gff_danRer10 -c gbrowse_gff_danRer10.gff3 (danRer10.fa)
```

Loading the genome sequences is not mandatory. It loaded, you can examine the sequences and submit the regions of sequences into BLAST or BLAT.








