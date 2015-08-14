# GFF3 for ancora
This explans how to use GFF3 for gbrowse which can provide much better performance as claimed.

## GFF3 from UCSC
```{sh}
genePredToGtf danRer10 refGene refGene_danRer10.gtf
gtf2gff3 refGene_danRer10.gtf > refGene_danRer10.gff
```

Then we need to manually add the chromosoms lines for GFF3 file which is required by gbrowse.
Check the script ```addChrLines.R```.

Then you can load the GFF3 file to mysql by
```{sh}
bp_seqfeature_load.pl -u $USER -p $PASS -d gbrowse_gff_danRer10_refGene -c refGene_danRer10.gff danRer10.fasta
```





