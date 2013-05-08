1. MGI_Coordinate.rpt downloaded from ftp://ftp.informatics.jax.org/pub/reports/index.html
2. Converted to bed file by running
  perl mgi2bed.pl <MGI_Coordinate.rpt >MGI_Coordinate.bed
3. Bed-file lifted to mm9 using the liftover tool at genome.ucsc.edu. Output saved as MGI_Coordinate_mm9.bed
4. MGI file for mm9 constructed by running:
  perl bed2mgi.pl MGI_Coordinate_mm9.bed MGI_Coordinate.rpt
