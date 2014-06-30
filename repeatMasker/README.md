# [RepeatMasker pipeline](http://genomewiki.ucsc.edu/index.php/RepeatMasker)

This is the documentation of running RepeatMasker pipeline.<br>
The pipeline is designed for the usgae on olifant server with Slurm as job scheduler.

## Setting
The script `repeatMasker/scripts/HgAutomate.pm` should be in the `PERL5LIB` path.
The alternative is to add `use lib '/mnt/biggley/home/gtan/Repos/CSC/repeatMasker/scripts';` to the 
script `src_hg_utils_automation_simplePartition.pl`.

## Partition
The script `src_hg_utils_automation_simplePartition.pl` is used to partition the genome.
It merely partitions the unmasked two Bit file genome into pieces with 500000 bases
into the folder RMPart.
```perl
perl /mnt/biggley/home/gtan/Repos/CSC/repeatMasker/scripts/src_hg_utils_automation_simplePartition.pl /export/data/goldenpath/DHAB/DHAB.unmasked.2bit 500000 RMPart
```

## Cluster run
Follow the pipeline script `pipeline.R`.


