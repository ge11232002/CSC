# RepeatMasker pipeline
This is the documentation of running RepeatMasker pipeline.<br>

## Run without job scheduler
When the genome is not assembled well and containing tens of thousands of contigs/scaffolds,
it is not necessary to partition the genome and run each piece in a single job.
```perl
/usr/local/bin/RepeatMasker/RepeatMasker -engine crossmatch  -pa 8 -align -species danio DHAB.fa
```

However, at the stage of collecting results,
the `RpeatMasker` script encounters a error
```
Get the error Can't call method "getScore" on unblessed reference 
at /usr/local/bin/RepeatMasker/PRSearchResult.pm line 159. 
after the computation.
```

## Run with job scheduler
The pipeline is designed for the usgae on olifant server with Slurm as job scheduler.
When the genome is assembled well and contains large chromosomes, 
it is better to parition the genome into smaller piecies and submit each pieces as a job.
The following is based on the UCSC [RpeatMasker pipeline](http://genomewiki.ucsc.edu/index.php/RepeatMasker).

### Setting
The script `repeatMasker/scripts/HgAutomate.pm` should be in the `PERL5LIB` path.
The alternative is to add `use lib '/mnt/biggley/home/gtan/Repos/CSC/repeatMasker/scripts';` to the 
script `src_hg_utils_automation_simplePartition.pl`.

### Partition
The script `src_hg_utils_automation_simplePartition.pl` is used to partition the genome.
It merely partitions the unmasked two Bit file genome into pieces with 500000 bases
into the folder RMPart.
```perl
perl /mnt/biggley/home/gtan/Repos/CSC/repeatMasker/scripts/src_hg_utils_automation_simplePartition.pl /export/data/goldenpath/DHAB/DHAB.unmasked.2bit 500000 RMPart
```

### Cluster run
Follow the pipeline script `pipeline.R`.

###  Collect cluster run results
```sh
liftUp DHAB.fa.out /dev/null carry RMPart/*/*/*.out
liftUp DHAB.fa.align /dev/null carry RMPart/*/*/*.align
```

###  Masking sequence
```sh
twoBitMask DHAB.unmasked.2bit ../DHAB_ncbi/DHAB.sorted.fa.out DHAB.rmsk.2bit
```


