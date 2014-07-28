#Whole genome pairwise alignment
This documentation describes the whole pipeline of 
generating the pairwise alignemnt between two different species genomes.
This pipeline is based on the 
[genomewiki](http://genomewiki.ucsc.edu/index.php/Whole_genome_alignment_howto) 
from UCSC, one [newsletter](https://lists.soe.ucsc.edu/pipermail/genome/2005-October/008710.html).
The parameters of lastz are chosen by the suggestions of the automation scripts
[doBlastzChainNet.pl](http://hgwdev.cse.ucsc.edu/~kent/src/unzipped/hg/utils/automation/doBlastzChainNet.pl)
and 
[RunLastzChain sh.txt](http://genomewiki.ucsc.edu/images/9/93/RunLastzChain_sh.txt).
This pipeline is implemented in R and 
designed for CSC cluster "alpha" with a scheduler Sun Grid Engine (SGE).

The parameters for lastz should refer to [ucsc](http://genomewiki.ucsc.edu/index.php/UCSC_Multiple_Alignments).

## Installation and Configuration
List of softwares need to be installed:

* [LASTZ](http://www.bx.psu.edu/miller_lab/dist/README.lastz-1.02.00/README.lastz-1.02.00a.html)
* [UCSC utilities](http://hgdownload.cse.ucsc.edu/admin/jksrc.zip). 
  In this pipeline, lavToPsl, axtChain, chainMergeSort, chainPreNet, chainNet,
  netSyntenic, netToAxt, axtSort are essential. 
  netClass is optional.
* R and the packages: [BatchJobs](http://cran.r-project.org/web/packages/BatchJobs/index.html), 
[rtracklayer](http://www.bioconductor.org/packages/release/bioc/html/rtracklayer.html)
* ```git clone git@github.com:ge11232002/CSC.git```

Configurations:
* Put the installed lastz under ```~/bin/lastz``` 
  and UCSC utilities under  ```~/bin/x86_64```.
  If you have a different location, 
  please create a symbolic link to the path mentioned above.
* The script for running the pipeline is under ```CSC/WholeGenomeAlignment/pipelines/SGEPairwiseAlignment.r```.
  Other support scripts used by the pipeline are in ```CSC/WholeGenomeAlignment/scripts```.
* Put ```CSC/WholeGenomeAlignment/pipelines/SGETemplate.tmpl``` under ```~/```.
  This is the SGE configuration template used by the R package BatchJobs.
  Modify it to fit your needs on specific server.
* Put ```CSC/WholeGenomeAlignment/pipelines/BatchJobs.R``` to ```~/.BatchJobs.R```.
  Similarly, this R script is sourced when the package BatchJobs is loaded.
  It determines the job schedulder used on the server.

## Running the pipeline
The whole pipeline starts with two soft-masked 2bit files, 
one as target genome, another as query genome. 
The output is the axtNet file which can usually be downloaded from [UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/vsMm10/axtNet/).

All we need to do is to set the paths of the two 2bit files 
and the distance of two species.
The distance can be "near", "medium" and "far". 
It will determine the parameters used in later stages.
Then follow all the functions in the script ```SGEPairwiseAlignment.r```.

## Pipeline details
To reduce the number of chromosome comparison pairs 
and since we are interested in ordinary chromosomes (e.g. chr\*),
the chromosomes names contain the "_" symbol are skipped.

The distance between two compared assemblies can refer to 
the [UCSC Multiple Alignments](http://genomewiki.ucsc.edu/index.php/UCSC_Multiple_Alignments).
and [UCSC topologies](http://genomewiki.ucsc.edu/index.php/Phylogenetic_Tree#Vertebrate_topology_used_at_UCSC_genome_browser).

## Other useful outputs from this pipeline
### hg19.galGal4.all.chain
This is the allChai. file you can get from UCSC download.

Example: [hg19.mm10.all.chain.gz](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/vsMm10/hg19.mm10.all.chain.gz)

### hg19.galGal4.all.pre.chain
This is the filtered file from allChain. 
The intention behind them is that there are no duplicates.
This pre.chain file is actually the liftOver files which can be downloaded from UCSC too.

Exmaple: [hg19ToMm10.over.chain.gz](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToMm10.over.chain.gz).

To load the chain file into local database,
```sh
hgLoadChain hg19 chainGalGal4 hg19.galGal4.all.pre.chain
```

### hg19.galGal4.noClass.net and hg19.galGal4.net
During the pipeline, an optional step can be used to "Add classification information using the database tables".
This step is only necessary for generating UCSC Genome Browser track
and has no impact on the following steps for axt file, according to
this [post](http://blog.gmane.org/gmane.science.biology.ucscgenome.general/month=20130301)

To run this step,
```sh
 netClass -noAr hg19.galGal4.noClass.net UCSC_hg19 UCSC_galGal4 hg19.galGal4.net
```

The resulting hg19.galGal4.net file is in the UCSC pairwise alignment download section.

Example: [hg19.mm10.net.gz](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/vsMm10/hg19.mm10.net.gz).

To load the net file into local database,
```sh
hgLoadNet UCSC_hg19 netGalGal4 hg19.galGal4.net
```


