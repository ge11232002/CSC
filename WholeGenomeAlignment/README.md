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

## Installation and Configuration
List of softwares need to be installed:

* [LASTZ](http://www.bx.psu.edu/miller_lab/dist/README.lastz-1.02.00/README.lastz-1.02.00a.html)
* [UCSC utilities](http://hgdownload.cse.ucsc.edu/admin/jksrc.zip). 
  In this pipeline, lavToPsl, axtChain, chainMergeSort, chainPreNet, chainNet,
  netSyntenic, netToAxt, axtSort are essential.
* R and the packages: [BatchJobs](http://cran.r-project.org/web/packages/BatchJobs/index.html), 
[rtracklayer](http://www.bioconductor.org/packages/release/bioc/html/rtracklayer.html)
* ```git clone git@github.com:ge11232002/CSC.git```
* ```git clone git@github.com:ge11232002/genomics.git```

Configurations:
* Put the installed lastz under ```~/bin/lastz``` 
  and UCSC utilities under  ```~/bin/x86_64```.
  If you have a different location, 
  please create a symbolic link to the path mentioned above.
* The script for running the pipeline is under ```CSC/WholeGenomeAlignment/pipelines/genomePairwiseAlignment.r```.
  Other support scripts used by the pipeline are defined in the repository ```genomics```.
  All the support scripts are sourced by the ```genomePairwiseAlignment.r``` script first.
  Put a copy under ```~/Repos/genomics```. 
* Put ```CSC/WholeGenomeAlignment/pipelines/SGETemplate.tmpl``` under ```~/```.
  This is the SGE configuration template used by the R package BatchJobs.
  Modify it to fit your needs on specific server.
* Put ```CSC/WholeGenomeAlignment/pipelines/BatchJobs.R``` to ```~/.BatchJobs.R```.
  Similarly, this R script is sourced when the package BatchJobs is loaded.
  It determines the job schedulder used on the server.

## Running the pipeline
The whole pipeline starts with two soft-masked 2bit files, 
one as target genome, another as query genome. 
The output is the axtNet file.
Later, perhaps we will include the procedure of repeat maskers.

All we need to do is to set the paths of the two 2bit files 
and the distance of two species.
The distance can be "near", "medium" and "far". 
It will determine the parameters used in later stages.
Then follow all the functions in the script ```genomePairwiseAlignment.r```.

## Pipeline details
To reduce the number of chromosome comparison pairs 
and since we are interested in ordinary chromosomes (e.g. chr\*),
the chromosomes names contain the "_" symbol are skipped.



