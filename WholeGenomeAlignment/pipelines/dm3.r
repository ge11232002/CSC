selfDir = "~/Repos/CSC/WholeGenomeAlignment/scripts"
selfScripts = list.files(path=selfDir, pattern='.*\\.r', full.names=TRUE, recursive=TRUE)
for(rs in selfScripts){message(rs);source(rs)}


## first set the two genomes are near, medium or far.
distance = "far"
assemblyTarget= "/export/data/goldenpath/dp4/dp4.2bit"
assemblyQuery = "/export/data/goldenpath/dm3/dm3.2bit"

### step 1: Alignments with lastz 
library(rtracklayer)
validchrsTarget = seqnames(seqinfo(TwoBitFile(assemblyTarget)))
validchrsQuery = seqnames(seqinfo(TwoBitFile(assemblyQuery)))

## To check the command of lastz
#lavs = lastz(assemblyTarget=assemblyTarget, assemblyQuery=assemblyQuery, chrsTarget="chr1", chrsQuery="chr1", distance=distance, format="lav", echoCommand=TRUE)

#reg = makeRegistry(id="lastz", seed=123, file.dir=file.path(getwd(), paste(basename(assemblyQuery), "-lastz_batchjobs", sep="")), work.dir=getwd(), skip=FALSE)
#batchExpandGrid(reg, function(chrsTarget, chrsQuery, assemblyTarget, assemblyQuery, distance, format){
lavs = c()
combinations <- expand.grid(validchrsTarget, validchrsQuery, 
                            stringsAsFactors=FALSE)
library(BiocParallel)
parallelFunc <- function(chrsTarget, chrsQuery, assemblyTarget, assemblyQuery, distance, format="lav"){
  lastz(assemblyTarget, assemblyQuery, chrsTarget, chrsQuery, distance, format=format)
}
multicoreParam <- MulticoreParam(workers=12)
### Step 1: Do the Lastz
dir.create("lav")
setwd("lav")
lavs <- bpmapply(parallelFunc, combinations$Var1, combinations$Var2, 
                 MoreArgs=list(assemblyTarget=assemblyTarget,
                               assemblyQuery=assemblyQuery,
                               distance=distance, format="lav"),
                 BPPARAM=multicoreParam)
setwd("../")

removeFiles <- FALSE
## lastz might exit without any errors. check the files whether are complete.
lavs <- file.path("lav", lavs)
filesNotCompleted = validateLastz(lavs)
stopifnot(is.null(filesNotCompleted))

### step 2: lav to psl 
dir.create("psl")
psls <- file.path("psl", sub("\\.lav$", ".psl", basename(lavs), 
                             ignore.case=TRUE))
psls = lavToPsl(lavs, psls, removeLav=removeFiles)

### step 3: Chaining
dir.create("chain")
outputs <- file.path("chain", sub(paste("\\.", "psl", "$", sep=""), 
                                  ".chain", basename(psls), ignore.case=TRUE))
chains = axtChain(psls, assemblyTarget, assemblyQuery, format="psl", 
                  outputs=outputs, distance=distance, removePsl=removeFiles)
allChain = chainMergeSort(path="chain", assemblyTarget, assemblyQuery, removeChains=removeFiles)

# now , this allChain is what we get from UCSC download e.g. hg19.danRer7.all.chain.gz

### step 3: Netting
  allPreChain = chainPreNet(allChain, assemblyTarget, assemblyQuery, removeAllChain=removeFiles)
# combine the chains into nets, add the synteny information to the net:
  netSyntenicFile = chainNetSyntenic(allPreChain, assemblyTarget, assemblyQuery)

#  Add classification information using the database tables:
# actually this step is not necessary in this pipeline according to http://blog.gmane.org/gmane.science.biology.ucscgenome.general/month=20130301. The class information will only be used for Genome Browser.
# Since it needs some specific modification of the table names for certain species, we skip this step now.
# If this step is done, then the generated class.net is the gzipped net file that you see in UCSC Downloads area.

### step 4: axtNet
  axtFile = netToAxt(netSyntenicFile, allPreChain, assemblyTarget, assemblyQuery, removeFiles=removeFiles)
### gzip the axtNet file
  my.system(paste("gzip", axtFile))

