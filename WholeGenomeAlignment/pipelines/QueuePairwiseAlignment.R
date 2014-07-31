selfDir = "~/Repos/CSC/WholeGenomeAlignment/scripts"
selfScripts = list.files(path=selfDir, pattern='.*\\.r', full.names=TRUE, 
                         recursive=TRUE, ignore.case=TRUE)
for(rs in selfScripts){message(rs);source(rs)}


## first set the two genomes are near, medium or far.
distance = "medium"
assemblyTarget= "/export/data/goldenpath/galGal4/galGal4.2bit"
assemblyQuery = "/export/data/goldenpath/hg19/hg19.2bit"

### step 1: Alignments with lastz 
# This step might be slow. run in parallel.
library(BatchJobs)
library(rtracklayer)
validchrsTarget = grep("chr", seqnames(seqinfo(TwoBitFile(assemblyTarget))), value=TRUE)
validchrsTarget = grep("_", validchrsTarget, invert=TRUE, value=TRUE)
validchrsQuery = grep("chr", seqnames(seqinfo(TwoBitFile(assemblyQuery))), value=TRUE)
validchrsQuery = grep("_", validchrsQuery, invert=TRUE, value=TRUE)

## To check the command of lastz
#lavs = lastz(assemblyTarget=assemblyTarget, assemblyQuery=assemblyQuery, chrsTarget="chr1", chrsQuery="chr1", distance=distance, format="lav", echoCommand=TRUE)
### Step 1: Do the Lastz
dir.create("lav")
setwd("lav")
reg = makeRegistry(id="lastz", seed=123, file.dir=file.path(getwd(), paste(basename(assemblyQuery), "-lastz_batchjobs", sep="")), work.dir=getwd(), skip=FALSE)
batchExpandGrid(reg, function(chrsTarget, chrsQuery, assemblyTarget, assemblyQuery, distance, format){
                selfDir = "~/Repos/CSC/WholeGenomeAlignment/scripts"
                selfScripts = list.files(path=selfDir, pattern='.*\\.r', full.names=TRUE, recursive=TRUE, ignore.case=TRUE)
                for(rs in selfScripts){source(rs)}
                lastz(assemblyTarget, assemblyQuery, chrsTarget, chrsQuery, distance, format)
                },
                chrsTarget=validchrsTarget,
                chrsQuery=validchrsQuery,
                more.args=list(assemblyTarget=assemblyTarget, 
                               assemblyQuery=assemblyQuery, distance=distance,
                               format="lav")
                )

testJob(reg, 1)
ids <- getJobIds(reg)
idsChunk <- chunk(ids, n.chunks=10, shuffle=TRUE)
submitJobs(reg)
showStatus(reg)
setwd("../")

#axtNetWholePipeline(lavs, assemblyTarget, assemblyQuery, distance, removeFiles=TRUE)

removeFiles <- FALSE
## lastz might exit without any errors. check the files whether are complete.
#lavs <- file.path("lav", lavs)
lavs <- list.files("lav", pattern=".*\\.lav", full.names=TRUE)
#filesNotCompleted = validateLastz(lavs)
#stopifnot(is.null(filesNotCompleted))

### step 2: lav to psl
dir.create("psl")
psls <- file.path("psl", sub(".lav", ".psl", basename(lavs), fixed=TRUE))
lavToPsl(lavs, psls, removeLav=removeFiles)

### step 3: Chaining
dir.create("chain")
outputs <- file.path("chain", sub(paste("\\.", "psl", "$", sep=""),
                                  ".chain", basename(psls), ignore.case=TRUE))
chains = axtChain(psls, assemblyTarget, assemblyQuery, format="psl",
                  outputs=outputs, distance=distance, removePsl=removeFiles)
allChain = chainMergeSort(path="chain", assemblyTarget, assemblyQuery, 
                          removeChains=removeFiles)

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
axtFile = netToAxt(netSyntenicFile, allPreChain, assemblyTarget, 
                   assemblyQuery, removeFiles=removeFiles)


### gzip the axtNet file
my.system(paste("gzip", axtFile))
