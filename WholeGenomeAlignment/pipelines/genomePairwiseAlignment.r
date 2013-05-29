selfDir = "~/Repos/genomics"
selfScripts = list.files(path=selfDir, pattern='.*\\.r', full.names=TRUE, recursive=TRUE)
for(rs in selfScripts){message(rs);source(rs)}


## first set the two genomes are near, medium or far.
distance = "near"
assemblyTarget= "/home/gt09/work/assembly/ce10/ce10.2bit"
assemblyQuery = "/home/gt09/work/assembly/cb3/cb3.2bit"

### step 1: Alignments with lastz 
# This step might be slow. run in parallel.
library(BatchJobs)
library(rtracklayer)
validchrsTarget = grep("_", seqnames(seqinfo(TwoBitFile(assemblyTarget))), invert=TRUE, value=TRUE)
validchrsQuery = grep("_", seqnames(seqinfo(TwoBitFile(assemblyQuery))), invert=TRUE, value=TRUE)
reg = makeRegistry(id="lastz", seed=123, file.dir=file.path(getwd(), "lastz_batchjobs"), work.dir=getwd(), skip=FALSE)
batchExpandGrid(reg, function(chrsTarget, chrsQuery, assemblyTarget, assemblyQuery, distance, format){
                selfDir = "~/Repos/genomics"
                selfScripts = list.files(path=selfDir, pattern='.*\\.r', full.names=TRUE, recursive=TRUE)
                for(rs in selfScripts){source(rs)}
                lastz(assemblyTarget, assemblyQuery, chrsTarget, chrsQuery, distance, format)
                },
                chrsTarget=validchrsTarget,
                chrsQuery=validchrsQuery,
                more.args=list(assemblyTarget=assemblyTarget, 
                               assemblyQuery=assemblyQuery, distance=distance,
                               format="lav")
                )
submitJobs(reg)
showStatus(reg)
lavs = reduceResultsVector(reg, fun=function(job, res) res)
#lavs = lastz(assemblyTarget=assemblyTarget, assemblyQuery=assemblyQuery, chrsTarget="chr1", chrsQuery="chr1", distance=distance, format="lav", echoCommand=TRUE)
filesNotCompleted = validateLastz(lavs)
stopifnot(is.null(filesNotCompleted))

### step 2: Chaining 
psls = lavToPsl(lavs, removeLav=TRUE)
chains = axtChain(psls, assemblyTarget, assemblyQuery, format="psl", distance=distance, removePsl=TRUE)
allChain = chainMergeSort(chains, assemblyTarget, assemblyQuery, removeChains=TRUE)
# now , this allChain is what we get from UCSC download e.g. hg19.danRer7.all.chain.gz

### step 3: Netting 
allPreChain = chainPreNet(allChain, assemblyTarget, assemblyQuery, removeAllChain=FALSE)
# combine the chains into nets, add the synteny information to the net:
netSyntenicFile = chainNetSyntenic(allPreChain, assemblyTarget, assemblyQuery)

#  Add classification information using the database tables:
# actually this step is not necessary in this pipeline according to http://blog.gmane.org/gmane.science.biology.ucscgenome.general/month=20130301. The class information will only be used for Genome Browser.
# Since it needs some specific modification of the table names for certain species, we skip this step now.
# If this step is done, then the generated class.net is the gzipped net file that you see in UCSC Downloads area.


### step 4: axtNet
axtFile = netToAxt(netSyntenicFile, allPreChain, assemblyTarget, assemblyQuery)



