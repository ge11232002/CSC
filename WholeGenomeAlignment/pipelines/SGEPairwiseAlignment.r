selfDir = "~/Repos/genomics"
selfScripts = list.files(path=selfDir, pattern='.*\\.r', full.names=TRUE, recursive=TRUE)
for(rs in selfScripts){message(rs);source(rs)}


## first set the two genomes are near, medium or far.
distance = "far"
assemblyTarget= "/home/gt09/work/assembly/hg19/hg19.2bit"
assemblyQuery = "/home/gt09/work/assembly/fr3/fr3.2bit"

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

reg = makeRegistry(id="lastz", seed=123, file.dir=file.path(getwd(), paste(basename(assemblyQuery), "-lastz_batchjobs", sep="")), work.dir=getwd(), skip=FALSE)
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

## Rerun the jobs automatically which expire for some reasons...Usually it happens when it hits the wall clock. But seems not to be the case on alpha.
monitorBatchJobs(reg)

lavs = reduceResultsVector(reg, fun=function(job, res) res)
stopifnot(length(lavs) == showStatus(reg)[1,"started"])

axtNetWholePipeline(lavs, assemblyTarget, assemblyQuery, distance, removeFiles=TRUE)

