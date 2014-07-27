selfDir = "~/Repos/CSC/repeatMasker/scripts"
selfScripts = list.files(path=selfDir, pattern='.*\\.R', full.names=TRUE, 
                         recursive=TRUE)
for(rs in selfScripts){message(rs);source(rs)}


## Compute together
#/usr/local/bin/RepeatMasker/RepeatMasker -engine crossmatch  -pa 8 -align -species danio DHAB.fa
# Get the error Can't call method "getScore" on unblessed reference at /usr/local/bin/RepeatMasker/PRSearchResult.pm line 159. after the computation.

## Cluster computation
### Partition
# perl /mnt/biggley/home/gtan/Repos/CSC/repeatMasker/scripts/src_hg_utils_automation_simplePartition.pl /export/data/goldenpath/DHAB/DHAB.unmasked.2bit 500000 RMPart

### RepeatMasker
inLstFns <- readLines("RMPart/partitions.lst")
inLstFns <- file.path(getwd(), "RMPart", inLstFns)


library(BatchJobs)
## require 1 core for each job
reg <- makeRegistry(id="repeatMasker", seed=123,
                    file.dir=file.path(getwd(), 
                                       paste("repeatMasker", "-batchjobs",
                                             sep="")),
                    work.dir=getwd(), skip=FALSE)
batchMap(reg, repeatMasker, inLstFns,
         more.args=list(species="danio", engine="ncbi",
                        binary="/usr/local/bin/RepeatMasker/RepeatMasker")
         )
ids <- getJobIds(reg)
idsChunk <- chunk(ids, n.chunk=10)
#testJob(reg, 1)
submitJobs(reg, ids=idsChunk)
showStatus(reg)

###  Collect cluster run results
#### liftUp DHAB.fa.out /dev/null carry RMPart/*/*/*.out
#### liftUp DHAB.fa.align /dev/null carry RMPart/*/*/*.align

###  Masking sequence
#### twoBitMask DHAB.unmasked.2bit ../DHAB_ncbi/DHAB.sorted.fa.out DHAB.rmsk.2bit



