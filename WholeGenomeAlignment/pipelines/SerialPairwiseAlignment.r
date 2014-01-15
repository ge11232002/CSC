selfDir = "~/Repos/CSC/WholeGenomeAlignment/scripts"
selfScripts = list.files(path=selfDir, pattern='.*\\.r', full.names=TRUE, recursive=TRUE)
for(rs in selfScripts){message(rs);source(rs)}


## first set the two genomes are near, medium or far.
distance = "near"
assemblyTarget= "/export/data/goldenpath/SC5314A21/SC5314A21.2bit"
assemblyQuery = "/export/data/goldenpath/CD36/CD36.2bit"

### step 1: Alignments with lastz 
library(rtracklayer)
validchrsTarget = seqnames(seqinfo(TwoBitFile(assemblyTarget)))
validchrsQuery = seqnames(seqinfo(TwoBitFile(assemblyQuery)))

## To check the command of lastz
#lavs = lastz(assemblyTarget=assemblyTarget, assemblyQuery=assemblyQuery, chrsTarget="chr1", chrsQuery="chr1", distance=distance, format="lav", echoCommand=TRUE)

#reg = makeRegistry(id="lastz", seed=123, file.dir=file.path(getwd(), paste(basename(assemblyQuery), "-lastz_batchjobs", sep="")), work.dir=getwd(), skip=FALSE)
#batchExpandGrid(reg, function(chrsTarget, chrsQuery, assemblyTarget, assemblyQuery, distance, format){
lavs = c()
for(chrsTarget in validchrsTarget){
  for(chrsQuery in validchrsQuery){
    lav = lastz(assemblyTarget, assemblyQuery, chrsTarget, chrsQuery, distance, format="lav")
    lavs = c(lavs, lav)
  }
}


axtNetWholePipeline(lavs, assemblyTarget, assemblyQuery, distance, removeFiles=TRUE)

