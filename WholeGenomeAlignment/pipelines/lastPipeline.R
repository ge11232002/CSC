selfDir = "~/Repos/CSC/WholeGenomeAlignment/scripts"
selfScripts = list.files(path=selfDir, pattern='.*\\.r', full.names=TRUE, 
                         recursive=TRUE, ignore.case=TRUE)
for(rs in selfScripts){message(rs);source(rs)}


## first set the two genomes are near, medium or far.
distance = "near"
targetDB= "/export/data/goldenpath/DHAB/DHAB"
assemblyQuery = "/export/data/goldenpath/danRer7/danRer7.fa"

last(targetDB, assemblyQuery, outputFn="DHAB.danRer7.maf",
     distance="near", format="MAF", mc.cores=1L)

## Convert maf to psl
### maf-convert.py psl DHAB.danRer7.maf > DHAB.danRer7.psl


## split psl
### split --lines=9767336 DHAB.danRer7.psl

## psl to chain
assemblyTarget = "/export/data/goldenpath/DHAB/DHAB.2bit"
assemblyQuery = "/export/data/goldenpath/danRer7/danRer7.2bit"
outputs = "DHAB.danRer7.chain"
psls = list.files(path=".", pattern="x.*")
outputs = paste0(psls, ".chain")

chains = axtChain(psls, assemblyTarget, assemblyQuery, format="psl",
                  outputs=outputs, distance=distance, removePsl=FALSE)


