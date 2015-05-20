selfDir = "~/Repos/CSC/WholeGenomeAlignment/scripts"
selfScripts = list.files(path=selfDir, pattern='.*\\.r', full.names=TRUE, 
                         recursive=TRUE, ignore.case=TRUE)
for(rs in selfScripts){message(rs);source(rs)}


## first set the two genomes are near, medium or far.
distance = "medium"
targetDB= "/export/data/goldenpath/danRer7/lastdb/danRer7"
assemblyQuery = "/export/data/goldenpath/cteIde1/cteIde1.fa"
outputFn <- "danRer7.cteIde1.maf"

## last
last(targetDB, assemblyQuery, outputFn=outputFn,
     distance=distance, format="MAF", mc.cores=8L,
     echoCommand=FALSE)

## Convert maf to psl
### maf-convert DHAB.danRer7.maf > DHAB.danRer7.psl
psls <- sub("\\.maf$", ".psl", outputFn)
cmd <- paste("maf-convert", outputFn, ">", psls)
my.system(cmd)

## psl to chain
assemblyTarget = file.path(dirname(dirname(targetDB)), 
                           paste0(basename(targetDB), ".2bit"))
assemblyQuery = sub("\\.fa$", ".2bit", assemblyQuery)
chains <- sub("\\.psl$", ".chain", psls)

removeFiles = FALSE
chains = axtChain(psls, assemblyTarget, assemblyQuery, format="psl",
                  outputs=chains, distance=distance, removePsl=FALSE)
allChain = chainMergeSort(path="chain", assemblyTarget, assemblyQuery, removeChains=removeFiles)

### step 3: Netting
allPreChain = chainPreNet(allChain, assemblyTarget, assemblyQuery, removeAllChain=removeFiles)
# combine the chains into nets, add the synteny information to the net:
netSyntenicFile = chainNetSyntenic(allPreChain, assemblyTarget, assemblyQuery)

### step 4: axtNet
axtFile = netToAxt(netSyntenicFile, allPreChain, assemblyTarget, 
                   assemblyQuery, removeFiles=removeFiles)
### gzip the axtNet file
my.system(paste("gzip", axtFile))


