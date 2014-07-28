selfDir = "~/Repos/CSC/WholeGenomeAlignment/scripts"
selfScripts = list.files(path=selfDir, pattern='.*\\.r', full.names=TRUE, recursive=TRUE)
for(rs in selfScripts){message(rs);source(rs)}


## first set the two genomes are near, medium or far.
distance = "near"
assemblyTarget= "/export/data/goldenpath/DHAB/DHAB"
assemblyQuery = "/export/data/goldenpath/ASM15152v1/ASM15152v1.2bit"

