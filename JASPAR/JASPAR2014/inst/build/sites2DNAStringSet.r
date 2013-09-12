


# filepath = "~/work/projects/JASPAR/data/sites2014"
sites2DNAStringSet = function(filepath){
  require(Biostrings)
  files = list.files(path=filepath, pattern="*\\.sites$", full.names=TRUE)
  GRL = list()
  oneFile = files[1]
  for(oneFile in files){
    sites = readDNAStringSet(filepath=oneFile, format="fasta")
    GRL = c(GRL, sites)
  }
  names(GRL) = sub("\\.sites$", "", basename(files))
  return(GRL)
}


JASPAR2014SitesSeqs = sites2DNAStringSet(filepath)
save(JASPAR2014SitesSeqs, file="/mnt/biggley/home/gtan/work/projects/JASPAR/JASPAR2014/data/JASPAR2014SitesSeqs.rda")



