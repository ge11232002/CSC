


# filepath = "~/work/projects/JASPAR/dataPackage/sites"
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
