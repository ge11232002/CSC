


filepath = "~/work/projects/JASPAR/data/sites2014"
sites2DNAStringSet = function(filepath){
  require(Biostrings)
  files = list.files(path=filepath, pattern="*\\.sites$", full.names=TRUE)
  GRL = list()
  oneFile = files[1]
  for(oneFile in files){
    message(oneFile)
    #sites = readDNAStringSet(filepath=oneFile, format="fasta")
    # Since there are "X" code in sites sequences, it conflicts with DNA_ALPAHBET in Biostrings. Just use BString here.
    sites = readBStringSet(filepath=oneFile, format="fasta")
    # in the original sites sequences, some of them use Xx to represent A,C,G or T, rather than Nn. Xx is not valid in DNAStringSet, so translate it first.
    sites = chartr("Xx", "Nn", sites)
    sites = DNAStringSet(sites)
    GRL = c(GRL, sites)
  }
  names(GRL) = sub("\\.sites$", "", basename(files))
  return(GRL)
}


JASPAR2014SitesSeqs = sites2DNAStringSet(filepath)
save(JASPAR2014SitesSeqs, file="/mnt/biggley/home/gtan/work/projects/JASPAR/JASPAR2014/data/JASPAR2014SitesSeqs.rda")



