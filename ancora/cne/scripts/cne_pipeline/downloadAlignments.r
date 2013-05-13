

downloadUCSCPA = function(build1, build2, outputDir=NULL){
  goldenPath = "rsync://hgdownload.cse.ucsc.edu/goldenPath/"
  outputDirBuild1 = file.path(outputDir, build1)
  if(!file.exists(outputDirBuild1)){
    dir.create(outputDirBuild1)
  }
  if(build1 == build2){
    build1URL = paste(goldenPath, build1, "/vsSelf/axtNet/*", sep="")
  }else{
    upperBuild2 = paste(toupper(substring(build2,1,1)), substring(build2,2), sep="")
    build1URL = paste(goldenPath, build1, "/vs", upperBuild2, "/axtNet/*", sep="")
    build1URL2 = paste(goldenPath, build1, "/vs", upperBuild2, "/", build1, ".", build2, ".net.axt.gz", sep="")
  }
  cmd = paste("rsync -avzP", build1URL, outputDirBuild1)
  message(cmd)
  res1 = try(system(cmd))
  cmd = paste("rsync -avzP", build1URL2, outputDirBuild1)
  message(cmd)
  res2 = try(system(cmd))
  if(class(res) == "try-error" | class(res2) == "try-error"){
    return(paste(build1, build2))
  }else{
    return("Success")
  }
}

builds = list("hg19"=c("mm10", "xenTro3", "tetNig2", "canFam3", "galGal3", "danRer7", "fr3", "danRer6", "equCab2", "danRer5", "anoCar2"),
              "hg18"=c("mm9", "ornAna1", "danRer5", "canFam2", "galGal3", "tetNig1", "monDom4", "xenTro2"),
              "mm10"=c("hg19", "galGal4", "danRer7", "fr3", "danRer6"),
              "mm9"=c("hg18", "canFam3", "equCab2", "danRer7", "danRer6", "danRer5", "tetNig2"),
              "canFam2"=c("hg18", "equCab2"),
              "canFam3"=c("hg19", "mm9"),
              "equCab2"=c("hg19", "mm9", "galGal3", "canFam2"),
              "galGal3"=c("hg19", "hg18", "equCab2"),
              "galGal4"=c("mm10"),
              "danRer7"=c("hg19", "mm10", "mm9", "oryLat2", "fr3", "tetNig2", "gasAcu1"),
              "danRer6"=c("hg19", "mm10", "mm9", "oryLat2", "tetNig2", "gasAcu1"),
              "danRer5"=c("hg19", "hg18", "mm9", "oryLat2", "tetNig1", "fr2"),
              "tetNig1"=c("hg18", "danRer5"),
              "tetNig2"=c("hg19", "mm9", "danRer7", "oryLat2", "fr3", "gasAcu1", "danRer6"),
              "dm3"=c("droAna3", "dp4", "droVir3", "droMoj3"),
              #"ce10"=c("cb4", "caeRem4", "caePb3"),
              "ce4"=c("cb3", "caeRem2", "caePb1"),
              "xenTro3"=c("hg19"),
              "xenTro2"=c("hg18"),
              "fr2"=c("danRer5"),
              "fr3"=c("hg19", "mm10", "danRer7", "tetNig2"),
              "ornAna1"=c("hg18"),
              "monDom4"=c("hg18"),
              "oryLat2"=c("danRer7", "danRer6", "danRer5", "tetNig2"),
              "gasAcu1"=c("danRer7", "danRer6", "tetNig2"),
              "droAna3"=c("dm3"),
              "dp4"=c("dm3"),
              "droVir3"=c("dm3"),
              "droMoj3"=c("dm3"),
              #"cb4"=c("ce10"),
              "cb3"=c("ce4"),
              "caeRem2"=c("ce4"),
              "caePb1"=c("ce4"),
              "anoCar2"=c("hg19"),
              #"caeRem4"=c("ce10"),
              #"caePb3"=c("ce10"),
              NULL)
######## download the pairwise alignment from UCSC############
outputDir = "/export/downloads/ucsc/axtNet"
report = c()
build1 = names(builds)[1]
for(build1 in names(builds)){
  for(build2 in builds[[build1]]){
    res = downloadUCSCPA(build1, build2, outputDir=outputDir)
    report = c(report, res)
  }
}

######## download the repeats from UCSC######################
builds = list("old"=c("hg18", "mm9", "equCab2", "tetNig1", "dm3", "ce4"),
              "new"=c("hg19", "mm10", "canFam3", "galGal4", "danRer7", "danRer6", "danRer5", "tetNig2"))
suffix = list("old"=c("*_gap.*", "*_rmsk.*", "refGene.*", "refLink.*", "knownGene.*", "knownIsoforms.*", "kgXref.*", "blastKG*", "")






