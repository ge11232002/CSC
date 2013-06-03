preBlatFilter = "/export/data/CNEs/pre-blatFilter"
blatFilter = "/export/data/CNEs/blatFiltered"

nr_cores = 8


files = list.files(preBlatFilter, pattern="^cne2w_.*", full.names=TRUE)
filesExist = list.files(blatFilter, pattern="^cne2wBf_.*")
filesNotToCompute = sub("Bf", "", filesExist)
files = setdiff(files, file.path(preBlatFilter, filesNotToCompute))

library(multicore)

cmd = "perl /opt/www/cne/scripts/cne_pipeline/blat_filter.pl --tmp "
setwd(blatFilter)
report = mclapply(files, function(oneFile){
                  message("Doing ", oneFile);
                  res=try(system(paste(cmd, Sys.getpid(),"-tmp ", oneFile, sep=""))); 
                  if(res != 0){return(oneFile)}}, mc.cores=nr_cores)


message("The running report: ", report)

