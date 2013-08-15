preBlatFilter = "/export/data/CNEs/pre-blatFilter"
blatFilter = "."

nr_cores = 8


files = list.files(preBlatFilter, pattern="^cne2w_.*", full.names=TRUE)
filesExist = list.files(blatFilter, pattern="^cne2w.*")
filesNotToCompute = sub("Bf", "", filesExist)
filesNotToCompute = sub("-tmp$", "", filesNotToCompute)

files = setdiff(files, file.path(preBlatFilter, filesNotToCompute))

library(multicore)

cmd = "perl /opt/www/cne/scripts/cne_pipeline/blat_filter.pl --tmp "
setwd(blatFilter)
report = mclapply(files, function(oneFile){
                  message("Doing ", oneFile);
                  res=system(paste0(cmd, basename(oneFile),"-tmp ", oneFile)); 
                  if(res != 0){return(oneFile)}}, mc.cores=nr_cores)


message("The running report: ", report)

