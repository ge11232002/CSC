preBlatFilter = "/export/data/CNEs/pre-blatFilter"
blatFilter = "."

nr_cores = 8


files = list.files(preBlatFilter, pattern="^cne2w_.*", full.names=TRUE)
filesExist = list.files(blatFilter, pattern="^cne2w.*")
filesNotToCompute = sub("Bf", "", filesExist)
filesNotToCompute = sub("-tmp$", "", filesNotToCompute)

files = setdiff(files, file.path(preBlatFilter, filesNotToCompute))

library(parallel)
library(BatchJobs)

cmd = "perl /opt/www/cne/scripts/cne_pipeline/blat_filter.pl --tmp "
setwd(blatFilter)
reg <- makeRegistry(id="blatfilter", seed=123,
                    file.dir=file.path(getwd(), paste("blatfilter", 
                                                      "-batchjobs",
                                                      sep="")),
                    work.dir=getwd(), skip=FALSE)
batchMap(reg, function(oneFile, cmd){
         res=system(paste0(cmd, basename(oneFile)," ", oneFile))
         return(oneFile)
                    }, files, more.args=list(cmd=cmd)
         )
submitJobs(reg)
showStatus(reg)

#report = mclapply(files, function(oneFile){
#                  message("Doing ", oneFile);
#                  res=system(paste0(cmd, basename(oneFile),"-tmp ", oneFile)); 
#                  if(res != 0){return(oneFile)}}, mc.cores=nr_cores)



