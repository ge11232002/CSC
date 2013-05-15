preBlatFilter = "/export/data/CNEs/pre-blatFilter"
blatFilter = "/export/data/CNEs/blatFiltered"

nr_cores = 4

preBlatFilter = "/export/data/CNEs/hg19/annotation"

files = list.files(preBlatFilter, full.names=TRUE)

library(multicore)

cmd = "perl /opt/www/cne/scripts/cne_pipeline/blat_filter.pl --tmp "

report = mclapply(files, function(oneFile){
                  res=try(system(paste(cmd, Sys.getpid(),"-tmp ", oneFile, sep=""))); 
                  if(res != 0){return(oneFile)}}, mc.cores=nr_cores)


message("The running report: ", report)

