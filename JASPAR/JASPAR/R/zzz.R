

.onLoad <- function(libname, pkgname){
  ns <- asNamespace(pkgname)
  path <- system.file("extdata", package=pkgname)
  files <- dir(path)
  # files = "/mnt/biggley/home/gtan/Repos/CSC/JASPAR/dataPackage/inst/extdata/JASPAR_2010.sqlite"
  for(i in seq_len(length(files))){
    jasparDb <- JASPARDb(metadata_dirpath=files[i])
    objname <- sub(".sqlite$","",files[i])
    assign(objname, jasparDb, envir=ns)
    namespaceExport(ns, objname)
  }
}


