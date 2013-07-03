

.onLoad <- function(libname, pkgname){
  ns <- asNamespace(pkgname)
  path <- system.file("extdata", package=pkgname, lib.loc=libname)
  files <- list.files(path, pattern="\\.sqlite$", full.names=TRUE)
  for(i in seq_len(length(files))){
    jasparDb <- JASPARDb(release_name=sub("\\.sqlite$", "", basename(files[i])), db_dirpath=files[i])
    objname <- sub(".sqlite$","",basename(files[i]))
    assign(objname, jasparDb, envir=ns)
    namespaceExport(ns, objname)
  }
}


