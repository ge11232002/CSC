### --------------------------------------------------------------
### SiteList objects
###

setClass("SiteList",
         contains="SimpleList",
         representation(
                        ),
         prototype(
                   elementType="Site"
                   )
         )

### -------------------------------------------------------------
### SiteList() constructor
###
SiteList = function(..., use.names=TRUE){
  listData = list(...)
  if(is(listData[[1]], "list"))
    listData = listData[[1]]
  ok = sapply(listData, is, "Site")
  if(!all(ok))
    stop("SiteList() only accepts Site objects!")
  if(!use.names)
    names(listData) = NULL
  IRanges:::newList("SiteList", listData)
}

### ----------------------------------------------------------------
### Methods
###
setMethod("writeGFF3", "SiteList",
          function(x){
            ans = do.call(rbind, lapply(x, writeGFF3))
            return(ans)
          }
          )
setMethod("writeGFF2", "SiteList",
           function(x){
             ans = do.call(rbind, lapply(x, writeGFF2))
             return(ans)
           }
           )



