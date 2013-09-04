### ----------------------------------------------------------------
### SitePairList obejct: holds the list of SitePair.
###

setClass("SitePairList",
         contains="SimpleList",
         representation(
                        ),
         prototype(
                   elementType="SitePair"
                   )
         )
### -----------------------------------------------------------
### The accessor-like methods
###
setMethod("site1", "SitePairList",
          function(x){
            ans = SiteList(lapply(x, site1))
            return(ans)
          }
          )
setMethod("site2", "SitePairList",
          function(x){
            ans = SiteList(lapply(x, site2))
            return(ans)
          }
          )

### ------------------------------------------------------------
### SitePairList constructor
### 
SitePairList = function(..., use.names=TRUE){
  listData = list(...)
  if(is(listData[[1]], "list")) # This is pretty ugly. better solution?
    listData = listData[[1]]
  ok = sapply(listData, is, "SitePair")
  if(!all(ok))
    stop("SitePairList() only accepts SitePair objects!")
  if(!use.names)
    names(listData) = NULL
  IRanges:::newList("SitePairList", listData)
}



