
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



