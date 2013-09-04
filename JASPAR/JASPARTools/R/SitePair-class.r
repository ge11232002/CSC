### -------------------------------------------------------------------
### SitePair object: a nucleotide sequence feature object representing (possibly putative) transcription factor binding site from A alignment

setClass("SitePair",
         slots=c(site1="Site",
                 site2="Site"
                 )
         )

### -------------------------------------------------------------------
### The accessor-like method
### 
setMethod("site1", "SitePair", function(x) x@site1)

setMethod("site2", "SitePair", function(x) x@site2)

setMethod("length", "SitePair", function(x) length(site1(x)))
### -----------------------------------------------------------------
### The constructor
###
SitePair = function(site1, site2){
  new("SitePair", site1=site1, site2=site2)
}

### ------------------------------------------------------------------
### Method
###
setMethod("writeGFF3", "SitePair",
          function(x){
            gff1 = writeGFF3(site1(x))
            gff2 = writeGFF3(site2(x))
            ans = rbind(gff1, gff2)
            return(ans)
          }
          )
setMethod("writeGFF2", "SitePair",
          function(x){
            gff1 = writeGFF2(site1(x))
            gff2 = writeGFF2(site2(x))
            ans = rbind(gff1, gff2)
            return(ans)
          }
          )


### -----------------------------------------------------------------
### The "show" method
### put them in a extended gff. any good idea?
setMethod("show", "SitePair",
          function(object){
            gff1 = writeGFF3(site1(object))
            gff2 = writeGFF3(site2(object))
            ans = cbind(gff1, gff2)
            cat("An object of class", class(object), "with",
                length(object), "site pair",
                ifelse(length(object)==1, "sequence", "sequences"))
            cat("\n")
            print(ans)
          }
          )


