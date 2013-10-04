
### -----------------------------------------------------------------------
### The "show" method.
###

setMethod("show", "XMatrix",
          function(object){
            cat("An object of class ", class(object), "\n", sep="")
            printList = list(
                             ID=ID(object),
                             Name=name(object),
                             "Matrix class"=matrixClass(object),
                             Strand=strand(object)
                             )
            if(is(object, "PWMatrix"))
              printList = c(printList, list(Pseudocounts=pseudocounts(object)))
            if(is(object, "ICMatrix"))
              printList = c(printList, list("Schneider correction"=schneider(object)))
            printList = as.data.frame(c(printList, tags(object)))
            print(printList)
            cat("Background:", bg(object), "\n")
            cat("Matrix:", "\n")
  # add the tags later and print pretty
            print(Matrix(object))
          }
          )

