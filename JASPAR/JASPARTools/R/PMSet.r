### -------------------------------------------------------------------
### XMatrixList objects
###

setClass("PFMatrixList",
         contains="CompressedList",
         representation(
                        "VIRTUAL",
                        unlistData="PFMatrix"
                        ),
         prototype(
                   elementType="PFMatrix"
                   )
         )

setClass("PWMatrixList",
         contains="CompressedList",
         representation(
                        "VIRTUAL",
                        unlistData="PWMatrix"
                        ),
         prototype(
                   elementType="PWMatrix"
                   )
         )

setClass("ICMatrixList",
         contains="CompressedList",
         representation(
                        "VIRTUAL",
                        unlistData="ICMatrix"
                        ),
         prototype(
                   elementType="ICMatrix"
                   )
         )

setClassUnion("XMatrixList", c("PFMatrixList", "ICMatrixList", "PWMatrixList"))

### ----------------------------------------------------------------------
### The accessor-like method
###



