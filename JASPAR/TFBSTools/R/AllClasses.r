### -----------------------------------------------------------------
### The position frequency matrix (PWM) class
###
setClass("PFMatrix",
         slots=c(ID="character",
                 name="character",
                 matrixClass="character",
                 strand="character",
                 bg="numeric",
                 tags="list",
                 matrix="matrix")
         )
setClass("PWMatrix", contains="PFMatrix",
         slots=c(pseudocounts="numeric")
         )
setClass("PWMatrix", contains="PFMatrix",
         slots=c(pseudocounts="numeric")
         )
setClassUnion("XMatrix", c("PFMatrix", "ICMatrix", "PWMatrix"))

### -----------------------------------------------------------------
### The set of PWM class
###
setClass("PFMatrixList",
         contains="SimpleList",
         representation(
                        ),
         prototype(
                   elementType="PFMatrix"
                   )
         )

setClass("PWMatrixList",
         contains="SimpleList",
         representation(
                        ),
         prototype(
                   elementType="PWMatrix"
                   )
         )

setClass("ICMatrixList",
         contains="PWMatrixList",
         representation(
                        ),
         prototype(
                   elementType="ICMatrix"
                   )
         )

setClassUnion("XMatrixList", c("PFMatrixList", "ICMatrixList", "PWMatrixList"))


