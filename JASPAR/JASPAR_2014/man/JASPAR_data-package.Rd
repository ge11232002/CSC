\name{JASPARData-package}
\alias{JASPARData-package}
\alias{JASPAR}
\docType{package}
\title{
Data package for JASPAR
}
\description{
Data package for JASPAR. Currently, JASPAR 2010 is included.
}
\details{
\tabular{ll}{
Package: \tab JASPARData\cr
Version: \tab 0.0.1\cr
Date: \tab 2013-07-03\cr
Depends: \tab R (>= 2.9.0), Biostrings(>= 2.26.3), RSQLite(>= 0.11.4)\cr
Suggests: \tab \cr
License: \tab GPL (>= 2)\cr
URL: \tab http://jaspar.genereg.net/\cr
Type: \tab Package\cr
NeedsCompilation: \tab no\cr
LazyData: \tab TRUE\cr
}

Index:
\preformatted{
  JASPAR-class            JASPAR object
  JASPARData-package     Data package for JASPAR
  JASPARDb-class          JASPARDb object
}

Overview:
A typical usage of the JASPAR data package is:
\enumerate{
  \item Load the package:
  \preformatted{
    library(JASPARData)
  }
  \item Query the database:
  \preformatted{
  ans = searchDb(JASPAR_2010, class="Zipper-Type", type="SELEX")
  }
  \item Display the results:
  \preformatted{
    TFBSinfo(ans)
    FMatrix(ans)
    seqs(ans)
  }
}
}
\author{
Ge Tan <gtan@me.com> 

Maintainer: Ge Tan <gtan@me.com>
}
\references{
  See \url{http://jaspar.genereg.net/} for more details about JASPAR.
}
\keyword{ package }
\examples{
## load the library
library(JASPARData)
## list the contents that are loaded into memory
ls("package:JASPARData")
## show the db object that is loaded by calling it's name
JASPAR_2010
}
