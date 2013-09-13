\name{JASPAR2014-package}
\alias{JASPAR2014-package}
\alias{JASPAR}
\docType{package}
\title{
Data package for JASPAR 2014
}
\description{
Data package for JASPAR 2014.
}
\details{
\tabular{ll}{
Package: \tab JASPAR2014\cr
Version: \tab 0.99.0\cr
Date: \tab 2013-09-12\cr
Depends: \tab R (>= 2.9.0), Biostrings(>= 2.26.3), RSQLite(>= 0.11.4)\cr
Suggests: \tab \cr
License: \tab GPL (>= 2)\cr
URL: \tab http://jaspar.binf.ku.dk\cr
Type: \tab Package\cr
NeedsCompilation: \tab no\cr
LazyData: \tab TRUE\cr
}

Index:
\preformatted{
  JASPAR-class            JASPAR object
  JASPAR2014-package      Data package for JASPAR 2014
%  JASPARDb-class          JASPARDb object
}

Overview:
A typical usage of the JASPAR data package is:
\enumerate{
  \item Load the package:
  \preformatted{
    library(JASPAR2014)
  }
  \item Query the database:
  \preformatted{
  ans = searchDb(JASPAR2014, class="Zipper-Type", type="SELEX")
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
  See \url{http://jaspar.binf.ku.dk} for more details about JASPAR.
}
\keyword{ package }
\examples{
## load the library
library(JASPAR2014)
## list the contents that are loaded into memory
ls("package:JASPAR2014")
## show the db object that is loaded by calling it's name
JASPAR2014
}



