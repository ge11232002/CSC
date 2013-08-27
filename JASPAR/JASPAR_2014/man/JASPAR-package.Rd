\name{JASPAR-package}
\alias{JASPAR-package}
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
Package: \tab JASPAR\cr
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
  JASPAR-package          Data package for JASPAR
  JASPARDb-class          JASPARDb object
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Ge Tan <gtan@me.com> 

Maintainer: Ge Tan <gtan@me.com>
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
## load the library
library(JASPAR)
## list the contents that are loaded into memory
ls("package:JASPAR")
## show the db object that is loaded by calling it's name
JASPAR_2010
}
