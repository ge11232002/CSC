\name{JASPAR2016-class}
\docType{class}
\alias{JASPAR2016-class}
\alias{JASPAR2016}

\title{JASPAR2016 object}
\description{
  The JASPAR2016 object class is a thin class for storing the path of 
  JASPAR2016 style SQLite file. 
}

\section{Slots}{
  \describe{
    \item{\code{db}:}{Object of class \code{"character"}: 
      a character string of the path of SQLite file.}
  }
}

\author{
  Ge Tan
}

\examples{
  \dontrun{
    library(JASPAR2016)
    JASPAR2016

    library(TFBSTools)

    opts <- list()
    opts[["species"]] <- 9606
    opts[["type"]] <- "SELEX"
    opts[["all_versions"]] <- TRUE
    PFMatrixList <- getMatrixSet(JASPAR2016, opts)

    opts2 <- list()
    opts2[["type"]] <- "SELEX"
    PFMatrixList2 <- getMatrixSet(JASPAR2016, opts2)
  }
}
\keyword{classes}

