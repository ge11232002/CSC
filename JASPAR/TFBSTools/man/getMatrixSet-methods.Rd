\name{get_MatrixSet}
\docType{methods}
\alias{get_MatrixSet}
\alias{get_MatrixSet,character-method}
\alias{get_MatrixSet,SQLiteConnection-method}
\title{Advanced JASPAR database search functions \code{get_MatrixSet}}
\description{
  This function fetches matrix data for all matrices in the database 
  matching criteria defined by the named arguments 
  and returns a XMatrixList object
}

\usage{
  \S4method{get_MatrixSet}{character}(x, opts)
  \S4method{get_MatrixSet}{SQLiteConnection}(x, opts)
}

\arguments{
  \item{x}{
    a character vector of length 1 for the path of JASPAR SQLite file, or a \code{SQLiteConnection} object.
  }
  \item{opts}{
    a search options list. See more details below.
  }
}

\details{
  The search options include three categories:

  (1) Database basic criterias: 
  
  \code{all=c(TRUE, FALSE)} 
  
  \code{ID}, \code{name} 
  
  \code{collection=c("CORE", "CNE", "PHYLOFACTS", ...)}
  
  \code{all_versions=c(TRUE,FALSE)}
  
  \code{species}, \code{matrixtype=c("PFM", "PWM", "ICM")}

  (2) Tags based criterias: 
  
  \code{class}, \code{type}, \code{family}
  
  \code{medline}, \code{tax_group}

  (3) Further criterias: 
  
  \code{min_ic} (minimum total information content of the matrix)
  
  \code{length} (minimum sites length)
  
  \code{sites} (minimum average sites number per base)

  When \code{all} is \code{TRUE}, it will get all the matrices 
  and has higher priority over other options. 
  Then ID has the second highest priority, 
  and will ignore all the followiing options.
  The rest options are combined in search with AND, 
  while multiple elements under one options have the logical operator OR.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{get_Matrix_by_ID}},
  \code{\link{get_Matrix_by_name}}
}

\examples{
  \dontrun{
    db = "/Users/gtan/CSC/JASPAR/JASPAR2014/inst/extdata/JASPAR2014.sqlite"
    opts = list()
    opts[["species"]] = 9606
    opts[["name"]] = "RUNX1"
    #opts[["class"]] = "Ig-fold"
    opts[["type"]] = "SELEX"
    opts[["all_versions"]] = TRUE
    siteList = get_MatrixSet(db, opts)
  }

}

\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
