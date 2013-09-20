\name{get_Matrix_by_ID}
\docType{methods}
\alias{get_Matrix_by_ID}
\alias{get_Matrix_by_name}
\alias{get_Matrix_by_ID,character-method}
\alias{get_Matrix_by_name,character-method}
\alias{get_Matrix_by_ID,SQLiteConnection-method}
\alias{get_Matrix_by_name,SQLiteConnection-method}

\title{Basic JASPAR database search functionis \code{get_Matrix_by_ID}, 
  \code{get_Matrix_by_name}}
\description{
  This method fetches matrix data under the given ID or name from the database and returns a XMatrix object.
}

\usage{
  \S4method{get_Matrix_by_ID}{character}(x, ID, type="PFM")
  \S4method{get_Matrix_by_ID}{SQLiteConnection}(x, ID, type="PFM")
  \S4method{get_Matrix_by_name}{character}(x, name, type="PFM")
  \S4method{get_Matrix_by_name}{SQLiteConnection}(x, name, type="PFM")
}

\arguments{
  \item{x}{
    a character vector of length 1 for the path of JASPAR SQLite file, or a \code{SQLiteConnection} object.
  }
  \item{ID}{
    a character vector of length 1 for the JASPAR stable ID. 
    See more details below.
  } 
  \item{name}{
    a character vector of length 1 for the JASPAR profile name.
  }
  \item{type}{
    the type can be "PFM", "PWM" or "ICM". 
    By default, only the freuency matrices are stored in JASPAR database.
    When the type "PWM" or "ICM" is selected, a default conversion will be done.
  }
}

\details{
  For get_Matrix_by_ID, 
  ID is a string which refers to the stable JASPAR ID 
  (usually something like "MA0001") with or without version numbers. 
  "MA0001" will give the latest version on MA0001, while "MA0001.2" will give the second version, if existing. 

  For get_Matrix_by_name, according to the current JASPAR data model, 
  name is not necessarily a unique identifier. 
  Also, names change over time. 
  In the case where there are several matrices with the same name in the database, 
  the function fetches the first one and prints a warning. 
  You've been warned. 
  Some matrices have multiple versions. 
  The function will return the latest version. 
  For specific versions, use get_Matrix_by_ID(ID.version)
}

\value{
  A XMatrix object is returned. 
  The exact type of the object depending on the \code{type} argument.
  NA is returned if matrix with the given ID/name is not found.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{get_MatrixSet}}
}

\examples{
  \dontrun{
    db = "/Users/gtan/CSC/JASPAR/JASPAR2014/inst/extdata/JASPAR2014.sqlite"
    pfm = get_Matrix_by_ID(db, ID="MA0003", type="PFM")
    pfm2 = get_Matrix_by_name(db, name="TFAP2A", type="PFM")
  }
}

\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
