\name{delete_Matrix_having_ID}
\docType{methods}
\alias{delete_Matrix_having_ID}
\alias{store_Matrix}
\alias{initializeJASPARDB}
%\alias{delete_Matrix_having_ID,character-method}
%\alias{store_Matrix,character-method}
%\alias{initializeJASPARDB, character-method}
%\alias{delete_Matrix_having_ID,SQLiteConnection-method}
%\alias{store_Matrix, SQLiteConnection-method}
%\alias{initializeJASPARDB, SQLiteConnection-method}

\title{JASPAR database operations}
\description{
  The functions to initialize, store matrix or delete matrix in JASPAR database.
}

\usage{
  \S4method{delete_Matrix_having_ID}{character}(x, IDs)
  \S4method{delete_Matrix_having_ID}{SQLiteConnection}(x, IDs)
  \S4method{store_Matrix}{character, PFMatrixList}(x, pfmList)
  \S4method{store_Matrix}{SQLiteConnection, PFMatrixList}(x, pfmList)
  \S4method{store_Matrix}{character, PFMatrix}(x, pfmList)
  \S4method{store_Matrix}{SQLiteConnection, PFMatrix}(x, pfmList)
  \S4method{initializeJASPARDB}{SQLiteConnection}(x)
  \S4method{initializeJASPARDB}{character}(x)
}

\arguments{
  \item{x}{
    a character vector of length 1 for the path of JASPAR SQLite file, or a \code{SQLiteConnection} object.
  }
  \item{IDs}{
    JASPAR stable IDs.
  }
  \item{pfmList}{
    The PFMatrixList object.
  }
}


\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
