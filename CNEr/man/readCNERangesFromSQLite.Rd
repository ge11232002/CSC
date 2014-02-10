\name{readCNERangesFromSQLite}
\alias{readCNERangesFromSQLite}

\title{
readCNERangesFromSQLite function
}
\description{
Query the SQLite database based on chromosome, coordiantes and some other criterias.
Usually not to be used directly.
For the CNE density plot, \code{fetchCNEDensity} function should be used.
}
\usage{
readCNERangesFromSQLite(dbName, tableName, chr, start, end,
                        whichAssembly=c("L","R"), minLength=NULL)
}

\arguments{
  \item{dbName}{
    A object of \code{character}, the path of the local SQLite database.
  }
  \item{tableName}{
    A object of \code{character}, the name of table for this CNE data table.
  }
  \item{chr}{
    A object of \code{character}, the chromosome to query
  }
  \item{start, end}{
    A object of \code{integer}, the start and end coordiante to fetch the CNEs.
  }
  \item{whichAssembly}{
    A object of \code{character}, the genome to fetch is in the "Left" columns or "Right" columns of the table.
  }
  \item{minLength}{
    A object of \code{integer}, the minimal length for selected CNEs.
  }
}

\value{
A object of \code{IRanges} is retu
}

\author{
Ge Tan
}

\examples{
    dbName = file.path(system.file("extdata", package="CNEr"),
                       "cne.sqlite")
    chr = "chr6"
    start = 19900000
    end =  28000000
    minLength = 50
    tableName = "cne_twoWay_danRer7_hg19_len50_id960_v1"
    fetchedCNERanges = readCNERangesFromSQLite(dbName, tableName, chr, 
      start, end, whichAssembly="L", minLength=minLength)
}
