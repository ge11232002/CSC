\name{saveCNEToSQLite}
\alias{saveCNEToSQLite}

\title{
saveCNEToSQLite function
}
\description{
This function save the CNE table into a local SQLite database.
}
\usage{
saveCNEToSQLite(CNE, dbName, tableName, overwrite=FALSE)
}

\arguments{
  \item{CNE}{
    A object of \code{data.frame}, the CNE data table.
  }
  \item{dbName}{
    A object of \code{character}, the path of the local SQLite database.
  }
  \item{tableName}{
    A object of \code{character}, the name of table for this CNE data table.
  }
  \item{overwrite}{
    A object of \code{boolean}, whether or not to overwrite the table with same table name.
  }
}
\details{
  The input CNE table should have the colnames 
  "chr1", "start1", "end1", "chr2", "start2", "end2", "strand", "similarity", "cigar".
  After the bin indexing, two additional columns "bin1" and "bin2" will be
  added before the column "chr1" and "chr2", respectively.
}


\author{
Ge Tan
}

\examples{
  dbName = tempfile()
  for(i in 1:length(cneBlatedDanRer7Hg19)){
    tableName = paste0("cne_twoWay_danRer7_hg19_len50_id", 
    as.integer(sub("_\\d+$", "", names(cneBlatedDanRer7Hg19)[i]))*20, "_v1")
    saveCNEToSQLite(cneBlatedDanRer7Hg19[[i]], dbName, tableName, overwrite=TRUE)
  }
}

