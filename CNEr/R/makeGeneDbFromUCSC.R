## This is used to download the annotation from UCSC and incooperate it for the tracks display

#.SUPPORTED_UCSC_TABLES = c(
#  ## tablename (unique key)   track             subtrack    auxiliary tablename
#  "knownGene",              "UCSC Genes",       NA,         
#  "refGene",                "RefSeq Genes",     NA,
#  "ensGene",                "Ensembl Genes",    NA
#  )

.SUPPORTED_UCSC_TABLES = list(
    "UCSC Genes"    = c("knownGene", "knownIsoforms", "kgXref"),
    "RefSeq Genes"  = c("refGene", "refLink"),
    "Ensembl Genes" = c("ensGene")
    )

supportedUCSCtables = function(){
  .SUPPORTED_UCSC_TABLES
}

queryRefSeq = function(con, tableNames){
  query = "SELECT distinct locusLinkId, refGene.name, refGene.name2, 'RefSeq Gene', chrom, strand, exonStarts, exonEnds, cdsStart, cdsEnd
                   FROM refGene, refLink WHERE refGene.name = refLink.mrnaAcc
                   ORDER BY locusLinkId, name2, name"
  ans = dbGetQuery(con, query)

}

makeGeneDbFromUCSC = function(genome="hg19",
                              tablename="RefSeq Genes",
                              host="genome-mysql.cse.ucsc.edu",
                              user="genome",
                              password=NULL,
                              opts=NULL){
  require(RMySQL)
  if(!isSingleString(genome))
    stop("'genome' must be a single string")
  if(!isSingleString(tablename))
    stop("'tablename' must be a single string")
  if(!tablename %in% names(.SUPPORTED_UCSC_TABLES))
    stop("table \"", tablename, "\" is not supported")
  if(!isSingleString(host))
    stop("'url' must be a single string")
  con = dbConnect(MySQL(), user=user, password=password, dbname=genome, host=host)
  tableNames = .SUPPORTED_UCSC_TABLES[[tablename]] 
  message("Download the ", tablename, " table ... ")
  ans = switch(tablename,
               "RefSeq Genes"=queryRefSeq(con)
               )

}



