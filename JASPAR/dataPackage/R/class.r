library(RSQLite)
library(Biostrings)


setClass("JASPARDb",
         representation(
                        provider="character",
                        provider_version="character",
                        release_date="character",
                        release_name="character",
                        source_url="character",
                        metadata_dirpath="character",
                        conn="SQLiteConnection")
         )

setClass("JASPAR",
         representation(
                        ID="character",
                        collection="Rle",
                        version="Rle",
                        name="character",
                        species="Rle",
                        TF_class="Rle",
                        medline="integer",
                        family="Rle",
                        tax_group="Rle",
                        acc="character",
                        type="Rle",
                        pazar_tf_id="character",
                        comment="character",
                        FMatrix="list",
                        seqs="DNAStringSetList"
                        )
         )
                        



## Constructor
JASPARDb = function(provider=character(), provider_version=character(), 
                  release_date=character(), release_name=character(), 
                  source_url=character(), metadata_dirpath=character()){
  new("JASPARDb", provider=provider, provider_version=provider_version,
      release_date=release_date, release_name=release_name,
      source_url=source_url, metadata_dirpath=metadata_dirpath, 
      conn=dbConnect(SQLite(), metadata_dirpath))
}

JASPAR = function(ID=character(), collection=Rle(), version=Rle(),
                  name=character(), species=Rle(), TF_class=Rle(),
                  medline=integer(), family=Rle(), tax_group=Rle(),
                  acc=character(), type=Rle(), pazar_tf_id=character(),
                  comments=character(), FMatrix=list(), seqs=DNAStringSetList()){
  new("JASPAR", ID=ID, collection=collection, version=version, name=name, species=species,
      TF_class=TF_class, medline=medline, family=family, tax_group=tax_group,
      acc=acc, type=type, pazar_tf_id=pazar_tf_id, comments=comments, 
      FMatrix=FMatrix, seqs=seqs)
}

## Methods
setGeneric("sqliteDir", function(x) standardGeneric("sqliteDir"))
setMethod("sqliteDir", "JASPARDb", function(x) x@metadata_dirpath)
setGeneric("conn", function(x) standardGeneric("conn"))
setMethod("conn", "JASPARDb", function(x) x@conn)
setGeneric("openDb", function(x) standardGeneric("openDb"))
setMethod("openDb", "JASPARDb", function(x) {
          conn(x)=dbConnect(SQLite(), qliteDir(x))
          return(x)})
setGeneric("closeDb", function(x) standardGeneric("closeDb"))
setMethod("closeDb", "JASPARDb", function(x) dbDisconnect(conn(x)))

## 
setGeneric("searchDb", function(x, ID, name, species, class, type) standardGeneric("searchDb"))
setMethod("searchDb", "JASPARDb", function(x, ID=NULL, name=NULL, species=NULL, 
                                         class=NULL, type=NULL){
      conn = conn(x)
      dbListTables(conn)
      callName = c("ID", "name", "species", "class", "type")
      dbField = c("BASE_ID", "NAME", "SPECIES", "TAG", "TAG")
      sqlCMD = paste("select ID from MATRIX")
      dbIDs = dbGetQuery(conn, sqlCMD)[[1]]
      # ID = c("MA0001", "MA0002")
      if(!is.null(ID)){
        sqlCMD = paste("select ID from MATRIX where BASE_ID IN ", "(\"", paste(ID, collapse="\",\""), "\")", sep="")
        dbIDs = intersect(dbIDs, dbGetQuery(conn, sqlCMD)[[1]])
      }
      # name = c("AGL3", "RUNX1")
      if(!is.null(name)){
        sqlCMD = paste("select ID from MATRIX where NAME IN ", "(\"", paste(name, collapse="\",\""), "\")", sep="")
        dbIDs = intersect(dbIDs, dbGetQuery(conn, sqlCMD)[[1]])
      }
      # species = c("Bacteria", "Azorhizobium")
      if(!is.null(species)){
        sqlCMD = paste("select TAX_ID from TAX where SPECIES IN ", "(\"", paste(species, collapse="\",\""), "\")", sep="")
        taxIDs = dbGetQuery(conn, sqlCMD)[[1]]
        sqlCMD = paste("select ID from MATRIX_SPECIES where TAX_ID IN ", "(", paste(taxIDs, collapse=","), ")", sep="")
        dbIDs = intersect(dbIDs, dbGetQuery(conn, sqlCMD)[[1]])
      }
      # class = c("Ig-fold")
      if(!is.null(class)){
        sqlCMD = paste("select ID from MATRIX_ANNOTATION where TAG == \"class\" and VAL IN ", "(\"", paste(class, collapse="\",\""), "\")", sep="")
        dbIDs = intersect(dbIDs, dbGetQuery(conn, sqlCMD)[[1]])
      }
      # type = c("SELEX")
      if(!is.null(type)){
        sqlCMD = paste("select ID from MATRIX_ANNOTATION where TAG == \"type\" and VAL IN ", "(\"", paste(type, collapse="\",\""), "\")", sep="")
        dbIDs = intersect(dbIDs, dbGetQuery(conn, sqlCMD)[[1]])
      }
      if(length(dbIDs) == 0){
        stop("No match result!")
      }
      res = data.frame(row.names=dbIDs)
      # collect from MATRIX table
      sqlCMD = paste("select * from MATRIX where ID IN ", "(", paste(dbIDs, collapse=","), ")", sep="")
      res = cbind(res, dbGetQuery(conn, sqlCMD))
      # collect from MATRIX_ANNOTATION table

      }
      
)
