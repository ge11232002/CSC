

setClass("JASPARDb",
         representation(
                        provider="character",
                        provider_version="character",
                        release_date="character",
                        release_name="character",
                        source_url="character",
                        #sites_seqs="character",
                        db_dirpath="character",
                        conn="SQLiteConnection")
## perhaps should put a metadata table in JASPAR.sqlite and fill this object dynamically
         )

setClass("JASPAR",
         representation(
                        ID="character",
                        collection="Rle",
                        version="Rle",
                        name="character",
                        species="Rle",
                        TF_class="Rle",
                        medline="character",
                        family="Rle",
                        tax_group="Rle",
                        acc="character",
                        type="Rle",
                        pazar_tf_id="character",
                        comment="character",
                        FMatrix="list",
                        seqs="list"
                        )
         )
                        



## Constructor
JASPARDb = function(provider=character(), provider_version=character(), 
                  release_date=character(), release_name=character(), 
                  source_url=character(), #sites_seqs=character(), 
                  db_dirpath=character()){
  new("JASPARDb", provider=provider, provider_version=provider_version,
      release_date=release_date, release_name=release_name,
      source_url=source_url, #sites_seqs=sites_seqs, 
      db_dirpath=db_dirpath, 
      conn=dbConnect(SQLite(), db_dirpath))
}

JASPAR = function(ID=character(), collection=Rle(), version=Rle(),
                  name=character(), species=Rle(), TF_class=Rle(),
                  medline=character(), family=Rle(), tax_group=Rle(),
                  acc=character(), type=Rle(), pazar_tf_id=character(),
                  comment=character(), FMatrix=list(), seqs=list()){
  new("JASPAR", ID=ID, collection=collection, version=version, name=name, species=species,
      TF_class=TF_class, medline=medline, family=family, tax_group=tax_group,
      acc=acc, type=type, pazar_tf_id=pazar_tf_id, comment=comment, 
      FMatrix=FMatrix, seqs=seqs)
}

### Updating and cloning.
setMethod("update", "JASPAR",
          function(object, ..., check=TRUE){
            initialize(object, ...)
          }
)

setGeneric("clone", function(x, ...) standardGeneric("clone"))  # not exported
setMethod("clone", "ANY",  # not exported
    function(x, ...)
    {
        if (nargs() > 1L)
            initialize(x, ...)
        else
            x
    }
)


## Slot getters and setters.
setGeneric("sqliteDir", function(x) standardGeneric("sqliteDir"))
setMethod("sqliteDir", "JASPARDb", function(x) x@db_dirpath)
setGeneric("conn", function(x) standardGeneric("conn"))
setMethod("conn", "JASPARDb", function(x) x@conn)
setGeneric("ID", function(x) standardGeneric("ID"))
setMethod("ID", "JASPAR", function(x) x@ID)
setGeneric("collection", function(x) standardGeneric("collection"))
setMethod("collection", "JASPAR", function(x) unique(x@collection))
setMethod("names", "JASPAR", function(x) x@name)
setGeneric("TFBSinfo", function(x) standardGeneric("TFBSinfo"))
setMethod("TFBSinfo", "JASPAR", function(x){
                                res = cbind(x@ID, x@name, as.vector(x@species), as.vector(x@TF_class), as.vector(x@family), as.vector(x@tax_group), x@acc, as.vector(x@type), x@medline, x@pazar_tf_id, x@comment)
                                colnames(res) = c("ID", "name", "species", "class", "family", "tax_group", "acc", "type", "medline", " Pazar ID", "comment")
                                rownames(res) = paste(x@ID, x@version, sep=".")
                                return(res)
                                })
setGeneric("FMatrix", function(x) standardGeneric("FMatrix"))
setMethod("FMatrix", "JASPAR", function(x){
                                          res = x@FMatrix
                                          names(res) = paste(x@ID, x@version, sep=".")
                                          return(res)})
setGeneric("seqs", function(x) standardGeneric("seqs"))
setMethod("seqs", "JASPAR", function(x) x@seqs)

#### setMethods
setGeneric("openDb", function(x) standardGeneric("openDb"))
setMethod("openDb", "JASPARDb", function(x) {
          conn(x)=dbConnect(SQLite(), qliteDir(x))
          return(x)})
setGeneric("closeDb", function(x) standardGeneric("closeDb"))
setMethod("closeDb", "JASPARDb", function(x) dbDisconnect(conn(x)))
setMethod("length", "JASPAR", function(x) length(x@ID))


# The main searchDb method
setGeneric("searchDb", function(x, ID=NULL, name=NULL, species=NULL, class=NULL, type=NULL) standardGeneric("searchDb"))
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
      res = list()
      # collect from MATRIX table
      sqlCMD = paste("select * from MATRIX where ID IN ", "(", paste(dbIDs, collapse=","), ")", sep="")
      matrix_table = dbGetQuery(conn, sqlCMD)
      res$ID = matrix_table[match(dbIDs, matrix_table$ID), "BASE_ID"]
      res$collection = matrix_table[match(dbIDs, matrix_table$ID), "COLLECTION"]
      res$version = matrix_table[match(dbIDs, matrix_table$ID), "VERSION"]
      res$name = matrix_table[match(dbIDs, matrix_table$ID), "NAME"]
      # collect from MATRIX_ANNOTATION table
      sqlCMD = paste("select * from MATRIX_ANNOTATION where ID IN ", "(", paste(dbIDs, collapse=","), ")", sep="")
      matrix_annotation_table = dbGetQuery(conn, sqlCMD)
      for(tag in c("class", "medline", "family", "tax_group", "type", "pazar_tf_id", "comment")){
        temp_table = matrix_annotation_table[matrix_annotation_table$TAG == tag, ]
        res[[tag]] = temp_table[match(dbIDs, temp_table$ID), "VAL"]
      }
      # collect from MATRIX_PROTEIN table
      sqlCMD = paste("select * from MATRIX_PROTEIN where ID IN ", "(", paste(dbIDs, collapse=","), ")", sep="")
      matrix_protein_table = dbGetQuery(conn, sqlCMD)
      res$acc = matrix_protein_table[match(dbIDs, matrix_protein_table$ID), "ACC"]
      # collect from MATRIX_SPECIES and TAX table
      sqlCMD = paste("select * from MATRIX_SPECIES where ID IN ", "(", paste(dbIDs, collapse=","), ")", sep="")
      matrix_species_table = dbGetQuery(conn, sqlCMD)
      sqlCMD = paste("select * from TAX where TAX_ID IN ", "(\"", paste(unique(matrix_species_table$TAX_ID), collapse="\",\""), "\")", sep="")
      tax_table = dbGetQuery(conn, sqlCMD)
      res$species = tax_table[match(matrix_species_table[match(dbIDs, matrix_species_table$ID), "TAX_ID"], tax_table$TAX_ID), "SPECIES"]
      # collect from MATRIX_DATA table
      sqlCMD = paste("select * from MATRIX_DATA where ID IN ", "(", paste(dbIDs, collapse=","), ")", sep="")
      matrix_data_table = dbGetQuery(conn, sqlCMD)
      res$FMatrix = list()
      for(dbID in dbIDs){
        res$FMatrix[[as.character(dbID)]] = matrix(matrix_data_table[matrix_data_table$ID == dbID, "val"], nrow=4, byrow=TRUE, dimnames=list(c("A", "C", "G", "T")))
      }
      # collect DNA Seq from sites.rda 
      res$seqs = switch(x@release_name,
                        "JASPAR_2010"=JASPAR_2010_SitesSeqs[paste(res$ID, res$version, sep=".")]
                        )
      #if(!exists(x@sites_seqs)){
      #  data(x@sites_seqs)
      #}
      #res$seqs = sitesSeqs[paste(res$ID, res$version, sep=".")]
      names(res$seqs) = paste(res$ID, res$version, sep=".")
      return(JASPAR(ID=res$ID, collection=Rle(res$collection), version=Rle(res$version),
             name=res$name, species=Rle(res$species), TF_class=Rle(res$class),
             medline=res$medline, family=Rle(res$family), tax_group=Rle(res$tax_group),
             acc=res$acc, type=Rle(res$type), pazar_tf_id=res$pazar_tf_id, 
             comment=res$comment, FMatrix=res$FMatrix, seqs=res$seqs))
      }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting and combining.
###

setMethod("[", "JASPAR",
          function(x, i, ..., drop){
            if(length(list(...)) > 0L)
              stop("invalid subsetting")
            if(missing(i))
              return(x)
            #i = IRanges:::normalizeSingleBracketSubscript(i, x)
            ans_ID = x@ID[i]
            ans_collection = x@collection[i]
            ans_version = x@version[i]
            ans_name = x@name[i]
            ans_species = x@species[i]
            ans_TF_class = x@TF_class[i]
            ans_medline = x@medline[i]
            ans_family = x@family[i]
            ans_tax_group = x@tax_group[i]
            ans_acc = x@acc[i]
            ans_type = x@type[i]
            ans_pazar_tf_id = x@pazar_tf_id[i]
            ans_comment = x@comment[i]
            ans_FMatrix = x@FMatrix[i]
            ans_seqs = x@seqs[i]
            clone(x, ID=ans_ID, collection=ans_collection, version=ans_version,
                  name=ans_name, species=ans_species, TF_class=ans_TF_class,
                  medline=ans_medline, family=ans_family, tax_group=ans_tax_group,
                  acc=ans_acc, type=ans_type, pazar_tf_id=ans_pazar_tf_id, 
                  comment=ans_comment, FMatrix=ans_FMatrix, seqs=ans_seqs)
          }
          )

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "show" method.
###
showJASPAR = function(x, half_nrow=5L){
  lx = length(x)
  if(is.null((head_nrow = getOption("showHeadLines"))))
    head_nrow = half_nrow
  if(is.null((tail_nrow = getOption("showTailLines"))))
    tail_nrow = half_nrow
  iW = nchar(as.character(lx))
  if(lx < (2*half_nrow+1L) | (lx < (head_nrow+tail_nrow+1L))){
    IDW = max(nchar(ID(x)))
    collectionW = max(nchar(collection(x)))
  }
}

setMethod("show", "JASPAR",
          function(object){
            lx = length(object)
            cat(" ", lx, " ", 
                ifelse(lx == 1L, "matched result", "matched results"), 
                " ", "in", class(object))
            cat("\n")
            if(lx != 0)
              print(TFBSinfo(object))
          }
          )


