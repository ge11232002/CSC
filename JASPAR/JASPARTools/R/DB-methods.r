.fillDBOptsWithDefaults = function(opts=list()){
  if(!"collection" %in% names(opts))
    opts[["collection"]] = "CORE"
  if(!"all_versions" %in% names(opts))
    opts[["all_versions"]] = FALSE
  if(!"matrixtype" %in% names(opts))
    opts[["matrixtype"]] = "PFM"
  return(opts)
}


.get_latest_version = function(con, baseID){
  sqlCMD = paste0("SELECT VERSION FROM MATRIX WHERE BASE_ID='", baseID, 
                  "' ORDER BY VERSION DESC LIMIT 1")
  latest = dbGetQuery(con, sqlCMD)[["VERSION"]]
  return(latest)
}

.get_internal_id = function(con, baseID, version){
  # picks out the internal id for a stable id + version. Also checks if this cobo exists or not
  sqlCMD = paste0("SELECT ID FROM MATRIX WHERE BASE_ID='", baseID,"' AND VERSION='", version, "'")
  ini_id = dbGetQuery(con, sqlCMD)[["ID"]]
  if(length(ini_id) != 1)
    stop("There are ", length(ini_id), " records with this based id and version combination!")
  return(ini_id)
}

.get_Matrix_by_int_id = function(con, int_id, type){
  # Get the pfm matrix
  # bases orders in ("A", "C", "G", "T")
  sqlCMD = paste0("SELECT val FROM MATRIX_DATA WHERE ID='", int_id, "' ORDER BY col, row")
  matrixVector = dbGetQuery(con, sqlCMD)[["val"]]
  if(length(matrixVector) %% 4 != 0)
    stop("The number of retrived elements ", length(matrixVector), " is incomplete!")
  FMatrix = matrix(as.integer(matrixVector), nrow=4, dimnames=list(c("A", "C", "G", "T")))
  
  # get remaining data in the matrix table: name, collection
  sqlCMD = paste0("SELECT BASE_ID,VERSION,COLLECTION,NAME FROM MATRIX WHERE ID='",
                  int_id, "'")
  tempTable = dbGetQuery(con, sqlCMD)
  baseID = tempTable[["BASE_ID"]]
  version = tempTable[["VERSION"]]
  collection = tempTable[["COLLECTION"]]
  name = tempTable[["NAME"]]

  # get species
  sqlCMD = paste0("SELECT TAX_ID FROM MATRIX_SPECIES WHERE ID='", int_id, "'")
  tempTable = dbGetQuery(con, sqlCMD)
  tax_ids = tempTable[["TAX_ID"]] ## need to convert to taxs, fix this here or some place.
  if(length(tax_ids) == 0)
    tax_ids = ""

  # get acc
  sqlCMD = paste0("SELECT ACC FROM MATRIX_PROTEIN WHERE ID='", int_id, "'")
  tempTable = dbGetQuery(con, sqlCMD)
  accs = tempTable[["ACC"]]
  if(length(accs) == 0)
    accs = ""

  # get remaining annotation as tags, form ANNOTATION table
  sqlCMD = paste0("SELECT TAG,VAL FROM MATRIX_ANNOTATION WHERE ID='", int_id, "'")
  tempTable = dbGetQuery(con, sqlCMD)
  tags = list()
  tags = mapply(function(x,y){tags[[x]]=y}, tempTable[["TAG"]], tempTable[["VAL"]], SIMPLIFY=FALSE)
  tags[["collection"]] = collection
  tags[["species"]] = tax_ids
  tags[["acc"]] = accs
  matrixClass = tags[["class"]]
  tags["class"] = NULL
  
  ans_pfm = PFMatrix(ID=paste0(baseID, ".", version),
                     name=name,
                     matrixClass=matrixClass,
                     tags=tags,
                     matrix=FMatrix
                     )
  if(type == "PFM")
    return(ans_pfm)
  else if(type == "PWM")
    return(toPWM(ans_pfm))
  else if(type == "ICM")
    return(toICM(ans_pfm))
  else
    stop("This should never happen")
}

### get_Matrix_by_ID fetches matrix data under the given ID from the database and returns a XMatrix object.
# Returns : a XMatrix object; the exact type of the object depending on the second argument (allowed values are 'PFM', 'ICM', and 'PWM'); returns NA if matrix with the given ID is not found.
# Args: 
    #ID: is a string which refers to the stable JASPAR ID (usually something like "MA0001") with or without version numbers. "MA0001" will give the latest version on MA0001, while "MA0001.2" will give the second version, if existing. Warnings will be given for non-existing matrices.
setMethod("get_Matrix_by_ID", "character",
          function(x, ID, type="PFM"){
            # here x is the path of SQLite db file.
            type = match.arg(type, c("PWM", "PFM", "ICM"))
            if(missing(ID))
              stop("ID needs to be specified!")
            con = dbConnect(SQLite(), x)
            on.exit(dbDisconnect(con))
            # separate stable ID and version number
            baseID = strsplit(ID, "\\.")[[1]][1]
            version = strsplit(ID, "\\.")[[1]][2]
            if(is.na(version))
              version = as.character(.get_latest_version(con, baseID))
            if(length(version) == 0)
              return(NA)
            # get internal ID - also a check for validity
            int_id = as.character(.get_internal_id(con, baseID, version))
            # get matrix using internal ID
            ans = .get_Matrix_by_int_id(con, int_id, type)
            return(ans)
          }
          )

### get_Matrix_by_name fetches matrix data under the given name from the database and returns a XMatrix object.
# Returns : a XMatrix object; the exact type of the object depending on the second argument (allowed values are 'PFM', 'ICM', and 'PWM'); returns NA if matrix with the given name is not found.
# Notes: According to the current JASPAR5 data model, name is not necessarily a unique identifier. Also, names change over time. In the case where there are several matrices with the same name in the database, the function fetches the first one and prints a warning on STDERR. You've been warned. Some matrices have multiple versions. The function will return the latest version. For specific versions, use get_Matrix_by_ID($ID.$version)
setMethod("get_Matrix_by_name", "character",
          function(x, name, type="PFM"){
            # here x is the path of SQLite db file
            type = match.arg(type, c("PWM", "PFM", "ICM"))
            if(missing(name))
              stop("name needs to be specified!")
            con = dbConnect(SQLite(), x)
            on.exit(dbDisconnect(con))
            sqlCMD = paste0("SELECT distinct BASE_ID  FROM MATRIX WHERE NAME='", name, "'")
            tempTable = dbGetQuery(con, sqlCMD)
            baseID = tempTable[["BASE_ID"]]
            if(length(baseID) == 0)
              return(NA)
            if(length(baseID) > 1)
              warning("There are ", length(baseID), " distinct stable IDs with name ", name, ": ", baseID)
            get_Matrix_by_ID(x, baseID[1], type=type)
          }
          )


### get_MatrixSet fetches matrix data under for all matrices in the database matching criteria defined by the named arguments and returns a XMatrixList object
# Returns : a XMatrixList object
# Notes: This method accepts named arguments, corresponding to arbitrary tags, and also some utility functions. Note that this is different from JASPAR2 and to some extent JASPAR4. As any tag is supported for database storage, any tag can be used for information retrieval. Additionally, arguments as 'name','class','collection' can be used (even though they are not tags). By default, only the last version of the matrix is given. The only way to get older matrices out of this to use an array of IDs with actual versions like MA0001.1, or set the argyment -all_versions=>1, in which  case you get all versions for each stable ID.
# Args: 
  # -all: gives absolutely all matrix entry, regardless of versin and collection. Only useful for backup situations and sanity checks. Takes precedence over everything else
  # -ID: a reference to an array of stable IDs (strings), with or without version, as above. tyically something like "MA0001.2" . Takes precedence over everything salve -all
  # -name: a reference to an array of transcription factor names (string). Will only take latest version. NOT a preferred way to access since names change over time.
  # -collection: a string corresponding to a JASPAR collection. Per default CORE
  # -all_versions: gives all matrix versions that fit with rest of criteria, including obsolete ones. Is off per default. Typical usage is in combiation with a stable IDs without versions to get all versinos of a particular matrix.
  ## typical tags:
  # -class: structural class names (strings)
  # -species: NCBI Taxonomy IDs (integers)
  # -taxgroup: higher taxonomic categories (string)
  ## Computed features of the matrices
  # -min_ic: float, minimum total information content of the matrix.
  # -matrixtype: string describing type of matrix to retrieve. If left out, the format will revert to the database format, which is PFM.
  # The arguments that expect list references are used in database query formulation: elements within lists are combined with 'OR' operators, and the lists of different types with 'AND'.
  # For example,
    # my $matrixset = $db->(-class => ['TRP_CLUSTER', 'FORKHEAD'],
    #                       -species => ['Homo sapiens', 'Mus musculus'],
    #                      );
    # gives a set of TFBS::Matrix::PFM objects (given that the matrix models are stored as such) whose (structural clas is 'TRP_CLUSTER' OR'FORKHEAD') AND (the species they are derived from is 'Homo sapiens'OR 'Mus musculus').
  # As above, unless IDs with version numbers are used, only one matrix per stable ID wil be returned: the matrix with the highest version number
  #The -min_ic filter is applied after the query in the sense that the matrices profiles with total information content less than specified are not included in the set.
setMethod("get_MatrixSet", "character",
         function(x, opts){
           opts = .fillDBOptsWithDefaults(opts)
           
         }
         )


