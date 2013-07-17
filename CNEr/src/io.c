#include "CNEr.h"


/*############################################*/

SEXP myReadBed(SEXP filepath){
  // load a filter file into R, and to be a GRanges in 1-based coordinates
  // This is tested and without memory leak!
  filepath = AS_CHARACTER(filepath);
  if(!IS_CHARACTER(filepath) || LENGTH(filepath) != 1)
    error("'filepath' must be a single string");
  if(STRING_ELT(filepath, 0) == NA_STRING)
    error("'filepath' is NA");
  // If filepath_elt is defined this way, the memory will be reclaimed by the end of .Call by R, do not need the free()
  char *filepath_elt = R_alloc(strlen(CHAR(STRING_ELT(filepath, 0))) + 1, sizeof(char));
  strcpy(filepath_elt, CHAR(STRING_ELT(filepath, 0)));
  Rprintf("Reading %s \n", filepath_elt);
  struct lineFile *lf = lineFileOpen(filepath_elt, TRUE);
  char *row[3];
  int nRanges = 0;
  while(lineFileRow(lf, row)){
    if(sameString(row[0], "track") || sameString(row[0], "browser")) continue;
    nRanges++;
  }
  lineFileClose(&lf);
  SEXP chromNames, starts, ends, returnList;
  PROTECT(chromNames = NEW_CHARACTER(nRanges));
  PROTECT(starts = NEW_INTEGER(nRanges));
  PROTECT(ends = NEW_INTEGER(nRanges));
  PROTECT(returnList = NEW_LIST(3));
  int *p_starts, *p_ends;
  int j = 0;
  p_starts = INTEGER_POINTER(starts);
  p_ends = INTEGER_POINTER(ends);
  lf = lineFileOpen(filepath_elt, TRUE);
  while(lineFileRow(lf, row)){
    if(sameString(row[0], "track") || sameString(row[0], "browser")) continue;
    p_starts[j] = lineFileNeedNum(lf, row, 1) + 1;
    p_ends[j] = lineFileNeedNum(lf, row, 2);
    if(p_starts[j] > p_ends[j])
      errAbort("start after end line %d of %s", lf->lineIx, lf->fileName);
    SET_STRING_ELT(chromNames, j, mkChar(row[0]));
    j++;
  }
  lineFileClose(&lf);
  SET_VECTOR_ELT(returnList, 0, chromNames);
  SET_VECTOR_ELT(returnList, 1, starts);
  SET_VECTOR_ELT(returnList, 2, ends);
  UNPROTECT(4);
  return(returnList);
}


/* ------------------------- Use the DNAStringSet way--------------------------*/
SEXP axt_info(SEXP filepath){
  // read a axt file and get the alignment length
  filepath = AS_CHARACTER(filepath);
  int nrAxtFiles, i;
      //nrec0, skip0;
  nrAxtFiles = GET_LENGTH(filepath);
  Rprintf("The number of axt files %d\n", nrAxtFiles);
  struct axt *curAxt;
  struct lineFile *lf;
  IntAE width_buf;
  width_buf = new_IntAE(0, 0, 0);
  char *filepath_elt;
  //nrec = AS_INTEGER(nrec);
  //skip = AS_INTEGER(skip);
  //nrec0 = INTEGER_POINTER(nrec)[0];
  //skip0 = INTEGER_POINTER(skip)[0];
  for(i = 0; i < nrAxtFiles; i++){
    filepath_elt = (char *) R_alloc(strlen(CHAR(STRING_ELT(filepath, i)))+1, sizeof(char));
    strcpy(filepath_elt, CHAR(STRING_ELT(filepath, i)));
    lf = lineFileOpen(filepath_elt, TRUE);
    while((curAxt = axtRead(lf)) != NULL){
      //if(i < skip0){continue;}
      //if(nrec != -1 || i >= nrec){break;}
      IntAE_insert_at(&width_buf, IntAE_get_nelt(&width_buf), curAxt->symCount);
      axtFree(&curAxt);
    }
    lineFileClose(&lf);
  }
  SEXP width;
  PROTECT(width = new_INTEGER_from_IntAE(&width_buf));
  Rprintf("The number of axt alignments is %d\n", GET_LENGTH(width));
  UNPROTECT(1);
  return(width);
}

SEXP readAxt(SEXP filepath){
  // load a axt file into R, and to be axt object
  // This is tested and without memory leak!
  // The jim kent's axt struct holds the starts in 0-based. Here we put it into 1-based.
  filepath = AS_CHARACTER(filepath);
  int nrAxtFiles, i, nrAxts;
  nrAxtFiles = GET_LENGTH(filepath);
  struct axt *axt=NULL, *curAxt;
  struct lineFile *lf;
  IntAE width_buf;
  SEXP ans_qSym, ans_tSym, width;
  PROTECT(width = axt_info(filepath));
  nrAxts = GET_LENGTH(width);
  PROTECT(ans_qSym = alloc_XRawList("BStringSet", "BString", width));
  PROTECT(ans_tSym = alloc_XRawList("BStringSet", "BString", width));
  cachedXVectorList cached_ans_qSym, cached_ans_tSym;
  cachedCharSeq cached_ans_elt;
  cached_ans_qSym = cache_XVectorList(ans_qSym);
  cached_ans_tSym = cache_XVectorList(ans_tSym);
  SEXP qNames, qStart, qEnd, qStrand, qSym, tNames, tStart, tEnd, tStrand, tSym, score, symCount, returnList;
  PROTECT(qNames = NEW_CHARACTER(nrAxts));
  PROTECT(qStart = NEW_INTEGER(nrAxts));
  PROTECT(qEnd = NEW_INTEGER(nrAxts));
  PROTECT(qStrand = NEW_CHARACTER(nrAxts));
  PROTECT(qSym = NEW_CHARACTER(nrAxts));
  PROTECT(tNames = NEW_CHARACTER(nrAxts));
  PROTECT(tStart = NEW_INTEGER(nrAxts));
  PROTECT(tEnd = NEW_INTEGER(nrAxts));
  PROTECT(tStrand = NEW_CHARACTER(nrAxts));
  PROTECT(tSym = NEW_CHARACTER(nrAxts));
  PROTECT(score = NEW_INTEGER(nrAxts));
  PROTECT(symCount = NEW_INTEGER(nrAxts));
  PROTECT(returnList = NEW_LIST(12));
  int *p_qStart, *p_qEnd, *p_tStart, *p_tEnd, *p_score, *p_symCount;
  p_qStart = INTEGER_POINTER(qStart);
  p_qEnd = INTEGER_POINTER(qEnd);
  p_tStart = INTEGER_POINTER(tStart);
  p_tEnd = INTEGER_POINTER(tEnd);
  p_score = INTEGER_POINTER(score);
  p_symCount = INTEGER_POINTER(symCount);
  int j = 0;
  i = 0;
  for(j = 0; j < nrAxtFiles; j++){
    char *filepath_elt = (char *) R_alloc(strlen(CHAR(STRING_ELT(filepath, j)))+1, sizeof(char));
    strcpy(filepath_elt, CHAR(STRING_ELT(filepath, j)));
    lf = lineFileOpen(filepath_elt, TRUE);
    while((axt = axtRead(lf)) != NULL){
      SET_STRING_ELT(qNames, i, mkChar(axt->qName));
      p_qStart[i] = axt->qStart + 1;
      p_qEnd[i] = axt->qEnd;
      if(axt->qStrand == '+')
        SET_STRING_ELT(qStrand, i, mkChar("+"));
      else
        SET_STRING_ELT(qStrand, i, mkChar("-"));
      cached_ans_elt = get_cachedXRawList_elt(&cached_ans_qSym, i);
      memcpy((char *) (&cached_ans_elt)->seq, axt->qSym, axt->symCount * sizeof(char));

      SET_STRING_ELT(tNames, i, mkChar(axt->tName));
      p_tStart[i] = axt->tStart + 1;
      p_tEnd[i] = axt->tEnd;
      if(axt->tStrand == '+')
        SET_STRING_ELT(tStrand, i, mkChar("+"));
      else
        SET_STRING_ELT(tStrand, i, mkChar("-"));
      cached_ans_elt = get_cachedXRawList_elt(&cached_ans_tSym, i);
      memcpy((char *) (&cached_ans_elt)->seq, axt->tSym, axt->symCount * sizeof(char));
      p_score[i] = axt->score;
      p_symCount[i] = axt->symCount;
      i++;
      axtFree(&axt);
    }
    lineFileClose(&lf);
  }
  SET_VECTOR_ELT(returnList, 0, tNames);
  SET_VECTOR_ELT(returnList, 1, tStart);
  SET_VECTOR_ELT(returnList, 2, tEnd);
  SET_VECTOR_ELT(returnList, 3, tStrand);
  SET_VECTOR_ELT(returnList, 4, ans_tSym);
  SET_VECTOR_ELT(returnList, 5, qNames);
  SET_VECTOR_ELT(returnList, 6, qStart);
  SET_VECTOR_ELT(returnList, 7, qEnd);
  SET_VECTOR_ELT(returnList, 8, qStrand);
  SET_VECTOR_ELT(returnList, 9, ans_qSym);
  SET_VECTOR_ELT(returnList, 10, score);
  SET_VECTOR_ELT(returnList, 11, symCount);
  UNPROTECT(16);
  //return R_NilValue;
  return returnList;
}

SEXP axt_info_scratch(SEXP filepath){
  filepath = AS_CHARACTER(filepath);
  if(!IS_CHARACTER(filepath) || LENGTH(filepath) != 1)
    error("'filepath' must be a single string");
  if(STRING_ELT(filepath, 0) == NA_STRING)
    error("'filepath' is NA");
  char *filepath_elt = R_alloc(strlen(CHAR(STRING_ELT(filepath, 0))) + 1, sizeof(char));
  strcpy(filepath_elt, CHAR(STRING_ELT(filepath, 0)));
  Rprintf("Reading %s \n", filepath_elt);
  struct lineFile *lf = lineFileOpen(filepath_elt, TRUE);
  char *row[1], *seq[1];
  int nRanges = 0;
  while(lineFileRow(lf, row)){
    //if(sameString(row[0], "track") || sameString(row[0], "browser")) continue;
    nRanges++;
    lineFileRow(lf, row);
    lineFileRow(lf, row);
  }
  lineFileClose(&lf);
  SEXP starts;
  PROTECT(starts = NEW_INTEGER(nRanges));
  int *p_starts;
  int j = 0;
  p_starts = INTEGER_POINTER(starts);
  lf = lineFileOpen(filepath_elt, TRUE);
  while(lineFileRow(lf, row)){
    //if(sameString(row[0], "track") || sameString(row[0], "browser")) continue;
    //p_starts[j] = lineFileNeedNum(lf, row, 1) + 1;
    p_starts[j] = lineFileNeedNum(lf, row, 0);
    //lineFileRow(lf, row);
    lineFileRow(lf, row);
    lineFileRow(lf, row);
    j++;
  }
  lineFileClose(&lf);
  UNPROTECT(1);
  return(starts);
  /*int nrAxtFiles, i;
  nrAxtFiles = GET_LENGTH(filepath);
  Rprintf("The number of axt files %d\n", nrAxtFiles);
  char str[20000];
  FILE *fp;
  char *filepath_elt;
  SEXP width;
  //PROTECT(width = NEW_INTEGER(246786));
  PROTECT(width = NEW_INTEGER(5344189));
  int *p_width = INTEGER_POINTER(width);
  int j = 0;
  //for(i = 0; i < nrAxtFiles; i++){
    //filepath_elt = (char *) R_alloc(strlen(CHAR(STRING_ELT(filepath, i)))+1, sizeof(char));
    //strcpy(filepath_elt, CHAR(STRING_ELT(filepath, i)));
    //fp = fopen(filepath_elt, "r");
    fp = fopen(CHAR(STRING_ELT(filepath, 0)), "r");
    if(fp == NULL){
      perror("Error opening file");
    }
    int k = 0;
    while(fgets(str, 20000, fp) != NULL){
      //if(k > 50)
      //  break;
      //k++;
      if(str[0] == '#')
        continue;
      //fgets(str, 20000, fp);
      p_width[j] = strlen(str);
      j++;
      //Rprintf("The read length is %d\n", strlen(str));
      //fgets(str, 20000, fp);
      //fgets(str, 20000, fp);
    }
    fclose(fp);
  //}
  UNPROTECT(1);
  //return R_NilValue;
  return(width);*/
}

SEXP axt_info_memory(SEXP filepath){
  // do not use axt read
  filepath = AS_CHARACTER(filepath);
  int nrAxtFiles, i;
  nrAxtFiles = GET_LENGTH(filepath);
  Rprintf("The number of axt files %d\n", nrAxtFiles);
  //IntAE width_buf;
  struct lineFile *lf;
  struct axt *curAxt;
  //width_buf = new_IntAE(0, 0, 0);
  char *words[10], *line;
  char *seqs[1];
  int wordCount;
  int nrAxts = 0;
  char *filepath_elt;
  for(i = 0; i < nrAxtFiles; i++){
    filepath_elt = (char *) R_alloc(strlen(CHAR(STRING_ELT(filepath, i)))+1, sizeof(char));
    strcpy(filepath_elt, CHAR(STRING_ELT(filepath, i)));
    struct lineFile *lf = lineFileOpen(filepath_elt, TRUE);
    while(lineFileChop(lf, words) > 0){
      nrAxts++;
      lineFileNeedNext(lf, &line, NULL);
      //lineFileRow(lf, seqs);
      //Rprintf("The seq is %s\n", seqs[0]);
    //  Rprintf("The seq is %s\n", line);
      //IntAE_insert_at(&width_buf, IntAE_get_nelt(&width_buf),  strlen(line));
      lineFileNeedNext(lf, &line, NULL);
      lineFileNeedNext(lf, &line, NULL);
    }
    /*while((curAxt = axtRead(lf)) != NULL){
      axtFree(&curAxt);
    }*/
    lineFileClose(&lf);
    //free(filepath_elt);
  }
  SEXP width;
  PROTECT(width = NEW_INTEGER(nrAxts));
  int *p_widths;
  p_widths = INTEGER_POINTER(width);
  int j = 0;
  for(i = 0; i < nrAxtFiles; i++){
    filepath_elt = (char *) R_alloc(strlen(CHAR(STRING_ELT(filepath, i)))+1, sizeof(char));
    strcpy(filepath_elt, CHAR(STRING_ELT(filepath, i)));
    lf = lineFileOpen(filepath_elt, TRUE);
    while(lineFileChop(lf, words) > 0){
      p_widths[j] = lineFileNeedNum(lf, words, 2);
      lineFileNeedNext(lf, &line, NULL);
      lineFileNeedNext(lf, &line, NULL);
      lineFileNeedNext(lf, &line, NULL);
      j++;
    }
    lineFileClose(&lf);
    //free(filepath_elt);
  }
  //Rprintf("The number of axt alignments is %d\n", GET_LENGTH(width));
  UNPROTECT(1);
  return width;
  //return R_NilValue;
}

/*void axt_info_dotC(char **filepath, int *width){
  //struct lineFile *lf;
  //struct axt *curAxt;
  //lf = lineFileOpen(filepath[0], TRUE);
  //Rprintf("The filepath is %s\n", filepath[0]);
  int i = 0;
  while((curAxt = axtRead(lf)) != NULL){
    width[i] = curAxt->symCount;
    axtFree(&curAxt);
    i++;
  }
  //lineFileClose(&lf);
  FILE *fp;
  char str[20000];
  fp = fopen(filepath[0], "r");
  while(fgets(str, 20000, fp) != NULL){
    if(str[0] == '#')
      continue;
    fgets(str, 20000, fp);
    width[i] = strlen(str);
    fgets(str, 20000, fp);
    fgets(str, 20000, fp);
    i++;
  }
  fclose(fp);
}*/

/*#define IOBUF_SIZE 20002
static char errmsg_buf[200];

static const char *AXT_comment_markup = "#", AXT_desc_markup = "0";

typedef struct axt_loader {
  const int *lkup;
  int lkup_length;
  void (*load_desc_line)(struct axt_loader *loader,
             const cachedCharSeq *desc_line);
  void (*load_empty_seq)(struct axt_loader *loader);
  void (*load_seq_data)(struct axt_loader *loader,
            const cachedCharSeq *seq_data);
  int nrec;
  void *ext;  // loader extension (optional)
} AXTloader;*/

/*
 * The AXTINFO loader only loads the lengths (and optionally the names)
 * of the sequences.
 */
/*typedef struct axtinfo_loader_ext {
  CharAEAE ans_names_buf;
  IntAE seqlengths_buf;
} AXTINFO_loaderExt;

static void AXTINFO_load_desc_line(AXTloader *loader,
    const cachedCharSeq *desc_line)
{
  AXTINFO_loaderExt *loader_ext;
  CharAEAE *ans_names_buf;

  loader_ext = loader->ext;
  ans_names_buf = &(loader_ext->ans_names_buf);
  // This works only because desc_line->seq is nul-terminated!
  append_string_to_CharAEAE(ans_names_buf, desc_line->seq);
  return;
}

static void AXTINFO_load_empty_seq(AXTloader *loader)
{
  AXTINFO_loaderExt *loader_ext;
  IntAE *seqlengths_buf;

  loader_ext = loader->ext;
  seqlengths_buf = &(loader_ext->seqlengths_buf);
  IntAE_insert_at(seqlengths_buf, IntAE_get_nelt(seqlengths_buf), 0);
  return;
}

static void AXTINFO_load_seq_data(AXTloader *loader,
    const cachedCharSeq *seq_data)
{
  AXTINFO_loaderExt *loader_ext;
  IntAE *seqlengths_buf;

  loader_ext = loader->ext;
  seqlengths_buf = &(loader_ext->seqlengths_buf);
  seqlengths_buf->elts[IntAE_get_nelt(seqlengths_buf) - 1] +=
    seq_data->length;
  return;
}

static AXTINFO_loaderExt new_AXTINFO_loaderExt()
{
  AXTINFO_loaderExt loader_ext;

  loader_ext.ans_names_buf = new_CharAEAE(0, 0);
  loader_ext.seqlengths_buf = new_IntAE(0, 0, 0);
  return loader_ext;
}

static AXTloader new_AXTINFO_loader(SEXP lkup, int load_descs,
    AXTINFO_loaderExt *loader_ext)
{
  AXTloader loader;

  if (lkup == R_NilValue) {
    loader.lkup = NULL;
  } else {
    loader.lkup = INTEGER(lkup);
    loader.lkup_length = LENGTH(lkup);
  }
  loader.load_desc_line = load_descs ? &AXTINFO_load_desc_line : NULL;
  loader.load_empty_seq = &AXTINFO_load_empty_seq;
  loader.load_seq_data = &AXTINFO_load_seq_data;
  loader.nrec = 0;
  loader.ext = loader_ext;
  return loader;
}*/

/*
 * Ignore empty lines and lines starting with 'FASTA_comment_markup' like in
 * the original Pearson FASTA format.
 */
/*static const char *parse_AXT_file(FILE *stream, int *recno, int *ninvalid,
    int nrec, int skip, AXTloader *loader)
{
  
}*/

/* --- .Call ENTRY POINT --- */
/*SEXP axt_info(SEXP efp_list, SEXP nrec, SEXP skip, SEXP use_names, SEXP lkup){
  int nrec0, skip0, load_descs, i, recno, ninvalid;
  AXTINFO_loaderExt loader_ext;
  AXTloader loader;
  FILE *stream;
  SEXP ans, ans_names;
  const char *errmsg;

  nrec0 = INTEGER(nrec)[0];
  skip0 = INTEGER(skip)[0];
  load_descs = LOGICAL(use_names)[0];
  loader_ext = new_AXTINFO_loaderExt();
  loader = new_AXTINFO_loader(lkup, load_descs, &loader_ext);
  recno = 0;
  for(i = 0; i < LENGTH(efp_list); i++){
    stream = R_ExternalPtrAddr(VECTOR_ELT(efp_list, i));
    ninvalid = 0;
    errmsg = parse_AXT_file(stream, &recno, &ninvalid,
        nrec0, skip0, &loader);
  }
}*/



