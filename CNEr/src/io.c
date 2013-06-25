#include "R.h"
#include  <ctype.h>  //for the tolower()
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "obscure.h"
#include "options.h"
#include "axt.h"
#include "Rdefines.h"

struct range
/* Start and end coordinate pair */
{
  int start;    /* Start (0 based) */
  int end;    /* End (non-inclusive) */
};

struct rangeArray
/* Array of start and end coordinate pairs */
{
  int n;
  struct range *ranges;
};

struct slRange
/* Start and end coordinate pair as linked list item */
{
  struct slRange *next;
  int start;    /* Start (0 based) */
  int end;    /* End (non-inclusive) */
};


/* Data structure used to represent different thresholds and intermediate results for each */

struct slThreshold
{
  struct slThreshold *next;
  int minScore;
  int winSize;
  int ceStart;
  int ceEnd;
  FILE *outFile;
};


struct hash *readBed(char *fileName)
/* Read a 3-column bed file into a hash, where keys are sequence names
 * and values are linked lists of coordinate ranges (slRange structures). */
{
  struct lineFile *lf = lineFileOpen(fileName, TRUE);
  struct hash *hash = newHash(0);
  struct hashEl *hel;
  struct slRange *range;
  char *row[3];

  while (lineFileRow(lf, row)) {
    if(sameString(row[0], "track") || sameString(row[0], "browser")) continue;
    AllocVar(range);
    range->next = NULL;
    range->start = lineFileNeedNum(lf, row, 1);
    range->end = lineFileNeedNum(lf, row, 2);
    if (range->start > range->end)
      errAbort("start after end line %d of %s", lf->lineIx, lf->fileName);

    hel = hashLookup(hash, row[0]);
    if (hel == NULL)
      hel = hashAdd(hash, row[0], range);
    else {
      slSafeAddHead(&hel->val, range);
    }
  }

  lineFileClose(&lf);

  return hash;
}

int slRangeCmpStart(const void *va, const void *vb)
/* Comparison function to sort linked list of ranges by start coordinate. */
{
  const struct slRange *a = *((struct slRange **)va);
  const struct slRange *b = *((struct slRange **)vb);
  return a->start - b->start;
}

void collapseRangeList(struct hashEl *hel)
/* Collapse a range list to a set of sorted, non-overlapping ranges. */
{
  struct slRange *a, *b;
  slSort(&hel->val, slRangeCmpStart); /* sort by start coord */
  a = hel->val;
  while((b = a->next)) {
    if(b->start <= a->end) {
      if(a->end < b->end) a->end = b->end;
      a->next = b->next;
      freez(&b);
    }
    else a = b;
  }
  /*for(a = hel->val; a; a=a->next) {
    printf("%d\t%d\n", a->start, a->end);
    }*/
}

void convertRangeListToArray(struct hashEl *hel)
/* Convert a linked list of ranges to an array.
 * The reason for doing this is that we can do a fast binary search on the array. */
{
  struct slRange *list, *slEl;
  struct range *arrayEl;
  struct rangeArray *arrayInfo;
  int n;

  list = hel->val;
  n = slCount(list)+1;
  AllocVar(arrayInfo);
  arrayInfo->n = n;
  arrayInfo->ranges = arrayEl = needMem(n * sizeof(*arrayEl));
  hel->val = arrayInfo;

  while((slEl = slPopHead(&list))) {
    arrayEl->start = slEl->start;
    arrayEl->end = slEl->end;
    free(slEl);
    arrayEl++;
  }

  /* The last array element is a "dummy" element that contains a coordinate pair
   * beyond any chromosome size. The presence of this element simplifies going
   * through the array in scanAxt() as it removes the need for an out-of-bounds check. */
  arrayEl->start = 1e9;
  arrayEl->end = 1e9+1;
}

SEXP readFilter(SEXP filepath)
/* Load a filter file. */
{
  //SEXP filepath_elt;
  PROTECT(filepath = AS_CHARACTER(filepath));
  if (!IS_CHARACTER(filepath) || LENGTH(filepath) != 1)
    error("'filepath' must be a single string");
  if(STRING_ELT(filepath, 0) == NA_STRING)
    error("'filepath' is NA");
  char *filepath_elt = R_alloc(strlen(CHAR(STRING_ELT(filepath, 0))), sizeof(char));
  strcpy(filepath_elt, CHAR(STRING_ELT(filepath, 0)));
  Rprintf(" %s \n", filepath_elt);
  struct hash *hash = readBed(filepath_elt);
  int i, nChroms, nRanges;
  struct hashEl *hel;
  struct slRange *list, *slEl;
  Rprintf("The hash size is %d\n", hash->size);
  Rprintf("The count of elements is %d\n", hash->elCount);
  
  nChroms = 0, nRanges = 0;
  for(i = 0; i < hash->size; ++i){
    for(hel = hash->table[i]; hel != NULL; hel = hel->next){
      list = hel->val;
      Rprintf("The is is %d\n", i);
      Rprintf("The name is %s\n", hel->name);
      while((slPopHead(&list))){
        nRanges++;
      }
    }
  }
  SEXP chromNames, starts, ends, returnList;
  PROTECT(chromNames = NEW_CHARACTER(nRanges));
  PROTECT(starts = NEW_INTEGER(nRanges));
  PROTECT(ends = NEW_INTEGER(nRanges));
  PROTECT(returnList = NEW_LIST(3));
  int *p_starts, *p_ends;
  p_starts = INTEGER_POINTER(starts);
  p_ends = INTEGER_POINTER(ends);
  int j = 0;
  for(i = 0; i < hash->size; ++i){
    for(hel = hash->table[i]; hel != NULL; hel = hel->next){
      list = hel->val;
      while((slEl = slPopHead(&list))){
        SET_STRING_ELT(chromNames, j, mkChar(hel->name));
        p_starts[j] = slEl->start;
        p_ends[j] = slEl->end;
        free(slEl);
        j++;
      }
    }
  }
  //SET_STRING_ELT(chromNames, 0, mkChar("test chrom"));
  Rprintf("The number of Ranges is %d\n", nRanges);
  Rprintf("The chrom name is %s\n", CHAR(STRING_ELT(chromNames, 0)));
  Rprintf("The start is %d\n", p_starts[0]);
  SET_VECTOR_ELT(returnList, 0, chromNames);
  SET_VECTOR_ELT(returnList, 1, starts);
  SET_VECTOR_ELT(returnList, 2, ends);
  UNPROTECT(5);
  return(returnList);
}

/*############################################*/

SEXP myReadBed(SEXP filepath){
  // load a filter file into R, and to be a GRanges
  // SEXP filepath_elt;
  PROTECT(filepath = AS_CHARACTER(filepath));
  if(!IS_CHARACTER(filepath) || LENGTH(filepath) != 1)
    error("'filepath' must be a single string");
  if(STRING_ELT(filepath, 0) == NA_STRING)
    error("'filepath' is NA");
  // If filepath_elt is defined this way, the memory will be reclaimed by the end of .Call by R, do not need the free()
  char *filepath_elt = R_alloc(strlen(CHAR(STRING_ELT(filepath, 0))), sizeof(char));
  strcpy(filepath_elt, CHAR(STRING_ELT(filepath, 0)));
  //filepath_elt = STRING_ELT(filepath, 0);
  //if(filepath_elt == NA_STRING)
  //  error("'filepath' is NA");
  Rprintf(" %s \n", filepath_elt);
  struct lineFile *lf = lineFileOpen(filepath_elt, TRUE);
  char *row[3];
  int nRanges = 0;
  while(lineFileRow(lf, row)){
    if(sameString(row[0], "track") || sameString(row[0], "browser")) continue;
    nRanges++;
  }
  lineFileClose(&lf);
  //Rprintf("The line number is %d\n", nRanges);
  SEXP chromNames, starts, ends, returnList;
  PROTECT(chromNames = NEW_CHARACTER(nRanges));
  PROTECT(starts = NEW_INTEGER(nRanges));
  PROTECT(ends = NEW_INTEGER(nRanges));
  PROTECT(returnList = NEW_LIST(3));
  int *p_starts, *p_ends;
  int j = 0;
  p_starts = INTEGER_POINTER(starts);
  p_ends = INTEGER_POINTER(ends);
  //lf = lineFileOpen(CHAR(filepath_elt), TRUE);
  lf = lineFileOpen(filepath_elt, TRUE);
  while(lineFileRow(lf, row)){
    if(sameString(row[0], "track") || sameString(row[0], "browser")) continue;
    p_starts[j] = lineFileNeedNum(lf, row, 1);
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
  UNPROTECT(5);
  return(returnList);
}

SEXP myReadAxt(SEXP filepath){
  // load a axt file into R, and to be axt2
  //SEXP filepath_elt;
  PROTECT(filepath = AS_CHARACTER(filepath));
  int nrAxtFiles, i, nrAxts;
  nrAxtFiles = GET_LENGTH(filepath);
  Rprintf("The number of axt files %d\n", nrAxtFiles);
  struct axt *axt=NULL, *curAxt;
  struct lineFile *lf;
  nrAxts = 0;
  for(i = 0; i < nrAxtFiles; i++){
    Rprintf("reading the axt file %s\n", CHAR(STRING_ELT(filepath, i)));
    char *filepath_elt = (char *) malloc(sizeof(char) * strlen(CHAR(STRING_ELT(filepath, i))));
    strcpy(filepath_elt, CHAR(STRING_ELT(filepath, i)));
    lf = lineFileOpen(filepath_elt, TRUE);
    //lf = lineFileOpen(CHAR(STRING_ELT(filepath, i)), TRUE);
    //Rprintf("Before reading axt\n");
    //Rprintf("The number of axt is %d\n", nrAxts);
    while((curAxt = axtRead(lf)) != NULL){
      //Rprintf("The name of query sequence is %s\n", curAxt->qName);
      //curAxt->next = axt;
      //axt = curAxt;
      slAddHead(&axt, curAxt);
      nrAxts++;
    }
    lineFileClose(&lf);
    free(filepath_elt);
  }
  axtFree(&curAxt);
  UNPROTECT(1);
  Rprintf("The total number of axt is %d\n", nrAxts);
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
  i = 0;
  char strand;
  while(axt){
    //Rprintf("The name of query seq is %s\n", axt->qName);
    SET_STRING_ELT(qNames, i, mkChar(axt->qName));
    // In Kent's axt struct, they use half open zero=based coordinates. It is different from the original coordinates in axt files. To present it in R, we still make it into 1-based coordinates.
    p_qStart[i] = axt->qStart + 1;
    p_qEnd[i] = axt->qEnd;
    if(axt->qStrand == '+')
      SET_STRING_ELT(qStrand, i, mkChar("+"));
    else
      SET_STRING_ELT(qStrand, i, mkChar("-"));
    //SET_STRING_ELT(qStrand, i, AS_CHARACTER(axt->qStrand));
    SET_STRING_ELT(qSym, i, mkChar(axt->qSym));
    SET_STRING_ELT(tNames, i, mkChar(axt->tName));
    p_tStart[i] = axt->tStart + 1;
    p_tEnd[i] = axt->tEnd;
    if(axt->tStrand == '+')
      SET_STRING_ELT(tStrand, i, mkChar("+"));
    else
      SET_STRING_ELT(tStrand, i, mkChar("-"));
    //SET_STRING_ELT(tStrand, i, mkChar(axt->tStrand));
    SET_STRING_ELT(tSym, i, mkChar(axt->tSym));
    p_score[i] = axt->score;
    p_symCount[i] = axt->symCount;
    i++;
    axt = axt->next;
  }
  SET_VECTOR_ELT(returnList, 0, tNames);
  SET_VECTOR_ELT(returnList, 1, tStart);
  SET_VECTOR_ELT(returnList, 2, tEnd);
  SET_VECTOR_ELT(returnList, 3, tStrand);
  SET_VECTOR_ELT(returnList, 4, tSym);
  SET_VECTOR_ELT(returnList, 5, qNames);
  SET_VECTOR_ELT(returnList, 6, qStart);
  SET_VECTOR_ELT(returnList, 7, qEnd);
  SET_VECTOR_ELT(returnList, 8, qStrand);
  SET_VECTOR_ELT(returnList, 9, qSym);
  SET_VECTOR_ELT(returnList, 10, score);
  SET_VECTOR_ELT(returnList, 11, symCount);
  UNPROTECT(13);
  return(returnList);
}



