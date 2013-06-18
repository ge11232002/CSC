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
  SEXP filepath_elt;
  PROTECT(filepath = AS_CHARACTER(filepath));
  if (!IS_CHARACTER(filepath) || LENGTH(filepath) != 1)
    error("'filepath' must be a single string");
  filepath_elt = STRING_ELT(filepath, 0);
  if(filepath_elt == NA_STRING)
    error("'filepath' is NA");
  Rprintf(" %s \n", CHAR(filepath_elt));
  struct hash *hash = readBed(CHAR(filepath_elt));
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




