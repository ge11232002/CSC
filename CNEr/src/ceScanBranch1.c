/* ceScanBranch1.c - scan axt alignment for conserved elements */

#include "R.h"
#include  <ctype.h>  //for the tolower()
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "obscure.h"
#include "options.h"
#include "axt.h"
/********************************************
 *  *** DATA STRUCTURES AND GLOBAL VARIABLES ***
 *   ********************************************/

/* Scoring matrix.
 *  * This will be set by setBpScores() to 1 for matches and 0 for mismatches and gaps. */

#define NR_CHARS 128
typedef int bpScores_t[NR_CHARS][NR_CHARS];
static bpScores_t bpScores;

/* Data structures to represent start and end coordinate pairs.
 * Used to store filters in memory. */

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

/*****************
 *** FUNCTIONS ***
 *****************/
void setBpScores(bpScores_t ss)
/* Set scoring matrix to 1 for matches and 0 for mismatches and gaps. */
{
  unsigned int i, j;
  int a, A;
  char validChars[] = "ACGT";

  // printf("%d\n", (int) sizeof(bpScores_t));

  for (i = 0; i < NR_CHARS; ++i)
    for (j = 0; j < NR_CHARS; ++j)
      ss[i][j] = 0;
  for (i = 0; i < sizeof(validChars)-1; ++i) {
    A = validChars[i];
    a = tolower(A);
    ss[A][A] = ss[a][A] = ss[A][a] = ss[a][a] = 1;
  }
}

struct hash *loadIntHash(char *fileName)
/* Read in a file full of name/number lines into a hash keyed
 * by name with number values. Adapted from axtToMaf.c. */
{
  struct lineFile *lf = lineFileOpen(fileName, TRUE);
  char *row[2];
  struct hash *hash = newHash(0);

  while (lineFileRow(lf, row)) {
    int num = lineFileNeedNum(lf, row, 1);
    hashAddInt(hash, row[0], num);
  }

  lineFileClose(&lf);
  return hash;
}

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


void ceScan(int *a, char **tFilterFile, char **qFilterFile, char **qSizeFile, int *winSize, int *minScore){
  int i = 0, n, minScore, winSize;
  struct slThreshold *trList = NULL, *tr;
  char rest, path[PATH_LEN];
  setBpScores(bpScores);
  for(i; i<*a;i++)
  {
    Rprintf("%d\n", bpScores['A']['A']);
    Rprintf("%s\n", *tFilterFile);
    Rprintf("Hello world!\n");
  }
}



