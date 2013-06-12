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

void printRangeArray(struct hashEl *hel)
/* Print a range array. For debugging purposes only. */
{
  struct rangeArray *arrayInfo = hel->val;
  struct range *ranges = arrayInfo->ranges;
  int i;
  printf("%s n=%d\n", hel->name, arrayInfo->n);
  for(i = 0; i < arrayInfo->n; i++) {
    printf("%02d: %d - %d\n", i, ranges[i].start, ranges[i].end);
  }
}

struct range *searchRangeArray(struct rangeArray *arrayInfo, int key)
/* Binary search range array. */
{
  struct range *array = arrayInfo->ranges;
  int low = 0;
  int high = arrayInfo->n - 1;
  int mid;

  while(low <= high) {
    mid = (low+high)/2;
    if(key <= array[mid].start) high = mid - 1;
    else if(key > array[mid].end) low = mid + 1;
    else return array+mid; /* return pointer to range that contains key */
  }

  /* key not found: return pointer to nearest higher range or abort if there is no higher range
   * (there should be one because we have added a dummy range with very high values) */
  if(low >= arrayInfo->n) errAbort("searchRangeArray: key %d out of bounds\n", key);
  return array+low;
}

struct hash *readFilter(char *fileName)
/* Load a filter file. */
{
  struct hash *hash = readBed(fileName);
  hashTraverseEls(hash, collapseRangeList);
  hashTraverseEls(hash, convertRangeListToArray);
  /* hashTraverseEls(hash, printRangeArray); */
  return hash;
}

struct hash *makeReversedFilter(struct hash *f1, struct hash *chrSizes)
/* Given a filter, create a reversed filter where coordinates increase in the opposite direction.
 * We use this for filtering alignments that have qStrand == '-'. */
{
  struct hash *f2 = newHash(0);
  struct hashCookie cookie = hashFirst(f1);
  struct hashEl *hel;
  struct rangeArray *fwd, *rev;
  struct range *arrayEl;
  int i, j, n, chrSize;

  /* Iterate over all sequences (chromosomes) in filter */
  while((hel = hashNext(&cookie))) {

    /* get sequence size */
    chrSize = hashIntVal(chrSizes, hel->name);

    /* get forward range array */
    fwd = hel->val;

    /* allocate memory for reversed range array */
    AllocVar(rev);
    n = rev->n = fwd->n; /* set nr of elements in range */
    rev->ranges = arrayEl = needMem(n * sizeof(struct range));

    /* copy dummy range */
    rev->ranges[n-1] = fwd->ranges[n-1];

    /* reverse other ranges */
    for(i = 0, j = n-2; j >= 0; i++, j--) {
      rev->ranges[i].start = chrSize - fwd->ranges[j].end;
      rev->ranges[i].end = chrSize - fwd->ranges[j].start;
    }

    /* add range array to hash keyed by sequence name */
    hashAdd(f2, hel->name, rev);
  }

  /* return reverse filter */
  return f2;
}

struct range *searchFilter(struct hash *filter, char *chrom, int pos)
/* Find the first filter at or following a given position */
{
  struct hashEl *hel;

  hel = hashLookup(filter, chrom);   /* find range array for sequence (chromosome) */
  if(hel) return searchRangeArray(hel->val, pos); /* search range array by position */
  else return NULL;
}

void printCigarString(FILE *fh, struct axt *axt, int i, int j)
/* Print CIGAR string that summarizes alignment */
{
  char type = 'M'; /* in our case first column is always match */
  char newType;
  int count = 0;

  for(; i <= j; i++) {
    /* Determine column type */
    if(axt->tSym[i] == '-') newType = 'D';
    else if(axt->qSym[i] == '-') newType = 'I';
    else newType = 'M';
    /* If same type as previous, just increase count, otherwise output previous */
    if(type == newType) count++;
    else {
      fprintf(fh, "%d%c", count, type);
      type = newType;
      count = 1;
    }
  }

  if(count) fprintf(fh, "%d%c", count, type);
}

void printElement(struct slThreshold *tr, struct axt *axt, struct hash *qSizes, int *profile, int *tPosList, int *qPosList)
/* Print one conserved element on stdout.
 * Arguments:
 * tr - contains threshold-specific information:
 *      parameters used to find CE, CE location in alignment, and filehandle to print to
 * axt - alignment
 * qSizes - query assembly chromosome sizes
 * profile - cumulative conservation profile for alignment
 * tPosList, qPosList - target and query position arrays for alignment
 */
{
  int score, qStart, qEnd, qSize;
  int i = tr->ceStart; /* start column of conserved element in alignment */
  int j = tr->ceEnd; /* end column of conserved element in alignment */

  /* truncate edges (mismatches and gaps) */
  while(bpScores[ (int) axt->qSym[i] ][ (int) axt->tSym[i] ] <= 0) i++;
  while(bpScores[ (int) axt->qSym[j] ][ (int) axt->tSym[j] ] <= 0) j--;

  /* compute score */
  score = profile[j] - profile[i] + bpScores[ (int) axt->qSym[i] ][ (int) axt->tSym[i] ];

  /* recompute query positions if query strand is - */
  if(axt->qStrand == '+') {
    qStart = qPosList[i];
    qEnd = qPosList[j];
  }
  else {
    qSize = hashIntVal(qSizes, axt->qName);
    qStart = qSize - qPosList[j] + 1;
    qEnd = qSize - qPosList[i] + 1;
  }

  /* output */
  fprintf(tr->outFile, "%s\t%d\t%d\t%s\t%d\t%d\t%c\t%.2f\t",
    axt->tName, tPosList[i]-1, tPosList[j],
    axt->qName, qStart-1, qEnd,
    axt->qStrand, 100.0 * score / (j-i+1));
  printCigarString(tr->outFile, axt, i, j);
  fputs("\n", tr->outFile);
}



void ceScan1(char *tFilterFile, char *qFilterFile, char *qSizeFile, struct slThreshold *thresholds, int nrAxtFiles, char *axtFiles[])
/* ceScan - Find conserved elements. */
{
  struct lineFile *lf;
  struct axt *axt;
  struct hash *tFilter, *qFilter, *qFilterRev, *qSizes;
  int i;
  
  setBpScores(bpScores);
  qSizes = loadIntHash(qSizeFile);
  Rprintf("Hello world!\n");
  tFilter = tFilterFile ? readFilter(tFilterFile) : NULL;
  qFilter = qFilterFile ? readFilter(qFilterFile) : NULL;
  qFilterRev = qFilter ? makeReversedFilter(qFilter, qSizes) : NULL;

  i = 0;
  lf = lineFileOpen(axtFiles[i], TRUE);
  axt = axtRead(lf);
  lineFileClose(&lf);
}

void ceScan(char **tFilterFile, char **qFilterFile, char **qSizeFile, int *winSize, int *minScore, int *nThresholds, char **axtFiles, int *nrAxtFiles, char **outFilePrefix){
  int i, n;
  struct slThreshold *trList = NULL, *tr;
  char rest, path[PATH_LEN];
  for(i=1; i<=*nThresholds;i++)
  {
    tr = needMem(sizeof(*tr));
    tr->minScore = *minScore++;
    tr->winSize = *winSize++;
    safef(path, sizeof(path), "%s_%d_%d", outFilePrefix, minScore, winSize);
    tr->outFile = mustOpen(path, "w");
    slAddHead(&trList, tr);
    Rprintf("The winsiez %d\n", tr->winSize);
  }
  //Rprintf("The filter1 is %s\n", *tFilterFile++);
  //Rprintf("The filter1 is %s\n", *tFilterFile);
  /* Call function ceScan with the arguments */
  ceScan1(*tFilterFile, *qFilterFile, *qSizeFile, trList, *nrAxtFiles, axtFiles);
  /* Close all output files */
  //for(tr = trList; tr != NULL; tr = tr->next)
  //  fclose(tr->outFile);

}



