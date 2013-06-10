/* ceScanBranch1.c - scan axt alignment for conserved elements */

#include "R.h"
#include "axt.h"
#include  <ctype.h>  //for the tolower()
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


void ceScan(int *a, char **tFilterFile, char **qFilterFile, char **qSizeFile, int *winSize, int *minScore){
  int i = 0;
  struct slThreshold *trList = NULL, *tr;
  setBpScores(bpScores);
  for(i; i<*a;i++)
  {
    Rprintf("%d\n", bpScores['A']['A']);
    Rprintf("%s\n", *tFilterFile);
    Rprintf("Hello world!\n");
  }
}



