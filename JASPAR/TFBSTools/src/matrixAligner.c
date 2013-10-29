#include "R.h"
#include "Rdefines.h"
#include "Rinternals.h"
#include <string.h>

struct entry
// box in score-table
{
  float score;          // the dynamically best-score-yet
  float cell_score;     // real score in this position
  struct entry *father; //  pointer to father box
  int insert;           // insertion, 1=yes, 0 =no
  int deletion;         // deletion,  1=yes, 0 =no
  int align[2];         // first is matrix_1 pos(i) ,second is second matrix (j)
  int aln_length;       // dynamically extended length of alignment so far
  char kind;
};

struct alignment
{
  float best_score; // the final score
  int length;       // the alignment length
  int gaps;         // number of gaps
  int over_string[30];  // string matrix 1 in alignment (gap represented with -1)
  int under_string[30]; // string matrix2 in alignment
};




/* ----------------.Call() Entry points: the main matrixAligner function ------------- */
SEXP matrixAligner(SEXP matrixQuery, SEXP matrixSubject, SEXP open_penalty, SEXP ext_penalty){
  int vidd1; // column number of matrixQuery
  int vidd2; // column number of matrixSubject

  vidd1 = INTEGER(GET_DIM(matrixQuery))[1];
  vidd2 = INTEGER(GET_DIM(matrixSubject))[1];
  Rprintf("the dim is %d\n", vidd1);
  
  // for simplicity with old code, but at the cost of assignment
  // make another matrix of profile with additional column so make pos 1 is index 1 in the matrix.
  float matris1[vidd1+1][4];
  int position_weights[vidd1+1]; 
  Rprintf("the position weight is %d\n", position_weights[0]); 
  Rprintf("the element is %f\n", matris1[0][0]);
  return R_NilValue;

}



