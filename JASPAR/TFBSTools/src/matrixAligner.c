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
  // matrixQuery and matrixSubject are matrix of integers.
  // open_penalty and ext_penalty are numerics.
  int vidd1; // column number of matrixQuery
  int vidd2; // column number of matrixSubject
  int i, j;

  vidd1 = INTEGER(GET_DIM(matrixQuery))[1];
  vidd2 = INTEGER(GET_DIM(matrixSubject))[1];
  Rprintf("the dim is %d\n", vidd1);
  
  // for simplicity with old code, but at the cost of assignment
  // make another matrix of profile with additional column so make pos 1 is index 1 in the matrix.
  float matris1[vidd1+1][4];
  float matris2[vidd2+1][4];
  float position_weights[vidd1+1]; // stores the number of sequences in each position
  float position_weights2[vidd2+1]; 
  
  // fill the first row of matris1 with 0 and position_weights with 0
  for(j=0; j<=3; j++){
    matris1[0][j] = 0;
    matris2[0][j] = 0;
  }
  for(i=0; i<=vidd1+1; i++){
    position_weights[i] = 0;
    position_weights2[i] = 0;
  }
  
  // fill in the matrix with these data:
  for(i=1; i<=vidd1; i++){
    for(j=0; j<=3; j++){
      matris1[i][j] = (float)INTEGER(matrixQuery)[(i-1)*4+j];
      position_weights[i] += matris1[i][j];
    }
  }
  for(i=1; i<=vidd2; i++){
    for(j=0; j<=3; j++){
      matris2[i][j] = (float)INTEGER(matrixSubject)[(i-1)*4+j];
      position_weights2[i] += matris2[i][j];
    }
  }

  // normalize the data
  for(j=0; j<=3; j++){
    for(i=1; i<=vidd1; i++){
      matris1[i][j] = matris1[i][j] / position_weights[i];
    }
  }
  Rprintf("the position weight is %f\n", matris1[1][0]);
  Rprintf("the position weight is %f\n", matris1[1][3]);  

  return R_NilValue;

}



