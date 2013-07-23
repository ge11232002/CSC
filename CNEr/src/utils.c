#include "CNEr.h"

SEXP bin_from_coord_range(SEXP starts, SEXP ends){
//Return the bin numbers that should be assigned to a feature spanning the given range.
//Here the inputs of starts and ends are 1-based
  starts = AS_INTEGER(starts);
  ends = AS_INTEGER(ends);
  int i, n;
  int *p_start, *p_end, *p_bin;
  n = GET_LENGTH(starts);
  SEXP bins;
  PROTECT(bins = NEW_INTEGER(n));
  p_start = INTEGER_POINTER(starts);
  p_end = INTEGER_POINTER(ends);
  p_bin = INTEGER_POINTER(bins);
  for(i=0; i<n; i++){
    p_bin[i] = binFromRange(p_start[i]-1, p_end[i]-1);
  }
  UNPROTECT(1);
  return bins;
}
