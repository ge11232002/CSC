#include "R.h"
#include "axt.h"
#include "common.h"
#include "linefile.h"
#include "obscure.h"
#include "options.h"
#include "hash.h"

void ceScan(int *a){
  int i;
  for(i=0;i<*a;i++){
    Rprintf("Hello world!\n");
  }
}
