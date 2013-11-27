#include <R_ext/Rdynload.h>

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {
  /* io.c */
  CALLMETHOD_DEF(myReadBed, 1),

  {NULL, NULL, 0}
};

void R_init_CNEr(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  return;
}

