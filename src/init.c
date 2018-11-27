#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void lasso(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void mult_lasso(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"lasso",      (DL_FUNC) &lasso,      12},
    {"mult_lasso", (DL_FUNC) &mult_lasso, 12},
    {NULL, NULL, 0}
};

void R_init_lasso2(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
