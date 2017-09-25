#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _bcSeq_CRISPR_matching(SEXP, SEXP, SEXP, SEXP, SEXP,
                                   SEXP, SEXP, SEXP, SEXP, SEXP,
                                   SEXP, SEXP, SEXP, SEXP);

extern SEXP _bcSeq_CRISPR_user_matching(SEXP, SEXP, SEXP, SEXP, SEXP,
                                        SEXP, SEXP, SEXP, SEXP, SEXP,
                                        SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_bcSeq_CRISPR_matching", (DL_FUNC) &_bcSeq_CRISPR_matching, 14},
    {"_bcSeq_CRISPR_user_matching", (DL_FUNC) &_bcSeq_CRISPR_user_matching, 14},
    {NULL, NULL, 0}
};

void R_init_bcSeq(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
