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
                                   SEXP, SEXP, SEXP, SEXP,SEXP);

extern SEXP _bcSeq_CRISPR_user_matching(SEXP, SEXP, SEXP, SEXP, SEXP,
                                        SEXP, SEXP, SEXP, SEXP, SEXP,
                                        SEXP, SEXP, SEXP, SEXP);

extern SEXP _bcSeq_CRISPR_matching_DNAString(SEXP, SEXP, SEXP, SEXP, SEXP,
                                   SEXP, SEXP, SEXP, SEXP, SEXP,
                                   SEXP, SEXP, SEXP, SEXP,
                                   SEXP, SEXP,SEXP,SEXP);

extern SEXP _bcSeq_CRISPR_user_matching_DNAString(SEXP, SEXP, SEXP, SEXP, SEXP,
                                        SEXP, SEXP, SEXP, SEXP, SEXP,
                                        SEXP, SEXP, SEXP, SEXP, 
                                        SEXP, SEXP,SEXP);
static const R_CallMethodDef CallEntries[] = {
    {"_bcSeq_CRISPR_matching", (DL_FUNC) &_bcSeq_CRISPR_matching, 15},
    {"_bcSeq_CRISPR_matching_DNAString", (DL_FUNC) &_bcSeq_CRISPR_matching_DNAString, 18},
    {"_bcSeq_CRISPR_user_matching", (DL_FUNC) &_bcSeq_CRISPR_user_matching, 14},
    {"_bcSeq_CRISPR_user_matching_DNAString", (DL_FUNC) &_bcSeq_CRISPR_user_matching_DNAString, 17},
    {NULL, NULL, 0}
};

void R_init_bcSeq(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
