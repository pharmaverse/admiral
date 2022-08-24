#ifndef VCTRS_H
#define VCTRS_H

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdbool.h>

extern SEXP (*vec_proxy)(SEXP);
extern SEXP (*vec_restore)(SEXP, SEXP, SEXP);
extern SEXP (*vec_assign_impl)(SEXP, SEXP, SEXP, bool);
extern SEXP (*vec_slice_impl)(SEXP, SEXP);
extern SEXP (*vec_names)(SEXP);
extern SEXP (*vec_set_names)(SEXP, SEXP);
extern SEXP (*vec_chop)(SEXP, SEXP);

void vctrs_init_api();

#endif
