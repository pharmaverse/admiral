#include "vctrs.h"

SEXP (*vec_proxy)(SEXP) = NULL;
SEXP (*vec_restore)(SEXP, SEXP, SEXP) = NULL;
SEXP (*vec_assign_impl)(SEXP, SEXP, SEXP, bool) = NULL;
SEXP (*vec_slice_impl)(SEXP, SEXP) = NULL;
SEXP (*vec_names)(SEXP) = NULL;
SEXP (*vec_set_names)(SEXP, SEXP) = NULL;
SEXP (*vec_chop)(SEXP, SEXP) = NULL;

void vctrs_init_api() {
  vec_proxy = (SEXP (*)(SEXP)) R_GetCCallable("vctrs", "vec_proxy");
  vec_restore = (SEXP (*)(SEXP, SEXP, SEXP)) R_GetCCallable("vctrs", "vec_restore");
  vec_assign_impl = (SEXP (*)(SEXP, SEXP, SEXP, bool)) R_GetCCallable("vctrs", "vec_assign_impl");
  vec_slice_impl = (SEXP (*)(SEXP, SEXP)) R_GetCCallable("vctrs", "vec_slice_impl");
  vec_names = (SEXP (*)(SEXP)) R_GetCCallable("vctrs", "vec_names");
  vec_set_names = (SEXP (*)(SEXP, SEXP)) R_GetCCallable("vctrs", "vec_set_names");
  vec_chop = (SEXP (*)(SEXP, SEXP)) R_GetCCallable("vctrs", "vec_chop");
}
