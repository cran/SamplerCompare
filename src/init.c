#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

#include "init.h"

void R_init_SamplerCompare(DllInfo *info) {
  R_CallMethodDef call_routines[] = {
      {"raw_symbol", (DL_FUNC)raw_symbol, 1},
      {"R_invoked_C_glue", (DL_FUNC)R_invoked_C_glue, 4},
      {"sampler_glue_R_dist", (DL_FUNC)sampler_glue_R_dist, 7},
      {"sampler_glue_C_dist", (DL_FUNC)sampler_glue_C_dist, 7},
      {"sr_draw", (DL_FUNC)sr_draw, 5},
      {"arms_sample", (DL_FUNC)arms_sample, 6},
      {"transition_sample", (DL_FUNC)transition_sample, 6},
      {"Gauss2_log_dens", (DL_FUNC)Gauss2_log_dens, 4},
      {"cone_log_dens", (DL_FUNC)cone_log_dens, 4},
      {NULL, NULL, 0}};
  R_NativePrimitiveArgType dchud_arg_types[] = {
      REALSXP, INTSXP,  INTSXP,  REALSXP, REALSXP, INTSXP,
      INTSXP,  REALSXP, REALSXP, REALSXP, REALSXP};
  R_NativePrimitiveArgType dchdd_arg_types[] = {
      REALSXP, INTSXP,  INTSXP,  REALSXP, REALSXP, INTSXP,
      INTSXP,  REALSXP, REALSXP, REALSXP, REALSXP, INTSXP};
  R_FortranMethodDef fortran_routines[] = {
      {"dchud", (DL_FUNC)dchud_, 11, dchud_arg_types},
      {"dchdd", (DL_FUNC)dchdd_, 12, dchdd_arg_types},
      {NULL, NULL, 0}};
  R_registerRoutines(info, NULL /* C routines */, call_routines,
                     fortran_routines, NULL /* external routines */);

  // Allow string lookups because we use them with `raw_symbol` for C
  // distributions.
  R_useDynamicSymbols(info, FALSE);
}
